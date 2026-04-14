import argparse
import json
import math
import random
import re
from collections import OrderedDict

import numpy as np
import pandas as pd
import scipy.stats
import torch
import wandb
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from torch import nn
from torch.utils.data import Dataset
from tqdm import tqdm

import mint
from mint.model.esm2 import ESM2


def standardize(x):
    return (x - x.mean(axis=0)) / (x.std(axis=0))


class DesautelsDataset(Dataset):

    score_cols = [
        "FoldX_Average_Whole_Model_DDG",
        "FoldX_Average_Interface_Only_DDG",
        "Statium",
        "Sum_of_Rosetta_Flex_single_point_mutations",
        "Sum_of_Rosetta_Total_Energy_single_point_mutations",
    ]

    def __init__(
        self, csv_path, train_test_split=0.05, split="train", split_seed=2023,
    ):
        super().__init__()
        self.csv_path = csv_path
        self.train_test_split = train_test_split
        assert split in ("train", "test")
        self.split = split
        self.split_seed = split_seed

        full_df = pd.read_csv(self.csv_path)
        full_df_filtered = full_df[~full_df["Statium"].isnull()]

        train_df = full_df_filtered.sample(frac=train_test_split, random_state=split_seed)
        test_df = full_df_filtered.drop(train_df.index)

        if split == "train":
            self.sequences = train_df["Antibody_Sequence"].tolist()
            self.labels = standardize(train_df[self.score_cols].to_numpy())
        elif split == "test":
            test_df = test_df.sample(frac=1.0)
            self.sequences = test_df["Antibody_Sequence"].tolist()
            self.labels = standardize(test_df[self.score_cols].to_numpy())

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, index):

        sequence = self.sequences[index]
        labels = self.labels[index]

        heavy_chain = sequence[:245]
        light_chain = sequence[245:]

        return heavy_chain, light_chain, labels


class DesautelsCollateFn:
    def __init__(self, wt_fasta_file, truncation_seq_length=None):
        self.alphabet = mint.data.Alphabet.from_architecture("ESM-1b")
        self.truncation_seq_length = truncation_seq_length

        with open(wt_fasta_file) as file:
            for i, line in enumerate(file):
                if i == 1:
                    wt_heavy_chain = line.strip()
                elif i == 3:
                    wt_light_chain = line.strip()

        self.wt_chains_single = [self.convert(c) for c in [[wt_heavy_chain], [wt_light_chain]]]
        self.wt_chain_ids_single = [
            torch.ones(c.shape, dtype=torch.int32) * i for i, c in enumerate(self.wt_chains_single)
        ]

    def __call__(self, batches):
        batch_size = len(batches)
        heavy_chain, light_chain, labels = zip(*batches)

        chains = [self.convert(c) for c in [heavy_chain, light_chain]]
        chain_ids = [torch.ones(c.shape, dtype=torch.int32) * i for i, c in enumerate(chains)]
        chains = torch.cat(chains, -1)
        chain_ids = torch.cat(chain_ids, -1)
        labels = torch.from_numpy(np.stack(labels, 0))

        wt_chains = torch.cat(self.wt_chains_single, -1).repeat(batch_size, 1)
        wt_chain_ids = torch.cat(self.wt_chain_ids_single, -1).repeat(batch_size, 1)

        return chains, chain_ids, wt_chains, wt_chain_ids, labels

    def convert(self, seq_str_list):
        batch_size = len(seq_str_list)
        seq_encoded_list = [
            self.alphabet.encode("<cls>" + seq_str.replace("J", "L") + "<eos>")
            for seq_str in seq_str_list
        ]
        if self.truncation_seq_length:
            for i in range(batch_size):
                seq = seq_encoded_list[i]
                if len(seq) > self.truncation_seq_length:
                    start = random.randint(0, len(seq) - self.truncation_seq_length + 1)
                    seq_encoded_list[i] = seq[start : start + self.truncation_seq_length]
        max_len = max(len(seq_encoded) for seq_encoded in seq_encoded_list)
        if self.truncation_seq_length:
            assert max_len <= self.truncation_seq_length
        tokens = torch.empty((batch_size, max_len), dtype=torch.int64)
        tokens.fill_(self.alphabet.padding_idx)

        for i, seq_encoded in enumerate(seq_encoded_list):
            seq = torch.tensor(seq_encoded, dtype=torch.int64)
            tokens[i, : len(seq_encoded)] = seq
        return tokens


def upgrade_state_dict(state_dict):
    """Removes prefixes 'model.encoder.sentence_encoder.' and 'model.encoder.'."""
    prefixes = ["encoder.sentence_encoder.", "encoder."]
    pattern = re.compile("^" + "|".join(prefixes))
    state_dict = {pattern.sub("", name): param for name, param in state_dict.items()}
    return state_dict


class DesautelsWrapper(nn.Module):
    def __init__(
        self,
        cfg,
        checkpoint_path,
        freeze_percent=0.0,
        use_multimer=True,
        device="cuda:0",
        finetune=True,
    ):
        super().__init__()
        self.cfg = cfg
        self.finetune = finetune
        self.model = ESM2(
            num_layers=cfg.encoder_layers,
            embed_dim=cfg.encoder_embed_dim,
            attention_heads=cfg.encoder_attention_heads,
            token_dropout=cfg.token_dropout,
            use_multimer=use_multimer,
        )
        checkpoint = torch.load(checkpoint_path, map_location=device)
        if use_multimer:
            # remove 'model.' in keys
            new_checkpoint = OrderedDict(
                (key.replace("model.", ""), value)
                for key, value in checkpoint["state_dict"].items()
            )
            self.model.load_state_dict(new_checkpoint)
        else:
            new_checkpoint = upgrade_state_dict(checkpoint["model"])
            self.model.load_state_dict(new_checkpoint)
        total_layers = 33
        for name, param in self.model.named_parameters():
            if "embed_tokens.weight" in name or "_norm_after" in name or "lm_head" in name:
                param.requires_grad = False
            else:
                layer_num = name.split(".")[1]
                if int(layer_num) <= math.floor(total_layers * freeze_percent):
                    param.requires_grad = False

        if self.finetune:
            self.project = nn.Sequential(
                nn.Linear(1280, 64), nn.SiLU(), nn.Dropout(0.2), nn.Linear(64, 5),
            )

    def forward_one(self, chains, chain_ids):
        mask = (
            (~chains.eq(self.model.cls_idx))
            & (~chains.eq(self.model.eos_idx))
            & (~chains.eq(self.model.padding_idx))
        )
        chain_out = self.model(chains, chain_ids, repr_layers=[33])["representations"][33]
        mask_expanded = mask.unsqueeze(-1).expand_as(chain_out)
        masked_chain_out = chain_out * mask_expanded
        sum_masked = masked_chain_out.sum(dim=1)
        mask_counts = mask.sum(dim=1, keepdim=True).float()  # Convert to float for division
        # mask_counts = mask_counts.where(mask_counts != 0, torch.ones_like(mask_counts))
        mean_chain_out = sum_masked / mask_counts
        return mean_chain_out

    def forward(self, chains, chain_ids, wt_chains, wt_chain_ids):
        chain_out = self.forward_one(chains, chain_ids)
        wt_chain_out = self.forward_one(wt_chains, wt_chain_ids)
        if self.finetune:
            out = self.project(wt_chain_out - chain_out)
        else:
            out = wt_chain_out - chain_out
        return out

    # def forward(self, chains, chain_ids, wt_chains, wt_chain_ids):
    #     chain_out = self.model(chains, chain_ids, repr_layers=[33])["representations"][33].mean(1)
    #     wt_chain_out = self.model(wt_chains, wt_chain_ids, repr_layers=[33])["representations"][33].mean(1)
    #     out = self.project(wt_chain_out - chain_out)
    #     return out


@torch.no_grad()
def evaluate(model, loader, args):

    device = args.device

    pred = []
    targets = []

    for step, eval_batch in enumerate(tqdm(loader)):

        chains, chain_ids, wt_chains, wt_chain_ids, target = eval_batch
        chains = chains.to(device)
        wt_chains = wt_chains.to(device)
        chain_ids = chain_ids.to(device)
        wt_chain_ids = wt_chain_ids.to(device)
        target = target.to(device).float()

        pred_ddg = model(chains, chain_ids, wt_chains, wt_chain_ids)

        pred.append(pred_ddg.detach().cpu().numpy())
        targets.append(target.cpu().numpy())

    pred = np.concatenate(pred)
    targets = np.concatenate(targets)

    mses = []
    pearsons = []
    spearmans = []
    for i in range(pred.shape[1]):
        mse = mean_squared_error(pred[:, i], targets[:, i], squared=False)
        pearson = scipy.stats.pearsonr(pred[:, i], targets[:, i])[0]
        spearman = scipy.stats.spearmanr(pred[:, i], targets[:, i])[0]

        mses.append(mse)
        pearsons.append(pearson)
        spearmans.append(spearman)

    return np.mean(mses), np.mean(pearsons), np.mean(spearmans)


@torch.no_grad()
def get_embeddings(model, loader, device="cuda"):
    model.to(device)
    embeddings = []
    targets = []
    for step, eval_batch in enumerate(tqdm(loader)):
        chains, chain_ids, wt_chains, wt_chain_ids, target = eval_batch
        chains = chains.to(device)
        wt_chains = wt_chains.to(device)
        chain_ids = chain_ids.to(device)
        wt_chain_ids = wt_chain_ids.to(device)
        target = target.to(device).float()

        embedding = model(chains, chain_ids, wt_chains, wt_chain_ids)
        embeddings.append(embedding.detach().cpu().numpy())
        targets.append(target.cpu().numpy())
    embeddings = np.concatenate(embeddings)
    targets = np.concatenate(targets)
    return embeddings, targets


def train(model, train_loader, val_loader, cfg, args):
    device = args.device

    optimizer = torch.optim.AdamW(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=cfg.lr[0],
        betas=json.loads(cfg.adam_betas),
        eps=cfg.adam_eps,
        weight_decay=cfg.weight_decay,
    )
    loss_fn = torch.nn.MSELoss()
    model.to(device)

    for epoch in range(args.num_epochs):
        print(f"Training at epoch {epoch}")
        loss_accum = 0
        for step, train_batch in enumerate(tqdm(train_loader)):

            model.train()
            optimizer.zero_grad()

            chains, chain_ids, wt_chains, wt_chain_ids, target = train_batch
            chains = chains.to(device)
            wt_chains = wt_chains.to(device)
            chain_ids = chain_ids.to(device)
            wt_chain_ids = wt_chain_ids.to(device)
            target = target.to(device).float()

            pred_ddg = model(chains, chain_ids, wt_chains, wt_chain_ids)
            loss = loss_fn(pred_ddg, target)

            loss.backward()
            optimizer.step()
            loss_accum += loss.detach().cpu().item()
        print(f"Loss at end of epoch {epoch}: {loss_accum/(step+1)}")

        if epoch == args.num_epochs - 1:
            print(f"Evaluating at epoch {epoch}")
            mse, pearson, spearman = evaluate(model, val_loader, args)
            print(mse, pearson, spearman)
            wandb.log(
                {
                    "train_loss": loss_accum / (step + 1),
                    "mse": mse,
                    "pearson": pearson,
                    "spearman": spearman,
                },
                step=epoch,
            )
        else:
            if args.wandb:
                wandb.log({"train_loss": loss_accum / (step + 1)}, step=epoch)


def main(args):

    for train_test_split in [0.005, 0.02, 0.2]:

        train_dataset = DesautelsDataset(
            "Desautels_insilico_data.csv", train_test_split=train_test_split
        )
        test_dataset = DesautelsDataset(
            "Desautels_insilico_data.csv", train_test_split=train_test_split, split="test"
        )

        train_loader = torch.utils.data.DataLoader(
            train_dataset,
            batch_size=args.batch_size,
            collate_fn=DesautelsCollateFn("rcsb_pdb_2G75.fasta"),
            shuffle=True,
        )
        test_loader = torch.utils.data.DataLoader(
            test_dataset,
            batch_size=args.batch_size * 2,
            collate_fn=DesautelsCollateFn("rcsb_pdb_2G75.fasta"),
            shuffle=False,
        )

        cfg = argparse.Namespace()
        with open(
            f"/data/cb/scratch/varun/mint/mint/models/esm2_t33_650M_UR50D.json"
        ) as f:
            cfg.__dict__.update(json.load(f))

        if args.use_mlp:
            model = DesautelsWrapper(
                cfg, args.checkpoint_path, args.freeze_percent, args.use_multimer, args.device
            )
            train(model, train_loader, test_loader, cfg, args)
        else:
            model = DesautelsWrapper(
                cfg,
                args.checkpoint_path,
                args.freeze_percent,
                args.use_multimer,
                args.device,
                finetune=False,
            )
            train_embeddings, train_targets = get_embeddings(
                model, train_loader, device=args.device
            )
            test_embeddings, test_targets = get_embeddings(model, test_loader, device=args.device)

            model = Ridge(alpha=0.01)
            model.fit(train_embeddings, train_targets)
            Y_pred = model.predict(test_embeddings)

            mses = []
            pearsons = []
            spearmans = []
            for i in range(Y_pred.shape[1]):
                mse = mean_squared_error(test_targets[:, i], Y_pred[:, i], squared=False)
                pearson = scipy.stats.pearsonr(test_targets[:, i], Y_pred[:, i])[0]
                spearman = scipy.stats.spearmanr(test_targets[:, i], Y_pred[:, i])[0]

                mses.append(mse)
                pearsons.append(pearson)
                spearmans.append(spearman)

            if args.wandb:
                wandb.log(
                    {
                        "pearson": np.mean(pearsons),
                        "mse": np.mean(mses),
                        "spearman": np.mean(spearmans),
                    }
                )

        if args.wandb:
            wandb.finish()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Finetuning on desautels dataset")

    # General args
    parser.add_argument("--batch_size", type=int, default=12)
    parser.add_argument("--num_epochs", type=int, default=20)
    parser.add_argument("--wandb", action="store_true", default=False)
    parser.add_argument("--wandb_key", type=str, default=None)
    parser.add_argument(
        "--checkpoint_path",
        type=str
    )
    parser.add_argument("--freeze_percent", type=float, default=1.0)
    parser.add_argument("--use_multimer", action="store_true", default=False)
    parser.add_argument("--device", type=str, default="cuda:0")
    parser.add_argument("--use_mlp", action="store_true", default=False)

    args = parser.parse_args()
    main(args)
