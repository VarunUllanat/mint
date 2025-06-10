import argparse
import json
import math
import os
import random
import re
from collections import OrderedDict

import torch
from tasks import get_task_datasets
from torch import nn
from tqdm import tqdm

import mint
from mint.model.esm2 import ESM2

CONFIG_DICT_PATH = (
    "../../data/esm2_t33_650M_UR50D.json"
)


class PPICollateFn:
    def __init__(self, truncation_seq_length=None):
        self.alphabet = mint.data.Alphabet.from_architecture("ESM-1b")
        self.truncation_seq_length = truncation_seq_length

    def __call__(self, batches):
        len(batches)
        heavy_chain, light_chain, labels = zip(*batches)
        chains = [self.convert(c) for c in [heavy_chain, light_chain]]
        chain_ids = [torch.ones(c.shape, dtype=torch.int32) * i for i, c in enumerate(chains)]
        chains = torch.cat(chains, -1)
        chain_ids = torch.cat(chain_ids, -1)
        # labels = torch.from_numpy(np.stack(labels, 0))
        return chains, chain_ids, 0

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


class MutationalPPICollateFn:
    def __init__(self, truncation_seq_length=None):
        self.alphabet = mint.data.Alphabet.from_architecture("ESM-1b")
        self.truncation_seq_length = truncation_seq_length

    def __call__(self, batches):
        len(batches)
        wt_ab, wt_ag, mut_ab, mut_ag, labels = zip(*batches)

        wt_chains = [self.convert(c) for c in [wt_ab, wt_ag]]
        wt_chain_ids = [
            torch.ones(c.shape, dtype=torch.int32) * i for i, c in enumerate(wt_chains)
        ]
        wt_chains = torch.cat(wt_chains, -1)
        wt_chain_ids = torch.cat(wt_chain_ids, -1)

        mut_chains = [self.convert(c) for c in [mut_ab, mut_ag]]
        mut_chain_ids = [
            torch.ones(c.shape, dtype=torch.int32) * i for i, c in enumerate(mut_chains)
        ]
        mut_chains = torch.cat(mut_chains, -1)
        mut_chain_ids = torch.cat(mut_chain_ids, -1)

        # labels = torch.from_numpy(np.stack(labels, 0))
        return wt_chains, wt_chain_ids, mut_chains, mut_chain_ids, 0

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


def get_sequences_by_chain(chain_ids, amino_acids):
    # Initialize list of lists for chains
    chains = []
    for i in range(torch.max(chain_ids) + 1):
        chains.append([])

    # Iterate over both chain ids and amino acids
    for chain_id, amino_acid in zip(chain_ids, amino_acids):
        # Append amino acid to the corresponding chain
        chains[chain_id].append(amino_acid)

    # Convert each sublist into a string to get the sequence of each chain
    return ["".join(chain) for chain in chains]


class PDBBindCollateFn:
    def __init__(self, max_length=1024):
        self.alphabet = mint.data.Alphabet.from_architecture("ESM-1b")
        self.max_length = max_length

    def __call__(self, batches):
        batch_size = len(batches)

        chains, chain_ids_old, labels = zip(*batches)

        batch_aa_chain = [
            get_sequences_by_chain(chain_ids_old[i], chains[i]) for i in range(batch_size)
        ]
        token_chains = torch.empty((batch_size, self.max_length), dtype=torch.int64)
        chain_ids = torch.empty((batch_size, self.max_length), dtype=torch.int64)

        for i in range(batch_size):
            tokens, chain_id = self.convert_single_batch(batch_aa_chain[i])
            token_chains[i, :] = tokens
            chain_ids[i, :] = chain_id

        targets = torch.tensor([labels])
        return token_chains, chain_ids, targets

    def convert_single_batch(self, aa_list):
        num_chains = len(aa_list)
        length_per_chain = self.max_length // num_chains
        seq_encoded_list = [
            self.alphabet.encode("<cls>" + seq_str.replace("J", "L") + "<eos>")
            for seq_str in aa_list
        ]

        tokens = torch.empty((1, self.max_length), dtype=torch.int64)
        tokens.fill_(self.alphabet.padding_idx)

        chain_ids = torch.empty((1, self.max_length), dtype=torch.int64)
        chain_ids.fill_(num_chains - 1)

        for i in range(num_chains):
            seq = seq_encoded_list[i]
            if len(seq) > length_per_chain:
                start = random.randint(0, len(seq) - length_per_chain + 1)
                seq_encoded_list[i] = seq[start : start + length_per_chain]
            start = i * length_per_chain

            tokens[0, start : start + len(seq_encoded_list[i])] = torch.tensor(
                seq_encoded_list[i], dtype=torch.int64
            )
            chain_ids[0, start : start + length_per_chain] = i

        return tokens, chain_ids


def upgrade_state_dict(state_dict):
    prefixes = ["encoder.sentence_encoder.", "encoder."]
    pattern = re.compile("^" + "|".join(prefixes))
    state_dict = {pattern.sub("", name): param for name, param in state_dict.items()}
    return state_dict


class ESMMultimerWrapper(nn.Module):
    def __init__(
        self,
        cfg,
        checkpoint_path,
        freeze_percent=0.0,
        use_multimer=True,
        sep_chains=True,
        device="cuda:0",
    ):
        super().__init__()
        self.cfg = cfg
        self.sep_chains = sep_chains
        self.model = ESM2(
            num_layers=cfg.encoder_layers,
            embed_dim=cfg.encoder_embed_dim,
            attention_heads=cfg.encoder_attention_heads,
            token_dropout=cfg.token_dropout,
            use_multimer=use_multimer,
        )
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)

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
        total_layers = cfg.encoder_layers
        for name, param in self.model.named_parameters():
            if "embed_tokens.weight" in name or "_norm_after" in name or "lm_head" in name:
                param.requires_grad = False
            else:
                layer_num = name.split(".")[1]
                if int(layer_num) <= math.floor(total_layers * freeze_percent):
                    param.requires_grad = False

    def get_one_chain(self, chain_out, mask_expanded, mask):
        masked_chain_out = chain_out * mask_expanded
        sum_masked = masked_chain_out.sum(dim=1)
        mask_counts = mask.sum(dim=1, keepdim=True).float()  # Convert to float for division
        mean_chain_out = sum_masked / mask_counts
        return mean_chain_out

    def forward(self, chains, chain_ids):
        mask = (
            (~chains.eq(self.model.cls_idx))
            & (~chains.eq(self.model.eos_idx))
            & (~chains.eq(self.model.padding_idx))
        )
        chain_out = self.model(chains, chain_ids, repr_layers=[33])["representations"][33]
        if self.sep_chains:
            mask_chain_0 = (chain_ids.eq(0) & mask).unsqueeze(-1).expand_as(chain_out)
            mask_chain_1 = (chain_ids.eq(1) & mask).unsqueeze(-1).expand_as(chain_out)
            mean_chain_out_0 = self.get_one_chain(
                chain_out, mask_chain_0, (chain_ids.eq(0) & mask)
            )
            mean_chain_out_1 = self.get_one_chain(
                chain_out, mask_chain_1, (chain_ids.eq(1) & mask)
            )
            return torch.cat((mean_chain_out_0, mean_chain_out_1), -1)
        else:
            mask_expanded = mask.unsqueeze(-1).expand_as(chain_out)
            masked_chain_out = chain_out * mask_expanded
            sum_masked = masked_chain_out.sum(dim=1)
            mask_counts = mask.sum(dim=1, keepdim=True).float()  # Convert to float for division
            mean_chain_out = sum_masked / mask_counts
            return mean_chain_out


class SKEMPI_ESMMultimerWrapper(ESMMultimerWrapper):
    def forward_one(self, chains, chain_ids):
        mask = (
            (~chains.eq(self.model.cls_idx))
            & (~chains.eq(self.model.eos_idx))
            & (~chains.eq(self.model.padding_idx))
        )
        chain_out = self.model(chains, chain_ids, repr_layers=[33])["representations"][33]
        if self.sep_chains:
            mask_chain_0 = (chain_ids.eq(0) & mask).unsqueeze(-1).expand_as(chain_out)
            mask_chain_1 = (chain_ids.eq(1) & mask).unsqueeze(-1).expand_as(chain_out)
            mean_chain_out_0 = self.get_one_chain(
                chain_out, mask_chain_0, (chain_ids.eq(0) & mask)
            )
            mean_chain_out_1 = self.get_one_chain(
                chain_out, mask_chain_1, (chain_ids.eq(1) & mask)
            )
            return torch.cat((mean_chain_out_0, mean_chain_out_1), -1)
        else:
            mask_expanded = mask.unsqueeze(-1).expand_as(chain_out)
            masked_chain_out = chain_out * mask_expanded
            sum_masked = masked_chain_out.sum(dim=1)
            mask_counts = mask.sum(dim=1, keepdim=True).float()  # Convert to float for division
            mean_chain_out = sum_masked / mask_counts
            return mean_chain_out

    def forward(self, wt_chains, wt_chain_ids, mut_chains, mut_chain_ids, cat=False):
        wt_chains_out = self.forward_one(wt_chains, wt_chain_ids)
        mut_chains_out = self.forward_one(mut_chains, mut_chain_ids)
        if cat:
            return torch.cat((wt_chains_out, mut_chains_out), -1)
        else:
            return wt_chains_out - mut_chains_out


@torch.no_grad()
def get_embeddings(model, loader, device, cat):
    embeddings = []
    for step, eval_batch in enumerate(tqdm(loader)):
        chains, chain_ids, target = eval_batch
        chains = chains.to(device)
        chain_ids = chain_ids.to(device)
        embedding = model(chains, chain_ids)
        embeddings.append(embedding.detach().cpu())
    embeddings = torch.cat(embeddings)
    return embeddings


@torch.no_grad()
def get_embeddings_two(model, loader, device, cat):
    embeddings = []
    for step, eval_batch in enumerate(tqdm(loader)):
        wt_chains, wt_chain_ids, mut_chains, mut_chain_ids, target = eval_batch
        wt_chains = wt_chains.to(device)
        wt_chain_ids = wt_chain_ids.to(device)
        mut_chains = mut_chains.to(device)
        mut_chain_ids = mut_chain_ids.to(device)
        embedding = model(wt_chains, wt_chain_ids, mut_chains, mut_chain_ids, cat)
        embeddings.append(embedding.detach().cpu())
    embeddings = torch.cat(embeddings)
    return embeddings


cfg = argparse.Namespace()
with open(CONFIG_DICT_PATH) as f:
    cfg.__dict__.update(json.load(f))


def main(args):

    devices = [int(s) for s in args.devices.split(",")]
    default_device = f"cuda:{str(devices[0])}"
    gpu_count = len(devices)

    if args.sep_chains:
        save_dir = f"./embeddings/{args.task}/{args.model_name}_sep"
    else:
        save_dir = f"./embeddings/{args.task}/{args.model_name}"

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    print(f"Evaluating on {args.task}")

    if args.task == "SKEMPI" or args.task == "MutationalPPI_cs":
        model = SKEMPI_ESMMultimerWrapper(
            cfg, args.checkpoint_path, 1.0, True, args.sep_chains, default_device
        )
    else:
        model = ESMMultimerWrapper(
            cfg, args.checkpoint_path, 1.0, True, args.sep_chains, default_device
        )

    model = torch.nn.DataParallel(model, device_ids=devices)
    model.to(devices[0])
    model.eval()

    train_dataset, val_dataset, test_dataset = get_task_datasets(args.task, args.test_run)

    if args.task == "SKEMPI" or args.task == "MutationalPPI_cs":
        emb_fn = get_embeddings_two
        train_loader = torch.utils.data.DataLoader(
            train_dataset,
            args.bs * gpu_count,
            collate_fn=MutationalPPICollateFn(args.max_seq_length),
            shuffle=False,
        )
        val_loader = torch.utils.data.DataLoader(
            val_dataset,
            args.bs * gpu_count,
            collate_fn=MutationalPPICollateFn(args.max_seq_length),
            shuffle=False,
        )
        test_loader = torch.utils.data.DataLoader(
            test_dataset,
            args.bs * gpu_count,
            collate_fn=MutationalPPICollateFn(args.max_seq_length),
            shuffle=False,
        )
    elif args.task == "Pdb-bind":
        emb_fn = get_embeddings
        train_loader = torch.utils.data.DataLoader(
            train_dataset,
            args.bs * gpu_count,
            collate_fn=PDBBindCollateFn(args.max_seq_length),
            shuffle=False,
        )
        if test_dataset is not None:
            test_loader = torch.utils.data.DataLoader(
                test_dataset,
                args.bs * gpu_count,
                collate_fn=PDBBindCollateFn(args.max_seq_length),
                shuffle=False,
            )
    else:
        emb_fn = get_embeddings
        train_loader = torch.utils.data.DataLoader(
            train_dataset,
            args.bs * gpu_count,
            collate_fn=PPICollateFn(args.max_seq_length),
            shuffle=False,
        )
        if val_dataset is not None:
            val_loader = torch.utils.data.DataLoader(
                val_dataset,
                args.bs * gpu_count,
                collate_fn=PPICollateFn(args.max_seq_length),
                shuffle=False,
            )
        if test_dataset is not None:
            test_loader = torch.utils.data.DataLoader(
                test_dataset,
                args.bs * gpu_count,
                collate_fn=PPICollateFn(args.max_seq_length),
                shuffle=False,
            )

    train_emb_file_name = f"{save_dir}/train.pt"
    if not os.path.isfile(train_emb_file_name):
        train_inputs = emb_fn(model, train_loader, default_device, args.cat)
        torch.save(train_inputs, train_emb_file_name)

    val_emb_file_name = train_emb_file_name.replace('train', 'val')
    if val_dataset is not None:
        if not os.path.isfile(val_emb_file_name):
            val_inputs = emb_fn(model, val_loader, default_device, args.cat)
            torch.save(val_inputs, val_emb_file_name)

    test_emb_file_name = train_emb_file_name.replace("train", "test")
    if test_dataset is not None:
        if not os.path.isfile(test_emb_file_name):
            test_inputs = emb_fn(model, test_loader, default_device, args.cat)
            torch.save(test_inputs, test_emb_file_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # General args
    parser.add_argument("--task", type=str, default="HumanPPI")
    parser.add_argument("--model_name", type=str, default="mint")
    parser.add_argument("--bs", type=int, default=4)
    parser.add_argument("--max_seq_length", type=int, default=2048)
    parser.add_argument("--devices", type=str, default="0")
    parser.add_argument("--checkpoint_path", type=str, default="")
    parser.add_argument("--sep_chains", action="store_true", default=False)
    parser.add_argument("--test_run", action="store_true", default=False)
    parser.add_argument("--cat", action="store_true", default=False)

    args = parser.parse_args()
    main(args)
