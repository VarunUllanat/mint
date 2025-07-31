import argparse
import json
import math
from collections import OrderedDict

import pandas as pd
import torch
from torch import nn
from torch.utils.data import Dataset

from ..data import Alphabet
from ..model.esm2 import ESM2


def load_config(path):
    cfg = argparse.Namespace()
    with open(path) as f:
        cfg.__dict__.update(json.load(f))
    return cfg


class CSVDataset(Dataset):
    def __init__(
        self, df_path, col1, col2,
    ):
        super().__init__()
        self.df = pd.read_csv(df_path)
        self.seqs1 = self.df[col1].tolist()
        self.seqs2 = self.df[col2].tolist()

    def __len__(self):
        return len(self.seqs1)

    def __getitem__(self, index):
        return self.seqs1[index], self.seqs2[index]


class CollateFn:
    def __init__(self, truncation_seq_length=None):
        self.alphabet = Alphabet.from_architecture("ESM-1b")
        self.truncation_seq_length = truncation_seq_length
        # self.batch_converter = alphabet.get_batch_converter(truncation_seq_length)

    def __call__(self, batches):
        chains = zip(*batches)
        chains = [self.convert(c) for c in chains]
        chain_ids = [torch.ones(c.shape, dtype=torch.int32) * i for i, c in enumerate(chains)]
        chains = torch.cat(chains, -1)
        chain_ids = torch.cat(chain_ids, -1)
        return chains, chain_ids

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

class MINTWrapper(nn.Module):
    def __init__(
        self,
        cfg,
        checkpoint_path,
        freeze_percent=1.0,
        use_multimer=True,
        sep_chains=False,
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
        total_layers = cfg.encoder_layers
        for name, param in self.model.named_parameters():
            if "embed_tokens.weight" in name or "_norm_after" in name or "lm_head" in name:
                param.requires_grad = False
            else:
                layer_num = name.split(".")[1]
                if int(layer_num) <= math.floor(total_layers * freeze_percent):
                    param.requires_grad = False
        self.model.to(device)

    def get_one_chain(self, chain_out, mask_expanded, mask):
        masked_chain_out = chain_out * mask_expanded
        sum_masked = masked_chain_out.sum(dim=1)
        mask_counts = mask.sum(dim=1, keepdim=True).float()  # Convert to float for division
        mean_chain_out = sum_masked / mask_counts
        return mean_chain_out

    def forward(self, chains, chain_ids):
        # Create a mask to exclude special tokens
        mask = (
            (~chains.eq(self.model.cls_idx)) &
            (~chains.eq(self.model.eos_idx)) &
            (~chains.eq(self.model.padding_idx))
        )
    
        # Get model representations
        chain_out = self.model(chains, chain_ids, repr_layers=[33])["representations"][33]
    
        if self.sep_chains:
            max_chain_id = chain_ids.max().item()
            mean_chain_outs = []
    
            for chain_id in range(max_chain_id + 1):
                chain_mask = (chain_ids == chain_id) & mask
                chain_mask_exp = chain_mask.unsqueeze(-1).expand_as(chain_out)
                mean_chain_out = self.get_one_chain(chain_out, chain_mask_exp, chain_mask)
                mean_chain_outs.append(mean_chain_out)
    
            # Concatenate outputs from each chain along the last dimension
            return torch.cat(mean_chain_outs, dim=-1)
        else:
            mask_expanded = mask.unsqueeze(-1).expand_as(chain_out)
            masked_chain_out = chain_out * mask_expanded
            sum_masked = masked_chain_out.sum(dim=1)
            mask_counts = mask.sum(dim=1, keepdim=True).float()
            mean_chain_out = sum_masked / mask_counts
            return mean_chain_out