import gzip
import random

import torch
from Bio import SeqIO

import mint


class CollateFn:
    def __init__(self, truncation_seq_length=None):
        self.alphabet = mint.data.Alphabet.from_architecture("ESM-1b")
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


class STRINGDataset(torch.utils.data.IterableDataset):
    def __init__(
        self,
        links_path,
        seqs_path,
        global_rank=0,
        world_size=1,
        concat=False,
        max_examples=None,
        max_len=None,
        overfit=False,
        seek=None,
    ):
        super().__init__()
        self.links_path = links_path
        self.seqs_path = seqs_path
        self.global_rank = global_rank
        self.world_size = world_size

        if max_examples:
            self.max_iters = int(max_examples // world_size)
        else:
            self.max_iters = None
        self.concat = concat
        self.max_len = max_len
        self.overfit = overfit
        self.seek = seek

    def __len__(self):
        return self.max_iters

    def __iter__(self):
        it = self.__iter_helper__()
        for i, n in enumerate(it):
            if self.seek and i < self.seek:
                pass
            else:
                yield n

    def __iter_helper__(self):
        self.seqs = {}
        links_f = iter(gzip.open(self.links_path, "rt"))
        seqs_f = iter(gzip.open(self.seqs_path, "rt"))
        i, j = 0, 0
        while True:
            try:
                name1, name2 = next(links_f).strip().split()[:2]
                if name1 not in self.seqs:
                    name, seq = next(seqs_f).strip().split()
                    self.seqs[name] = seq
                if name2 not in self.seqs:
                    name, seq = next(seqs_f).strip().split()
                    self.seqs[name] = seq
                if self.max_len and not (
                    len(self.seqs[name1]) <= self.max_len and len(self.seqs[name2]) <= self.max_len
                ):
                    continue  # don't increment i
                if i % self.world_size == self.global_rank:
                    if self.concat:
                        if self.overfit:
                            while True:
                                yield (self.seqs[name1] + "G" * 25 + self.seqs[name2],)
                        else:
                            yield (self.seqs[name1] + "G" * 25 + self.seqs[name2],)
                    else:
                        if self.overfit:
                            while True:
                                yield self.seqs[name1], self.seqs[name2]
                        else:
                            yield self.seqs[name1], self.seqs[name2]
                    j += 1
                i += 1
                if j == self.max_iters:
                    break
            except StopIteration:
                links_f = iter(gzip.open(self.links_path, "rt"))


class ESMDataset(torch.utils.data.IterableDataset):
    def __init__(self, path):
        super().__init__()
        self.seqs = SeqIO.parse(open(path), "fasta")

    def __iter__(self):
        seqs = iter(self.seqs)
        while True:
            seq = next(seqs)
            yield seq.name, seq.seq


class DummyDataset(torch.utils.data.Dataset):
    def __init__(self):
        super().__init__()
        self.alphabet = esm.data.Alphabet.from_architecture("ESM-1b")

    def __len__(self):
        return 1024

    def __getitem__(self, idx):
        return torch.randint(4, 24, size=(1024,))
        seq = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        return self.alphabet.get_batch_converter()([("P68871", seq)])[2].squeeze(0)
