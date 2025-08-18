#!/usr/bin/env python3
"""
Split sequence pairs into train/val sets so no sequences ≥40% similar appear in both.
- Prompts for input paths with sensible defaults.
- Fixed 90/10 split by pair-similarity components (derived from MMseqs2 clusters).
- Auto-installs: fire (pip) and mmseqs2 (conda).
"""

import os, re, sys, shutil, tempfile, subprocess as sp
from pathlib import Path
from random import Random

# ---------------- Helpers ----------------

def run(cmd, **kwargs):
    sp.run(cmd, check=True, **kwargs)

def ensure_package(pkg, conda_channels=None):
    """Ensure a Python/conda package is available."""
    try:
        __import__(pkg)
    except ImportError:
        if conda_channels:  # conda package
            print(f"[i] Installing {pkg} via conda...")
            run(["conda", "install", "-y", *sum([["-c", ch] for ch in conda_channels], []), pkg])
        else:
            print(f"[i] Installing {pkg} via pip...")
            run([sys.executable, "-m", "pip", "install", pkg])

# Ensure Fire is available before import
ensure_package("fire")
import fire  # type: ignore

def ensure_mmseqs():
    """Ensure mmseqs2 binary is available on PATH (conda install if missing)."""
    try:
        run(["mmseqs", "--version"], stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        return "mmseqs"
    except Exception:
        ensure_package("mmseqs2", conda_channels=["conda-forge", "bioconda"])
        return "mmseqs"

def detect_real_cpu_count():
    try:
        return int(sp.check_output(["nproc"]).strip())
    except Exception:
        return os.cpu_count() or 1

def normalize_to_oneline_fasta(src: Path, dst: Path):
    """
    Convert input (FASTA or 2-col 'ID SEQ') to one-line-per-seq FASTA.
    Headers are truncated to the first token to match MMseqs2 ID handling.
    """
    # Peek first non-empty line
    first = ""
    with src.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip():
                first = line
                break
    is_fasta = first.startswith(">")

    with src.open("r", encoding="utf-8", errors="ignore") as fin, dst.open("w", encoding="utf-8") as fout:
        if is_fasta:
            header, seq = None, []
            for raw in fin:
                s = raw.rstrip("\n\r")
                if not s: continue
                if s.startswith(">"):
                    if header is not None:
                        fout.write(f">{header}\n{''.join(seq)}\n")
                    header = s[1:].strip().split()[0]  # first token only
                    seq = []
                else:
                    seq.append(re.sub(r"[ \t]", "", s))
            if header is not None:
                fout.write(f">{header}\n{''.join(seq)}\n")
        else:
            # Two-column: each non-empty, non-comment line is "ID SEQ"
            for raw in fin:
                s = raw.strip()
                if not s or s.startswith("#"): continue
                parts = re.split(r"[ \t]+", s, maxsplit=1)
                if len(parts) != 2:
                    sys.exit(f"Bad line (need 'ID SEQ'): {s}")
                pid, seq = parts[0], re.sub(r"[ \t]", "", parts[1])
                fout.write(f">{pid}\n{seq}\n")

def fasta_to_dict(fasta_path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    with fasta_path.open() as f:
        cur = None
        for line in f:
            line = line.rstrip("\n\r")
            if not line: continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = ""
            else:
                seqs[cur] += re.sub(r"[ \t]", "", line)
    return seqs

class DSU:
    def __init__(self, n: int):
        self.p = list(range(n))
    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x
    def union(self, x: int, y: int):
        rx, ry = self.find(x), self.find(y)
        if rx != ry: self.p[ry] = rx

# ---------------- CLI ----------------

class PairSplitter:
    def run(
        self,
        seqs: str | None = None,
        links: str | None = None,
        out_train_seqs: str = "train.seqs.txt",
        out_val_seqs: str   = "val.seqs.txt",
        out_train_links: str = "train.links.txt",
        out_val_links: str   = "val.links.txt",
        seed: int | None = None,
        keep_temp: bool = False,
    ):
        """
        Split sequence pairs into train/val using pair-similarity components built from
        MMseqs2 clusters at ≥40% identity (no leakage across splits).

        Args:
          seqs: path to sequences file; if not provided, you will be prompted
          links: path to links (pairs) file; if not provided, you will be prompted
          out_train_seqs/out_val_seqs: outputs (ID<TAB>SEQ)
          out_train_links/out_val_links: outputs (pairs kept within split)
          seed: optional random seed (component shuffling)
          keep_temp: if True, preserve temp directory
        """
        # Prompt with defaults if not provided
        default_seqs = "training.seqs.txt"
        default_links = "training.links.txt"
        if seqs is None:
            inp = input(f"Path to sequences file [{default_seqs}]: ").strip()
            seqs = inp or default_seqs
        if links is None:
            inp = input(f"Path to links file [{default_links}]: ").strip()
            links = inp or default_links

        seqs_path = Path(seqs)
        links_path = Path(links)
        if not seqs_path.exists():
            sys.exit(f"Missing sequences file: {seqs_path}")
        if not links_path.exists():
            sys.exit(f"Missing links file: {links_path}")

        # Fixed split: 90% pairs → train (by pair count), 10% → val
        split_ratio = 0.9

        # Ensure deps
        mmseqs_bin = ensure_mmseqs()
        threads = detect_real_cpu_count()

        workdir = Path(tempfile.mkdtemp(prefix="mint_split_tmp_"))
        tmpdir = workdir / "mmseqs_tmp"
        tmpdir.mkdir(exist_ok=True)

        try:
            # 0) Normalize to one-line FASTA
            print("[0/9] Normalizing sequences → one-line FASTA")
            fa = workdir / "seqs.fasta"
            normalize_to_oneline_fasta(seqs_path, fa)

            # 1) Build DB
            print(f"[1/9] mmseqs createdb (threads={threads})")
            run([mmseqs_bin, "createdb", str(fa), str(workdir / "seqDB")])

            # 2) Cluster at 40% identity
            print("[2/9] mmseqs linclust (≥40% id)")
            run([
                mmseqs_bin, "linclust", str(workdir / "seqDB"), str(workdir / "clu40"), str(tmpdir),
                "--min-seq-id", "0.4", "--cov-mode", "1", "-c", "0.8",
                "--threads", str(threads)
            ])

            # 3) Export rep→member mapping
            print("[3/9] mmseqs createtsv (rep→member)")
            clu_tsv = workdir / "clu40.tsv"
            run([mmseqs_bin, "createtsv", str(workdir / "seqDB"), str(workdir / "seqDB"),
                 str(workdir / "clu40"), str(clu_tsv)])

            seq2clu: dict[str, str] = {}
            with clu_tsv.open() as f:
                for line in f:
                    s = line.strip()
                    if not s: continue
                    rep, member = s.split("\t")[:2]
                    seq2clu[member] = rep
            num_clusters = len(set(seq2clu.values()))

            # 4) Load pairs; keep only those with both IDs clustered
            print("[4/9] Loading pairs & indexing by clusters")
            pairs: list[tuple[str, str]] = []
            with links_path.open() as fin:
                for line in fin:
                    s = line.strip()
                    if not s or s.startswith("#"): continue
                    a, b, *_ = re.split(r"[ \t]+", s)
                    if a in seq2clu and b in seq2clu:
                        pairs.append((a, b))
            m = len(pairs)
            print(f"     Usable pairs: {m}")

            # 5) Build pair-similarity components via DSU (pairs linked if share a cluster)
            print("[5/9] Building pair-similarity components")
            dsu = DSU(m)
            clu2pairs: dict[str, list[int]] = {}
            for i, (a, b) in enumerate(pairs):
                for sid in (a, b):
                    clu = seq2clu[sid]
                    clu2pairs.setdefault(clu, []).append(i)
            for plist in clu2pairs.values():
                if len(plist) > 1:
                    base = plist[0]
                    for idx in plist[1:]:
                        dsu.union(base, idx)

            comp_map: dict[int, list[int]] = {}
            for i in range(m):
                root = dsu.find(i)
                comp_map.setdefault(root, []).append(i)
            comps = list(comp_map.values())
            print(f"     Components (pair-sim groups): {len(comps)}")

            # 6) Split components to reach ~90% pairs in train
            print("[6/9] Splitting components into train/val (target 90% pairs in train)")
            rng = Random(seed)
            rng.shuffle(comps)
            total_pairs = sum(len(c) for c in comps)
            train_target = int(total_pairs * split_ratio)

            train_idx, val_idx, acc = set(), set(), 0
            for comp in comps:
                if acc + len(comp) <= train_target:
                    train_idx.update(comp)
                    acc += len(comp)
                else:
                    val_idx.update(comp)

            # 7) Write split links
            print("[7/9] Writing links")
            out_train_links_p = Path(out_train_links)
            out_val_links_p   = Path(out_val_links)
            with out_train_links_p.open("w") as ft, out_val_links_p.open("w") as fv:
                for i, (a, b) in enumerate(pairs):
                    if i in train_idx: ft.write(f"{a}\t{b}\n")
                    elif i in val_idx: fv.write(f"{a}\t{b}\n")

            # 8) Write split sequences (only those appearing in each split's pairs)
            print("[8/9] Writing sequences (ID\\tSEQ)")
            seq_dict = fasta_to_dict(fa)
            train_seqs, val_seqs = set(), set()
            for i, (a, b) in enumerate(pairs):
                if i in train_idx:
                    train_seqs.update((a, b))
                elif i in val_idx:
                    val_seqs.update((a, b))
            out_train_seqs_p = Path(out_train_seqs)
            out_val_seqs_p   = Path(out_val_seqs)
            with out_train_seqs_p.open("w") as f:
                for sid in sorted(train_seqs):
                    if sid in seq_dict: f.write(f"{sid}\t{seq_dict[sid]}\n")
            with out_val_seqs_p.open("w") as f:
                for sid in sorted(val_seqs):
                    if sid in seq_dict: f.write(f"{sid}\t{seq_dict[sid]}\n")

            # 9) Stats
            train_pairs = len(train_idx)
            val_pairs   = len(val_idx)
            train_seq_count = sum(1 for _ in out_train_seqs_p.open())
            val_seq_count   = sum(1 for _ in out_val_seqs_p.open())
            print("[9/9] Done ✅")
            print(f"  Total clusters (≥40% id): {num_clusters}")
            print(f"  Components: {len(comps)}")
            print(f"  Pairs  → train: {train_pairs}, val: {val_pairs} (train ratio ~ {train_pairs/(train_pairs+val_pairs):.3f})")
            print(f"  Seqs   → train: {train_seq_count}, val: {val_seq_count}")

        finally:
            if keep_temp:
                print(f"[i] Temp preserved at: {workdir}")
            else:
                shutil.rmtree(workdir, ignore_errors=True)

if __name__ == "__main__":
    fire.Fire(PairSplitter().run)
