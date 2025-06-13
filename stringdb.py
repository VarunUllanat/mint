import gzip
import random

from Bio import SeqIO
from mint.utils.logging import get_logger

logger = get_logger(__name__)

logger.info("===Reading seqs====")
seqs = {}
for seq in SeqIO.parse(open("../protein.sequences.v12.0.fa"), "fasta"):
    seqs[seq.name] = str(seq.seq)
    if len(seqs) % 1000000 == 0:
        logger.info(f"{len(seqs) / 1e6} million seqs read")
logger.info(f"Done, {len(seqs)} seqs total")

logger.info("===Reading reps====")
reps = {}
for line in open("../clu50.tsv"):
    rep, seq = line.strip().split()
    reps[seq] = rep
    if len(reps) % 1000000 == 0:
        logger.info(f"{len(reps) / 1e6} million reps read")
logger.info(f"Done, {len(reps)} reps total, {len(set(reps.values()))} clusters")

logger.info("===Reading links====")
f = gzip.open("../protein.physical.links.full.v12.0.txt.gz", "rt")
f = iter(f)
next(f)  # skip first line
links = []
i = 0
while True:
    try:
        line = next(f).strip()
    except StopIteration:
        break
    i += 1
    links.append(line)
    if i % 1000000 == 0:
        logger.info(f"{i / 1e6} million links read")
        # if i / 1e6 == 10: break
logger.info(f"Done, {len(links)} links total")

logger.info("===Shuffling links===")
random.seed(137)
random.shuffle(links)
logger.info("Done shuffling links")


logger.info("===Filtering links===")
linked_clusters = set()
filtered_links = []
i = 0
for link in links:
    i += 1
    name1, name2 = link.split()[:2]
    clu1, clu2 = reps[name1], reps[name2]
    clu1, clu2 = tuple(sorted((clu1, clu2)))
    if (clu1, clu2) not in linked_clusters:
        linked_clusters.add((clu1, clu2))
        filtered_links.append(link)
    if i % 1000000 == 0:
        logger.info(f"{i / 1e6} million links filtered, {len(filtered_links) / 1e6} million kept")

links = filtered_links
logger.info(f"Done, {i} links filtered, {len(links)} kept")

logger.info("===Shuffling links===")
random.seed(731)
random.shuffle(links)
logger.info("Done shuffling links")

with gzip.open("../filtered.links.txt.gz", "wt") as links_file:
    for link in links:
        links_file.write(link + "\n")

num_val = 250000
validation = links[:num_val]
training = links[num_val:]

logger.info("===Writing validation links===")
written_seqs = set()
with gzip.open("validation.links.txt.gz", "wt") as links_file:
    with gzip.open("validation.seqs.txt.gz", "wt") as seqs_file:
        for link in validation:
            links_file.write(link + "\n")
            name1, name2 = link.split()[:2]
            if name1 not in written_seqs:
                seqs_file.write(name1 + " " + seqs[name1] + "\n")
                written_seqs.add(name1)
            if name2 not in written_seqs:
                seqs_file.write(name2 + " " + seqs[name2] + "\n")
                written_seqs.add(name2)
logger.info(f"Done, {num_val} validation links written, {len(written_seqs)} seqs")

logger.info("===Writing training links===")
i = 0
written_seqs = set()
with gzip.open("training.links.txt.gz", "wt") as links_file:
    with gzip.open("training.seqs.txt.gz", "wt") as seqs_file:
        for link in training:
            i += 1
            links_file.write(link + "\n")
            name1, name2 = link.split()[:2]
            if name1 not in written_seqs:
                seqs_file.write(name1 + " " + seqs[name1] + "\n")
                written_seqs.add(name1)
            if name2 not in written_seqs:
                seqs_file.write(name2 + " " + seqs[name2] + "\n")
                written_seqs.add(name2)
            if i % 1000000 == 0:
                logger.info(
                    f"{i / 1e6} million training links written, {len(written_seqs) / 1e6} million seqs"
                )
logger.info(f"Done, {i} training links written, {len(written_seqs)} seqs")

logger.info("===Extracting validation clusters===")
val_clus = []
for link in validation:
    name1, name2 = link.split()[:2]
    val_clus.append(reps[name1])
    val_clus.append(reps[name2])
val_clus = set(val_clus)
logger.info(f"Done, {len(val_clus)} validation clusters")

i, j = 0, 0
logger.info("===Writing filtered training links===")
written_seqs = set()
with gzip.open("training_filtered.links.txt.gz", "wt") as links_file:
    with gzip.open("training_filtered.seqs.txt.gz", "wt") as seqs_file:
        for link in training:
            i += 1
            name1, name2 = link.split()[:2]
            clu1, clu2 = reps[name1], reps[name2]
            if clu1 not in val_clus and clu2 not in val_clus:
                j += 1
                links_file.write(link + "\n")
                if name1 not in written_seqs:
                    seqs_file.write(name1 + " " + seqs[name1] + "\n")
                    written_seqs.add(name1)
                if name2 not in written_seqs:
                    seqs_file.write(name2 + " " + seqs[name2] + "\n")
                    written_seqs.add(name2)
            if i % 1000000 == 0:
                logger.info(
                    f"{i / 1e6} million training links filtered, {j / 1e6} million written, {len(written_seqs) / 1e6} million seqs"
                )
logger.info(f"{i} training links filtered, {j} kept, {len(written_seqs)} seqs")
