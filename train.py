from mint.utils.parsing import parse_train_args

args = parse_train_args()

import argparse
import json
import re

import lightning as pl
import torch
from lightning.pytorch.callbacks import ModelCheckpoint
from lightning.pytorch.strategies import DDPStrategy

from mint.utils.dataset import CollateFn, STRINGDataset
from mint.utils.wrapper import ESMWrapper

torch.set_float32_matmul_precision("medium")


def upgrade_state_dict(state_dict):
    """Removes prefixes 'model.encoder.sentence_encoder.' and 'model.encoder.'."""
    prefixes = ["encoder.sentence_encoder.", "encoder."]
    pattern = re.compile("^" + "|".join(prefixes))
    state_dict = {pattern.sub("", name): param for name, param in state_dict.items()}
    return state_dict


trainer = pl.Trainer(
    default_root_dir=f"./checkpoints/{args.run_name}",
    accelerator="gpu",
    devices=[0, 1, 2, 3, 4, 5, 6, 7],
    max_steps=args.max_steps,
    num_sanity_val_steps=2,
    enable_progress_bar=not args.wandb,
    gradient_clip_val=args.grad_clip,
    enable_checkpointing=True,
    callbacks=[ModelCheckpoint(dirpath=f"./checkpoints/{args.run_name}", save_top_k=-1,)],
    accumulate_grad_batches=args.accumulate_grad,
    val_check_interval=args.val_check_interval,
    strategy=DDPStrategy(find_unused_parameters=True)
    if args.freeze_self_attn
    else "ddp_find_unused_parameters_false",
)

if args.dataset_split in ["filtered", "full"]:
    val_links_file = "../validation.links.txt.gz"
    val_seqs_file = "../validation.seqs.txt.gz"
    if args.dataset_split == "filtered":
        train_links_file = "../training_filtered.links.txt.gz"
        train_seqs_file = "../training_filtered.seqs.txt.gz"
    else:
        pass
elif args.dataset_split == "filtered_50":
    val_links_file = "../validation.links.50.txt.gz"
    val_seqs_file = "../validation.seqs.50.txt.gz"
    train_links_file = "../training.links.50.txt.gz"
    train_seqs_file = "../training.seqs.50.txt.gz"

val_ds = STRINGDataset(
    val_links_file,
    val_seqs_file,
    global_rank=trainer.global_rank,
    world_size=trainer.world_size,
    max_examples=args.val_examples,
    concat=args.concat,
    max_len=args.val_max_len,
)

val_loader = torch.utils.data.DataLoader(
    val_ds, batch_size=args.batch_size, collate_fn=CollateFn(args.crop_length)
)

train_ds = STRINGDataset(
    train_links_file,
    train_seqs_file,
    global_rank=trainer.global_rank,
    world_size=trainer.world_size,
    concat=args.concat,
    overfit=args.overfit,
    seek=args.dataset_seek,
)
train_loader = torch.utils.data.DataLoader(
    train_ds, batch_size=args.batch_size, collate_fn=CollateFn(args.crop_length)
)

model_name = {
    "8M": "esm2_t6_8M_UR50D",
    "35M": "esm2_t12_35M_UR50D",
    "150M": "esm2_t30_150M_UR50D",
    "650M": "esm2_t33_650M_UR50D",
    "3B": "esm2_t36_3B_UR50D",
    "15B": "esm2_t48_15B_UR50D",
}[args.model]

cfg = argparse.Namespace()
with open(f"models/{model_name}.json") as f:
    cfg.__dict__.update(json.load(f))

model = ESMWrapper(cfg, args)

if not (args.ckpt or args.reinitialize):
    state_dict = torch.load(f"models/{model_name}.pt")["model"]
    model.model.load_state_dict(upgrade_state_dict(state_dict), strict=False)

if (not args.no_multimer) and args.copy_weights:
    for layer in model.model.layers:
        layer.multimer_attn.load_state_dict(layer.self_attn.state_dict(), strict=False)

if args.validate:
    if args.ckpt:
        ckpt = torch.load(args.ckpt, weights_only=False)
        model.load_state_dict(ckpt["state_dict"], strict=False)
    trainer.validate(model, val_loader)
else:
    trainer.fit(model, train_loader, val_loader, ckpt_path=args.ckpt)

# reinit
# python train.py --batch_size 2 --crop_len 512 --model 650M --val_check_interval 320000 --reinitialize --accumulate_grad 32 --run_name 650M_reinit_filtered --wandb --dataset_split filtered

# freeze
# python train.py --batch_size 2 --crop_len 512 --model 650M --val_check_interval 320000 --copy_weights --accumulate_grad 32 --freeze_self_attn --run_name 650M_freeze_filtered --wandb --dataset_split filtered

# nofreeze
# python train.py --batch_size 2 --crop_len 512 --model 650M --val_check_interval 320000 --accumulate_grad 32 --run_name 650M_nofreeze_filtered --copy_weights --wandb --dataset_split filtered
