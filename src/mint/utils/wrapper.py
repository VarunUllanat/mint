import json
import time
from collections import defaultdict

import lightning as pl
import numpy as np
import torch
import wandb

from ..model.esm2 import ESM2
from .utils.logging import get_logger

logger = get_logger(__name__)


def gather_log(log, world_size):
    if world_size == 1:
        return log
    log_list = [None] * world_size
    torch.distributed.all_gather_object(log_list, log)
    log = {key: sum([l[key] for l in log_list], []) for key in log}
    return log


def get_log_mean(log):
    out = {}
    for key in log:
        try:
            out[key] = np.mean(log[key])
        except:
            pass
    if log:
        out["entries"] = len(log[key])
    return out


class ESMWrapper(pl.LightningModule):
    def __init__(self, cfg, args):
        super().__init__()
        self.save_hyperparameters()
        self.cfg = cfg
        self.args = args
        self.model = ESM2(
            num_layers=cfg.encoder_layers,
            embed_dim=cfg.encoder_embed_dim,
            attention_heads=cfg.encoder_attention_heads,
            token_dropout=cfg.token_dropout,
            use_multimer=not args.no_multimer,
        )
        self.iter_step = -1
        self._log = defaultdict(list)
        self.last_log_time = time.time()

    def training_step(self, batch, batch_idx):
        self.stage = "train"
        loss = self.forward(batch)
        # Manual checkpointing
        # if self.iter_step % 15000 == 0:
        # if self.trainer.is_global_zero:
        # torch.save(self.model.state_dict(), f'./workdir/3B_nofreeze/checkpoint_iter_{self.iter_step}.pt')
        return loss

    def validation_step(self, batch, batch_idx):
        self.stage = "val"
        self.forward(batch)
        if self.args.validate:
            self.try_print_log()

    def forward(self, batch):
        self.iter_step += 1
        # 15% of tokens randomly sampled from the sequence. For those 15% of tokens, we change the input token to a special “masking”
        # token with 80% probability, a randomly-chosen alternate amino acid token with 10% probability, and the original input token
        # (i.e. no change) with 10% probability. We take the loss to be the whole batch average cross entropy loss between the model’s
        # predictions and the true token for these 15% of amino acid tokens.

        tokens, chain_ids = batch
        mask = (
            (~tokens.eq(self.model.cls_idx))
            & (~tokens.eq(self.model.eos_idx))
            & (~tokens.eq(self.model.padding_idx))
        )
        mask = (torch.rand(tokens.shape, device=tokens.device) < 0.15) & mask

        rand = torch.rand(tokens.shape, device=tokens.device)
        randaa = torch.randint(4, 24, tokens.shape, device=tokens.device)

        inp = tokens
        inp = torch.where((rand < 0.8) & mask, self.model.mask_idx, inp)
        inp = torch.where((rand > 0.9) & mask, randaa, inp)

        out = self.model(inp, chain_ids)["logits"]
        loss = torch.nn.functional.cross_entropy(out.transpose(1, 2), tokens, reduction="none")
        loss = (loss * mask).sum() / mask.sum()

        self.log("tokens", mask.sum())
        self.log("loss", loss)
        self.log("perplexity", torch.exp(loss))
        self.log("dur", time.time() - self.last_log_time)
        self.last_log_time = time.time()
        return loss

    def try_print_log(self):
        step = self.iter_step if self.args.validate else self.trainer.global_step
        if (step + 1) % self.args.print_freq == 0:
            log = self._log
            log = {key: log[key] for key in log if "iter_" in key}

            log = gather_log(log, self.trainer.world_size)
            mean_log = get_log_mean(log)
            mean_log.update(
                {
                    "epoch": self.trainer.current_epoch,
                    "step": self.trainer.global_step,
                    "iter_step": self.iter_step,
                }
            )
            if self.trainer.is_global_zero:
                logger.info(str(mean_log))
                if self.args.wandb:
                    wandb.log(mean_log)
            for key in list(log.keys()):
                if "iter_" in key:
                    del self._log[key]

    def log(self, key, data):
        if isinstance(data, torch.Tensor):
            data = data.detach().cpu().item()
        log = self._log
        log["iter_" + key].append(data)
        log[self.stage + "_" + key].append(data)

    def on_train_epoch_end(self):
        log = self._log
        log = {key: log[key] for key in log if "train_" in key}
        log = gather_log(log, self.trainer.world_size)
        mean_log = get_log_mean(log)
        mean_log.update(
            {
                "epoch": self.trainer.current_epoch,
                "step": self.trainer.global_step,
                "iter_step": self.iter_step,
            }
        )

        if self.trainer.is_global_zero:
            logger.info(str(mean_log))
            if self.args.wandb:
                wandb.log(mean_log)

            # path = os.path.join(
            #     os.environ["MODEL_DIR"], f"train_{self.trainer.current_epoch}.csv"
            # )
            # pd.DataFrame(log).to_csv(path)
        for key in list(log.keys()):
            if "train_" in key:
                del self._log[key]

    def on_validation_epoch_end(self):
        log = self._log
        log = {key: log[key] for key in log if "val_" in key}
        log = gather_log(log, self.trainer.world_size)
        if self.trainer.is_global_zero:
            logger.info(str(get_log_mean(log)))
            if self.args.wandb:
                wandb.log(get_log_mean(log))

        for key in list(log.keys()):
            if "val_" in key:
                del self._log[key]

    def on_before_optimizer_step(self, optimizer, x):
        self.try_print_log()
        if self.args.check_grad:
            for name, p in self.model.named_parameters():
                if p.requires_grad and p.grad is None:
                    print(name)

    def configure_optimizers(self):
        # For model training optimization, we used Adam with 𝛽𝛽1 = 0.9, 𝛽𝛽2 = 0.98, 𝜖𝜖 = 10−8 and 𝐿𝐿2 weight decay of
        # 0.01 for all models except the 15 billion parameter model, where we used a weight decay of 0.1. The learning rate is
        # warmed up over the first 2,000 steps to a peak value of 4e-4 (1.6e-4 for the 15B parameter model), and then linearly
        # decayed to one tenth of its peak value over the 90% of training duration
        if self.args.freeze_self_attn:
            self.model.requires_grad_(False)
            for name, p in self.model.named_parameters():
                if "multimer_attn" in name:
                    p.requires_grad = True

        optimizer = torch.optim.AdamW(
            filter(lambda p: p.requires_grad, self.model.parameters()),
            lr=self.cfg.lr[0],
            betas=json.loads(self.cfg.adam_betas),
            eps=self.cfg.adam_eps,
            weight_decay=self.cfg.weight_decay,
        )

        warmup = torch.optim.lr_scheduler.LinearLR(
            optimizer, start_factor=1e-12, end_factor=1.0, total_iters=self.cfg.warmup_updates
        )
        decay = torch.optim.lr_scheduler.LinearLR(
            optimizer,
            start_factor=1.0,
            end_factor=self.cfg.end_learning_rate / self.cfg.lr[0],
            total_iters=int(0.9 * int(self.cfg.total_num_update)),
        )
        scheduler = torch.optim.lr_scheduler.SequentialLR(
            optimizer, schedulers=[warmup, decay], milestones=[self.cfg.warmup_updates]
        )

        return {
            "optimizer": optimizer,
            "lr_scheduler": {"scheduler": scheduler, "interval": "step"},
        }
