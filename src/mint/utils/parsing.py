import os
from argparse import ArgumentParser


def parse_train_args():
    parser = ArgumentParser()
    parser.add_argument("--ckpt", type=str, default=None)
    parser.add_argument("--max_steps", type=int, default=450000)
    parser.add_argument("--val_examples", type=int, default=250000)
    parser.add_argument("--no_multimer", action="store_true")
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--concat", action="store_true")
    parser.add_argument("--reinitialize", action="store_true")
    parser.add_argument("--overfit", action="store_true")
    parser.add_argument("--copy_weights", action="store_true")
    parser.add_argument("--dataset_seek", type=int, default=None)
    parser.add_argument("--batch_size", type=int, default=2)
    parser.add_argument("--crop_length", type=int, default=512)
    parser.add_argument("--val_max_len", type=int, default=None)
    parser.add_argument(
        "--model", choices=["8M", "35M", "150M", "650M", "3B", "15B"], default="650M"
    )
    parser.add_argument("--check_grad", action="store_true")
    parser.add_argument("--freeze_self_attn", action="store_true")
    parser.add_argument("--accumulate_grad", type=int, default=1)
    parser.add_argument("--grad_clip", type=float, default=1.0)
    parser.add_argument("--val_check_interval", type=int, default=10000)
    parser.add_argument("--print_freq", type=int, default=100)
    parser.add_argument("--wandb", action="store_true")
    parser.add_argument("--run_name", type=str, default="default")
    parser.add_argument("--dataset_split", type=str, choices=["full", "filtered", "filtered_50"])

    args = parser.parse_args()
    os.environ["MODEL_DIR"] = os.path.join("workdir", args.run_name)
    os.environ["WANDB_LOGGING"] = str(int(args.wandb))
    # if args.wandb:
    #     if subprocess.check_output(["git", "status", "-s"]):
    #         exit()
    # args.commit = (
    #     subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()
    # )
    return args
