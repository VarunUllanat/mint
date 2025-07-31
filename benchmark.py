
import torch
import time
import pandas as pd
import matplotlib.pyplot as plt
from mint.helpers.extract import load_config, MINTWrapper

# === Config ===
device = 'cuda:0'
cfg = load_config("data/esm2_t33_650M_UR50D.json")
checkpoint_path = 'mint.ckpt'

# === Initialize model ===
wrapper = MINTWrapper(cfg, checkpoint_path, sep_chains=True, device=device)
wrapper.train()
for param in wrapper.model.parameters():
    param.requires_grad = True

# === Parameters ===
sequence_lengths = [64, 128, 256, 512]
chain_counts = [2, 3, 4]
max_batch_try = 128

# === Helper: Dummy input batch generator ===
def get_dummy_batch(batch_size, seq_len, chains):
    x = []
    ids = []
    for c in range(chains):
        x.append(torch.randint(1, 32, (batch_size, seq_len), dtype=torch.long))
        ids.append(torch.full((batch_size, seq_len), c, dtype=torch.long))
    x = torch.cat(x, dim=1).to(device)
    ids = torch.cat(ids, dim=1).to(device)
    return x, ids

# === Forward Pass Benchmark ===
forward_results = []

wrapper.eval()
for num_chains in chain_counts:
    for seq_len in sequence_lengths:
        input_tokens = []
        chain_ids = []

        for c in range(num_chains):
            tokens = torch.randint(low=1, high=32, size=(seq_len,), dtype=torch.long)
            input_tokens.append(tokens)
            chain_ids.append(torch.full((seq_len,), c, dtype=torch.long))

        input_tokens = torch.cat(input_tokens).unsqueeze(0).to(device)
        chain_ids = torch.cat(chain_ids).unsqueeze(0).to(device)

        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()
        torch.cuda.synchronize()
        start = time.perf_counter()

        with torch.no_grad():
            _ = wrapper(input_tokens, chain_ids)

        torch.cuda.synchronize()
        elapsed_time = time.perf_counter() - start
        peak_memory = torch.cuda.max_memory_allocated(device) / 1e6

        forward_results.append({
            "chains": num_chains,
            "seq_len": seq_len,
            "time": elapsed_time,
            "memory": peak_memory,
        })


# === Plot Forward Pass Results ===
df_forward = pd.DataFrame(forward_results)

plt.figure()
for chains in sorted(df_forward['chains'].unique()):
    subset = df_forward[df_forward['chains'] == chains]
    plt.plot(subset['seq_len'], subset['memory'], label=f'{chains} chains', marker='o')
plt.xlabel('Sequence Length per Chain')
plt.ylabel('Peak Memory (MB)')
plt.title('Forward Pass: Memory Usage vs Sequence Length')
plt.legend()
plt.tight_layout()
plt.savefig("mint_memory_plot.png")

plt.figure()
for chains in sorted(df_forward['chains'].unique()):
    subset = df_forward[df_forward['chains'] == chains]
    plt.plot(subset['seq_len'], subset['time'], label=f'{chains} chains', marker='o')
plt.xlabel('Sequence Length per Chain')
plt.ylabel('Forward Pass Time (s)')
plt.title('Forward Pass: Time vs Sequence Length')
plt.legend()
plt.tight_layout()
plt.savefig("mint_time_plot.png")

print("Forward pass plots saved: mint_memory_plot.png and mint_time_plot.png")

# === Backward Pass Max Batch Size Benchmark ===
wrapper.train()
backward_results = []

for chains in chain_counts:
    for seq_len in sequence_lengths:
        batch_size = 1
        max_ok = 0

        while batch_size <= max_batch_try:
            try:
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats()

                x, ids = get_dummy_batch(batch_size, seq_len, chains)
                output = wrapper(x, ids)
                loss = output.sum()
                loss.backward()

                max_ok = batch_size
                batch_size *= 2

            except RuntimeError as e:
                if 'out of memory' in str(e):
                    break
                else:
                    raise e

        backward_results.append({
            'sequence_length': seq_len,
            'chain_count': chains,
            'max_batch_size': max_ok
        })

df_backward = pd.DataFrame(backward_results)
df_backward.to_csv("mint_backward_batchsize.csv", index=False)
print("\n Backward pass results saved: mint_backward_batchsize.csv")
print("\n")
print(df_backward)




