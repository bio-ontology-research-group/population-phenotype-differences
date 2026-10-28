#!/bin/bash
# Qwen3.8-27B-FP8 (vLLM 0.28.0 cu129 build, same model and version as the BORG unimatrix server) serving a
# local client that extracts patients/phenotypes/origin from PubMed case-report abstracts.
# One H200 per array task; the server is bound to localhost with a per-job random key.
#SBATCH --partition=pi-hohndor
#SBATCH --qos=pi-hohndor
#SBATCH --job-name=qwen-sel
#SBATCH --array=0-2
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --gres=gpu:nvidia_h200:1
set -euo pipefail
cd /mnt/data/u/hohndor/qwen-serve
export PATH=$HOME/.pixi/bin:$PATH HF_HOME=/mnt/data/u/hohndor/hf HF_HUB_OFFLINE=1
T=${SLURM_ARRAY_TASK_ID}; PORT=$((8100 + T))
echo "[$(date -Is)] node=$(hostname) job=${SLURM_JOB_ID} task=$T"
nvidia-smi --query-gpu=name,memory.total --format=csv
umask 077; mkdir -p keys out; KEYF=keys/key.${SLURM_JOB_ID}.$T
python3 -c "import secrets; print(secrets.token_hex(24))" > "$KEYF"
pixi run --manifest-path /mnt/data/u/hohndor/qwen-serve-cu129/pixi.toml vllm serve Qwen/Qwen3.8-27B-FP8 --served-model-name qwen3.8-27b --host 127.0.0.1 --port $PORT \
  --api-key "$(cat $KEYF)" --max-model-len 32768 --kv-cache-dtype fp8 --gpu-memory-utilization 0.92 \
  --reasoning-parser qwen3 --speculative-config '{"method":"mtp","num_speculative_tokens":2}' > logs/vllm-${SLURM_JOB_ID}_$T.log 2>&1 &
VPID=$!
for i in $(seq 1 120); do curl -sf http://127.0.0.1:$PORT/health >/dev/null && break; sleep 10; done
curl -sf http://127.0.0.1:$PORT/health >/dev/null || { echo "vllm did not start"; tail -50 logs/vllm-${SLURM_JOB_ID}_$T.log; kill $VPID; exit 1; }
echo "[$(date -Is)] vllm up"
pixi run --manifest-path /mnt/data/u/hohndor/qwen-serve-cu129/pixi.toml python -I qwen_tasks.py select --cands shards/sel_$T.jsonl --out out/sel_$T.jsonl \
  --key-file "$KEYF" --endpoint http://127.0.0.1:$PORT --jobs 128 --batch 10
# retry abstracts that failed to parse (resumable: only non-ok entries are redone)
pixi run --manifest-path /mnt/data/u/hohndor/qwen-serve-cu129/pixi.toml python -I qwen_tasks.py abstract --cands shards/cands_$T.jsonl --out out/abs_$T.jsonl \
  --key-file "$KEYF" --endpoint http://127.0.0.1:$PORT --jobs 64
kill $VPID; rm -f "$KEYF"
echo "[$(date -Is)] DONE"
