# Phase G3 — Fine-tune Geneformer on CSC-high vs CSC-low
#
# Base model: Geneformer-V2-104M_CLcancer (pre-trained on cancer single cells)
# Task:       Binary classification — CSC-high (1) vs CSC-low (0)
# Method:     Freeze first N transformer layers, fine-tune top layers + classifier head
#
# Input:  data/geneformer/brca_csc.dataset
# Output: results/geneformer_finetuned/   ← fine-tuned model weights

import os
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("results/geneformer_finetuned", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

MODEL_DIR    = "geneformer_repo/Geneformer-V2-104M_CLcancer"
DATASET_PATH = "data/geneformer/brca_csc.dataset"
OUT_DIR      = "results/geneformer_finetuned"
OUT_PREFIX   = "brca_csc"

# ─────────────────────────────────────────────
# STEP 1 — Load and inspect tokenized dataset
# ─────────────────────────────────────────────
from datasets import load_from_disk

print("Loading tokenized dataset...")
dataset = load_from_disk(DATASET_PATH)
print(f"Dataset: {len(dataset)} cells, columns: {dataset.column_names}")

import collections
label_counts = collections.Counter(dataset["csc_label"])
print(f"Label distribution: {dict(label_counts)}")

# ─────────────────────────────────────────────
# STEP 2 — Prepare data splits with Classifier
# ─────────────────────────────────────────────
from geneformer import Classifier

cc = Classifier(
    classifier="cell",
    cell_state_dict={"state_key": "csc_label", "states": "all"},
    training_args={
        "num_train_epochs":              5,
        "learning_rate":                 5e-5,
        "lr_scheduler_type":             "cosine",
        "warmup_ratio":                  0.05,
        "per_device_train_batch_size":   1,    # keep MPS/CPU buffer under ~750 MB
        "per_device_eval_batch_size":    1,
        "gradient_accumulation_steps":   12,   # effective batch = 12
        "gradient_checkpointing":        True, # trade compute for memory
        "weight_decay":                  0.01,
        "output_dir":                    OUT_DIR,
        "evaluation_strategy":           "epoch",
        "save_strategy":                 "epoch",
        "load_best_model_at_end":        True,
        "metric_for_best_model":         "eval_macro_f1",
        "logging_steps":                 100,
        "fp16":                          False,
        "use_mps_device":                False, # force CPU; MPS has 9 GiB buffer limit on attn
        "no_cuda":                       True,
    },
    freeze_layers=6,
    num_crossval_splits=1,
    split_sizes={"train": 0.8, "valid": 0.1, "test": 0.1},
    forward_batch_size=12,
    nproc=4,
    ngpu=0,
    model_version="V2",
)

print("\nPreparing train/valid/test splits...")
# Do NOT pass attr_to_split — prepare_data renames the state_key column
# before splitting, which breaks the attr_to_split lookup. Use random splits.
cc.prepare_data(
    input_data_file=DATASET_PATH,
    output_directory=OUT_DIR,
    output_prefix=OUT_PREFIX,
)
print("Data preparation complete.")

# prepare_data creates _labeled_train and _labeled_test (no separate valid)
train_dataset = load_from_disk(os.path.join(OUT_DIR, f"{OUT_PREFIX}_labeled_train.dataset"))
test_dataset  = load_from_disk(os.path.join(OUT_DIR, f"{OUT_PREFIX}_labeled_test.dataset"))

print(f"\nSplit sizes — train: {len(train_dataset)}, test: {len(test_dataset)}")

# Load id_class_dict (maps integer label IDs back to string names)
id_class_path = os.path.join(OUT_DIR, f"{OUT_PREFIX}_id_class_dict.pkl")
with open(id_class_path, "rb") as f:
    id_class_dict = pickle.load(f)
print(f"Label mapping: {id_class_dict}")
num_classes = len(id_class_dict)

# ─────────────────────────────────────────────
# STEP 3 — Fine-tune
# ─────────────────────────────────────────────
print(f"\nFine-tuning Geneformer ({MODEL_DIR})...")
print(f"  Base model: Geneformer-V2-104M_CLcancer")
print(f"  Frozen layers: 6 (fine-tuning top 6 + classifier head)")
print(f"  Epochs: 10  |  LR: 5e-5  |  Batch: 12")

all_metrics = cc.train_classifier(
    model_directory=MODEL_DIR,
    num_classes=num_classes,
    train_data=train_dataset,
    eval_data=test_dataset,
    output_directory=OUT_DIR,
    predict=True,
)
print("\nTraining complete.")
print(f"Final metrics: {all_metrics}")

# ─────────────────────────────────────────────
# STEP 4 — Evaluate on held-out test set
# ─────────────────────────────────────────────
print("\nEvaluating on test set...")
test_metrics = cc.evaluate_saved_model(
    model_directory=OUT_DIR,
    id_class_dict_file=id_class_path,
    test_data_file=os.path.join(OUT_DIR, f"{OUT_PREFIX}_labeled_test.dataset"),
    output_directory=OUT_DIR,
    output_prefix=f"{OUT_PREFIX}_test_eval",
)
print(f"Test metrics: {test_metrics}")

# Save metrics summary
metrics_df = pd.DataFrame([{
    "phase": "train/valid",
    **{k: v for k, v in all_metrics.items() if isinstance(v, (int, float))},
}])
metrics_df.to_csv(f"results/tables/G3_training_metrics.csv", index=False)

print(f"\n✓ Phase G3 complete.")
print(f"  Fine-tuned model saved: {OUT_DIR}")
print(f"  Next: Phase G4 — Extract attention weights")
