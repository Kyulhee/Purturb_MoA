"""
Cross-Perturbation MoA Framework v2 -- Training & Evaluation Pipeline

Key changes from v1:
- Uses rFCFP fingerprint input instead of drug index embedding
- Evaluates both classification head and metric (prototype) head
- Supports zero-shot MoA evaluation via nearest-neighbor in latent space
"""
import os
import json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.metrics import (
    accuracy_score, top_k_accuracy_score, f1_score,
    confusion_matrix, classification_report
)
from pathlib import Path
from config import ModelConfig, RESULTS_DIR, MODEL_DIR, SEED
from model import CrossPerturbMoAModel


def set_seed(seed: int = SEED):
    """Set random seeds for reproducibility."""
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def train_epoch(model, dataloader, optimizer, disc_optimizer, device, config):
    """Train one epoch."""
    model.train()
    epoch_losses = {}

    for batch in dataloader:
        x, moa_labels, drug_labels, cl_labels, dose, fingerprints = [b.to(device) for b in batch]

        # Main model step
        optimizer.zero_grad()
        total_loss, loss_dict = model.compute_loss(
            x, fingerprints, dose, cl_labels, moa_labels,
            n_covariates=config.covariate_dim, perturb_idx=drug_labels
        )
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

        # Discriminator step (separate optimizer)
        disc_optimizer.zero_grad()
        with torch.no_grad():
            _, _, _, z_perturb, _, _ = model(x, fingerprints, dose, cl_labels, drug_labels)
        cov_pred = model.discriminator(z_perturb.detach())
        loss_disc = nn.functional.cross_entropy(cov_pred, cl_labels)
        loss_disc.backward()
        disc_optimizer.step()

        # Accumulate losses
        for k, v in loss_dict.items():
            if k not in epoch_losses:
                epoch_losses[k] = []
            epoch_losses[k].append(v.item() if torch.is_tensor(v) else v)

    # Average losses
    return {k: np.mean(v) for k, v in epoch_losses.items()}


def evaluate(model, dataloader, device, n_moa_classes: int, top_k: list = None):
    """Evaluate model on test set with both classification and metric heads."""
    model.eval()
    top_k = top_k or [1, 3, 5]

    all_moa_true = []
    all_moa_pred_cls = []
    all_moa_pred_metric = []
    all_moa_probs_cls = []
    all_moa_probs_metric = []
    all_z_perturb = []
    all_recon_losses = []
    all_fingerprints = []

    with torch.no_grad():
        for batch in dataloader:
            x, moa_labels, drug_labels, cl_labels, dose, fingerprints = [b.to(device) for b in batch]

            x_recon, moa_logits, metric_logits, z_perturb, _, _ = model(
                x, fingerprints, dose, cl_labels, drug_labels
            )
            recon_loss = nn.functional.mse_loss(x_recon, x, reduction='none').mean(dim=-1)

            probs_cls = torch.softmax(moa_logits, dim=-1)
            probs_metric = torch.softmax(metric_logits, dim=-1)

            all_moa_true.append(moa_labels.cpu().numpy())
            all_moa_pred_cls.append(moa_logits.argmax(dim=-1).cpu().numpy())
            all_moa_pred_metric.append(metric_logits.argmax(dim=-1).cpu().numpy())
            all_moa_probs_cls.append(probs_cls.cpu().numpy())
            all_moa_probs_metric.append(probs_metric.cpu().numpy())
            all_z_perturb.append(z_perturb.cpu().numpy())
            all_recon_losses.append(recon_loss.cpu().numpy())
            all_fingerprints.append(fingerprints.cpu().numpy())

    moa_true = np.concatenate(all_moa_true)
    moa_pred_cls = np.concatenate(all_moa_pred_cls)
    moa_pred_metric = np.concatenate(all_moa_pred_metric)
    moa_probs_cls = np.concatenate(all_moa_probs_cls)
    moa_probs_metric = np.concatenate(all_moa_probs_metric)
    z_perturb = np.concatenate(all_z_perturb)
    recon_losses = np.concatenate(all_recon_losses)

    # Classification head metrics
    results_cls = {
        'cls_accuracy_top1': accuracy_score(moa_true, moa_pred_cls),
        'cls_f1_macro': f1_score(moa_true, moa_pred_cls, average='macro', zero_division=0),
        'cls_f1_weighted': f1_score(moa_true, moa_pred_cls, average='weighted', zero_division=0),
    }

    # Metric (prototype) head metrics
    results_metric = {
        'metric_accuracy_top1': accuracy_score(moa_true, moa_pred_metric),
        'metric_f1_macro': f1_score(moa_true, moa_pred_metric, average='macro', zero_division=0),
        'metric_f1_weighted': f1_score(moa_true, moa_pred_metric, average='weighted', zero_division=0),
    }

    # Top-k accuracy for both heads
    for k in top_k:
        if n_moa_classes >= k:
            try:
                results_cls[f'cls_accuracy_top{k}'] = top_k_accuracy_score(
                    moa_true, moa_probs_cls, labels=np.arange(n_moa_classes), k=k
                )
            except Exception:
                results_cls[f'cls_accuracy_top{k}'] = 0.0
            try:
                results_metric[f'metric_accuracy_top{k}'] = top_k_accuracy_score(
                    moa_true, moa_probs_metric, labels=np.arange(n_moa_classes), k=k
                )
            except Exception:
                results_metric[f'metric_accuracy_top{k}'] = 0.0

    # Reconstruction MSE
    results_recon = {'recon_mse': float(np.mean(recon_losses))}

    # Combined results
    results = {**results_cls, **results_metric, **results_recon}

    # Per-class classification report (metric head)
    cls_report = classification_report(moa_true, moa_pred_metric, output_dict=True,
                                        zero_division=0)
    results['classification_report'] = cls_report

    return results, z_perturb, moa_true


def evaluate_zero_shot(model, train_z: np.ndarray, train_moa: np.ndarray,
                       test_z: np.ndarray, test_moa: np.ndarray,
                       device, top_k: list = None):
    """
    Evaluate zero-shot MoA prediction via nearest-neighbor in latent space.

    For leave-MoA-out: test MoA classes are unseen during training.
    The metric head's classification output cannot predict unseen classes,
    but nearest-neighbor matching in latent space can find the most similar
    TRAINING sample and report its MoA.
    """
    top_k = top_k or [1, 3, 5]
    model.eval()

    z_ref = torch.tensor(train_z, dtype=torch.float32).to(device)
    ref_labels = torch.tensor(train_moa, dtype=torch.long).to(device)
    z_query = torch.tensor(test_z, dtype=torch.float32).to(device)

    results = {}
    for k in top_k:
        predicted_moa, similarities = model.predict_zero_shot(
            z_query, z_ref, ref_labels, top_k=k
        )
        # Top-1 accuracy: does the nearest training sample have the same MoA?
        top1_correct = (predicted_moa[:, 0].cpu().numpy() == test_moa).mean()
        results[f'zero_shot_top{k}_acc'] = float(top1_correct)

    return results


def pretrain_fingerprint_encoder(model, train_ds, metadata: dict,
                                  device: torch.device, n_epochs: int = 200,
                                  lr: float = 1e-3):
    """
    Pretrain the rFCFP drug encoder to predict MoA from fingerprint alone.

    This is critical for leave-compound-out: the fingerprint encoder must learn
    a meaningful fingerprint->MoA mapping before the full model training, because
    during full model training, gradients from 20K cells dilute the fingerprint signal.

    After pretraining, the drug encoder's weights serve as a warm start that already
    captures fingerprint->MoA structure, making it much easier for the full model
    to leverage fingerprint information for unseen compounds.
    """
    print("Pretraining fingerprint encoder (fingerprint -> MoA)...")

    # Extract unique drug fingerprints and their MoA labels from training data
    # Use the first occurrence of each drug to get unique (fingerprint, MoA) pairs
    drug_labels_all = train_ds.tensors[2].numpy()  # drug indices
    moa_labels_all = train_ds.tensors[1].numpy()    # MoA labels
    fp_all = train_ds.tensors[5].numpy()            # fingerprints

    # Get unique drugs and their MoA (majority vote)
    unique_drugs = np.unique(drug_labels_all)
    drug_fps = []
    drug_moas = []
    for d in unique_drugs:
        mask = drug_labels_all == d
        fp = fp_all[mask][0]  # all same for same drug
        moa = np.bincount(moa_labels_all[mask]).argmax()
        drug_fps.append(fp)
        drug_moas.append(moa)

    X_fp = torch.tensor(np.array(drug_fps), dtype=torch.float32).to(device)
    y_moa = torch.tensor(np.array(drug_moas), dtype=torch.long).to(device)

    print(f"  Pretraining on {len(unique_drugs)} unique drugs, {metadata['n_moa']} MoA classes")

    # Train just the drug_encoder + a temporary classification head
    encoder = model.drug_encoder
    perturb_dim = model.config.perturbation_dim
    temp_head = nn.Sequential(
        nn.Linear(perturb_dim, 64),
        nn.ReLU(),
        nn.Linear(64, metadata['n_moa']),
    ).to(device)

    optimizer = optim.Adam(
        list(encoder.parameters()) + list(temp_head.parameters()),
        lr=lr
    )

    for ep in range(1, n_epochs + 1):
        encoder.train()
        temp_head.train()
        z = encoder(X_fp, torch.ones(len(X_fp), 1, device=device))  # dose=1 for pretraining
        logits = temp_head(z)
        loss = nn.functional.cross_entropy(logits, y_moa)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if ep % 50 == 0:
            with torch.no_grad():
                pred = logits.argmax(dim=-1)
                acc = (pred == y_moa).float().mean().item()
            print(f"  Pretrain epoch {ep}: loss={loss.item():.4f}, acc={acc:.4f}")

    # Final accuracy
    with torch.no_grad():
        encoder.eval()
        temp_head.eval()
        z = encoder(X_fp, torch.ones(len(X_fp), 1, device=device))
        acc = (temp_head(z).argmax(dim=-1) == y_moa).float().mean().item()
    print(f"  Pretrain final accuracy: {acc:.4f}")

    return model


def train_pipeline(train_ds, test_ds, metadata: dict, config: ModelConfig = None,
                   save_dir: Path = RESULTS_DIR, pretrain_encoder: bool = True):
    """
    Full training pipeline for Cross-Perturbation MoA Framework v2.
    """
    config = config or ModelConfig()
    set_seed(SEED)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # DataLoaders
    train_loader = DataLoader(train_ds, batch_size=config.batch_size,
                               shuffle=True, num_workers=0)
    test_loader = DataLoader(test_ds, batch_size=config.batch_size,
                              shuffle=False, num_workers=0)

    # Model
    fp_bits = metadata.get('fingerprint_bits', 1024)
    model = CrossPerturbMoAModel(
        n_genes=metadata['n_genes'],
        n_perturbations=metadata['n_drugs'] + 1,
        n_covariates=metadata['n_cell_lines'],
        n_moa_classes=metadata['n_moa'],
        config=config,
        fingerprint_bits=fp_bits,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    print(f"Model parameters: {n_params:,}")

    # Pretrain fingerprint encoder
    if pretrain_encoder:
        model = pretrain_fingerprint_encoder(model, train_ds, metadata, device)
        # Freeze drug_encoder after pretraining so MoA structure is preserved.
        # Without freezing, full training corrupts the fingerprint->MoA mapping
        # by overfitting to seen drugs, destroying generalization to unseen compounds.
        for param in model.drug_encoder.parameters():
            param.requires_grad = False
        n_frozen = sum(p.numel() for p in model.drug_encoder.parameters())
        print(f"Frozen drug_encoder: {n_frozen:,} parameters")

    # Optimizers (exclude frozen drug_encoder)
    trainable_params = [p for p in model.parameters() if p.requires_grad]
    optimizer = optim.Adam(trainable_params, lr=config.learning_rate, weight_decay=1e-5)
    disc_optimizer = optim.Adam(model.discriminator.parameters(),
                                 lr=config.learning_rate * 2, weight_decay=1e-5)

    # Training loop
    best_acc = 0.0
    patience_counter = 0
    history = []

    print(f"\nTraining for {config.n_epochs} epochs...")
    print(f"Train: {metadata['train_size']}, Test: {metadata['test_size']}")
    print(f"MoA classes: {metadata['n_moa']}, Drugs: {metadata['n_drugs']}, Genes: {metadata['n_genes']}")
    print(f"Fingerprint bits: {fp_bits}")
    print("-" * 80)

    for epoch in range(1, config.n_epochs + 1):
        # Train
        train_losses = train_epoch(model, train_loader, optimizer, disc_optimizer,
                                    device, config)

        # Evaluate
        if epoch % config.eval_every == 0 or epoch == config.n_epochs:
            test_results, z_test, moa_true = evaluate(
                model, test_loader, device, metadata['n_moa'], config.top_k
            )

            # Log
            log_str = (f"Epoch {epoch:3d} | "
                       f"Train Loss: {train_losses.get('total', 0):.4f} | "
                       f"Cls Top-1: {test_results['cls_accuracy_top1']:.4f} | "
                       f"Metric Top-1: {test_results['metric_accuracy_top1']:.4f} | "
                       f"Recon MSE: {test_results['recon_mse']:.4f}")
            print(log_str)

            # Save history
            record = {'epoch': epoch, 'train_losses': train_losses, **test_results}
            history.append(record)

            # Early stopping (use metric head accuracy)
            metric_acc = test_results['metric_accuracy_top1']
            if metric_acc > best_acc:
                best_acc = metric_acc
                patience_counter = 0
                save_dir.mkdir(parents=True, exist_ok=True)
                torch.save(model.state_dict(), save_dir / "best_model.pt")
            else:
                patience_counter += 1
                if patience_counter >= config.early_stop_patience:
                    print(f"Early stopping at epoch {epoch}")
                    break

    # Final evaluation
    best_model_path = save_dir / "best_model.pt"
    if best_model_path.exists():
        model.load_state_dict(torch.load(best_model_path,
                                          map_location=device, weights_only=True))
    final_results, z_test, moa_true = evaluate(
        model, test_loader, device, metadata['n_moa'], config.top_k
    )

    # Also get training set embeddings for zero-shot evaluation
    train_loader_eval = DataLoader(train_ds, batch_size=config.batch_size,
                                    shuffle=False, num_workers=0)
    _, z_train, moa_train = evaluate(
        model, train_loader_eval, device, metadata['n_moa'], config.top_k
    )

    # Zero-shot evaluation (nearest-neighbor in latent space)
    zero_shot_results = evaluate_zero_shot(
        model, z_train, moa_train, z_test, moa_true, device, config.top_k
    )
    final_results.update(zero_shot_results)

    print(f"\n--- Final Results ({metadata['split_strategy']}) ---")
    print(f"  Classification Head  - Top-1: {final_results['cls_accuracy_top1']:.4f}")
    print(f"  Metric (Prototype)   - Top-1: {final_results['metric_accuracy_top1']:.4f}")
    print(f"  Zero-shot NN Top-1:            {final_results.get('zero_shot_top1_acc', 0):.4f}")
    print(f"  Recon MSE:                     {final_results['recon_mse']:.4f}")

    # Save results
    save_results(final_results, history, metadata, save_dir)

    return model, final_results, z_test, moa_true


def save_results(results: dict, history: list, metadata: dict, save_dir: Path):
    """Save evaluation results to JSON."""
    save_dir.mkdir(parents=True, exist_ok=True)

    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    # Save final results
    final_save = {}
    for k, v in results.items():
        if k != 'classification_report':
            final_save[k] = convert(v)
    final_save['metadata'] = {k: convert(v) for k, v in metadata.items()}

    with open(save_dir / "final_results.json", 'w') as f:
        json.dump(final_save, f, indent=2)

    # Save training history
    history_save = []
    for record in history:
        r = {}
        for k, v in record.items():
            if k == 'train_losses':
                r[k] = {kk: convert(vv) for kk, vv in v.items()}
            elif k == 'classification_report':
                continue
            else:
                r[k] = convert(v)
        history_save.append(r)

    with open(save_dir / "training_history.json", 'w') as f:
        json.dump(history_save, f, indent=2)

    print(f"\nResults saved to {save_dir}")


def evaluate_cross_perturbation(model, drug_z: np.ndarray, drug_moa: np.ndarray,
                                 crispr_z: np.ndarray = None, crispr_moa: np.ndarray = None,
                                 top_k: int = 5):
    """Evaluate cross-perturbation matching."""
    device = next(model.parameters()).device

    if crispr_z is not None:
        z_query = torch.tensor(crispr_z, dtype=torch.float32).to(device)
        z_ref = torch.tensor(drug_z, dtype=torch.float32).to(device)
        ref_labels = torch.tensor(drug_moa, dtype=torch.long).to(device)

        topk_idx, topk_sim, matched_moa = model.cross_perturbation_match(
            z_query, z_ref, ref_labels, top_k=top_k
        )

        query_moa = torch.tensor(crispr_moa, dtype=torch.long).to(device)
        top1_correct = (matched_moa[:, 0] == query_moa).float().mean().item()

        return {
            'cross_perturb_top1_acc': top1_correct,
            'n_queries': len(crispr_z),
        }
    else:
        z = torch.tensor(drug_z, dtype=torch.float32).to(device)
        labels = torch.tensor(drug_moa, dtype=torch.long).to(device)

        topk_idx, topk_sim, matched_moa = model.cross_perturbation_match(
            z, z, labels, top_k=top_k
        )

        n = len(drug_z)
        top1_correct = sum(
            1 for i in range(n)
            if matched_moa[i, 0].item() == labels[i].item()
        )

        return {
            'within_drug_top1_match': top1_correct / max(n, 1),
            'n_samples': n,
        }
