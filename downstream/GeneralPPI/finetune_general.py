import argparse

import numpy as np
import scipy
import torch
import wandb
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    accuracy_score,
    auc,
    f1_score,
    mean_squared_error,
    precision_recall_curve,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.preprocessing import PowerTransformer, StandardScaler
from tasks import get_task_datasets
from torch import nn
from torch.utils.data import Dataset
from tqdm import tqdm

MIN_EMB_DIM = 640


def clean_tn(input_tensor):
    if type(input_tensor) == tuple:
        input_tensor = input_tensor[0]
    input_tensor = input_tensor.squeeze()
    return input_tensor


def regression_metrics(test_targets, Y_pred):
    r2_score(test_targets, Y_pred)
    p_corr = scipy.stats.pearsonr(test_targets, Y_pred)
    s_corr = scipy.stats.spearmanr(test_targets, Y_pred)
    mse = mean_squared_error(test_targets, Y_pred)
    return {"pearson": p_corr[0], "spearman": s_corr[0], "rmse": mse ** (1.0 / 2.0)}


def classification_metrics(targets, predictions, threshold=0.5):
    binary_predictions = (predictions >= threshold).astype(int)
    accuracy = accuracy_score(targets, binary_predictions)
    f1 = f1_score(targets, binary_predictions)
    auc_score = roc_auc_score(targets, predictions)
    precision_vals, recall_vals, _ = precision_recall_curve(targets, predictions)
    auprc = auc(recall_vals, precision_vals)
    return {
        "Accuracy": accuracy,
        "AUPRC": auprc,
        "F1 Score": f1,
        "AUROC": auc_score,
    }


def multilabel_metrics(label, output):
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    pre_y = (output > 0.5).astype(int)
    truth_y = label
    N, C = pre_y.shape
    for i in range(N):
        for j in range(C):
            if pre_y[i][j] == truth_y[i][j]:
                if truth_y[i][j] == 1:
                    TP += 1
                else:
                    TN += 1
            elif truth_y[i][j] == 1:
                FN += 1
            elif truth_y[i][j] == 0:
                FP += 1

        # Accuracy = (TP + TN) / (N*C + 1e-10)
        Precision = TP / (TP + FP + 1e-10)
        Recall = TP / (TP + FN + 1e-10)
        F1_score = 2 * Precision * Recall / (Precision + Recall + 1e-10)

    return {"f1": F1_score, "precision": Precision, "recall": Recall}


class MLPDataset(Dataset):
    def __init__(self, x, y):
        super().__init__()
        self.x = x
        self.y = y

    def __len__(self):
        return len(self.y)

    def __getitem__(self, index):
        return self.x[index], self.y[index]


def return_ridge_model(use_mlp, embed_dim, seed=0):
    if use_mlp:
        model = MLPRegressor()
        param_grid = {
            "activation": ["relu"],
            "alpha": [0.0001],
            "learning_rate": ["adaptive"],
            "solver": ["adam"],
            "learning_rate_init": [0.001],
            "max_iter": [100],
            "hidden_layer_sizes": [(640,)],
            "random_state": [seed],
            "early_stopping": [True],
            "validation_fraction": [0.1],
            "tol": [1e-4],
        }
    else:
        model = Ridge()
        lambda_grid = np.logspace(0, -6, num=7).tolist()
        lambda_grid.append(0)
        param_grid = {"alpha": lambda_grid, "max_iter": [15000]}
    return model, param_grid, regression_metrics


def return_logistic_model(use_mlp, embed_dim, seed=0):
    if use_mlp:
        model = MLPClassifier()
        param_grid = {
            "activation": ["relu"],
            "alpha": [0.0001],
            "learning_rate": ["adaptive"],
            "solver": ["adam"],
            "learning_rate_init": [0.001],
            "max_iter": [100],
            "hidden_layer_sizes": [(640,),],
            "early_stopping": [True],
            "random_state": [seed],
            "validation_fraction": [0.1],
            "tol": [1e-4],
        }
    else:
        model = LogisticRegression()
        param_grid = {
            "penalty": ["l2"],
            "C": [0.0001, 0.0005, 0.001, 0.0025, 0.005],
            "max_iter": [15000],
        }
    return model, param_grid, classification_metrics


def calculate_mean_std(metrics):
    aggregated_metrics = {}
    for metric_dict in metrics:
        for key, value in metric_dict.items():
            if key not in aggregated_metrics:
                aggregated_metrics[key] = []
            aggregated_metrics[key].append(value)
    mean_std_metrics = {}
    for key, values in aggregated_metrics.items():
        mean_std_metrics[key + "_mean"] = np.mean(values)
        mean_std_metrics[key + "_std"] = np.std(values)
    return mean_std_metrics


def convert_train_test_labels(train, test):
    t = PowerTransformer()
    train = t.fit_transform(train.reshape(-1, 1))
    test = t.transform(test.reshape(-1, 1))
    return train[:, 0], test[:, 0]


def convert_train_test_features(train, test):
    t = StandardScaler()
    train = t.fit_transform(train)
    test = t.transform(test)
    return train, test


def pre_defined_cv(embeddings, targets, input_df, task_type, use_mlp, reps):

    X_scaled = embeddings.numpy()
    embed_dim = X_scaled.shape[-1]

    cv_cols = input_df.filter(like="split").columns

    for cv_col in cv_cols:
        metrics_rep = []
        for rep in range(reps):
            if task_type == "reg":
                model, param_grid, metric_fn = return_ridge_model(use_mlp, embed_dim, seed=rep)
                scoring = "neg_mean_squared_error"
            elif task_type == "bc":
                model, param_grid, metric_fn = return_logistic_model(use_mlp, embed_dim, seed=rep)
                scoring = "roc_auc"

            train_idx = input_df[cv_col] == "train"
            test_idx = input_df[cv_col] == "test"

            train_embeddings = X_scaled[train_idx]
            test_embeddings = X_scaled[test_idx]

            train_embeddings, test_embeddings = convert_train_test_features(
                train_embeddings, test_embeddings
            )

            train_targets = targets[train_idx]
            test_targets = targets[test_idx]

            if task_type == "reg":
                train_targets, test_targets = convert_train_test_labels(
                    train_targets, test_targets
                )

            grid_search = GridSearchCV(
                estimator=model, param_grid=param_grid, cv=3, scoring=scoring, verbose=10
            )
            grid_search.fit(train_embeddings, train_targets)
            best_model = grid_search.best_estimator_
            if task_type == "reg":
                Y_pred = best_model.predict(test_embeddings)
            elif task_type == "bc":
                Y_pred = best_model.predict_proba(X_test)[:, 1]

            metrics_dict = metric_fn(test_targets, Y_pred)
            metrics_dict["test_size"] = len(test_targets)
            metrics_dict = {k + "_" + cv_col: v for k, v in metrics_dict.items()}

            metrics_rep.append(metrics_dict)

        all_metrics_rep = calculate_mean_std(metrics_rep)
        # wandb.log(all_metrics_rep)

    return all_metrics_rep


def cross_validation(embeddings, targets, task_type, use_mlp, normalize):

    X_scaled = embeddings.numpy()
    embed_dim = X_scaled.shape[-1]

    if task_type == "reg":
        model, param_grid, metric_fn = return_ridge_model(use_mlp, embed_dim)
        scoring = "neg_mean_squared_error"
    elif task_type == "bc":
        model, param_grid, metric_fn = return_logistic_model(use_mlp, embed_dim)
        scoring = "roc_auc"

    outer_cv = KFold(n_splits=10, shuffle=True, random_state=0)
    inner_cv = KFold(n_splits=5, shuffle=True, random_state=0)
    clf = GridSearchCV(
        estimator=model, param_grid=param_grid, cv=inner_cv, scoring=scoring, verbose=1
    )

    metrics = []

    for train_idx, test_idx in tqdm(outer_cv.split(X_scaled), total=10):
        X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
        Y_train, Y_test = targets[train_idx], targets[test_idx]

        if normalize:
            X_train, X_test = convert_train_test_features(X_train, X_test)

        if task_type == "reg":
            Y_train, Y_test = convert_train_test_labels(Y_train, Y_test)

        # convert to dataloader
        # train()

        clf.fit(X_train, Y_train)
        best_model = clf.best_estimator_
        if task_type == "reg":
            Y_pred = best_model.predict(X_test)
        elif task_type == "bc":
            Y_pred = best_model.predict_proba(X_test)[:, 1]

        metrics.append(metric_fn(Y_test, Y_pred))

    all_metrics = calculate_mean_std(metrics)

    # wandb.log(all_metrics)

    return all_metrics


class SimpleMLP(nn.Module):
    def __init__(self, input_size, output_size):
        super(SimpleMLP, self).__init__()

        self.project = torch.nn.Sequential(
            nn.Linear(input_size, input_size),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(input_size, output_size),
        )

    def forward(self, x):
        return self.project(x)


@torch.no_grad()
def evaluate(model, loader, name, task_type, device):
    preds = []
    targets = []
    for step, eval_batch in enumerate(loader):
        embs, target = eval_batch
        embs = embs.to(device)
        target = target.to(device)
        pred = model(embs).squeeze(-1)

        if "c" in task_type:
            pred = torch.sigmoid(pred)

        preds.append(pred.detach().cpu().numpy())
        targets.append(target.cpu().numpy())
    preds = np.concatenate(preds)
    targets = np.concatenate(targets)

    if task_type == "reg":
        metrics = regression_metrics(targets, preds)
    elif task_type == "bc":
        metrics = classification_metrics(targets, preds)
    elif task_type == "mc":
        metrics = multilabel_metrics(targets, preds)

    return {f"{name}_{k}": i for k, i in metrics.items()}


def train_mlp(
    train_loader, val_loader, test_loader, metadata, input_size, monitor_metric, device, rep
):
    torch.manual_seed(rep)
    model = SimpleMLP(input_size, metadata["output_size"]).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=metadata["lr"])
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=10
    )

    if metadata["task_type"] == "reg":
        loss_fn = torch.nn.MSELoss()
    else:
        loss_fn = torch.nn.BCEWithLogitsLoss()

    best_val_metric = 0

    for epoch in tqdm(range(metadata["num_epochs"])):
        loss_accum = 0
        for step, train_batch in enumerate(train_loader):

            model.train()
            optimizer.zero_grad()

            embs, target = train_batch
            embs = embs.to(device)
            target = target.to(device)

            pred = model(embs)
            loss = loss_fn(pred.squeeze(), target.float())

            loss.backward()
            optimizer.step()

            loss_accum += loss.detach().cpu().item()

        test_metrics_dict = evaluate(model, test_loader, "test", metadata["task_type"], device)
        val_metrics_dict = evaluate(model, val_loader, "val", metadata["task_type"], device)
        all_metrics_dict = {**test_metrics_dict, **val_metrics_dict}
        all_metrics_dict["train_loss"] = loss_accum / (step + 1)

        if val_metrics_dict[f"val_{monitor_metric}"] >= best_val_metric:
            best_test_metrics = test_metrics_dict
            best_val_metric = val_metrics_dict[f"val_{monitor_metric}"]

        scheduler.step(loss_accum / (step + 1))

    best_test_metrics = {f"best_{k}": i for k, i in best_test_metrics.items()}
    return best_test_metrics


def clean_targets(t):
    if type(t[0]) == str:
        t = [i.split(";") for i in t]
        return np.array(t).astype(int)
    else:
        return np.array(t)


def main(args):
    # wandb.login(key="03c57b40de2b2f02c0d3c1357868f3ef656696ce")
    # wandb.init()

    train_dataset, val_dataset, test_dataset, metadata = get_task_datasets(
        args.task, return_metadata=True
    )
    train_targets = clean_targets(train_dataset.targets)
    val_targets = clean_targets(val_dataset.targets) if val_dataset is not None else None
    test_targets = clean_targets(test_dataset.targets) if test_dataset is not None else None

    if args.model == "all":
        model_list = [
            "esm2_t30_150M_UR50D",
            "esm2_t33_650M_UR50D",
            "esm1b_t33_650M_UR50S",
            "esm2_t36_3B_UR50D",
            "progen2-large",
            "prot_t5_xl_uniref50",
            "prot_t5_xl_bfd",
            "mint",
        ]
    else:
        model_list = [args.model]

    for model_name in model_list:

        if args.sep:
            model_name = model_name + "_sep"

        try:
            train_embs = clean_tn(torch.load(f"embeddings/{args.task}/{model_name}/train.pt"))
            val_embs = (
                clean_tn(torch.load(f"embeddings/{args.task}/{model_name}/val.pt"))
                if val_dataset is not None
                else None
            )
            test_embs = (
                clean_tn(torch.load(f"embeddings/{args.task}/{model_name}/test.pt"))
                if test_dataset is not None
                else None
            )
        except:
            print(f"No {args.task} embeddings file found for {model_name}")
            continue

        input_size = train_embs.shape[-1]

        assert len(train_embs) == len(train_targets)

        if metadata["method"] == "cv":
            metrics = cross_validation(
                train_embs,
                train_targets,
                metadata["task_type"],
                args.use_mlp_for_cv,
                args.normalize,
            )
        elif metadata["method"] == "pcv":
            metrics = pre_defined_cv(
                train_embs,
                train_targets,
                metadata["train_df"],
                metadata["task_type"],
                args.use_mlp_for_cv,
                reps=args.rep,
            )
        else:
            all_metrics = []
            for rep in range(args.rep):
                train_dataset = MLPDataset(train_embs, train_targets)
                # test = train if no test provided
                test_dataset = (
                    MLPDataset(test_embs, test_targets)
                    if test_dataset is not None
                    else train_dataset
                )
                # val = test if no val provided
                val_dataset = (
                    MLPDataset(val_embs, val_targets) if val_dataset is not None else test_dataset
                )

                train_loader = torch.utils.data.DataLoader(
                    train_dataset, batch_size=args.batch_size, shuffle=True, drop_last=True
                )
                val_loader = torch.utils.data.DataLoader(
                    val_dataset, batch_size=args.batch_size, shuffle=False
                )
                test_loader = torch.utils.data.DataLoader(
                    test_dataset, batch_size=args.batch_size, shuffle=False
                )

                metrics = train_mlp(
                    train_loader,
                    val_loader,
                    test_loader,
                    metadata,
                    input_size,
                    metadata["monitor_metric"],
                    args.device,
                    rep,
                )
                all_metrics.append(metrics)
            # wandb.log(calculate_mean_std(all_metrics))
            print(calculate_mean_std(all_metrics))

        # wandb.finish()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # General args
    parser.add_argument("--task", type=str, default="HumanPPI")
    parser.add_argument("--model", type=str, default="all")
    parser.add_argument("--batch_size", type=int, default=16)
    parser.add_argument("--rep", type=int, default=3)
    parser.add_argument("--device", type=str, default="cuda:0")
    parser.add_argument("--use_mlp_for_cv", action="store_true", default=False)
    parser.add_argument("--normalize", action="store_true", default=False)
    parser.add_argument("--sep", action="store_true", default=False)

    args = parser.parse_args()
    main(args)
