import copy
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn import svm
from sklearn.model_selection import cross_validate
from sklearn.pipeline import make_pipeline
from sklearn import metrics

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.metrics import RocCurveDisplay, roc_curve, auc, precision_recall_curve, average_precision_score, classification_report
from sklearn.model_selection import train_test_split, RandomizedSearchCV, StratifiedKFold
from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.base import clone

import matplotlib.pyplot as plt


class MLP(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, 1)
        self.sigmoid = nn.Sigmoid() 
        
    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.sigmoid(x)
        return x
    
    
def evaluate_model(criterion, model, validation_data_x, validation_data_y,
                   plot=False, plot_folder=None, experiment_name=None):
    # prepare data
    X_val_tensor = torch.tensor(validation_data_x, dtype=torch.float32)
    y_val_tensor = torch.tensor(validation_data_y, dtype=torch.float32)
    val_dataset = TensorDataset(X_val_tensor, y_val_tensor)
    val_loader = DataLoader(val_dataset, batch_size=64, shuffle=False)

    # Evaluate the model
    model.eval()
    y_pred = []
    y_true = []
    evaluation_results = []
    with torch.no_grad():
        val_loss = 0.0
        for inputs, labels in val_loader:
            outputs = model(inputs)
            y_pred.extend(outputs.cpu().numpy())
            y_true.extend(labels.cpu().numpy())
            loss = criterion(outputs, labels.view(-1, 1))
            val_loss += loss.item()
    roc_auc = metrics.roc_auc_score(y_true, y_pred)
    pr_auc = metrics.average_precision_score(y_true, y_pred)
    val_loss /= len(val_loader)
    print(f"Validation Loss: {val_loss}\nValidation ROC-AUC: {roc_auc}\nValidation PR-AUC: {pr_auc}")
    
    if plot:
        fpr, tpr, _ = metrics.roc_curve(y_true, y_pred)
        plt.plot(fpr, tpr)
        plt.xlabel('false positive rate')
        plt.ylabel('true positive rate')
        plt.title(f"{experiment_name} - roc_auc {roc_auc}")
        plt.savefig(f"{plot_folder}/{experiment_name}_roc_auc.png")
        plt.close()
    
        precision, recall, _ = metrics.precision_recall_curve(y_true, y_pred, drop_intermediate=True)
        plt.plot(precision, recall)
        plt.xlabel('precision')
        plt.ylabel('recall')
        plt.title(f"{experiment_name} - pr_auc {pr_auc}")
        plt.savefig(f"{plot_folder}/{experiment_name}_pr_auc.png")
        plt.close()
    return val_loss, roc_auc, pr_auc

    
def train_model(training_data_X, training_data_y,
                validation_data_x, validation_data_y):
    # define model
    input_size = len(training_data_X[0])
    print(f"input_size: {input_size}")
    hidden_size = 64
    output_size = 1
    model = MLP(input_size, hidden_size)

    # Define loss function and optimizer
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    # prepare data
    X_train_tensor = torch.tensor(training_data_X, dtype=torch.float32)
    y_train_tensor = torch.tensor(training_data_y, dtype=torch.float32)
    print(X_train_tensor.shape, y_train_tensor.shape)
    train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
    train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)

    # train
    num_epochs = 100
    current_best = None
    best_loss = float('inf')
    patience = 10
    wait = 0
    best_model = None
    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels.view(-1, 1))
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
        print(f"Epoch {epoch+1}, Training Loss: {running_loss / len(train_loader)}")
        
        # eval
        val_loss, _, _ = evaluate_model(criterion, model, validation_data_x, validation_data_y)
        if val_loss < best_loss:
            best_loss = val_loss
            wait = 0
            # Save the model if needed
            best_model = copy.deepcopy(model)
        else:
            wait += 1
            if wait >= patience:
                print("Model training completes")
                break
                
        print("==========================")
    return best_model


def predict(model, test_data_X):
    # prepare data    
    X_test_tensor = torch.tensor(test_data_X, dtype=torch.float32)
    test_dataset = TensorDataset(X_test_tensor)
    test_loader = DataLoader(test_dataset, batch_size=64, shuffle=False)

    # Evaluate the model
    model.eval()
    y_pred = []
    with torch.no_grad():
        for inputs in test_loader:
            outputs = model(inputs[0])
            y_pred.extend(outputs.cpu().numpy())
    return y_pred
    

def run(X_groups, y_groups, test_X):
    roc_auc_list, pr_auc_list = [], []
    test_y_pred_groups = [-np.inf for _ in range(len(test_X))]
    for group_idx in range(5):
        test_data_X, test_data_y = X_groups[group_idx], y_groups[group_idx]        
        training_data_X, training_data_y = [], []
        for current_group_idx in range(4):
            if current_group_idx != group_idx:
                training_data_X += X_groups[current_group_idx]
                training_data_y += y_groups[current_group_idx]
        
        model = train_model(training_data_X, training_data_y)
        roc_auc, pr_auc = evaluate_model(model, test_data_X, test_data_y)
        roc_auc_list.append(roc_auc)
        pr_auc_list.append(pr_auc)
        
        test_y_pred = predict(model, test_X)
        for idx, test_y_pred_val in enumerate(predict(model, test_X)):
            test_y_pred_groups[idx] = max(test_y_pred_val[0], test_y_pred_groups[idx])
        
    return roc_auc_list, pr_auc_list, test_y_pred_groups


def plot_roc_pr(X, y, model, fpr_cutoff=0.1, n_splits=10):
    cv = StratifiedKFold(n_splits=n_splits)
    
    feature_importances = []
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    
    precisions = []
    aps = []
    mean_recall = np.linspace(0, 1, 100)
    
    closest_model = None
    closest_fpr_diff = np.inf
    
    fig, ax = plt.subplots(1, 2, figsize=(10, 8))
    
    for i, (train, test) in enumerate(cv.split(X, y)):
        X_train, X_test = X.iloc[train], X.iloc[test]
        y_train, y_test = y[train], y[test]
        model_clone = clone(model)
        model_clone.fit(X_train, y_train)
        
        y_score = model_clone.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_score)
        roc_auc = auc(fpr, tpr)
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        aucs.append(roc_auc)
        feature_importances.append(model_clone.feature_importances_)

        
        precision, recall, _ = precision_recall_curve(y_test, y_score)
        ap = average_precision_score(y_test, y_score)
        precisions.append(np.interp(mean_recall, recall[::-1], precision[::-1]))
        aps.append(ap)
        
        # Plot individual ROC and PR curves
        ax[0].plot(fpr, tpr, alpha=0.3, label=f'ROC fold {i+1} (AUC = {roc_auc:.2f})')
        ax[1].plot(recall, precision, alpha=0.3, label=f'PR fold {i+1} (AP = {ap:.2f})')
        
        # Check if this model is the closest to the desired FPR cutoff
        fpr_diff = np.abs(fpr - fpr_cutoff)
        min_fpr_diff_index = np.argmin(fpr_diff)
        if fpr_diff[min_fpr_diff_index] < closest_fpr_diff:
            closest_fpr_diff = fpr_diff[min_fpr_diff_index]
            closest_model = model_clone
    
    # Plot mean ROC and PR curves with error bars
    mean_tpr = np.mean(tprs, axis=0)
    mean_precision = np.mean(precisions, axis=0)
    mean_tpr[-1] = 1.0
    std_tpr = np.std(tprs, axis=0)
    std_precision = np.std(precisions, axis=0)
    
    ax[0].plot(mean_fpr, mean_tpr, color='blue',
               label=f'Mean ROC (AUC = {np.mean(aucs):.2f} $\pm$ {np.std(aucs):.2f})', lw=2, alpha=0.8)
    ax[0].fill_between(mean_fpr, mean_tpr-std_tpr, mean_tpr+std_tpr, color='grey', alpha=0.2)
    ax[0].set_xlabel('False Positive Rate')
    ax[0].set_ylabel('True Positive Rate')
    ax[0].set_title('ROC Curves Across K-Fold Validation')
    ax[0].legend(loc="lower right")
    
    ax[1].plot(mean_recall, mean_precision, color='blue',
               label=f'Mean PR (AP = {np.mean(aps):.2f} $\pm$ {np.std(aps):.2f})', lw=2, alpha=0.8)
    ax[1].fill_between(mean_recall, mean_precision-std_precision, mean_precision+std_precision, color='grey', alpha=0.2)
    ax[1].set_xlabel('Recall')
    ax[1].set_ylabel('Precision')
    ax[1].set_title('Precision-Recall Curves Across K-Fold Validation')
    ax[1].legend(loc="lower left")
    
    #plt.tight_layout()
    #plt.show()
    
    return closest_model, feature_importances, fig
