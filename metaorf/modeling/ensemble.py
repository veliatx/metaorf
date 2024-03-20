import copy
import numpy as np


from sklearn.pipeline import make_pipeline
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.metrics import RocCurveDisplay, roc_curve, auc, precision_recall_curve, average_precision_score, classification_report
from sklearn.model_selection import train_test_split, RandomizedSearchCV, RepeatedStratifiedKFold, StratifiedKFold, cross_validate
from sklearn.base import clone

import matplotlib.pyplot as plt


class Dataset:
    def __init__(self, X, y, name, orf_ids, model=None, feature_importance_df=None):
        self.X = X
        self.y = y
        self.name = name
        self.orf_ids = orf_ids
        self.model = model
        self.feature_importance_df = feature_importance_df


def plot_roc_pr(X, y, model, fpr_cutoff=0.1, n_splits=10):
    """
    """
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
    
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    
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


def plot_pr(X, y, model, n_splits=5):
    """
    """
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=2)

    feature_importances = []

    precisions = []
    aps = []
    thresholds = []
    mean_recall = np.linspace(0, 1, 100)

    fig, ax = plt.subplots(figsize=(10, 8))

    for i, (train, test) in enumerate(cv.split(X, y)):
        X_train, X_test = X.iloc[train], X.iloc[test]
        y_train, y_test = y[train], y[test]
        model_clone = clone(model)
        model_clone.fit(X_train, y_train)
        feature_importances.append(model_clone.feature_importances_)

        y_score = model_clone.predict_proba(X_test)[:, 1]
        precision, recall, threshold = precision_recall_curve(y_test, y_score)
        threshold = np.append(threshold, 1)
        ap = average_precision_score(y_test, y_score)
        precisions.append(np.interp(mean_recall, recall[::-1], precision[::-1]))
        aps.append(ap)
        thresholds.append(np.interp(mean_recall, recall[::-1], threshold[::-1]))
        
        ax.plot(recall, precision, alpha=0.3, label=f'PR fold {i+1} (AP = {ap:.2f})')

    mean_precision = np.mean(precisions, axis=0)
    std_precision = np.std(precisions, axis=0)

    mean_thresh = np.mean(thresholds, axis=0)
    std_thresh = np.std(thresholds, axis=0)

    ax.plot(mean_recall, mean_thresh, color='black',
            label=f'Mean Decision Probability', lw=2, alpha=0.8)
    ax.fill_between(mean_recall, mean_thresh-(std_thresh/2), mean_thresh+(std_thresh/2), color='grey', alpha=0.2)

    ax.plot(mean_recall, mean_precision, color='blue',
            label=f'Mean PR (AP = {np.mean(aps):.2f} $\pm$ {np.std(aps):.2f})', lw=2, alpha=0.8)
    ax.fill_between(mean_recall, mean_precision-std_precision, mean_precision+std_precision, color='grey', alpha=0.2)

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('Precision-Recall Curves Across K-Fold Validation')
    ax.legend(loc="lower left")

    return feature_importances, fig


def plot_roc(ds, classifier, n_splits=5, fpr_cutoff=.05):
    """
    """
    cv = StratifiedKFold(n_splits=5)

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots()

    fpr_cutoff=0.05
    closest_fpr_diff = float('inf')
    feature_importances = []

    for i, (train, test) in enumerate(cv.split(ds.X, ds.y)):
        
        cloned_classifier = clone(classifier)

        cloned_classifier.fit(ds.X.iloc[train], ds.y[train])
        y_test = cloned_classifier.predict_proba(ds.X.iloc[test])[:, 1]
        y = ds.y[test]
        viz = roc_curve(y, y_test)
        fpr, tpr, thresholds = viz
        feature_importances.append(cloned_classifier.feature_importances_)

        tprs.append(np.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        ax.plot(fpr, tpr, lw=1, alpha=0.3, label=f'ROC fold {i+1} (AUC = {roc_auc:.2f})')
        

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')
    ax.set_xlabel('False Positive Rate (FPR)')
    ax.set_ylabel('True Positive Rate (TPR)')
    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
        title="Receiver Operating Characteristic")
    ax.legend(loc="lower right")
    
    return feature_importances, fig

