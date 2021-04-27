# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 21:23:27 2021

@author: Dell
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import KFold, train_test_split
from sklearn.svm import SVC
from sklearn.metrics import auc, roc_curve, roc_auc_score, plot_roc_curve, accuracy_score
from sklearn.naive_bayes import GaussianNB

# #############################################################################
# Load data

df = pd.read_csv('C:/Users/Dell/Documents/lfp-svm-data-bla-dmpfc.csv')
df = df.sample(n = 1000, replace = True)

# Data preprocessing
X = df.dropna(axis = 1)        # get rid of NA columns in the data frame
X = X.drop('stress', axis = 1) # get rid of labels
y = df['stress']               # designate labels

# #############################################################################
# Classification and ROC analysis

# Run classifier with cross-validation and plot ROC curves
acc_score = [];
fpr = [];
tpr = [];
thres = [];
ra_score = [];
num_iter = 1000;
fig, ax = plt.subplots()
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

for i in range(num_iter):
    val = i + 1
    print(val)
    x_train, x_test, y_train, y_test = train_test_split(X, y, test_size = 0.30, random_state = None)
    
    #Train the model
    model = GaussianNB()
    model.fit(x_train, y_train)
    viz = plot_roc_curve(model, x_test, y_test,
                           alpha=0.3, lw=1, ax=ax)
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)
    y_score = model.predict_proba(x_test)[:, 1]
    false_positive_rate, true_positive_rate, threshold = roc_curve(y_test, y_score)
    y_predict = model.predict(x_test)
    acc = accuracy_score(y_predict, y_test)
    acc_score.append(acc)
    fpr.append(false_positive_rate)
    tpr.append(true_positive_rate)
    thres.append(threshold)
    ra_score.append(roc_auc_score(y_test, y_score))
    
avg_acc_score = sum(acc_score)/num_iter
print('accuracy of each fold -{}'.format(acc_score))
print('Avg accuracy : {}'.format(avg_acc_score))

fpr = pd.DataFrame(fpr);
tpr = pd.DataFrame(tpr);
thres = pd.DataFrame(thres)
ra_score = pd.DataFrame(ra_score)

# tprs = []
# aucs = []
# mean_fpr = np.linspace(0, 1, 100)


# for i, (train, test) in enumerate(cv.split(X, y)):
#     classifier.fit(X[train], y[train])
#     viz = plot_roc_curve(classifier, X[test], y[test],
#                          name='ROC fold {}'.format(i),
#                          alpha=0.3, lw=1, ax=ax)
#     interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
#     interp_tpr[0] = 0.0
#     tprs.append(interp_tpr)
#     aucs.append(viz.roc_auc)

ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
        label='Chance', alpha=.8)

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

ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
        title="Receiver operating characteristic example")
ax.legend(loc="lower right")
plt.show()