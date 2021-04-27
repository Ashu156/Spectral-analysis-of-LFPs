# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 08:38:20 2021

@author: Dell
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
#from sklearn.model_selection import KFold
#from sklearn.model_selection import cross_val_score 
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt
#import seaborn as sns

#Load data
df = pd.read_csv('C:/Users/Dell/Documents/lfp-svm-data.csv')
df = df.sample(n = 1000, replace = True)

#Data preprocessing
x = df.dropna(axis = 1)        # get rid of NA columns in the data frame
x = x.drop('stress', axis = 1) # get rid of labels
y = df['stress']               # designate labels

# Alternative way to evaluate model accuracy

acc_score = [];
fpr = [];
tpr = [];
thres = [];
ra_score = [];
num_iter = 1000;

for i in range(num_iter):
    val = i + 1
    print(val)
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.30, random_state = None)
    
    #Train the model
    model = GaussianNB()
    model.fit(x_train, y_train)
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

  
   
#Calculate ROC for the model

y_score = model.predict_proba(x_test)[:,1]
false_positive_rate, true_positive_rate, threshold = roc_curve(y_test, y_score)
print('roc_auc_score for Naive Bayes Classifier: ', roc_auc_score(y_test, y_score))
    
