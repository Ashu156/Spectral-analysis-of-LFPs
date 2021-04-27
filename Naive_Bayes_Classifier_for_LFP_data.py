# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:39:05 2021

@author: Dell
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score 
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns

#Load data
df = pd.read_csv('C:/Users/Dell/Documents/lfp-svm-data.csv')
df = df.sample(n = 1000, replace = True)

#Data preprocessing
x = df.dropna(axis = 1)        # get rid of NA columns in the data frame
x = x.drop('stress', axis = 1) # get rid of labels
y = df['stress']               # designate labels


    

# x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.3, random_state = 42)
 
# Implementing cross validation

k = 5
kf = KFold(n_splits = k, random_state = None)

#Train the model
model = GaussianNB()

result = cross_val_score(model, x, y, cv = kf)
print("Avg accuracy: {}".format(result.mean()))

# Alternative way to evaluate model accuracy

#acc_score = [];
#
#for train_index , test_index in kf.split(x):
#    x_train , x_test = x.iloc[train_index,:],x.iloc[test_index,:]
#    y_train , y_test = y[train_index] , y[test_index]
#    
#    model.fit(x_train,y_train)
#    pred_values = model.predict(x_test)
#    
#    acc = accuracy_score(pred_values , y_test)
#    acc_score.append(acc)
#    
#    avg_acc_score = sum(acc_score)/k
# 
#print('accuracy of each fold - {}'.format(acc_score))
#print('Avg accuracy : {}'.format(avg_acc_score))


x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.2, random_state = 42)

# Train the model
model = GaussianNB()
model.fit(x_train, y_train)

#Calculate ROC for the model
from sklearn.metrics import roc_curve, roc_auc_score
y_score = model.predict_proba(x_test)[:,1]
false_positive_rate, true_positive_rate, threshold = roc_curve(y_test, y_score)
print('roc_auc_score for Naive Bayes Classifier: ', roc_auc_score(y_test, y_score))

## Plotting ROC for the classifier
plt.subplots(1, figsize=(10,10))
plt.title('Receiver Operating Characteristic - Naive Bayes Classifer')
plt.plot(false_positive_rate, true_positive_rate)
plt.plot([0, 1], ls="--")
plt.plot([0, 0], [1, 0] , c="0.7"), plt.plot([1, 1] , c="0.7")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()

###########################################################################
# Test against a shuffled model version
###########################################################################

# Load data
df = pd.read_csv('C:/Users/Dell/Documents/lfp-svm-data.csv')
#df = df.sample(n = 1000, replace = True)

# Data preprocessing
x = df.dropna(axis = 1)        # get rid of NA columns in the data frame
x = x.drop('stress', axis = 1) # get rid of labels
y = df['stress']              # designate labels


# Shuffle labels
np.random.shuffle(y)
x = x.sample(n = 1000, replace = True)
y = y.sample(n = 1000, replace = True)

# Implementing cross validation (k-fold)
k = 5
kf = KFold(n_splits = k, random_state = 42)

#Train the model
model = GaussianNB()

result = cross_val_score(model, x, y, cv = kf)
print("Avg accuracy: {}".format(result.mean()))

# Alternative way to evaluate model accuracy

#acc_score = [];
#
#for train_index , test_index in kf.split(x):
#    x_train , x_test = x.iloc[train_index,:],x.iloc[test_index,:]
#    y_train , y_test = y[train_index] , y[test_index]
#    
#    model.fit(x_train,y_train)
#    pred_values = model.predict(x_test)
#    
#    acc = accuracy_score(pred_values , y_test)
#    acc_score.append(acc)
#    
#    avg_acc_score = sum(acc_score)/k
# 
#print('accuracy of each fold - {}'.format(acc_score))
#print('Avg accuracy : {}'.format(avg_acc_score))

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.2, random_state = 42)

# Train the model
model = GaussianNB()
model.fit(x_train, y_train)

#Calculate ROC for the model
from sklearn.metrics import roc_curve, roc_auc_score
y_score = model.predict_proba(x_test)[:,1]
false_positive_rate, true_positive_rate, threshold = roc_curve(y_test, y_score)
print('roc_auc_score for Naive Bayes Classifier: ', roc_auc_score(y_test, y_score))

## Plotting ROC for the classifier
plt.subplots(1, figsize=(10,10))
plt.title('Receiver Operating Characteristic - Naive Bayes Classifer')
plt.plot(false_positive_rate, true_positive_rate)
plt.plot([0, 1], ls="--")
plt.plot([0, 0], [1, 0] , c="0.7"), plt.plot([1, 1] , c="0.7")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()

#plt.rcParams['text.usetex'] = True
#plt.rcParams['text.latex.unicode'] = True
#plt.savefig('NBC.eps', format='eps', dpi=600)