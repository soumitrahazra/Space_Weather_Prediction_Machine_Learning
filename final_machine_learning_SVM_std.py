import sys
import scipy
import numpy as np
import matplotlib 
import sklearn
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import Imputer
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn import svm
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from collections import Counter
from sklearn.svm import SVC


df=pd.read_csv('shazra_final_ready_another_loop24span12.csv')
#where_are_NaNs = isnan(df)
#df[where_are_NaNs] = 0
df=df.fillna(df.mean())
#df=df.dropna()

xs, ys = df.loc[:, 'col1':'col119'], df['col120'] #X <- 1st column name : last column name, y <- label column name

rs=33
Xa_train, Xa_test, ya_train, ya_test = train_test_split(xs, ys, test_size = 0.3, random_state=rs, stratify = ys)
#scalers = preprocessing.StandardScaler().fit(Xa_train)
scalers = preprocessing.MinMaxScaler().fit(Xa_train)
Xa_train = scalers.transform(Xa_train)
Xa_test = scalers.transform(Xa_test)
classifier=SVC(kernel='rbf') 
parameters=[{'C':[0.001, 0.01, 0.1, 1, 10], 'gamma':[0.001, 0.01, 0.1, 1]}] 
Grid_search=GridSearchCV(estimator=classifier, param_grid=parameters, cv=10) 
Grid_search.fit(Xa_train, ya_train)
best=Grid_search.best_params_

def classification_report_with_accuracy_score(y_true, y_pred):
    #print (classification_report(y_true, y_pred) )# print classification report
    #print( confusion_matrix(y_true, y_pred))
    TN, FP, FN, TP = confusion_matrix(y_true, y_pred).ravel()
    P= TP+FN
    N=TN+FP
    CH= ((TP+FP)*(TP+FN))/(P+N)
    HSS_1= (TP +TN -N)/P
    HSS_2= (2*((TP*TN) - (FN*FP)))/((P*(FN+TN))+(N*(TP+FP)))
    GS= (TP-CH)/(TP+FP+FN-CH)
    TSS = (TP/(TP+FN)) - (FP/(FP+TN))
    Pr = precision_score(y_true, y_pred, pos_label=1)
    Pr1 = precision_score(y_true, y_pred, pos_label=0)
    Re = recall_score(y_true, y_pred, pos_label=1)
    Re1= recall_score(y_true, y_pred, pos_label=0)
    F1 = f1_score(y_true, y_pred, pos_label=1)
    F11 = f1_score(y_true, y_pred, pos_label=0)
    a=accuracy_score(y_true, y_pred) # return accuracy score
    
    
    return a,Pr,Pr1,Re,Re1,F1,F11,HSS_1, HSS_2, GS, TSS, TP, TN, FP, FN


## For calculate average and standard deviation over 100 realization

a_f=[]
Pr_f=[]
Pr1_f=[]
Re_f=[]
Re1_f=[]
F1_f=[]
F11_f=[]
HSS_1_f=[]
HSS_2_f=[]
GS_f=[]
TSS_f=[]
for i in range(100):
    X_train, X_test, y_train, y_test = train_test_split(xs, ys, test_size = 0.3, stratify = ys)
#standardization (almost like normalizing)  
    #scaler = preprocessing.StandardScaler().fit(X_train)
    scaler = preprocessing.MinMaxScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)
    #classifier_svm= SVC(C=best['C'], gamma=best['gamma'], kernel='rbf', class_weight={1:10})
    classifier_svm= SVC(C=best['C'], gamma=best['gamma'], kernel='rbf', class_weight='balanced')
    #classifier_svm= SVC(C=best['C'], gamma=best['gamma'], kernel='rbf')
    classifier_svm.fit(X_train, y_train)
    y_pred = classifier_svm.predict(X_test)
    a,Pr,Pr1,Re,Re1,F1,F11,HSS_1, HSS_2, GS, TSS, TP, TN, FP, FN=classification_report_with_accuracy_score(y_test, y_pred)
    a_f=a_f+[a]
    Pr_f=Pr_f+[Pr]
    Pr1_f=Pr1_f+[Pr1]
    Re_f=Re_f+[Re]
    Re1_f=Re1_f+[Re1]
    F1_f=F1_f+[F1]
    F11_f=F11_f+[F11]
    HSS_1_f=HSS_1_f+[HSS_1]
    HSS_2_f=HSS_2_f+[HSS_2]
    GS_f=GS_f+[GS]
    TSS_f=TSS_f+[TSS]


print ("TP:", TP, "FP:", FP, "TN:", TN, "FN:", FN)
acc=np.array(a_f)
print("acc_avg: ", round(acc.mean(),3),"acc_std: ",round(acc.std(),3))

print()
Prec=np.array(Pr_f)
print("Prec(P)_avg: ", round(Prec.mean(),3),"Prec(P)_std: ",round(Prec.std(),3))


Prec1=np.array(Pr1_f)
print("Prec(N)_avg: ", round(Prec1.mean(),3),"Prec(N)_std: ",round(Prec1.std(),3))

print()
Recal=np.array(Re_f)
print("Recall(P)_avg: ", round(Recal.mean(),3),"Recall(P)_std: ",round(Recal.std(),3))
Recal1=np.array(Re1_f)
print("Recall(N)_avg: ", round(Recal1.mean(),3),"Recall(N)_std: ",round(Recal1.std(),3))

print()
F1_score=np.array(F1_f)
print("F1(P)_avg: ", round(F1_score.mean(),3),"F1(P)_std: ",round(F1_score.std(),3))
F1_score1=np.array(F11_f)
print("F1(N)_avg: ", round(F1_score1.mean(),3),"F1(N)_std: ",round(F1_score1.std(),3))

print()
HSS_score=np.array(HSS_1_f)
print("HSS(P)_avg: ", round(HSS_score.mean(),3),"HSS(P)_std: ",round(HSS_score.std(),3))


print()
HSS_score2=np.array(HSS_2_f)
print("HSS_2(P)_avg: ", round(HSS_score2.mean(),3),"HSS_2(P)_std: ",round(HSS_score2.std(),3))

print()
GS_score=np.array(GS_f)
print("GS(P)_avg: ", round(GS_score.mean(),3),"GS(P)_std: ",round(GS_score.std(),3))

print()
TSS_score=np.array(TSS_f)
print("TSS(P)_avg: ", round(TSS_score.mean(),3),"TSS(P)_std: ",round(TSS_score.std(),3))



