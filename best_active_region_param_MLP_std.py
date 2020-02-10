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

from keras.models import Model
from keras.models import Sequential
from keras.layers import Dense, Activation
from keras import regularizers

#Different parameters: col1:col7- USFLUX, col8:col14- MEANGAM
#col15:col21- MEANGBT, col22:col28- MEANGBZ, col29:col35- MEANGBH, col36:col42- MEANJZD
#col43:col49- TOTUSJZ, col50:col56- MEANALP, col57:col63- MEANJZH, col64:col70- TOTUSJH
#col71:col77- ABSNJZH, col78:col84- SVANCPP, col85:col91- MEANPOT, col92:col98- TOTPOT
#col99:col105- MEANSHR, col106:col112- SHRGT45, col113:col119- AREA_ACR


model1=Sequential()
#model1.add(Dense(units=3,input_dim=112, activation="relu",kernel_regularizer=regularizers.l2(0.01),name="first"))
#model1.add(Dense(units=28,activation="relu",name="second"))
#model1.add(Dense(units=14,activation="relu",name="third"))
#model1.add(Dense(units=7,activation="relu",name="fourth"))
#model1.add(Dense(units=1,activation="sigmoid"))
#model1.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

model1.add(Dense(units=3,input_dim=7, activation="relu",kernel_regularizer=regularizers.l2(0.01),name="first"))
#model1.add(Dense(units=35,activation="relu",name="second"))
#model1.add(Dense(units=15,activation="relu",name="third"))
#model1.add(Dense(units=7,activation="relu",name="fourth"))
model1.add(Dense(units=1,activation="sigmoid"))
model1.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

df=pd.read_csv('shazra_final_ready_another_loop24span12.csv')
#df=df.replace({0:np.nan})
#df=df.fillna(df.mean())
#df=df.ffill()
df=df.fillna(df.mean())

xs, ys = df.loc[:, 'col113':'col119'], df['col120']
class_weight={0:0.25, 1:0.75}

def classification_report_with_accuracy_score(y_true, y_pred):
    #print (classification_report(y_true, y_pred) )# print classification report
    #print( confusion_matrix(y_true, y_pred))
    A = confusion_matrix(y_true, y_pred)
    TP = A[0,0]
    #print(TP)
    FN = A[0,1]
    #print(FN)
    FP = A[1,0]
    #print(FP)
    TN = A[1,1]
    #print(TN)
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
    
    
    return a,Pr,Pr1,Re,Re1,F1,F11,HSS_1, HSS_2, GS, TSS


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
for i in range(30):
    X_train, X_test, y_train, y_test = train_test_split(xs, ys, test_size = 0.3, stratify = ys)
    #standardization   (almost like normalizing)  
    scaler = preprocessing.MinMaxScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)
    model1.fit(X_train, y_train,epochs=500,batch_size=100, class_weight=class_weight, verbose=0)
    #model1.fit(X_res,y_res,epochs=500,batch_size=100)
    predicted=np.where(model1.predict(X_train)>0.5,1,0)    #on training data
    predicted_test=np.where(model1.predict(X_test)>0.5,1,0)  #on test data
    a,Pr,Pr1,Re,Re1,F1,F11,HSS_1, HSS_2, GS, TSS=classification_report_with_accuracy_score(y_test, predicted_test)
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




