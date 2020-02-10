import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import itertools

df=pd.read_csv('tss_contribution_final_logist.csv')
col1=df['Col1']
data3=df['Svm-Tss'].sort_values()
err3=df['Svm-std'].sort_values()
data2=df['Logistic-tss'].sort_values()
err2=df['Logistic-std'].sort_values()


cols=df['Col1'][data2.index]
colsa=df['Col1'][data3.index]

x_label = cols
x_tick = np.arange(len(cols))
xs_label = colsa
xs_tick = np.arange(len(colsa))


#fig = plt.figure(figsize = (10,10))
#plt.bar(x_tick, data2, align = 'center', alpha = 0.5)
#fig.suptitle("CO2 Emissions by Electric Power Sector", fontsize= 25)
#plt.xticks(x_tick, x_label, rotation = 70, fontsize = 10)
#plt.yticks(fontsize = 20)
#plt.xlabel('Magnetic Parameters', fontsize = 15)
#plt.ylabel("TSS", fontsize = 15)
#plt.show()
#plt.savefig('TSS.png', dpi=500) 


# Plots
plt.ion()
fig, ((ax1, ax2)) = plt.subplots(2, 1, figsize=(20,20))
ax1.bar(x_tick, data2, yerr=err2, align='center', alpha=0.9, ecolor='black', capsize=5)
#ax1.set_xticks(x_tick, x_label, rotation = 70, fontsize = 10)
ax1.set_xticks(x_tick)
ax1.set_xticklabels(x_label, rotation=70, fontsize=14)
#ax1.set_yticks(labelsize=14)
#ax1.set_yticks(fontsize = 20)
ax1.tick_params(axis="y", labelsize=20)
ax1.set_xlabel('Magnetic Parameters', fontsize = 16)
ax1.set_ylabel("TSS", fontsize = 16)
ax = plt.gca();
ax.FontSize = 16;
ax2.bar(xs_tick, data3, yerr=err3, align='center', alpha=0.9, ecolor='black', capsize=5)
ax2.set_xticks(xs_tick)
ax2.set_xticklabels(xs_label, rotation=70, fontsize=14)
#ax2.set_yticklabels(fontsize = 20)
ax2.tick_params(axis="y", labelsize=20)
ax2.set_xlabel('Magnetic Parameters', fontsize = 16)
ax2.set_ylabel("TSS", fontsize = 16)
fig.savefig('fig10a.png', format='png', dpi=500)

fig = plt.figure(figsize = (10,10))
x_label = cols
x_tick = np.arange(len(cols))
plt.bar(x_tick, data2, align = 'center', alpha = 0.9)
#fig.suptitle("CO2 Emissions by Electric Power Sector", fontsize= 25)
plt.xticks(x_tick, x_label, rotation = 70, fontsize = 10)
plt.yticks(fontsize = 20)
plt.xlabel('Magnetic Parameters', fontsize = 15)
plt.ylabel("TSS", fontsize = 15)
#plt.show()
plt.savefig('TSS.pdf') 
