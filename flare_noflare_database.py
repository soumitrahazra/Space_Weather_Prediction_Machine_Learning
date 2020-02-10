import numpy as np, matplotlib.pylab as plt, matplotlib.mlab as mlab 
import pandas as pd, mpld3, requests, urllib, json 
import sunpy, sunpy.map 
from datetime import datetime as dt_obj 
from bs4 import BeautifulSoup 
from datetime import timedelta 
from mpld3 import plugins 
from sunpy.time import TimeRange 
import sunpy.instr.goes 
from scipy.stats import pearsonr as pearse 
   
   
search_param = {"eventName": "Solar Flare",  
                    "catalogName": "ALL", 
                    "startSearchDate": "2010-05-01", 
                    "endSearchDate": "2019-03-14"} 
     
t_start = search_param['startSearchDate'] 
t_end   = search_param['endSearchDate'] 


# Grab all the data from the GOES database 
#dpeak= listofresults[i]['goes_location'][0] (for Goes location checking)
time_range = TimeRange(t_start,t_end) 
listofresults = sunpy.instr.goes.get_goes_event_list(time_range,'M')   


df=pd.DataFrame()
ar_data=[]
xflare_data=[]
peaktime_data=[]
goescoord1=[]
goescoord2=[]

for i in range(len(listofresults)): 
    ds = listofresults[i]['noaa_active_region'] 
    dsa = listofresults[i]['goes_class']
    dpeak= listofresults[i]['peak_time']
    goeseventcord1= listofresults[i]['goes_location'][0]
    goeseventcord2= listofresults[i]['goes_location'][1]
    ar_data.append(ds)
    xflare_data.append(dsa) 
    goescoord1.append(goeseventcord1)
    goescoord2.append(goeseventcord2)
    peaktime_data.append(dpeak)
    

df['Class']= xflare_data
df['AR']= ar_data
df['Goes Cord1']= goescoord1
df['Goes Cord2']= goescoord2 
df['Peak Time']= peaktime_data
df.columns = ['Class', 'AR', 'Goes Cord1', 'Goes Cord2', 'Peak Time']

df=df[df['AR'] >10000]


listofactiveregions = list(df['AR'].values.flatten())
flareclass= list(df['Class'].values.flatten())
goescord1= list(df['Goes Cord1'].values.flatten())
goescord2= list(df['Goes Cord2'].values.flatten())
peaktime= list(df['Peak Time'].values.flatten())
answer = pd.read_csv('http://jsoc.stanford.edu/doc/data/hmi/harpnum_to_noaa/all_harps_with_noaa_ars.txt',sep=' ')
n_elements= len(listofactiveregions)
positive_harp_data=[]
positive_data=[]

for i in range(n_elements): 
    idx = answer[answer['NOAA_ARS'].str.contains(str(int(listofactiveregions[i])))] 
    if (idx.empty == True):
          positive_harp= np.nan
          positive_harp_data.append(positive_harp) 
          continue 
    positive_harp= idx.HARPNUM.values[0] 
    positive_harp_data.append(positive_harp) 
    

harpnum=pd.read_csv('total_harp_number.csv')
dfms=pd.DataFrame()
dfms['Positive']= positive_harp_data
dfm=dfms.dropna()
dfs=pd.DataFrame()
dfs['Positive']= harpnum['HARP Number']

dfd=pd.DataFrame()
dfd= pd.merge(dfs, dfm, how='outer', indicator=True)

rows_in_dfs_not_in_dfm = dfd[dfd['_merge']=='left_only'][dfs.columns]


df_positive=pd.DataFrame()
df_positive['Class']= flareclass
df_positive['NOAA'] = listofactiveregions
df_positive['HARP_Num'] = dfm['Positive']
df_positive['Goes Cord1']= goescord1
df_positive['Goes Cord2']= goescord2
df_positive['Peak Time']= peaktime
df_positive.columns=['Class', 'NOAA AR', 'Harp No', 'Goes Cord1', 'Goes Cord2', 'Peak Time']
df_positive=df_positive[abs(df_positive['Goes Cord2']) <72] 
df_positive=df_positive[abs(df_positive['Goes Cord1']) <72] 
df_positive=df_positive.dropna()
dfs_positive=pd.DataFrame()
dfs_positive['Class']= df_positive['Class']
dfs_positive['NOAA']= df_positive['NOAA AR']
dfs_positive['HARP_Num']= df_positive['Harp No']
dfs_positive['Peak Time']= df_positive['Peak Time']
dfs_positive.columns=['Class', 'NOAA AR', 'Harp No', 'Peak Time']

df_negative= pd.DataFrame()
df_negative['HARP_Num']= rows_in_dfs_not_in_dfm['Positive']
df_negative.columns= ['Harp No']


dfs_positive.to_csv('final_flare_positive_harp_data3.csv', sep='\t', index = False)
df_positive.to_csv('final_flare_positive_harp_check_data3.csv', sep='\t', index = False)
df_negative.to_csv('final_noflare_negative_harp_data3.csv', sep='\t', index = False)

