import os
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler

def normalize_save(Data):
    train_data_np = np.array(Data, dtype=float)

    mean = np.mean(train_data_np, axis=0, keepdims=True)
    std = np.std(train_data_np, axis=0, ddof=1, keepdims=True)
    index = np.where(std == 0) 
    std[index] = 1e-7
    train_data_np = (train_data_np - mean) / std
    return Data

DATA_SET = 'RPI7317' #choose any of the datasets here

script_dir, script_name = os.path.split(os.path.abspath(sys.argv[0]))
parent_dir = os.path.dirname(script_dir)
raw_dir = 'D:/CSU/RPI-MD-/RPI-MD/data/raw_data' 
inputfile = raw_dir + DATA_SET + "RPI7317.csv"

data_start = pd.read_csv(inputfile)

label_P = np.ones(int('243'))
label_N = np.zeros(int('245'))

label_start = np.hstack((label_P, label_N))
label = np.array(label_start)
data1 = np.array(data_start)
data = data1[:, 1:]
data_nor = normalize_save(data)

Zongshu = data_nor
RNA_shu = Zongshu[:,0:674]
pro_shu = Zongshu[:,674:]

def MDS_select(data,n_components=300):
    embedding = SelectFromModel(n_components=n_components)
    new_data = embedding.fit_transform(data)
    return new_data

new_RNA_data = MDS_select(RNA_shu,n_components=39)
new_pro_data = MDS_select(pro_shu,n_components=135)

'''savepath'''
[m1,n1] = np.shape(new_RNA_data)
[m2,n2] = np.shape(new_pro_data)
cc1 = str(n1)
cc2 = str(n2)
data_new = np.hstack((new_RNA_data,new_pro_data))
optimal_features = pd.DataFrame(data=data_new)
savep = '_MDS_'+cc1+'_'+cc2 +'.csv'
outpath = 'D:\CSU\RPI-MD-\RPI-MD'
optimal_features.to_csv(outpath + DATA_SET + savep)