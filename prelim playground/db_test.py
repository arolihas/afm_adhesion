import pandas
from os import listdir, getcwd
from os.path import isfile, join

path = getcwd()
files=[f for f in listdir(path) if isfile(join(path,f)) and '.csv' in f]
data_list = {}
for f in files:
	data_list[f] = pandas.read_csv(f)

data_list[files[0]] = data_list[files[0]].drop([' No'], axis=1)
data_list[files[0]] = data_list[files[0]].rename(columns={'ZIndex':'No'})
for f in files:
	data_list[f] = data_list[f].reindex(columns=['No', 'zSensr', 'defl'])

approach_curve = []
approach_constants = []
full_curve = []
full_constants = []
for k in data_list:
    if 'Apr' in k:
        approach_curve.append(data_list[k])
        approach_constants.append(float(k[-9:-4]))
    else:
        full_curve.append(data_list[k])
        full_constants.append(float(k[-9:-4]))
for i in range(len(approach_curve)):
    approach_curve[i].defl *= approach_constants[i]

import matplotlib.pyplot as plt
n = 0
d = approach_curve[n].defl * approach_constants[n]
z = approach_curve[n].zSensr
t = approach_curve[n].No

from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
import numpy as np

X = pandas.concat([d,z], axis=1).values
scaler = StandardScaler()
X = scaler.fit_transform(X)

db = DBSCAN()
db.fit(X)
plt.figure(figsize=(9, 8))
plt.suptitle('Approach Curve ' + str(n))
plt.subplot(211)
plt.scatter(t, d, c=db.labels_)
plt.ylabel('Force (N)')
plt.subplot(212)
plt.plot(t, z, c=db.labels_)
plt.ylabel('Z position')
plt.show()
db.fit()
