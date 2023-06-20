# %% Import zeoPy tool
from zeoPy import zeoPy

# %% Import required packages

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os 
from sklearn.model_selection import train_test_split
import time
from scipy.optimize import nnls
from sklearn.metrics import r2_score

# %% Create an instance of the class
f = zeoPy()

# %% Get the unique segments of the IZA database       
          
uniqueSegmentsIZA,bounds = f.getUniqueSegments("../MaterialsInfo",removeIsolated=True)
# uniqueSegmentsIZA,bounds = f.getUniqueSegments("../MaterialsInfoFiltered",removeIsolated=True)

# %% Get the unique segments big enough to capture H2 and compute the distance matrix
uniqueSegmentsIZABig = uniqueSegmentsIZA[uniqueSegmentsIZA["Volume"] >= 7.2]
uniqueSegmentsIZABig.reset_index(inplace=True,drop=True)

Distance_Mat_IZA = f.computeDistanceMatrix(uniqueSegmentsIZABig)


# %% Compute the clustering
#Error threshold above which segments are not clustered together
distance_threshold = 1.0/0.99

Sim_Mat, finalClusteringIZA, trained_model = f.spectralClustering(Distance_Mat_IZA,distance_threshold,repeatedValues = True) 


#Remove duplicates from the results
nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)
nonRepeatedClusters = nonRepeatedIZAClustering.Cluster.to_numpy()

# %% Get the total volume of the materials capable of absorbing H2

#Get the total number of materials
materials = np.unique(uniqueSegmentsIZABig.Material.to_numpy())
numberOfMaterials = materials.shape[0]

materialsTotalVolume = np.zeros(numberOfMaterials)

for i in range(numberOfMaterials):
    mat = materials[i]
    materialFound = uniqueSegmentsIZABig[uniqueSegmentsIZABig["Material"] == mat]
    volumes = materialFound.Volume.to_numpy()
    reps = materialFound.TimesRepeated.to_numpy()
    volume = 0
    for j in range(materialFound.shape[0]):
        vol = volumes[j]
        rep = reps[j]
        volume += vol*rep
    materialsTotalVolume[i] = volume
materialsVolume = pd.DataFrame(materialsTotalVolume, index=materials,columns=["Volume"])

# %% Create the A matrix with the volumes
numberOfClusters = finalClusteringIZA.Cluster.max() + 1

A_mat = np.zeros((materials.shape[0],numberOfClusters))
nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)


for i in range(A_mat.shape[0]):
    mat = materials[i]
    materialFound = uniqueSegmentsIZABig[uniqueSegmentsIZABig["Material"] == mat]
    materialFound.reset_index(inplace=True,drop=True)
    
    for j in range(materialFound.shape[0]):
        segmentID = materialFound.ID[j]
        #Find the segment in the clustering
        cluster = nonRepeatedIZAClustering[nonRepeatedIZAClustering["ID"] == segmentID].to_numpy()[0,1]
        timesRepeated = uniqueSegmentsIZABig[uniqueSegmentsIZABig["ID"] == segmentID].TimesRepeated.to_numpy()[0]
        A_mat[i,cluster] = A_mat[i,cluster] + timesRepeated
        
A_mat = pd.DataFrame(A_mat,index=materials)

# %% Check results 

x, ress = nnls(A_mat.to_numpy(),materialsVolume.to_numpy().flatten(),maxiter=100000000000000000)

#Test 
pred = np.matmul(A_mat.to_numpy(),x)
print("R2= ",r2_score(pred,materialsVolume.to_numpy().flatten()))

compare = pd.DataFrame(pred,columns=["PredTotal"],index=materialsVolume.index)
compare["SimTotal"] = materialsVolume.to_numpy()

# %% Checking local volumes
nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)
nonRepeatedIZAClustering.reset_index(inplace=True,drop=True)
results = np.zeros((nonRepeatedIZAClustering.shape[0],6))

segmentList = []
for i in range(nonRepeatedIZAClustering.shape[0]):
    segmentID = nonRepeatedIZAClustering.ID[i]
    mat = segmentID.split("_")[0]
    segmentList.append(segmentID)
    segmentCluster = nonRepeatedIZAClustering.Cluster[i]
    clusteredAlone = False
    checkCluster = A_mat.to_numpy()[:,segmentCluster]
    shared = checkCluster[checkCluster > 0]
    if shared.shape[0] <=1:
        clusteredAlone = True
    segmentSimVol = uniqueSegmentsIZABig[uniqueSegmentsIZABig["ID"] == segmentID].Volume.to_numpy()[0]
    segmentTimesRepeated = uniqueSegmentsIZABig[uniqueSegmentsIZABig["ID"] == segmentID].TimesRepeated.to_numpy()[0]
    segmentPredVol = x[segmentCluster]
    results[i,0] = segmentSimVol
    results[i,1] = segmentPredVol
    results[i,2] = clusteredAlone
    results[i,3] = segmentTimesRepeated
    results[i,4] = compare.loc[mat].SimTotal
    results[i,5] = compare.loc[mat].PredTotal
    
    
    

resData = pd.DataFrame(results,columns=["SimVolume","PredVolume","ClusteredAlone","TimesRepeated","MaterialTotalVolume","MaterialPredVolume"],index=segmentList)
print("R2 Score - Local Values: ",r2_score(results[:,0],results[:,1]))
print("R2 Score - Total Values: ",r2_score(results[:,4],results[:,5]))
# %% Checking outliers



        
    
    