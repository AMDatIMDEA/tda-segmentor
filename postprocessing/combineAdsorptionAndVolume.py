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
from sklearn.metrics import mean_squared_error

# %% Create an instance of the class
f = zeoPy()
f.removeMaterialsWithoutAdsorption("../MaterialsInfo/", "../IZA_Hydrogen_Adsorption_gL.csv", "../MaterialsInfoFiltered/")


# %% Get the unique segments of the IZA database       
          
# uniqueSegmentsIZA,bounds = f.getUniqueSegments("../MaterialsInfo",removeIsolated=True)
# Only use the materials with adsorption data
uniqueSegmentsIZA,bounds = f.getUniqueSegments("../MaterialsInfoFiltered",removeIsolated=True)

# %% Get the unique segments big enough to capture H2 and compute the distance matrix
uniqueSegmentsIZABig = uniqueSegmentsIZA[uniqueSegmentsIZA["Volume"] >= 7.2]
uniqueSegmentsIZABig.reset_index(inplace=True,drop=True)

Distance_Mat_IZA = f.computeDistanceMatrix(uniqueSegmentsIZABig)

# %% Compute the clustering
#Error threshold above which segments are not clustered together
distance_threshold = 1.0/0.94

Sim_Mat, finalClusteringIZA, trained_model = f.spectralClustering(Distance_Mat_IZA,distance_threshold,repeatedValues = True) 


#Remove duplicates from the results
nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)
nonRepeatedClusters = nonRepeatedIZAClustering.Cluster.to_numpy()

# %% Get the total accessible volume of the materials in a supercell
materialsVolume = f.computeSuperCellsAccessibleVolume(uniqueSegmentsIZABig, "../materialsVolume.csv",superCell = True)

# %% Now we need to compute the matrix A to solve the equations. First, we perform the slice
# of the matrix A corresponding to the direct adsorption equations.

Aadsorption = f.computeA(finalClusteringIZA, "../materialsVolume.csv")

# %% Get the A matrix for the volume data
# Avolume = f.computeAVolume(uniqueSegmentsIZABig,finalClusteringIZA, "../materialsVolume.csv",superCell=True)

Adensity = f.computeAVolumeDensity(uniqueSegmentsIZABig,finalClusteringIZA,"../materialsVolume.csv","../IZA_Hydrogen_Adsorption_gL.csv",materialsVolume,superCell=True)

# %% Now get the b vectors in order to concatenate them
bAdsorption, bVolume,x, r2,spearman,rel_error,rmse,predictedResult = f.combineAdsorptionVolume(Aadsorption, Adensity, materialsVolume, "../IZA_Hydrogen_Adsorption_gL.csv")

# %% Check if there is any relationship between the local adsorption and the volume/number of conexions
nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)
nonRepeatedIZAClustering.reset_index(inplace=True,drop=True)

bIZA = bAdsorption

segmentsInfo = np.zeros((uniqueSegmentsIZABig.shape[0],3))
materialNames = []
counter = 0
for i in range(uniqueSegmentsIZABig.shape[0]):
    currentSegment = uniqueSegmentsIZABig.loc[i]
    material = currentSegment.Material
    try:
        #Check if we have adsorption data of the material
        bIZA.loc[material]
        #Segment identifier
        materialID = currentSegment.ID
        #Assigned cluster to the segment
        assignedCluster =nonRepeatedIZAClustering[nonRepeatedIZAClustering["ID"] == materialID].Cluster.to_numpy()[0]
        #Local adsorption of the segment
        assignedAdsorption = x[assignedCluster]
        materialNames.append(materialID)
        segmentsInfo[counter,0] = assignedAdsorption
        segmentsInfo[counter,1] = currentSegment.Volume
        segmentsInfo[counter,2] = currentSegment.NumberOfConexions
        
        # print(materialID,assignedAdsorption,currentSegment.Volume)
        
        counter += 1
        
    except:
        print("Material without adsorption data: ",material)
        
segmentsInfo = segmentsInfo[~np.all(segmentsInfo == 0,axis=1)]
segmentsInfo = pd.DataFrame(segmentsInfo,columns=["Adsorption","Volume","Conexions"],index=materialNames)

plt.figure()
plt.plot(segmentsInfo.Volume,segmentsInfo.Adsorption,'.')
plt.xlabel("Segments Volume")
plt.ylabel("Local Adsorption")
plt.title("Segments Volume vs Local Adsorption")

plt.figure()
plt.plot(segmentsInfo.Conexions,segmentsInfo.Adsorption,'.')  
plt.xlabel("Segments Number of Connections")
plt.ylabel("Local Adsorption")
plt.title("Segments Number of Connections vs Local Adsorption")  

plt.figure()
plt.plot(segmentsInfo.Conexions,segmentsInfo.Volume,'.')  
plt.xlabel("Segments Number of Connections")
plt.ylabel("Segments Volume")
plt.title("Segments Number of Connections vs Segments Volume") 

#%% Check RMSE in adsorption

predAds = predictedResult[:bAdsorption.shape[0]]

rmseAds = mean_squared_error(bAdsorption, predAds,squared=False)


