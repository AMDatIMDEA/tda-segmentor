# %% Import zeoPy tool
from zeoPy import zeoPy

# %% Import required packages

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os 
from sklearn.model_selection import train_test_split
import time

# %% 1. Create a folder filtered with the materials with adsorption data

f = zeoPy()
f.removeMaterialsWithoutAdsorption("../MaterialsInfo/", "../IZA_Hydrogen_Adsorption_gL.csv", "../MaterialsInfoFiltered/")


# %% Get the unique segments of the IZA database(filtered)        
          
uniqueSegmentsIZA,bounds = f.getUniqueSegments("../MaterialsInfo",removeIsolated=True)
#uniqueSegmentsIZA,bounds = f.getUniqueSegments("../MaterialsInfoFiltered",removeIsolated=True)


#Plot the number of unique segments per material
materialNames = uniqueSegmentsIZA.Material.to_numpy()
materialNamesUnique = np.unique(materialNames)

materialCountsIZA = np.zeros_like(materialNamesUnique)
for i in range(materialNamesUnique.shape[0]):
    count = (materialNames == materialNamesUnique[i]).sum()
    materialCountsIZA[i] = count
bins = np.arange(np.max(materialCountsIZA)+1)
plt.figure()
plt.hist(materialCountsIZA, bins=bins,edgecolor="black")
plt.title("IZA: Number of unique segments per material")
plt.xlabel("Number of Unique Segments")
plt.ylabel("# Materials")


# %% Get the unique segments big enough to capture H2 and compute the distance matrix
uniqueSegmentsIZABig = uniqueSegmentsIZA[uniqueSegmentsIZA["Volume"] >= 7.2]
uniqueSegmentsIZABig.reset_index(inplace=True,drop=True)

Distance_Mat_IZA = f.computeDistanceMatrix(uniqueSegmentsIZA)

# %% Compute the spectral clustering based on the similarity matrix of the combined databases

#Error threshold above which segments are not clustered together
distance_threshold = 1.0/0.99

Sim_Mat, finalClusteringIZA, trained_model = f.spectralClustering(Distance_Mat_IZA,distance_threshold,repeatedValues = True) 


#Remove duplicates from the results
nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)
nonRepeatedClusters = nonRepeatedIZAClustering.Cluster.to_numpy()

#Plot the number of unique segments per cluster
numberOfClusters = nonRepeatedIZAClustering.Cluster.max()
clusters = np.arange(numberOfClusters+1)
segmentsPerCluster = np.zeros_like(clusters)
for i in range(clusters.shape[0]):
    segments = nonRepeatedIZAClustering[nonRepeatedIZAClustering["Cluster"] == clusters[i]]
    segmentsNumber = segments.shape[0]
    segmentsPerCluster[i] = segmentsNumber
    
bins = np.arange(np.max(segmentsPerCluster) + 1)
plt.figure()
hist = plt.hist(segmentsPerCluster, bins=bins,edgecolor="black")
plt.title("Number of unique segments per cluster")
plt.xlabel("Number of Unique Segments in the cluster")
plt.ylabel("Clusters") 



# %% Compute the matrix A required to solve the equation Ax = b
Atotal = f.computeA(finalClusteringIZA, "../materialsVolume.csv")

#Matrix A with the values changed to set bigger contrasts in the Images
A_normalized = np.zeros_like(Atotal)
Atotal_A = Atotal.to_numpy()
for i in range(Atotal.shape[0]):
    for j in range(Atotal.shape[1]):
        currentValue = Atotal_A[i,j]
        if currentValue > 0:
            A_normalized[i,j] = 1e6

plt.figure()
plt.imshow(A_normalized,interpolation='none',cmap='YlGn')
plt.title("A matrix")
plt.xlabel("Clusters of unique segments")
plt.ylabel("Materials")

# %% Find out how many materials have more than one segment clustered alone
AtotalA = Atotal.to_numpy()
markMat = []
for i in range(Atotal.shape[0]):
    counter = 0
    for j in range(Atotal.shape[1]):
        value = AtotalA[i,j]
        if value > 0:
            cluster = AtotalA[:,j]
            sharedMat = cluster[cluster[:] > 0]
            if sharedMat.shape[0] == 1:
                counter += 1
    markMat.append(counter)
    
markMat = np.array(markMat)

onlyOne = markMat[markMat[:] == 1]
moreThanOne = markMat[markMat[:] > 1]
    
print("Mats with one segment clustered alone: ",onlyOne.shape[0],"/",AtotalA.shape[0])
print("Mats with more than one segment clustered alone: ",moreThanOne.shape[0],"/",AtotalA.shape[0])
    

# %% Solve the adsorption equations system Ax = b
xIZA, rmseIZA, r2IZA, rel_errorIZA,spearmanIZA,bIZA,predictedResultIZA,matInfoIZA,A_Adsorption_IZA,b_Adsorption_IZA = f.solveSystem("../IZA_Hydrogen_Adsorption_gL.csv", Atotal)


plt.figure()
plt.plot(bIZA,predictedResultIZA,'.')
plt.plot(bIZA,bIZA)
plt.xlabel("True adsorption")
plt.ylabel("Predicted adsorption")
title = "H2 Adsorption: R2= " + str(round(r2IZA,3)) + " Spearman= " + str(round(spearmanIZA,3)) 
plt.title(title)

# %% Check if there is any relationship between the local adsorption and the volume/number of conexions
nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)
nonRepeatedIZAClustering.reset_index(inplace=True,drop=True)

bIZA = b_Adsorption_IZA

segmentsInfo = np.zeros((uniqueSegmentsIZA.shape[0],3))
materialNames = []
counter = 0
for i in range(uniqueSegmentsIZA.shape[0]):
    currentSegment = uniqueSegmentsIZA.loc[i]
    material = currentSegment.Material
    try:
        #Check if we have adsorption data of the material
        bIZA.loc[material]
        #Segment identifier
        materialID = currentSegment.ID
        #Assigned cluster to the segment
        assignedCluster =nonRepeatedIZAClustering[nonRepeatedIZAClustering["ID"] == materialID].Cluster.to_numpy()[0]
        #Local adsorption of the segment
        assignedAdsorption = xIZA[assignedCluster]
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

# %% Checking errors in the outliers
material,bins,store,distance_mat,groups,uniques,uniqueSegments,volumes,repeats,steps,numberOfConexions = f.getUnique("../MaterialsInfoFiltered/SBE.csv", removeIsolated=True,minValue=bounds[0],maxValue=(bounds[1]))


# %% Check how many materials share all their segments with other materials

Aa = Atotal.to_numpy()
# Aa = A_with_adsorption.to_numpy()

#Materials with at least one segment connected to other materials
SharedMaterials = []
#Materials with all of their segments connected to other materials
superSharedMaterials = []

for i in range(Aa.shape[0]):
    shared = False
    totallyShared = False
    segments = []
    for j in range(Aa.shape[1]):
        current = Aa[i,j]
        #If the current cluster appears in the material
        if current > 0 :
            #Check the column corresponding to that segment to check if any other
            #materials share that segment
            checkColumn = Aa[:,j]
            sharedMats = checkColumn[checkColumn[:] > 0]
            
            #If more than one material share that segment
            if sharedMats.shape[0] > 1:
                shared = True
                SharedMaterials.append(shared)
                segments.append(True)
            else:
                segments.append(False)
                
        if j == (Aa.shape[1] - 1):
            SharedMaterials.append(shared)
    segments = np.array(segments)
    check = segments[segments[:] == False]
    if check.shape[0] == 0:
        totallyShared = True
    superSharedMaterials.append(totallyShared)
    

SharedMaterials = np.array(SharedMaterials)  
#Materials not sharing any segment with other materials          
nonShared = SharedMaterials[SharedMaterials[:] != True]
print("Not connected materials: ",len(nonShared), "/",len(Aa))

print("Materials with at least one segment connected: ", len(Aa) - len(nonShared), "/",len(Aa))

       
superSharedMaterials = np.array(superSharedMaterials)
posSuperSharedMaterials = superSharedMaterials[superSharedMaterials[:] == True]
print("Materials with all their segments connected: ",len(posSuperSharedMaterials),"/",len(Aa))  
          

# %% Create a Test Set to check if we have enough data to predict other materials adsorption

#Take the super shared Materials
# superConnected = A_with_adsorption[superSharedMaterials]
superConnected = Atotal[superSharedMaterials]
# notSuperConnected = A_with_adsorption[superSharedMaterials == False]
notSuperConnected = Atotal[superSharedMaterials == False]

nonSelected,test_set = train_test_split(superConnected,test_size=0.15,random_state=123)
print("Test set size: ", test_set.shape[0])

trainSet = notSuperConnected.append(nonSelected)

# %% Check results
xTrain, rmseTrain, r2Train, rel_errorTrain,spearmanTrain,bTrain,predictedResultTrain,matInfoTrain,A_Adsorption_Train,b_Adsorption_Train= f.solveSystem("../IZA_Hydrogen_Adsorption_gL.csv", trainSet)

pred = xTrain
print("Test set")

rmseTest, r2Test, rel_errorTest,spearmanTest,bTest,predictedResultTest = f.checkTestSet(test_set, "../IZA_Hydrogen_Adsorption_gL.csv", "../PCOD_H2_Adsorption.csv", pred)

# %% Graph true vs predicted

plt.figure()
plt.plot(bTrain,predictedResultTrain,'.')
plt.plot(bTrain,bTrain)
plt.xlabel("True adsorption")
plt.ylabel("Predicted adsorption")
plt.title("Train Set")

plt.figure()
plt.plot(bTest,predictedResultTest,'.')
plt.plot(bTest,bTest)
plt.xlabel("True adsorption")
plt.ylabel("Predicted adsorption")
plt.title("Test Set")

 
# %% Try only with the IZA database

start = time.time()

f7 = zeoPy()

uniqueSegmentsIZA,bounds = f7.getUniqueSegments("../MaterialsInfo",removeIsolated=(True))
# %%
#---------------------------------------------
# Test bank

#What if we only select the segments with a volume bigh enough to capture H2?
uniqueSegmentsIZABig = uniqueSegmentsIZA[uniqueSegmentsIZA["Volume"] >= 7.0]
Distance_Mat_IZA = f7.computeDistanceMatrix(uniqueSegmentsIZABig)
dist_Threshold = 1.0/0.98 #Default:1.0/0.99
dist_Threshold = 1.0/0.99
Sim_Mat, finalClusteringIZA, trained_model = f7.spectralClustering(Distance_Mat_IZA,dist_Threshold,repeatedValues = True) 

A_IZA = f7.computeA(finalClusteringIZA, "../materialsVolume.csv")
xIZA, rmseIZA, r2IZA, rel_errorIZA,spearmanIZA,bIZA,predictedResultIZA,matInfoIZA,A_Adsorption_IZA,b_Adsorption_IZA = f7.solveSystem("../IZA_Hydrogen_Adsorption_gL.csv", A_IZA)

#Matrix A with the values changed to set bigger contrasts in the Images
A_normalized = np.zeros_like(A_Adsorption_IZA)
A_Adsorption_IZA_A = A_Adsorption_IZA.to_numpy()
for i in range(A_Adsorption_IZA.shape[0]):
    for j in range(A_Adsorption_IZA.shape[1]):
        currentValue = A_Adsorption_IZA_A[i,j]
        if currentValue > 0:
            A_normalized[i,j] = 1e6

plt.figure()
plt.imshow(A_normalized,interpolation='none',cmap='YlGn')
plt.title("A matrix")
plt.xlabel("Unique segments clusters")
plt.ylabel("Materials")

trueAdsorption = b_Adsorption_IZA.to_numpy()

plt.figure()
plt.plot(trueAdsorption,predictedResultIZA,'.')
plt.plot(trueAdsorption,trueAdsorption)
plt.xlabel("Simulated Adsorption")
plt.ylabel("Predicted Adsorption")
title = "H2 Adsorption: R2= " + str(r2IZA) + " Spearman= " + str(spearmanIZA) 
plt.title(title)

# Check how many segments we have per Material
materialNames = uniqueSegmentsIZABig.Material.to_numpy()
materialNamesUnique = np.unique(materialNames)

materialCountsIZA = np.zeros_like(materialNamesUnique)
for i in range(materialNamesUnique.shape[0]):
    count = (materialNames == materialNamesUnique[i]).sum()
    materialCountsIZA[i] = count
bins = np.arange(np.max(materialCountsIZA)+1)
plt.figure()
plt.hist(materialCountsIZA, bins=bins,edgecolor="black")
plt.title("Number of unique segments per material not removing isolated segments")
plt.xlabel("Number of Unique Segments")
plt.ylabel("Materials")

nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)
nonRepeatedClusters = nonRepeatedIZAClustering.Cluster.to_numpy()

numberOfClusters = nonRepeatedIZAClustering.Cluster.max()
clusters = np.arange(numberOfClusters+1)
segmentsPerCluster = np.zeros_like(clusters)
for i in range(clusters.shape[0]):
    segments = nonRepeatedIZAClustering[nonRepeatedIZAClustering["Cluster"] == clusters[i]]
    segmentsNumber = segments.shape[0]
    segmentsPerCluster[i] = segmentsNumber
    
bins = np.arange(np.max(segmentsPerCluster) + 1)
plt.figure()
hist = plt.hist(segmentsPerCluster, bins=bins,edgecolor="black")
plt.title("Number of unique segments per cluster")
plt.xlabel("Number of Unique Segments in the cluster")
plt.ylabel("Clusters") 


#----
# Local adsorption graphs
segWithAds = xIZA[xIZA > 0]
segWithAds =np.sort(segWithAds)

bins = np.linspace(0,np.max(xIZA),30)
plt.figure()
plt.hist(xIZA,bins=bins,edgecolor="black")
plt.hist(xIZA[xIZA==0],bins=bins,edgecolor="black")
labels = ["Cluster with Adsorption","Clusters with Zero adsorption"]
plt.legend(labels)
plt.title("Local adsorption of the clusters")
plt.xlabel("Local adsorption")
plt.ylabel("Number of clusters")


#---------------------------------------------
# %%
Distance_Mat_IZA = f7.computeDistanceMatrix(uniqueSegmentsIZA)
dist_Threshold = 1.0/0.98 #Default:1.0/0.99
dist_Threshold = 1.0/0.99
Sim_Mat, finalClusteringIZA, trained_model = f7.spectralClustering(Distance_Mat_IZA,dist_Threshold,repeatedValues = True) 
A_IZA = f7.computeA(finalClusteringIZA, "../materialsVolume.csv")
xIZA, rmseIZA, r2IZA, rel_errorIZA,spearmanIZA,bIZA,predictedResultIZA,matInfoIZA,A_Adsorption_IZA,b_Adsorption_IZA = f7.solveSystem("../IZA_Hydrogen_Adsorption_gL.csv", A_IZA)

end = time.time()

print("Time elapsed for the IZA computation ",end-start)


# %% Check how many segments we have per Material
materialNames = uniqueSegmentsIZA.Material.to_numpy()
materialNamesUnique = np.unique(materialNames)

materialCountsIZA = np.zeros_like(materialNamesUnique)
for i in range(materialNamesUnique.shape[0]):
    count = (materialNames == materialNamesUnique[i]).sum()
    materialCountsIZA[i] = count
bins = np.arange(np.max(materialCountsIZA)+1)
plt.figure()
plt.hist(materialCountsIZA, bins=bins,edgecolor="black")
plt.title("Number of unique segments per material not removing isolated segments")
plt.xlabel("Number of Unique Segments")
plt.ylabel("Materials")

# %% Check how many segments we have per cluster

nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)
nonRepeatedClusters = nonRepeatedIZAClustering.Cluster.to_numpy()

numberOfClusters = nonRepeatedIZAClustering.Cluster.max()
clusters = np.arange(numberOfClusters+1)
segmentsPerCluster = np.zeros_like(clusters)
for i in range(clusters.shape[0]):
    segments = nonRepeatedIZAClustering[nonRepeatedIZAClustering["Cluster"] == clusters[i]]
    segmentsNumber = segments.shape[0]
    segmentsPerCluster[i] = segmentsNumber
    
bins = np.arange(np.max(segmentsPerCluster) + 1)
plt.figure()
hist = plt.hist(segmentsPerCluster, bins=bins,edgecolor="black")
plt.title("Number of unique segments per cluster")
plt.xlabel("Number of Unique Segments in the cluster")
plt.ylabel("Clusters") 

# %% Create an Image with the Materials and Segments relations

#Matrix A with the values changed to set bigger contrasts in the Images
A_normalized = np.zeros_like(A_Adsorption_IZA)
A_Adsorption_IZA_A = A_Adsorption_IZA.to_numpy()
for i in range(A_Adsorption_IZA.shape[0]):
    for j in range(A_Adsorption_IZA.shape[1]):
        currentValue = A_Adsorption_IZA_A[i,j]
        if currentValue > 0:
            A_normalized[i,j] = 1e6

plt.figure()
plt.imshow(A_normalized,interpolation='none',cmap='YlGn')
plt.title("A matrix")
plt.xlabel("Unique segments clusters")
plt.ylabel("Materials")

# %% Check how many materials share each segment
A_A = A_IZA.to_numpy()
# A_A = A_Adsorption_IZA.to_numpy()

connectedMaterials = np.zeros(A_Adsorption_IZA_A.shape[1])

for j in range(A_Adsorption_IZA_A.shape[1]):
    column = A_A[:,j]
    pos_column = column[column[:] > 0]
    sharing = pos_column.shape[0]
    connectedMaterials[j] = sharing

bins = np.arange(np.max(connectedMaterials) + 1)
plt.figure()
plt.hist(connectedMaterials,bins=bins,edgecolor="black")
plt.title("Number of materials connected in each cluster")
plt.xlabel("Number of materials connected to the cluster")
plt.ylabel("Number of clusters")

# %% True vs Predicted Graph
trueAdsorption = b_Adsorption_IZA.to_numpy()

plt.figure()
plt.plot(trueAdsorption,predictedResultIZA,'.')
plt.plot(trueAdsorption,trueAdsorption)
plt.xlabel("Simulated Adsorption")
plt.ylabel("Predicted Adsorption")
title = "H2 Adsorption: R2= " + str(r2IZA) + " Spearman= " + str(spearmanIZA) 
plt.title(title)

# %% Check if there is any relationship between the local adsorption and the volume/number of conexions
nonRepeatedIZAClustering = finalClusteringIZA.drop_duplicates(keep='first',inplace=False)

segmentsInfo = np.zeros((uniqueSegmentsIZA.shape[0],4))
counter = 0
for i in range(uniqueSegmentsIZA.shape[0]):
    currentSegment = uniqueSegmentsIZA.loc[i]
    material = currentSegment.Material
    try:
        #Check if we have adsorption data of the material
        b_Adsorption_IZA.loc[material]
        #Segment identifier
        materialID = currentSegment.ID
        #Assigned cluster to the segment
        assignedCluster =nonRepeatedIZAClustering[nonRepeatedIZAClustering["ID"] == materialID].Cluster.to_numpy()[0]
        #Local adsorption of the segment
        assignedAdsorption = xIZA[assignedCluster]
        segmentsInfo[counter,0] = 
        segmentsInfo[counter,1] = assignedAdsorption
        segmentsInfo[counter,2] = currentSegment.Volume
        segmentsInfo[counter,3] = currentSegment.NumberOfConexions
        
        counter += 1
        
    except:
        print("Material without adsorption data: ",material)
        
segmentsInfo = segmentsInfo[~np.all(segmentsInfo == 0,axis=1)]

plt.figure()
plt.plot(segmentsInfo[:,1],segmentsInfo[:,0],'.')
plt.xlabel("Segments Volume")
plt.ylabel("Local Adsorption")
plt.title("Segments Volume vs Local Adsorption")
plt.figure()
plt.plot(segmentsInfo[:,2],segmentsInfo[:,0],'.')  
plt.xlabel("Segments Number of Connections")
plt.ylabel("Local Adsorption")
plt.title("Segments Number of Connections vs Local Adsorption")  

# %% Local adsorption graphs
segWithAds = xIZA[xIZA > 0]
segWithAds =np.sort(segWithAds)

bins = np.linspace(0,np.max(xIZA),30)
plt.figure()
plt.hist(xIZA,bins=bins,edgecolor="black")
plt.hist(xIZA[xIZA==0],bins=bins,edgecolor="black")
labels = ["Cluster with Adsorption","Clusters with Zero adsorption"]
plt.legend(labels)
plt.title("Local adsorption of the clusters")
plt.xlabel("Local adsorption")
plt.ylabel("Number of clusters")

            
            

# %% Check how many materials share each segment
Aa = A.to_numpy()

sharingMaterials = np.zeros((Aa.shape[1],2))
for j in range(Aa.shape[1]):
    currentSegment = Aa[:,j]
    sharedMaterials = currentSegment[currentSegment[:] > 0]
    numberOfSharedMaterials = sharedMaterials.shape[0]
    sharingMaterials[j,0] = j
    sharingMaterials[j,1] = numberOfSharedMaterials
    
plt.figure()
plt.hist(sharingMaterials[:,1])

notShared = sharingMaterials[sharingMaterials[:,1] > 1]
print("Not shared segments: ", notShared.shape[0],"/",sharingMaterials.shape[0])

# %% Try only with the PCOD database
start = time.time()

f8 = UniqueSegments()

uniqueSegmentsPCOD,bounds = f8.getUniqueSegments("../PCOD_MaterialsInfo")
Distance_Mat_PCOD = f8.computeDistanceMatrix(uniqueSegmentsPCOD)
dist_Threshold = 1.0/0.98 #Default:1.0/0.99
Sim_Mat, finalClusteringPCOD, trained_model = f8.spectralClustering(Distance_Mat_PCOD,dist_Threshold,repeatedValues = True) 
A = f8.computeA(finalClusteringPCOD, "../PCOD_MaterialsVolume.csv")
xPCOD, rmsePCOD, r2PCOD, rel_errorPCOD,spearmanPCOD,bPCOD,predictedResultPCOD,filtTotalPCOD,matInfoPCOD = f8.solveSystem("../PCOD_H2_Adsorption.csv", A)

end = time.time()

print("Time elapsed for the PCOD computation ",end-start)

# %% Check how many materials share each segment
Aa = A.to_numpy()

sharingMaterials = np.zeros((Aa.shape[1],2))
for j in range(Aa.shape[1]):
    currentSegment = Aa[:,j]
    sharedMaterials = currentSegment[currentSegment[:] > 0]
    numberOfSharedMaterials = sharedMaterials.shape[0]
    sharingMaterials[j,0] = j
    sharingMaterials[j,1] = numberOfSharedMaterials
    
plt.figure()
plt.hist(sharingMaterials[:,1])

notShared = sharingMaterials[sharingMaterials[:,1] > 1]
print("Not shared segments: ", notShared.shape[0],"/",sharingMaterials.shape[0])

# %% Ejercicios previos

# %% Compute the distance matrix with all the unique segments of the IZA  database
       
prueba = uniqueSegments[uniqueSegments["Volume"] >= 0.0]       
Distance_Mat = f.computeDistanceMatrix(prueba)
prueba = prueba.to_numpy()   

# %% Find the clusters of similar unique segments

dist_Threshold = 1.0/0.98 #Default:1.0/0.99
Sim_Mat, finalClustering, trained_model = f.spectralClustering(Distance_Mat,dist_Threshold,repeatedValues = True)  

# %%
A = f.computeA(finalClustering, "../materialsVolume.csv")
x, rmse, r2, rel_error,spearman,b,predictedResult,filtTotal,matInfo = f.solveSystem("../IZA_Hydrogen_Adsorption_gL.csv", A)

# %%
        
completed = f.getClusterHisto(finalClustering,uniqueSegments)   

# %%
'''
Do the same with the PCOD
'''

# %% Get the unique segments of the PCOD using the same thresholds in the histograms
f2 = UniqueSegments()    
uniqueSegments2,bounds2 = f.getUniqueSegments("../PCOD_MaterialsInfo",bounds=bounds)  
    
# %% Compare the unique segments in the IZA database with the PCOD database

maxVol = uniqueSegments.Volume.max()

f3 = UniqueSegments()
prevHist,currentHist,results, finalResults,compareResults,error = f3.comparePrevious(uniqueSegments2, completed,x,"../PCOD_MaterialsVolume.csv","../materialsVolume.csv","../PCOD_H2_Adsorption.csv")

        


# %%
        
# %% 2. Get the unique segments of the materials in the folder

f = UniqueSegments()
# material,bins, histSto,distance_mat,groups,unique,uniquehisto,volume = f.getUnique("MaterialsInfo/SFV.csv", True)
# material,bins, histSto,distance_mat,groups,unique,uniquehisto,vols,repeated,steps = f.getUnique("../MaterialsInfo/MFI.csv", True,minValue=-7.9447,maxValue=(-1.2))
# uniqueSegments = f.getUniqueSegments("MaterialsInfo")
uniqueSegments = f.getUniqueSegments("../MaterialsInfo")
# uniqueSegments = f.getUniqueSegments("../PCOD_MaterialsInfo")

# %%
# material,bins, histSto,distance_mat,groups,unique,uniquehisto,vols,repeated,steps = f.getUnique("../PCOD_MaterialsInfo/8289716.csv", False,minValue=-9.06688,maxValue=(-1.2))


# %% 3. Get segments with  volume greater than a threshold and compute the 
# distance matrix

# prueba = uniqueSegments[uniqueSegments["Volume"] >= 7.2]
# prueba = uniqueSegments[uniqueSegments["Volume"] >= 0.2]
prueba = uniqueSegments[uniqueSegments["Volume"] >= 0.0]

#Compute the distance matrix
Distance_Mat = f.computeDistanceMatrix(prueba)
prueba = prueba.to_numpy()

# %%
dist_Threshold = 1.0/0.98 #Default:1.0/0.99
Sim_Mat, finalClustering, trained_model = f.spectralClustering(Distance_Mat,dist_Threshold,repeatedValues = True)

# %%


# %%
# f = UniqueSegments()
A = f.computeA(finalClustering, "../PCODmaterialsVolume.csv")

PCOD_ads = pd.read_csv("../PCOD_H2_Adsorption.csv",index_col=(0))
temperature = str(275.9)
pressure = str(403.4)
col = temperature + "_K_" + pressure + "_bar"
my_col = PCOD_ads[col]


# x, rmse, r2, rel_error,spearman,b,predictedResult,filtTotal,matInfo = f.solveSystem("../IZA_Hydrogen_Adsorption_gL.csv", A)
x, rmse, r2, rel_error,spearman,b,predictedResult,filtTotal,matInfo = f.solveSystem("../PCOD_H2_Adsorption.csv", A)
print("FILTERED")
# x2, rmse2, r22, rel_error2,spearman2,b2,predictedResult2,filtTotal2,matInfo2 = f.solveSystemFiltered("../IZA_Hydrogen_Adsorption_gL.csv", A)
x2, rmse2, r22, rel_error2,spearman2,b2,predictedResult2,filtTotal2,matInfo2 = f.solveSystemFiltered("../PCOD_H2_Adsorption.csv", A)

# %%
plt.figure()
plt.plot(b,predictedResult,'.')
title = "H2 Adsorption: R2= " + str(r2) + " \nSpearman= " + str(spearman) + " \nRel Error= " + str(rel_error)
plt.title(title)

plt.figure()
plt.plot(b2,predictedResult2,'.')
title = "H2 Adsorption: R2= " + str(r22) + " \nSpearman= " + str(spearman2) + " \nRel Error= " + str(rel_error2)
plt.title(title)




 
