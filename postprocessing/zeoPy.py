#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class to cluster the segments of a material in similar groups. There are functions
to cluster them inside the same material and between materials.

@author: Jorge Zorrilla Prieto
"""

# %% 0. Import required Python packages

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os 
import shutil
from sklearn.metrics import r2_score
from scipy import stats
from sklearn.cluster import AgglomerativeClustering
from scipy.optimize import nnls
from sklearn.metrics import mean_squared_error
import time


# %% 1. Class


class zeoPy:
    
    def __init__(self):
        
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("Welcome to the Unique Segments Class")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
    
    def computeMaterialsHistograms(self, inputFolder,outputFolder, histogramSteps = 150, removeIsolated= False, bounds = None):
        #Files in the input folder
        inputFiles = os.listdir(inputFolder)
        #Number of files in the folder
        numberOfFiles = len(inputFiles)
        
        if os.path.exists(outputFolder):
            #To avoid a folder with old/unnecessary files
            shutil.rmtree(outputFolder,ignore_errors=True)
            print("Removed")
            os.mkdir(outputFolder)
        else:
            os.mkdir(outputFolder)
        
        print("Getting scalar bounds of the materials in the folder. Wait")
        #If not predefined bound values
        if bounds == None:
            #Find the bound values by reading all the files in the folder
            maxValue,minValue = self.getBounds(inputFolder)
        else:
            minValue= bounds[0]
            maxValue = bounds[1]
        #Number of steps of the histograms. A high number of steps could lead 
        # to higher errors
        steps = histogramSteps
        
        #Bin edges of the histograms
        bins = np.linspace(minValue,maxValue,steps)
        
        for i in range(numberOfFiles):
            #Current file
            inputFile = inputFolder + "/" + inputFiles[i]
            #File without extension to make it more readable
            file_without_extension = inputFiles[i].split(".")[0]
            
            #Read the input file
            material = pd.read_csv(inputFile)
            
            #Remove isolated regions if needed
            if(removeIsolated == True):
                material = material[material["NumberOfConexions"] >= 1]
                
            
            
            #Get the IDs of the segments of the material
            segmentsID = material.regionID.to_numpy()
            segmentsID = np.unique(segmentsID) #Remove repeated names
            
            #Store the histograms of the segments in order to compare them
            store = np.zeros((segmentsID.shape[0],steps))
            for j in range(segmentsID.shape[0]):
                segment = material[material["regionID"] == segmentsID[j]]
                
                histo = np.histogram(segment.Scalar,bins=bins)
                
                store[j,0] = segmentsID[j]
                store[j,1:] = histo[0]
            
            df = pd.DataFrame(data=store)
            df=df.rename(columns = {0:'Segment'})
            df.to_csv(outputFolder + file_without_extension + ".csv",index = False)
            print("Computed histograms of: ", i+1,"/",numberOfFiles," files")
            
            
            
            
            
            
            
            

    def getUnique(self,inputFile,removeIsolated,minValue = None,maxValue=None):
        '''
        From a .csv file, remove(or not) the isolated segments. Then compute
        the histogram of each segment based on their scalar values and compare
        them in order to find similar segments

        Parameters
        ----------
        inputFile : String(Path to the .csv file)
            File that stores the information of the segment obtained from the 
            C++ zeoTDA class \n
            
        removeIsolated : bool
            Remove segments with a number of conexiones equal to 0 \n
            
        minValue : double, optional
            Min value of the bins of the histograms. Default value is None so 
            the function will look for these bounds
        maxValue : double, optional
            Max value of the bins of the histograms. Default value is None so 
            the function will look for these bounds

        Returns
        -------
        material : DataFrame
            The readed material
        bins : TYPE
            DESCRIPTION.
        store : TYPE
            DESCRIPTION.
        distance_mat : DataFrame
            Distance Matrix of the segments of the material
        groups : TYPE
            DESCRIPTION.
        uniques : List of lists
            Unique groups with the segments ID they gather
        uniqueSegments : TYPE
            DESCRIPTION.
        volumes : TYPE
            DESCRIPTION.
        repeated : TYPE
            DESCRIPTION.
        steps : TYPE
            DESCRIPTION.

        '''
        
        #Read the input file
        material = pd.read_csv(inputFile)
        
        #Remove isolated regions if needed
        if(removeIsolated == True):
            material = material[material["NumberOfConexions"] >= 1]
        
        
        #Get the histogram bin edges for the material if needed
        if minValue == None:
            scalarMin = material.Scalar.min()
        else:
            scalarMin = minValue
        if maxValue == None:
            scalarMax = material.Scalar.max()
        else:
            scalarMax = maxValue
        
        #Number of steps of the histograms. A high number of steps could lead 
        # to higher errors
        steps = 150
        
        #Bin edges of the histograms
        bins = np.linspace(scalarMin,scalarMax,steps)
        
        #Get the IDs of the segments of the material
        segmentsID = material.regionID.to_numpy()
        segmentsID = np.unique(segmentsID) #Remove repeated names
        
        #Store the histograms of the segments in order to compare them
        store = np.zeros((segmentsID.shape[0],steps))
        for i in range(segmentsID.shape[0]):
            segment = material[material["regionID"] == segmentsID[i]]
            
            histo = np.histogram(segment.Scalar,bins=bins)
            
            store[i,0] = segmentsID[i]
            store[i,1:] = histo[0]
        
        
        distance_matrix = np.zeros((segmentsID.shape[0],segmentsID.shape[0]))
        #Compute the distance matrix between the segments
        for i in range(segmentsID.shape[0]):
            material1 = store[i,1:]
            for j in range(segmentsID.shape[0]):
                material2 = store[j,1:]
                
                #Different metrics. You can try a couple of them
                #----------------------------------------
                
                # error = np.linalg.norm(material2-material1)
                # error = (np.linalg.norm(material2-material1) / np.linalg.norm(material1)) * 100
                # error = stats.spearmanr(material1,material2)[0]
                error = r2_score(material1,material2)
                # error = 1 - spatial.distance.cosine(material1,material2)
                
                #----------------------------------------
                distance_matrix[i,j] = error
        
        #Store the results
        distance_mat = pd.DataFrame(data = distance_matrix,columns=segmentsID,index=segmentsID)
        
        #Get the groups of similar segments
        groups = []
        for i in range(distance_matrix.shape[0]):
            group = [segmentsID[i]]
            for j in range(distance_matrix.shape[1]):
                error = distance_matrix[i,j]
                if (i != j) and (error > 0.98): #Threshold to set similar segments
                    group.append(segmentsID[j])
            groups.append(np.sort(group))
         
        #Get the unique groups   
        uniques = []
        for i in range(len(groups)):
            group = groups[i] #Current group
            
            if i == 0:
                uniques.append(group)
            else:
                
                signal = False #Is array already there
                for j in range(len(uniques)):
                    compare = uniques[j]
                    equal = np.array_equal(group, compare)
                    if equal == True:
                        signal = True
                        
                if signal == False:
                    uniques.append(group)
        
                    
        
        #Check if there are groups containing others
        signals = [True] #Suppose all the groups to be unique
        signals = signals * len(uniques)
        
        for i in range(len(uniques)):
            currentGroup = uniques[i]
            if signals[i] == True:
                for j in range(len(uniques)):
                    
                    if signals[j] == True:
                        #Avoid comparing the same group
                        if i != j:
                            compareGroup = uniques[j]
                            #If the current group is smaller than the one we are comparing to
                            if currentGroup.shape[0] < compareGroup.shape[0]:
                                signal = True #Delete this group
                                for element in currentGroup:
                                    if (element in compareGroup) == False:
                                        signal = False #Groups different -> don't delete
                                        break
                                if signal == True:
                                    signals[i] = False
                                    uniques[i] = False
                                    break
        
        #Unique groups with the segments ID they gather
        uniques = [x for x in uniques if type(x) != bool]
        
        
                    
                    
        #Save one sample histogram of each unique group
        uniqueSegments = np.zeros((len(uniques),steps-1))
        volumes = np.zeros(len(uniques))
        numberOfCon = np.zeros(len(uniques))
        repeated = np.zeros(len(uniques))
        for i in range(len(uniques)):
            group = uniques[i]
            
            numberOfTimesRepeated = len(group)
            repeated[i] = numberOfTimesRepeated
            #Get one of the segments in the group to save its histogram
            sampleSegment = group[0]
            
            sampleRegion = material[material["regionID"] == sampleSegment]
            sampleVolume = np.unique(sampleRegion.Volume.to_numpy())[0]
            sampleNumberOfConex = np.unique(sampleRegion.NumberOfConexions.to_numpy())[0]
            
            volumes[i] = sampleVolume
            numberOfCon[i] = sampleNumberOfConex
            
            sampleHisto = store[store[:,0] == sampleSegment]
            sampleHisto = sampleHisto[:,1:].flatten()
            uniqueSegments[i,:] = sampleHisto
        
        
            
        
        return material,bins,store,distance_mat,groups,uniques,uniqueSegments, volumes,repeated,steps,numberOfCon
    
    def getUniqueSegments(self,inputFolder,removeIsolated=False, bounds=None):
        '''
        Use the UniqueSegments function for all the files in a folder

        Parameters
        ----------
        inputFolder : string
            Path to the folder where the files with the materials info is stored
            
        removeIsolated : bool, optional
            (False) Remove isolated segments of a material
            
        bounds : list, optional
            (None) Predefined values of the bounds for the computation of the histograms

        Returns
        -------
        df : DataFrame
            Dataframe with the results of the UniqueSegments functions stored
            
        list
            Bounds used for the histograms computation in case they are needed

        '''
        
        print("zeoPy: getUniqueSegments function ")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        #Files in the input folder
        inputFiles = os.listdir(inputFolder)
        #Number of files in the folder
        numberOfFiles = len(inputFiles)
        
        print("Getting scalar bounds of the materials in the folder. Wait")
        #If not predefined bound values
        if bounds == None:
            #Find the bound values by reading all the files in the folder
            maxValue,minValue = self.getBounds(inputFolder)
        else:
            minValue= bounds[0]
            maxValue = bounds[1]
        
        
        UniqueHistos = []
        MaterialsNames = []
        GroupsID = []
        Volumes = []
        Conexions = []
        Repeats = []
        for i in range(numberOfFiles):
            #Current file
            inputMaterial = inputFolder + "/" + inputFiles[i]
            #File without extension to make it more readable
            file_without_extension = inputFiles[i].split(".")[0]
            
            #Find the unique segments of the material
            material,bins,store,distance_mat,groups,uniques,uniqueSegments,volumes,repeats,steps,numberOfConexions = self.getUnique(inputMaterial, removeIsolated,minValue=(minValue),maxValue=(maxValue))
            
            #Number of unique segments in the material
            numberOfUnique = uniqueSegments.shape[0]
            
            segmentIDs = np.arange(numberOfUnique)
            #Repeat the material name one time for each of the unique segments
            materialNames = [file_without_extension] * numberOfUnique
            
            UniqueHistos.append(uniqueSegments)
            MaterialsNames.append(materialNames)
            GroupsID.append(segmentIDs)
            Volumes.append(volumes)
            Repeats.append(repeats)
            Conexions.append(numberOfConexions)
            
            print("Got unique segments of: ", i+1, "/",numberOfFiles, " materials")
        
        
        #Variables to store the results
        Histos = np.zeros((1,steps-1))
        Vol = np.zeros(1)
        Rep = np.zeros(1)
        Con = np.zeros(1)
        namesColumn = []
        segIDs = np.zeros((1))
        
        for i in range(len(UniqueHistos)):
            current = UniqueHistos[i]
            Histos = np.concatenate((Histos,current), axis=0)
            namesColumn = namesColumn + MaterialsNames[i]
            segIDs =  np.concatenate((segIDs,GroupsID[i]),axis=0)
            Vol = np.concatenate((Vol,Volumes[i]), axis=0)
            Con = np.concatenate((Con,Conexions[i]),axis=0)
            Rep = np.concatenate((Rep,Repeats[i]), axis=0)
            
        Histos = Histos[1:,:]
        segIDs = segIDs[1:]
        Vol = Vol[1:]
        Rep = Rep[1:]
        Con = Con[1:]
        
        IDs = []
        
        for i in range(len(namesColumn)):
            current = namesColumn[i] + "_" + str(int(segIDs[i]))
            IDs.append(current)
        
        df = pd.DataFrame(data = Histos)
        df.insert(0,"Material",namesColumn)
        df.insert(1,"GroupID",segIDs)
        df.insert(2,"Volume",Vol)
        df.insert(3,"TimesRepeated",Rep)
        df.insert(4,"ID",IDs)
        df.insert(5,"NumberOfConexions",Con)
        
        self.MaterialsData = df
            
            
        return df,[minValue,maxValue]
            
            
    def groupBySize(self,df,sizeThreshold):
        df = df.copy()
        #Number of points of each segment
        values = df.Volume.to_numpy()
        #Sort the values in ascending order
        sorted_values = np.sort(values)
        #Store the limits between size clusters
        limits = []

        #Previous value to compare(initially set it to negative to ensure that the first
        # value is a limit)
        prev_value = sorted_values[0]
        limits.append(prev_value)
        for i in range(sorted_values.shape[0]):
            if i != 0:
                current = sorted_values[i] #Current value compared
                difference = current - prev_value #Difference between current and prev
                if difference > sizeThreshold: #If the difference is greater than threshold
                    limits.append(current)
                prev_value = current
                
        groups = []
        
        for i in range(df.shape[0]):
            currentValue = df.Volume[i]
            
            group = 0
            for j in range(1,len(limits)):
                currentMax = limits[j]
                if currentValue < currentMax:
                    groups.append(group)
                    break
                else:
                    if currentValue == limits[-1]:
                        groups.append(len(limits)-1)
                        break
                    else:
                        
                        group += 1
                    
            
                
        df.insert(4,"SizeGroup",groups)
                
        
        return limits, df
            
            
        
        
        
    def getBounds(self,inputFolder):
        '''
        Get the bound values of all the files in a folder

        Parameters
        ----------
        inputFolder : string
            File Path to the input folder

        Returns
        -------
        maximumValue : float
            Files max value
        minimumValue : float
            Files min value

        '''
        
        
        #Input files in the folder
        inputFiles = os.listdir(inputFolder)
        #Number of files in the folder
        numberOfFiles = len(inputFiles)
        
        # Max/min values of the files
        maxValues = np.zeros((numberOfFiles))
        minValues = np.zeros((numberOfFiles))
        for i in range(numberOfFiles):
            #Read the input file
            material = pd.read_csv(inputFolder + "/" + inputFiles[i])
            #Find max/min values in the file
            maxValues[i] = material.Scalar.max()
            minValues[i] = material.Scalar.min()
            
            # print("Computed: ", i, "/",numberOfFiles)
        #Find max/min values of the folder
        maximumValue = np.max(maxValues)
        minimumValue = np.min(minValues)
        
        return maximumValue,minimumValue
    
    
    def computeDistanceMatrix(self,df):
        '''
        Compute the distance matrix between all the unique segments of 
        the materials. It uses their histograms to compute them

        Parameters
        ----------
        df : DataFrame
            Data from the getUniqueSegments() function with the information of 
            the unique segments of the materials

        Returns
        -------
        distance_mat : DataFrame
            Distance Matrix of all the unique segments of the materials in the 
            database

        '''
        print("zeoPy: Distance Matrix Function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        
        materials = df.Material.to_numpy()
        groupsid = df.GroupID.to_numpy()
        full_names = []
        for i in range(df.shape[0]):
            name = materials[i] + "_" + str(int(groupsid[i]))
            full_names.append(name)
        
        dfA = df.to_numpy()
        dfA = dfA[:,5:]
        
        distance_matrix = np.zeros((dfA.shape[0],dfA.shape[0]))
        #Compute the distance matrix between the segments
        for i in range(dfA.shape[0]):
            material1 = dfA[i,:]
            for j in range(dfA.shape[0]):
                material2 = dfA[j,:]
                # error = np.linalg.norm(material2-material1)
                # error = (np.linalg.norm(material2-material1) / np.linalg.norm(material1)) * 100
                # error = stats.spearmanr(material1,material2)[0]
                error = r2_score(material1,material2)
                distance_matrix[i,j] = error
            print("Distance Matrix: Row Computed:",i+1,"/",dfA.shape[0])
        
        distance_mat = pd.DataFrame(data = distance_matrix,columns=full_names,index=full_names)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        return distance_mat
    
    def computeDistanceMatrixCombined(self,df1,df2):
        '''
        Compute the distance matrix between all the unique segments of 
        the materials. It uses their histograms to compute them

        Parameters
        ----------
        df1 : DataFrame
            DataFrame with the unique segments information of the first database.
            You can obtain it from the getUniqueSegments function
            
        df2 : DataFrame
            DataFrame with the unique segments information of the second database.
            You can obtain it from the getUniqueSegments function

        Returns
        -------
        distance_mat : DataFrame
            Distance Matrix of all the segments with the R2 values between 
            their histograms

        '''
        #Computation start time
        start = time.time()
        
        print("zeoPy: Distance Matrix Combined Function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        #Combine both dataframes
        df = pd.concat([df1,df2])
        #Reset index of the new dataframe to avoid bugs
        df.reset_index(inplace=True, drop=True)
        
        
        #Get the histograms of the combined dataframe
        dfA = df.to_numpy()
        dfA = dfA[:,5:]
        
        #Variable to store the results
        distance_matrix = np.zeros((dfA.shape[0],dfA.shape[0]))
        
        #Compute the distance matrix between the segments
        for i in range(dfA.shape[0]):
            material1 = dfA[i,:]
            for j in range(dfA.shape[0]):
                material2 = dfA[j,:]
                
                #We use the R2 score to compare the histograms
                error = r2_score(material1,material2)
                #Store the results
                distance_matrix[i,j] = error
            print("Distance Matrix: Row Computed:",i+1,"/",dfA.shape[0])
        
        #Create the dataframe of the distance matrix to attach each result to its material
        distance_mat = pd.DataFrame(data = distance_matrix,columns=df.ID,index=df.ID)
        
        end = time.time()
        print("Time elapsed to compute the distance matrix: ", (end - start))
        
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        
        return distance_mat
        
      
    
    def spectralClustering(self,distanceMatrix,distanceThreshold,repeatedValues = False):
        '''
        Compute the agglomerative clustering(spectral clustering) of the previously computed
        distance matrix(DistanceMatrixCombined function).

        Parameters
        ----------
        distanceMatrix : DataFrame
            Distance Matrix of the unique segments of the  dataset
            (from the DistanceMatrix function)
            
        distanceThreshold : float
            Similarity value above which the segments are NOT clustered
            
        repeatedValues : bool, optional
            Take into account how many times each segment is repeated
            in a unit cell. The default is False.

        Returns
        -------
        similarityMatrix : Numpy Array
            Similarity Matrix of the unique segments of the dataset
        clustering : DataFrame
            Clustering information of the unique segments of the dataset
           
        '''
        
        
        
        print("zeoPy: Spectral Clustering Function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        Sim_Mat = np.zeros_like(distanceMatrix)
        
        columns = distanceMatrix.columns.to_numpy()
        Distance_A = distanceMatrix.to_numpy()

        for i in range(Sim_Mat.shape[0]):
            for j in range(Sim_Mat.shape[1]):
                if i == j:
                    Sim_Mat[i,j] = 0
                else:
                    if Distance_A[i,j] <= 0:
                        Sim_Mat[i,j] = 100
                    else:
                        Sim_Mat[i,j] = 1.0 / float(Distance_A[i,j])

        for i in range(Sim_Mat.shape[0]):
            for j in range(Sim_Mat.shape[1]):
                if i > j:
                    if Sim_Mat[i,j] >= Sim_Mat[j,i]:
                        
                        Sim_Mat[i,j] = Sim_Mat[j,i]
                    elif Sim_Mat[i,j] < Sim_Mat[j,i]:
                        Sim_Mat[j,i] = Sim_Mat[i,j]
                    
    

        clustering = AgglomerativeClustering(n_clusters=None,affinity="precomputed",distance_threshold=(distanceThreshold),linkage="average")
        labels = clustering.fit_predict(Sim_Mat)
        print(np.max(labels))
        
        if repeatedValues == False:

            prueba = np.vstack((columns,labels))
            prueba = prueba.T
    
            final = pd.DataFrame(data=prueba, columns=["ID","Cluster"])
            
            return Sim_Mat,final
        elif repeatedValues == True:
            
            newColumn = []
            for i in range(self.MaterialsData.shape[0]):
                material = self.MaterialsData.Material[i]
                segID = int(self.MaterialsData.GroupID[i])
                complete = material + "_" + str(segID)
                newColumn.append(complete)
            self.MaterialsData["ID"] = newColumn
            
            repeatedColumns = []
            repeatedLabels = []
            for i in range(columns.shape[0]):
                material = columns[i]
                label = labels[i]
                timesRepeated = self.MaterialsData[self.MaterialsData["ID"] == material].TimesRepeated.to_numpy()
                timesRepeated = int(timesRepeated[0])
                rep = [material]
                rep = rep * timesRepeated
                
                repLabel = [label]
                repLabel = repLabel * timesRepeated
                
                repeatedColumns = repeatedColumns + rep
                repeatedLabels = repeatedLabels + repLabel
            
            prueba = np.vstack((repeatedColumns,repeatedLabels))
            prueba = prueba.T
    
            final = pd.DataFrame(data=prueba, columns=["ID","Cluster"])
            final = final.astype({'Cluster':'int32'})
            
            return Sim_Mat,final
        
        
    def spectralClusteringCombined(self,distanceMatrix,distanceThreshold,repeatedValues = False,uniqueSegments1 = None, uniqueSegments2 = None):
        '''
        Compute the agglomerative clustering(spectral clustering) of the previously computed
        distance matrix(DistanceMatrixCombined function).

        Parameters
        ----------
        distanceMatrix : DataFrame
            Distance Matrix of the unique segments of the combined dataset
            (from the DistanceMatrixCombined function)
            
        distanceThreshold : float
            Similarity value above which the segments are NOT clustered
            
        repeatedValues : bool, optional
            Take into account how many times each segment is repeated
            in a unit cell. The default is False.
            
        uniqueSegments1 : DataFrame, optional
            First database DataFrame with the required information to 
            take into account the repeated values in this function. It is the
            result dataframe from the getUniqueFunction. The default is None.
            CAUTION: Required if "repeatedValues" is set to True
            
        uniqueSegments2 : DataFrame, optional
            Second database DataFrame with the required information to 
            take into account the repeated values in this function. It is the
            result dataframe from the getUniqueFunction. The default is None.
            CAUTION: Required if "repeatedValues" is set to True

        Returns
        -------
        list
            Similarity matrix and final clustering of the segments/Dataframe with
            the final clustering information

        '''
        
        
        
        print("zeoPy: Spectral Clustering Combined Function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        
        #Create the similarity matrix from the distance matrix
        #----------------------------------------------------------------------
        Sim_Mat = np.zeros_like(distanceMatrix)
        
        #Segment identifiers
        columns = distanceMatrix.columns.to_numpy()
        
        Distance_A = distanceMatrix.to_numpy()

        for i in range(Sim_Mat.shape[0]):
            for j in range(Sim_Mat.shape[1]):
                #Values equal to zero in the diagonal
                if i == j:
                    
                    Sim_Mat[i,j] = 0
                else:
                    
                    #Set a big error for negative R2 scores
                    if Distance_A[i,j] <= 0:
                        
                        Sim_Mat[i,j] = 100
                    #Set the inverse for the remaining values
                    else:
                        
                        Sim_Mat[i,j] = 1.0 / float(Distance_A[i,j])
        
        #Make the similarity matrix symmetric
        for i in range(Sim_Mat.shape[0]):
            for j in range(Sim_Mat.shape[1]):
                if i > j:
                    
                    if Sim_Mat[i,j] >= Sim_Mat[j,i]:
                        
                        Sim_Mat[i,j] = Sim_Mat[j,i]
                    elif Sim_Mat[i,j] < Sim_Mat[j,i]:
                        Sim_Mat[j,i] = Sim_Mat[i,j]
        
        #----------------------------------------------------------------------
                    
        #Compute the Agglomerative Clustering
        #----------------------------------------------------------------------
        
        clustering = AgglomerativeClustering(n_clusters=None,affinity="precomputed",distance_threshold=(distanceThreshold),linkage="average")
        #Labels of the clusters
        labels = clustering.fit_predict(Sim_Mat)
        print("Clustering. Number of Clusters: ", np.max(labels))
        
        if repeatedValues == False:
            
            #Link the segments identifiers to their cluster labels
            merge = np.vstack((columns,labels))
            merge = merge.T
    
            final = pd.DataFrame(data=merge, columns=["ID","Cluster"])
            
            return Sim_Mat,final
        
        elif repeatedValues == True:
            
            if (type(uniqueSegments1) == pd.core.frame.DataFrame) and  (type(uniqueSegments2) == pd.core.frame.DataFrame):
                
                #Link the segments identifiers to their cluster labels
                merge = np.vstack((columns,labels))
                merge = merge.T
        
                final = pd.DataFrame(data=merge, columns=["ID","Cluster"])
                #Variable to store the results
                timesRepeateds = np.zeros((final.shape[0]))
                
                #For each of the unique segments of the total clustering
                for i in range(final.shape[0]):
                    #Current unique segment
                    currentID = final.ID[i]
                    
                    #Find it in the first database
                    segmentInfo = uniqueSegments1[uniqueSegments1["ID"] == currentID]
                    #If not found in the first database, look for it on the second database
                    if segmentInfo.shape[0] == 0:
                        segmentInfo = uniqueSegments2[uniqueSegments2["ID"] == currentID]
                        
                    timesRepeated = segmentInfo.TimesRepeated.to_numpy()[0]
                    timesRepeateds[i] = int(timesRepeated)
                
                final["TimesRepeated"] = timesRepeateds
                
                return Sim_Mat,final
                
            else:
                
                 print("Repeated values selected, please include the data of the unique segments for both datasets")
                 return None   
        
        #----------------------------------------------------------------------
                    
        
        
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                
    
    def computeA(self,inputDataFrame,volumeFile):
        '''
        Given the results from the Spectral Clustering in the zeoPy.Spectral 
        Clustering function, compute the matrix A of the equations system
        Ax= b for the adsorption equations

        Parameters
        ----------
        inputDataFrame : DataFrame
            Results from the zeoPy.spectralClustering() function.
        volumeFile : String
            Path to the file storing adsorption data of the materials of the database

        Returns
        -------
        df : DataFrame
            Matrix A

        '''
        print("zeoPy: Compute A function")
        
        #Read the data file with the volume information of the materials in
        #database
        dfVolume = pd.read_csv(volumeFile)
        #Get the max volume value in order to compute the supercells
        maxVol = dfVolume.Volume.max()
        
        #Data with the information of the clusters of unique segments from the
        #zeoPy.spectralClustering function
        df = inputDataFrame.copy()
        
        #Get the clusters labels and sort them 
        clustersID = []
        for i in range(df.shape[0]):
            clustersID.append(df.Cluster[i])
        
        clustersID = list(dict.fromkeys(clustersID))
        clustersID.sort()
        
        #Get the unique segments identifiers
        materials = []
        for i in range(df.shape[0]):
            region  = df.ID[i]
            material = region.split("_")[0]
            materials.append(material)
            
        materials = list(dict.fromkeys(materials))

        #Store the results
        storage = np.zeros((len(materials),len(clustersID)))
        
        #For each of the unique segments' identifiers
        for i in range(df.shape[0]):
            region  = df.ID[i]
            material = region.split("_")[0]
            
            
            if type(dfVolume.Material.to_numpy()[0]) != str:
                material = int(material)
            
            volume = dfVolume[dfVolume["Material"] == material].to_numpy()
            
            
            vol = volume[0,1]
            #Compute the ratio between a supercell and the current unit cell
            ratio = maxVol/vol
            # ratio = 1
            
            if type(dfVolume.Volume.to_numpy()[0]) != str:
                material = str(material)
            mIndex = materials.index(material)
            
            nCluster = df.Cluster[i]
            nIndex = clustersID.index(nCluster)
            storage[mIndex,nIndex] = storage[mIndex,nIndex] + 1*ratio
            

        #Store the results in a dataframe   
        dfA = pd.DataFrame(data=storage,index=materials)    
        dfA.to_csv('A_matrix.csv')
        
        
        return dfA
    
    def computeACombined(self,inputDataFrame,volumeFile1,volumeFile2):
        '''
        Compute the matrix A required for solving the adsorption equation Ax = b
        where A is the matrix storing the unique segments clusters, x is the 
        adsorption of each cluster(to be found) and b is the total adsorption 
        of the materials

        Parameters
        ----------
        inputDataFrame : DataFrame
            DataFrame with the information from the spectralClusteringCombined 
            function. Includes information about the segments and the clusters 
            they are attached to
            
        volumeFile1 : string
            Path to the file with the volume information of the materials in the first dataset
            
        volumeFile2 : string
            Path to the file with the volume information of the materials in the second dataset

        Returns
        -------
        dfA : DataFrame
            Matrix A

        '''
        print("zeoPy: ComputeACombined function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        #Read the volume file from the first database
        dfVolume1 = pd.read_csv(volumeFile1)
        
        #Maximum volume in the first database
        maxVol = dfVolume1.Volume.max()
        
        #Read the volume file from the second database
        dfVolume2 = pd.read_csv(volumeFile2)
        
        
        #Create a copy of the dataframe to avoid non required changes
        df = inputDataFrame.copy()
        
        #Get the clusters number from the input file and remove repeated values and sort them
        #in order to have a list with all the clusters ID
        clustersID = []
        for i in range(df.shape[0]):
            clustersID.append(df.Cluster[i])
        
        clustersID = list(dict.fromkeys(clustersID))
        clustersID.sort()
        
        
        #Get the materials' ID from the input file and remove repeated values and sort them
        # in order to have a list of all the segments ID
        materials = []
        for i in range(df.shape[0]):
            region  = df.ID[i]
            material = region.split("_")[0]
            materials.append(material)
            
        materials = list(dict.fromkeys(materials))

        #Store the results
        storage = np.zeros((len(materials),len(clustersID)))
        
        #For each of the unique segments in the file
        for i in range(df.shape[0]):
            region  = df.ID[i]
            #Get the material's name
            material = region.split("_")[0]
            
            #Get the volume of the material in the first database
            volume = dfVolume1[dfVolume1["Material"] == material].to_numpy()
            
            #If not found, look it  for in the second database
            if volume.shape[0] == 0:
                volume = dfVolume2[dfVolume2["Material"] == int(material)].to_numpy()
            
            
            vol = volume[0,1]
            #Compute the ration between the unit cell of the material and the maximum unit cell
            ratio = maxVol/vol
            
            #Find the material by its index
            mIndex = materials.index(material)
            #Find the cluster by its index
            nCluster = df.Cluster[i]
            nIndex = clustersID.index(nCluster)
            storage[mIndex,nIndex] = storage[mIndex,nIndex] + 1 * ratio * df.TimesRepeated[i]
            
            
        #Create a dataframe to store the matrix A   
        dfA = pd.DataFrame(data=storage,index=materials)    
    
        return dfA
    
    def solveSystem(self,adsorptionFile,ADataFrame):
        '''
        Solve the adsorption equations system Ax = b to find the local adsorption
        contribution of each type of segment to the total adsorption of the materials

        Parameters
        ----------
        adsorptionFile : string
            Path to the adsorption files of the  database.
            
        ADataFrame : DataFrame
            Matrix A of the equations system.

        Returns
        -------
        x : Numpy array
            Adsorption contribution of each type of segment.
            
        rmse : float
            Root Mean Squared Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments.
            
        r2 : float
            R2 Score Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments.
            
        meanRelativeError : float
            Mean Relative Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments.
            
        spearman : float
            Spearman Coefficient between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments.
            
        b : Numpy Array
            Adsorption Data of the materials.
            
        predictedResult : DataFrame
            Predicted adsorption of the materials to compare to their true
            adsorption.
        
        materialsInfo : DataFrame
            DataFrame with the info required of the segments for future tasks
        A_adsorption_filtered : DataFrame
            A matrix filtered with the material that have adsorption data
            
        b_adsorption_filtered : DataFrame
            b vector filtered with the materials that have adsorption data



        '''
        
        print("zeoPy: Solve System  function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        
        if type(adsorptionFile) == str:
            dfAdsorption = pd.read_csv(adsorptionFile,index_col=(0))
            dfAdsorption.sort_index(inplace=True)
        elif type(adsorptionFile) == pd.core.frame.DataFrame:
            dfAdsorption = adsorptionFile
        
        dfA = ADataFrame
        materials = dfA.index.to_numpy()
        
        # return dfA, dfAdsorption
        storage = np.zeros((dfA.shape[0],dfA.shape[1]+1))

        for i in range(dfA.index.shape[0]):
            try:
                material = dfA.index[i]
                aRow = dfA.loc[[material]].to_numpy()
                if dfAdsorption.index.dtype == int:
                    adsorptionRow = dfAdsorption.loc[[int(material)]].to_numpy()
                else:
                    adsorptionRow = dfAdsorption.loc[[(material)]].to_numpy()
                total = np.append(aRow,adsorptionRow[0,7])
                storage[i,:] = total
            except:
                print("Not adsorption data of ", material)
                
                
        filteredStorage = storage[~np.all(storage == 0, axis=1)]
        filteredMaterials = materials[~np.all(storage == 0, axis=1)]
        
        
        
        A = filteredStorage[:,:-1]
        b = filteredStorage[:,-1]
        
        A_adsorption_filtered = pd.DataFrame(data = A, index=filteredMaterials)
        b_adsorption_filtered = pd.DataFrame(data = b, index = filteredMaterials)
        
        
        
        Aa = A
        #Number of materials connected to each material
        NumberOfConexions = []
        #Number of connected segments of each material
        NumberOfConnectedSegments = []
        for i in range(Aa.shape[0]):
            numberOfConnectedMaterials = 0
            numberOfConSegments = 0
            for j in range(Aa.shape[1]):
                currentValue = Aa[i,j]
                if currentValue > 0:
                    group = Aa[:,j].flatten()
                    nonNegativegroup = group[group[:] != 0]
                    #If there is more than one material sharing this cluster
                    if nonNegativegroup.shape[0] > 1:
                        conexions = nonNegativegroup.shape[0] - 1
                        numberOfConnectedMaterials += conexions
                        numberOfConSegments += 1
            NumberOfConexions.append(numberOfConnectedMaterials)
            NumberOfConnectedSegments.append(numberOfConSegments)
        data = np.vstack([NumberOfConexions,NumberOfConnectedSegments])
        materialsInfo = pd.DataFrame(data = data.T, columns=["ConnectedMaterials","ConnectedSegments"],index=filteredMaterials)
        
        
        
        x, ress = nnls(A,b,maxiter=1000000000000000)
        
        
        
        
        predictedResult = np.matmul(A,x)
        
        spearman = stats.spearmanr(b,predictedResult)[0]
        
        rel_error = (b-predictedResult)/b
        rel_error = abs(rel_error)
        
        # return filteredStorage, A, b, x, predictedResult
        
        rmse = mean_squared_error( b, predictedResult,squared = False)
        r2 = r2_score(b, predictedResult)
        
        return x, rmse, r2, np.mean(rel_error),spearman,b,predictedResult,materialsInfo, A_adsorption_filtered,b_adsorption_filtered
    
    def solveSystemCombined(self,ADataFrame,adsorptionFile1,adsorptionFile2,selectedColumn=7):
        '''
        Solve the adsorption equations system Ax = b to find the local adsorption
        contribution of each type of segment to the total adsorption of the materials

        Parameters
        ----------
        ADataFrame : DataFrame
            Matrix A of the equations system
            
        adsorptionFile1 : string
            Path to the adsorption files of the first database
            
        adsorptionFile2 : string
            Path to the adsorptoin files of the second database
            
        selectedColumn : int (optional)
        Selected column in the adsortion dataset. It will determine what 
        temperature and pressure values were used to compute the adsorption data
        Default: 8 th column --> Temperature:275.9 K Pressure: 403.4 bar

        Returns
        -------
        x : Numpy array
            Adsorption contribution of each type of segment
            
        rmse : float
            Root Mean Squared Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments
            
        r2 : float
            R2 Score Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments
            
        meanRelativeError : float
            Mean Relative Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments
            
        spearman : float
            Mean Relative Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments
            
        b : Numpy Array
            Adsorption Data of the materials
            
        predictedResult : DataFrame
            Predicted adsorption of the materials to compare to their true
            adsorption
        A_with_adsorption : DataFrame
            A merged with b for future manipulations

        '''
        
        print("zeoPy: Solve System Combined function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        #Set the different options for the input adsorption files
        #----------------------------------------------------------------------
        if type(adsorptionFile1) == str:
            dfAdsorption1 = pd.read_csv(adsorptionFile1,index_col=(0))
            dfAdsorption1.sort_index(inplace=True)
        elif type(adsorptionFile1) == pd.core.frame.DataFrame:
            dfAdsorption1 = adsorptionFile1
            
        if type(adsorptionFile2) == str:
            dfAdsorption2 = pd.read_csv(adsorptionFile2,index_col=(0))
            dfAdsorption2.sort_index(inplace=True)
        elif type(adsorptionFile2) == pd.core.frame.DataFrame:
            dfAdsorption2 = adsorptionFile2
        
        #----------------------------------------------------------------------
        
        dfA = ADataFrame
        materials = dfA.index.to_numpy()
        
        storage = np.zeros((dfA.shape[0],dfA.shape[1]+1))

        for i in range(dfA.index.shape[0]):
            #Get the material name from the index
            material = dfA.index[i]
            aRow = dfA.loc[[material]].to_numpy()
            
            
            
            try:
                #Look for the material in the first database
                adsorptionRow = dfAdsorption1.loc[[(material)]].to_numpy()
                total = np.append(aRow,adsorptionRow[0,selectedColumn])
                storage[i,:] = total
            except:
                try:
                    #If not found, look for the material in the PCOD database
                    adsorptionRow = dfAdsorption2.loc[[(material)]].to_numpy()
                    storage[i,:] = total
                except:
                    try:
                        #If not found, look for the material in the PCOD database changing the index to int
                        adsorptionRow = dfAdsorption2.loc[[int(material)]].to_numpy()
                        storage[i,:] = total
                    except:
                        #If not found, then it does not exist
                        print("Not adsorption data of ", material)
            
        
                
        #Remove te materials without adsorption data       
        filteredStorage = storage[~np.all(storage == 0, axis=1)]
        filteredMaterials = materials[~np.all(storage == 0, axis=1)]
        
        
        # A term for the Ax = b
        A = filteredStorage[:,:-1]
        # b term fro the Ax = b
        b = filteredStorage[:,-1]
        
        #Scipy's Non Negative Leasts Squares Function to solve the system and
        #find x
        x, ress = nnls(A,b,maxiter=1000000000000000)
        
        #Matrix A filtered only with the materials with adsorption data
        A_with_adsorption = pd.DataFrame(data=A,index=filteredMaterials)
        
        #Materials' predicted total adsorption
        predictedResult = np.matmul(A,x)
        
        #Spearman coefficient
        spearman = stats.spearmanr(b,predictedResult)[0]
        
        #Mean relative error
        rel_error = (b-predictedResult)/b
        rel_error = abs(rel_error)
        
        #Root Mean Squared Error
        rmse = mean_squared_error( b, predictedResult,squared = False)
        #R2 score error
        r2 = r2_score(b, predictedResult)
        
        return x, rmse, r2, np.mean(rel_error),spearman,b,predictedResult, A_with_adsorption
        
        
        
    def checkTestSet(self,testSet,adsorptionFile1,adsorptionFile2,trainSolutions):
        '''
        Check if we have enough information to predict adsorption data in a test
        set

        Parameters
        ----------
        testSet : DataFrame
            Test Set to check
            
        adsorptionFile1 : string
            File path to the adsorption file of the first database.
            
        adsorptionFile2 : string
            File path to the adsorption file of the second database.
            
        trainSolutions : Numpy Array
            Predicted adsorptions of the clusters in the train dataset.

        Returns
        -------
        rmse : float
            Root Mean Squared Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments
            
        r2 : float
            R2 Score Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments
            
        meanRelativeError : float
            Mean Relative Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments
            
        spearman : float
            Mean Relative Error between the total adsorption of the materials
            in their adsorption files and summing up the local contributions 
            of their segments
            
        b : Numpy Array
            Adsorption Data of the materials
            
        predictedResult : DataFrame
            Predicted adsorption of the materials to compare to their true
            adsorption
        

        '''
        
        print("zeoPy: Check Test Function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        #Set the different options for the input adsorption files
        #----------------------------------------------------------------------
        
        if type(adsorptionFile1) == str:
            dfAdsorption1 = pd.read_csv(adsorptionFile1,index_col=(0))
            dfAdsorption1.sort_index(inplace=True)
        elif type(adsorptionFile1) == pd.core.frame.DataFrame:
            dfAdsorption1 = adsorptionFile1
            
        if type(adsorptionFile2) == str:
            dfAdsorption2 = pd.read_csv(adsorptionFile2,index_col=(0))
            dfAdsorption2.sort_index(inplace=True)
        elif type(adsorptionFile2) == pd.core.frame.DataFrame:
            dfAdsorption2 = adsorptionFile2
            
        #----------------------------------------------------------------------
        
        dfA = testSet
        materials = dfA.index.to_numpy()
        
        
        
        storage = np.zeros((dfA.shape[0],dfA.shape[1]+1))

        for i in range(dfA.index.shape[0]):
            #Get the material name from the index
            material = dfA.index[i]
            aRow = dfA.loc[[material]].to_numpy()
            
            try:
                #Look for the material in the IZA database
                adsorptionRow = dfAdsorption1.loc[[(material)]].to_numpy()
                total = np.append(aRow,adsorptionRow[0,7])
                storage[i,:] = total
            except:
                try:
                    #Look for the material in the PCOD database
                    adsorptionRow = dfAdsorption2.loc[[(material)]].to_numpy()
                    total = np.append(aRow,adsorptionRow[0,7])
                    storage[i,:] = total
                except:
                    try:
                        #Look for the material in the PCOD database changing the index to int
                        adsorptionRow = dfAdsorption2.loc[[int(material)]].to_numpy()
                        total = np.append(aRow,adsorptionRow[0,7])
                        storage[i,:] = total
                    except:
                        #If not found, then it does not exist
                        print("Not adsorption data of ", material)
            
        
                
        #Remove te materials without adsorption data       
        filteredStorage = storage[~np.all(storage == 0, axis=1)]
        filteredMaterials = materials[~np.all(storage == 0, axis=1)]
        
        
        # A term for the Ax = b
        A = filteredStorage[:,:-1]
        # b term fro the Ax = b
        b = filteredStorage[:,-1]
        
        predictedResult = np.matmul(A,trainSolutions)
        
        spearman = stats.spearmanr(b,predictedResult)[0]
        
        rel_error = (b-predictedResult)/b
        rel_error = abs(rel_error)
        
        # return filteredStorage, A, b, x, predictedResult
        
        rmse = mean_squared_error( b, predictedResult,squared = False)
        r2 = r2_score(b, predictedResult)
        
        return rmse, r2, np.mean(rel_error),spearman,b,predictedResult
        
        
    def getClusterHisto(self,finalClustering,uniqueHistos):
        uniqueHistos = uniqueHistos.copy()
        groups = []
        for i in range(uniqueHistos.shape[0]):
            try:
                mat = uniqueHistos.iloc[i]
                identifier = str(mat.Material) + "_" + str(int(mat.GroupID))
                group = float(finalClustering[finalClustering['ID'] == identifier].to_numpy()[0,1])
                groups.append(int(group))
            except:
                print(identifier)
                groups.append(-1)
           
        uniqueHistos["Cluster"] = groups
        uniqueHistos = uniqueHistos[uniqueHistos["Cluster"] != (-1)]
        uniqueHistos = uniqueHistos.iloc[:,4:]
        return uniqueHistos
    
    
        
    
    
    def comparePreviousTotal(self,histos1,previousHistos,clusterAdsorption,materialsVolume1,materialsVolume2,adsorptionFile):
        
        
        materialsVolume1 = pd.read_csv(materialsVolume1)
        materialsVolume2 = pd.read_csv(materialsVolume2)
        
        superCellVol = materialsVolume2.Volume.max()
        
        
        
        clusters = previousHistos.Cluster.to_numpy()
        
        prevHist = previousHistos.to_numpy()
        prevHist = prevHist[:,:-2]
        
        currentMaterials = histos1.Material.to_numpy()
        timesRepeated = histos1.TimesRepeated.copy().to_numpy()
 
        
        currentHistos = histos1.to_numpy()
        currentHistos = currentHistos[:,4:]
        
        results = np.zeros((currentHistos.shape[0],2))
        for i in range(currentHistos.shape[0]):
            histo1 = currentHistos[i]
            material = currentMaterials[i]
            
            regionsDifference = np.zeros((prevHist.shape[0],2))
            for j in range(prevHist.shape[0]):
                histo2 = prevHist[j]
                r2 = r2_score(histo1,histo2)
                regionsDifference[j,0] = clusters[j]
                regionsDifference[j,1] = r2
            closestMaterialIndex = np.argmax(regionsDifference[:,1])
            closestCluster = regionsDifference[closestMaterialIndex,0]
            
            if regionsDifference[closestMaterialIndex,1] < 0.98:
                print("Not that close:", material, " Error: ", regionsDifference[closestMaterialIndex,1])
            else:
                print("GOOD")
            results[i,0] = material
            results[i,1] = closestCluster
            
        res = pd.DataFrame(data=results, columns=["Material", "Cluster"])
        res["TimesRepeated"] = timesRepeated
        res["Volume"] = histos1.Volume.to_numpy()
        res = res.astype({"Material" : int})
        
        for i in range(res.shape[0]):
            mat = res.Material[i]
            vol = materialsVolume1[materialsVolume1["Material"] == mat]
            vol = vol.Volume
            
            ratio = superCellVol / vol
            res.TimesRepeated[i] = res.TimesRepeated[i] * ratio
            
        regionsAdsorption = []
        for i in range(res.shape[0]):
            cluster = int(res.Cluster[i])
            
            clusAds = clusterAdsorption[cluster]
            clusAds = clusAds * res.TimesRepeated[i]
            
            regionsAdsorption.append(clusAds)
        res["Adsorption"] = regionsAdsorption
            
        mats = res.Material.unique()
        
        finalRes = np.zeros((mats.shape[0],2))
        for i in range(mats.shape[0]):
            current = res[res["Material"] == mats[i]]
            adsorptions = current.Adsorption.to_numpy()
            totalAds = np.sum(adsorptions)
            finalRes[i,0] = int(mats[i])
            finalRes[i,1] = totalAds
        finalRes = pd.DataFrame(data=finalRes,columns=["Material","TotalAdsorption"])
            
        trueAdsorption = pd.read_csv(adsorptionFile, index_col=0)
        trueAds = trueAdsorption.iloc[:,7]
        
        compare = np.zeros((finalRes.shape[0],3))
        
        for i in range(finalRes.shape[0]):
            current = int(finalRes.Material[i])
            trueAd = trueAds[trueAds.index == current].to_numpy()[0]
            
            compare[i,0] = current
            compare[i,1] = trueAd
            compare[i,2] = finalRes.TotalAdsorption[i]
            
        compareData = pd.DataFrame(data=compare,columns=["Material","TrueAds","PredAds"])
        
        r2 = r2_score(compare[:,1],compare[:,2])
        
        return prevHist,currentHistos,res, finalRes,compareData,r2
        
            

    def removeMaterialsWithoutAdsorption(self,inputDirectory,inputAdsorptionData, outputDirectory):
        '''
        Read the material's files in a directory and copy the ones with adsorption data
        to a new directory

        Parameters
        ----------
        inputDirectory : string
            Input folder's directory
        inputAdsorptionData : string
            Input Adsorption file
        outputDirectory : string
            Ouput folder's directory'

        Returns
        -------
        None.

        '''
        print("zeoPy: Remove Materials Without Adsorption Function")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        #Read input directory
        inputFiles = os.listdir(inputDirectory)
        numberOfFiles = len(inputFiles)
        print("Number of files: " ,numberOfFiles)
        
        #Read adsorption file
        adsorption = pd.read_csv(inputAdsorptionData,index_col=(0))
        
        #Create output directory
        if os.path.exists(outputDirectory):
            #To avoid a folder with old/unnecessary files
            shutil.rmtree(outputDirectory,ignore_errors=True)
            print("Removed")
            os.mkdir(outputDirectory)
        else:
            os.mkdir(outputDirectory)
        
        for file in inputFiles:
            material = file.split(".")[0]
            
            #Look for the material in the adsorption file
            try:
                adsorption.loc[material]
                src_path = inputDirectory + file
                dst_path = outputDirectory + file
                shutil.copy(src_path,dst_path)
            #If not found, maybe it is because of the index type
            except:
                try:
                    adsorption.loc[int(material)]
                    
                    #Copy the material in the new folder
                    src_path = inputDirectory + file
                    dst_path = outputDirectory + file
                    shutil.copy(src_path,dst_path)
                except:
                    print("Not adsorption data of: ", material)
                    
    def computeSuperCellsAccessibleVolume(self,inputData,volumeFile,superCell=False):
            
            #Read the volume file in order to create superCells
            materialsVolume = pd.read_csv(volumeFile,index_col=(0))
            #Max volume in the dataset
            maxMaterialVolume = materialsVolume.Volume.max()
            print("Max material volume in the dataset: ", maxMaterialVolume)
            
           
            
            
            #Get the materials names
            materials = np.unique(inputData.Material.to_numpy())
            numberOfMaterials = materials.shape[0]
            
            #Save the materials Accessible Volume
            materialsTotalVolume = np.zeros(numberOfMaterials)
            
            #For each of the materials
            for i in range(numberOfMaterials):
                mat = materials[i]
                #Find the material in the dataframe
                materialFound = inputData[inputData["Material"] == mat]
                #Get the volumes of the segment types in the material
                volumes = materialFound.Volume.to_numpy()
                #Get the number of times each segment type is repeated
                reps = materialFound.TimesRepeated.to_numpy()
                volume = 0
                for j in range(materialFound.shape[0]):
                    vol = volumes[j]
                    rep = reps[j]
                    volume += vol*rep
                if superCell == False:
                    materialsTotalVolume[i] = volume
                if superCell == True:
                    materialTotalVolume = materialsVolume.loc[mat]
                    ratio = maxMaterialVolume / materialTotalVolume
                    materialsTotalVolume[i] = volume * ratio
                
            materialsVolumeData = pd.DataFrame(materialsTotalVolume, index=materials,columns=["Volume"])
            return materialsVolumeData
    def computeAVolume(self,inputData,clustering,volumeFile,superCell=False):
        
        
        #Get the number of possible clusters
        numberOfClusters = clustering.Cluster.max() + 1
        
        #Read the volume file in order to create superCells
        materialsVolume = pd.read_csv(volumeFile,index_col=(0))
        #Volume 
        maxMaterialVolume = materialsVolume.Volume.max()
        
        #Get the materials names
        materials = np.unique(inputData.Material.to_numpy())
        numberOfMaterials = materials.shape[0]
        
        A_mat = np.zeros((materials.shape[0],numberOfClusters))
        nonRepeatedClustering = clustering.drop_duplicates(keep="first",inplace=False)
        
        
        for i in range(A_mat.shape[0]):
            mat = materials[i]
            materialFound = inputData[inputData["Material"] == mat]
            materialFound.reset_index(inplace=True,drop=True)
            volumes = materialFound.Volume.to_numpy()
            reps = materialFound.TimesRepeated.to_numpy()
            
            for j in range(materialFound.shape[0]):
                segmentID = materialFound.ID[j]
                #Find the segment in the clustering
                cluster = nonRepeatedClustering[nonRepeatedClustering["ID"] == segmentID].to_numpy()[0,1]
                vol = volumes[j]
                rep = reps[j]
                volume = vol*rep
                if superCell == True:
                    materialTotalVolume = materialsVolume.loc[mat]
                    ratio = maxMaterialVolume / materialTotalVolume
                    volume = volume * ratio
                A_mat[i,cluster] = A_mat[i,cluster] + volume
                
        A_mat = pd.DataFrame(A_mat,index=materials)
        return A_mat
    
    def computeAVolumeDensity(self,inputData,clustering,volumeFile,adsorptionFile,accessibleVolume,superCell=False):
        
        #Read the adsorption file
        adsorption = pd.read_csv(adsorptionFile,index_col=(0))
        
        #Get the number of possible clusters
        numberOfClusters = clustering.Cluster.max() + 1
        
        #Read the volume file in order to create superCells
        materialsVolume = pd.read_csv(volumeFile,index_col=(0))
        #Volume 
        maxMaterialVolume = materialsVolume.Volume.max()
        
        #Get the materials names
        materials = np.unique(inputData.Material.to_numpy())
        numberOfMaterials = materials.shape[0]
        
        A_mat = np.zeros((materials.shape[0],numberOfClusters))
        nonRepeatedClustering = clustering.drop_duplicates(keep="first",inplace=False)
        
        #For each of the materials
        for i in range(A_mat.shape[0]):
            mat = materials[i]
            #Find the material in the info dataframe
            materialFound = inputData[inputData["Material"] == mat]
            materialFound.reset_index(inplace=True,drop=True)
           
            #Number of times repeated of each segment type
            reps = materialFound.TimesRepeated.to_numpy()
            
            #Material total Accessible Volume
            materialAccessibleVolume = accessibleVolume.loc[mat].Volume
            
            #Materal adsorption density
            materialTotalAdsorption = adsorption.loc[materials[i]].to_numpy()[7]
            materialGeneralDensity = materialTotalAdsorption/materialAccessibleVolume
            
            #For each of the segment types of the material
            for j in range(materialFound.shape[0]):
                segmentID = materialFound.ID[j]
                #Find the segment in the clustering
                cluster = nonRepeatedClustering[nonRepeatedClustering["ID"] == segmentID].to_numpy()[0,1]
               
                rep = reps[j]
                value = rep/materialGeneralDensity
                if superCell == True:
                    materialTotalVolume = materialsVolume.loc[mat]
                    ratio = maxMaterialVolume / materialTotalVolume
                    value = value * ratio
                A_mat[i,cluster] = A_mat[i,cluster] + value
                
        A_mat = pd.DataFrame(A_mat,index=materials)
        return A_mat
    
    def combineAdsorptionVolume(self, Aads,Avol,matVolume,adsorptionFile):
    
        materials = Aads.index.to_numpy()
        
        adsorption = pd.read_csv(adsorptionFile,index_col=(0))
        
        bAds = np.zeros_like(materials)
        for i in range(materials.shape[0]):
            materialAds = adsorption.loc[materials[i]].to_numpy()[7]
            bAds[i] = materialAds
        
        bAdsorption = pd.DataFrame(data=bAds,columns=["Adsorption[g/L]"],index=materials)
        materials = Avol.index.to_numpy()
        
        bVol = np.zeros_like(materials)
        for i in range(materials.shape[0]):
            materialAccessibleVolume = matVolume.loc[materials[i]].Volume
            # materialAds = adsorption.loc[materials[i]].to_numpy()[7]
            # materialVol = matVolume.loc[materials[i]].to_numpy()[0]
            bVol[i] = materialAccessibleVolume
        
        bVolume = pd.DataFrame(data=bVol,columns=["AccessibleVolume[Ang3]"],index=materials)
            
        Atotal = np.concatenate([Aads,Avol])
        btotal = np.concatenate([bAds,bVol])
        
        #Scipy's Non Negative Leasts Squares Function to solve the system and
        #find x
        x, ress = nnls(Atotal,btotal,maxiter=1000000000000000)
        

        #Materials' predicted total adsorption
        predictedResult = np.matmul(Atotal,x)
        
        #Spearman coefficient
        spearman = stats.spearmanr(btotal,predictedResult)[0]
        
        #Mean relative error
        rel_error = (btotal-predictedResult)/btotal
        rel_error = np.linalg.norm(rel_error)
        
        #Root Mean Squared Error
        rmse = mean_squared_error( btotal, predictedResult,squared = False)
        #R2 score error
        r2 = r2_score(btotal, predictedResult)
        
        return bAdsorption, bVolume,x, r2,spearman,rel_error,rmse, predictedResult
        
        
            
            
            
        
        
        
        
            