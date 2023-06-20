
<p align="right">
  <img width="300" height="100" src="https://user-images.githubusercontent.com/92568531/169989478-a156733e-fbaa-4f7a-bdcb-e974ddd02d9b.png">
</p>

<!---
<p align="center">
  <img width="400" height="200" src="https://user-images.githubusercontent.com/92568531/169777821-9411a6e3-45cb-4c84-a9c7-3100413bbf46.png">
</p>
-->

<!---
<p align="center">
  <img width="400" height="200" src="https://user-images.githubusercontent.com/92568531/169775759-1c646cf5-2265-4689-9e92-c76436bceeaf.png">
</p>
-->
<!---
<p align="center">
  <img width="500" height="250" src="https://user-images.githubusercontent.com/92568531/169784943-9c671061-5c0c-4cac-8b85-9c0c23a07002.png">
</p>
-->
<!---
<p align="center">
  <img width="500" height="250" src="https://user-images.githubusercontent.com/92568531/169898882-5699a4f0-18b5-4ae8-bdd1-517daae5058e.png">
</p>
-->

<p align="center">
  <img width="500" height="250" src="https://user-images.githubusercontent.com/92568531/169987520-814d98b1-c198-4606-941c-4ec38a83c6f3.png">
</p>



# **ZeoPy User's guide**
Post-proccessing Python tool used to analyze the results from the C++ zeoTDA tool. Here, we will use the information from the segments of the void space of the nanoporous materials obtained in the zeoTDA tool.

## 0. Initial steps
As mentioned above, in this tool we will analyze the information obtained in the zeoTDA tool. So first, we would need to compute the MSC analysis of the materials(**zeoTDA::MSC()**)and then get the information of the accessible segments of the void space for an invited molecule(**zeoTDA::accessibleVoidSpace()**).
Having this information from the materials of a dataset stored in a folder, we can proceed to use the zeoPy tool.

## 1. Installation of the tool
To use this tool, we only need to have Python installed on our computer and a IDE/text editor capable of reading Python files.

---

**TIP:** Although we can use whatever IDE/text editor we prefer, we recommend using **Anaconda's Spyder IDE** which we find extremly useful for data analysis tasks.

---

## 2. Tool's description

### 2.1 Tool's structure
The tool is organized in a Python class. This class contains all the required functions needed to compute the post-processing analysis of the Topological Data Analysis(TDA) results. Although it is neccessary to compute the required functions in a certain order, the structure is quite simple.

## 3. zeoPy Usage Guide

### 3.1 Functions Glossary
- [getUnique()](#getUnique)
- [getUniqueSegments()](#getUniqueSegments)
- [computeDistanceMatrix()](#computeDistanceMatrix)
- [spectralClustering()](#spectralClustering)
- [computeA()](#computeA)
- [solveSystem()](#solveSystem)
- [computeDistanceMatrixCombined()](#computeDistanceMatrixCombined)
- [spectralClusteringCombined()](#spectralClusteringCombined)
- [computeACombined()](#computeACombined)
- [solveSystemCombined()](#solveSystemCombined)

### 3.2 Functions Description

#### `zeoPy.getUnique(self,inputFile,removeIsolated,minValue=None,maxValue=None)`  <a name="getUnique"></a>
One of the main functions in the class. It reads the segment's information file of a material obtained in the zeoTDA tool. Once it is readed, it computes the histograms of the fields included in the file for all the segments of the material. Finally, it compares the histograms to find which segments are similar in the material. This way, we get the unique segments of the material.

##### Parameters:
- `inputFile`: (string) Path to the input .csv file. File that stores the information of the segments obtained in the C++ zeoTDA.
- `removeIsolated`:(bool) Remove isolated segments of the material.
- `minValue`: (float)(default:None) Field min value used to create the bins of the histograms. If None, the tool will find this value itself. Useful if we need to set the same limits for all the materials in order to compare them. 
- `minValue`: (float)(default:None) Field max value used to create the bins of the histograms. If None, the tool will find this value itself. Useful if we need to set the same limits for all the materials in order to compare them.

##### Return:
- `material`:(DataFrame) Input material's DataFrame.
- `bins`: (list) Bins edges used to compute the histograms.
- `distance_mat`:(DataFrame) Distance Matrix of the segments of the material.
- `groups`:(list) Segment members of each unique segment group.

#### `zeoPy.getUniqueSegments(self, inputFolder,removeIsolated=False,bounds=None)` <a name="getUniqueSegments"></a>
Use the getUnique function for all the files in a folder. First, it finds the field bounds(if not given) of the whole materials in the folder. Then, it gets the unique segments of each of the materials. It returns a DataFrame with the information of this unique segments for all the materials.

##### Parameters:
- `inputFolder`:(string) Path to the folder where the files with the segments info is stored.
- `removeIsolated`:(bool) Remove isolated regions of the materials.
- `bounds`: (list)(default:None) Histogram bounds in case we want to set a predefined values.

##### Return:
- `df`:(DataFrame) Data of all the unique segments of the materials
- `bounds`:(list) Bounds used for the histograms computation in case they are needed.

#### `zeoPy.computeDistanceMatrix(self,df)` <a name="computeDistanceMatrix"></a>
Compute the distance matrix between all the unique segments of the materials. It uses their histograms to compare them.

##### Parameters:
-`df`:(DataFrame) Data of the materials' unique segments obtained in the zeoPy.getUniqueSegments() function.

##### Return:
- `distance_mat`: Distance Matrix comparing all the unique segments of the materials.

#### zeoPy.spectralClustering(self,distanceMatrix,distanceThreshold,repeatedValues=True) <a name="spectralClustering"></a>
Compute the Agglomerative Clustering(spectral Clustering) of the previously computed distance matrix(zeoPy.computeDistanceMatrix()).

##### Parameters:
- `distanceMatrix`:(DataFrame) Distance matrix of the unique segments of the dataset computed in the zeoPy.computeDistanceMatrix() function.
- `distanceThreshold`:(float) Similarity value above which the unique segments are NOT clustered together.
- `repeatedValues`:(bool) Take into account how many times each unique segment is repeated in a material.

##### Return:
- `similarityMatrix`:(Numpy Array) Similarity Matrix(based on the Distance Matrix) of the unique segments of the dataset.
- `clustering`:(DataFrame) Clustering information of the unique segments of the dataset.

#### `zeoPy.computeA(self, inputDataFrame,volumeFile)` <a name="computeA"></a>
Given the results from the zeoPy.spectralClustering() function, compute the matrix A of the adsorption equations system Ax= b.

##### Parameters:
- `inputDataFrame`:(DataFrame) Results from the zeoPy.spectralClustering() function.
- `volumeFile`:(String) Path to the file storing adsorption data of the materials of the database.

##### Return:
- `matrixA`:(DataFrame) Matrix A of the adsorption equations system Ax = b.

#### `zeoPy.solveSystem(self, adsorptionFile, ADataFrame)` <a name="solveSystem"></a>
Solve the adsorption equations system Ax=b to find the local adsorption contribution of each type of segment to the total adsorption of the materials.

##### Parameters:
- `adsorptionFile`:(String) Path to the adsorption file of the database
- `ADataFrame`:(DataFrame) Matrix A of the equations system computed in the zeoPy.computeA() function.

###### Return:
- `x`:(Numpy Array) x vector of the adsorption equations. Adsorption contribution of each cluster of segments.
- `rmse`:(Float) Root Mean Squared Error between the total adsorption of the materials(from their adsorption files) and the predicted adsorption obtained summing up the local adsorption of the unique segments of that material.
- `R2`:(Float) R2 Score Error between the total adsorption of the materials(from their adsorption files) and the predicted adsorption obtained summing up the local adsorption of the unique segments of that material.
- `mean Relative Error`:(Float) Mean Relative Error between the total adsorption of the materials(from their adsorption files) and the predicted adsorption obtained summing up the local adsorption of the unique segments of that material.
- `spearman`:(Float) Spearman Coefficient between the total adsorption of the materials(from their adsorption files) and the predicted adsorption obtained summing up the local adsorption of the unique segments of that material.
- `b`:(Numpy Array) b vector from the adsorption equations system. Simulated Adsorption Data of the materials.
- `predictedResult`:(DataFrame) Predicted total adsorption of the materials to compare to their true(simulated) adsorption.
- `materialsInfo`:(DataFrame) DataFrame with info of the unique segments' clusters to use in future tasks.
- `A_adsorption_filtered`:(DataFrame) A matrix filtered with the materials that have adsorption data.
- `b_adsorptoin_filered`:(DataFrame) b vector filtered with the materials that have adsorption data.


#### `zeoPy.computeDistanceMatrixCombined(self,df1,df2)` <a name="computeDistanceMatrixCombined"></a>
Compute the distance matrix between all the unique segments from two databases. It uses the results from the zeoPy.getUniqueSegments() for each of the databases. The distance matrix is computed based on the histograms of the unique segments.

##### Parameters:
- `df1`:(DataFrame) Results datarframe from the zeoPy.getUniqueSegments() function for the first database.
- `df2`:(DataFrame) Results datarframe from the zeoPy.getUniqueSegments() function for the second database.

##### Return:
-`distance_mat`:(DataFrame) Distance matrix of all the unique segments of both datasets.

#### `zeoPy.spectralClusteringCombined(self,distanceMatrix,distanceThreshold,repeatedValues=False,uniqueSegments1=False,uniqueSegments2=False)` <a name="spectralClusteringCombined"></a>
Compute the Agglomerative Clustering(spectral Clustering) of the previously compute distance matrix(from the zeoPy.DistanceMatrixCombined function).

##### Parameters:
- `distanceMatrix`:(DataFrame) Distance Matrix of the unique segments of the combined dataset
- `distanceThreshold:`(Float) Similarity value above which the unique segments are NOT clustered together
- `repeatedValues=False`: Take into account how many times each unique segments is repeated in a material
- `uniqueSemgnents1=False`: First database DataFrame with the required information to take into account the repeated values in this function. It is the result dataframe from the zeoPy.getUniqueSegments function(). The default is None. **CAUTION:** Required if "repeatedValues" is set to True
- `uniqueSemgnents2=False`: Second database DataFrame with the required information to take into account the repeated values in this function. It is the result dataframe from the zeoPy.getUniqueSegments function(). The default is None. **CAUTION:** Required if "repeatedValues" is set to True

##### Return:
- `similarityMatrix`:(Numpy Array) Similarity Matrix(based on the Distance Matrix) of the unique segments of the dataset.
- `clustering`:(DataFrame) Clustering information of the unique segments of the dataset.

#### `zeoPy.computeACombined(self, inputDataFrame,volumeFile1,volumeFile2)` <a name="computeACombined"></a>
  Compute the matrix A required for solving the adsorption equation Ax = b where A is the matrix storing the unique segments clusters, x is the adsorption of each cluster(to be found) and b is the total adsorption of the materials.

##### Parameters:
- `inputDataFrame`:(DataFrame) Data with the information from the zeoPy.spectralClusteringCombined() function. Includes information about the segments and the clusters they are clustered into.
- `volumeFile1`:(String) Path to the file with the volume information of the materials in the first dataset.
- `volumeFile2`:(String) Path to the file with the volume information of the materials in the second dataset.

##### Return:
- `matrixA`:(DataFrame) Matrix A of the adsorption equations system Ax = b.

#### `zeoPy.solveSystemCombined(self, ADataFrame, adsorptionFile1, adsorptionFile2)` <a name="solveSystemCombined"></a>
Solve the adsorption equations system Ax=b to find the local adsorption contribution of each type of segment to the total adsorption of the materials.

##### Parameters:
- `ADataFrame`:(DataFrame) Matrix A of the equations system computed in the zeoPy.computeACombined() function.
- `adsorptionFile1`:(String) Path to the adsorption file of the first database
- `adsorptionFile2`:(String) Path to the adsorption file of the second database

###### Return:
- `x`:(Numpy Array) x vector of the adsorption equations. Adsorption contribution of each cluster of segments.
- `rmse`:(Float) Root Mean Squared Error between the total adsorption of the materials(from their adsorption files) and the predicted adsorption obtained summing up the local adsorption of the unique segments of that material.
- `R2`:(Float) R2 Score Error between the total adsorption of the materials(from their adsorption files) and the predicted adsorption obtained summing up the local adsorption of the unique segments of that material.
- `mean Relative Error`:(Float) Mean Relative Error between the total adsorption of the materials(from their adsorption files) and the predicted adsorption obtained summing up the local adsorption of the unique segments of that material.
- `spearman`:(Float) Spearman Coefficient between the total adsorption of the materials(from their adsorption files) and the predicted adsorption obtained summing up the local adsorption of the unique segments of that material.
- `b`:(Numpy Array) b vector from the adsorption equations system. Simulated Adsorption Data of the materials.
- `predictedResult`:(DataFrame) Predicted total adsorption of the materials to compare to their true(simulated) adsorption.
- `materialsInfo`:(DataFrame) DataFrame with info of the unique segments' clusters to use in future tasks.
- `A_adsorption_filtered`:(DataFrame) A matrix filtered with the materials that have adsorption data.


## 4. Usage Examples

### Example 1:



## Authors:

- Maciek Haranczyk(IMDEA Materials Institute): maciej.haranczyk@imdea.org
- Jorge Zorrilla Prieto(IMDEA Materials Institute): jorge.zorrilla.prieto@gmail.com




