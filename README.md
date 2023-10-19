<p align="center">
  <img src="https://github.com/AMDatIMDEA/tda-segmentor/blob/main/images/logo.png?raw=true" alt="Sublime's custom image"/>
</p>

# Overview

Welcome to tda-segmentor! This is a software package written in C++ that segments a nanoporous material into different regions, which can then be used to identify unique features and properties. This tool uses ideas from topological data analysis (TDA), to generate what are called as morse smale complexes, which decompose a given space into monotonically disjoint segments. This tool is built on the [Topology Toolkit (TTK)](https://topology-tool-kit.github.io/), an open source TDA library that does most of the TDA analysis. The package tda-segmentor, outputs in a high-throughput way the segment information of the void space, the void space accessible to a guest molecule, graph representations of the accessible void space, graph representation of the solid region etc. 

# Installation Instructions

Compiling the tda-segmentor code should be done in two steps :

* Install the ttk library from source.
* Compile the tda-segmentor code linking the ttk library.

## Installing the Topology Toolkit (TTK)

* TTK must be installed from the source code and detailed instructions
  can be found [here](https://topology-tool-kit.github.io/installation.html).
* To install TTK from source, first Paraview must be compiled from sources. The code tda-segmentor has been tested with  
  Paraview-v5.11.1 and ttk-1.2.0 on a linux machine and Paraview-v5.10.1 and ttk-1.2.0 on a Mac machine. 

**Important notes for multicore advantages**

While compiling Paraview and TTK, the following modules and flags
need to be set active so that all the algorithms are compiled 
with multicore capabilities: 

While compiling Paraview:

* PARAVIEW_USE_MPI ------------------------------- **ON**
* VTK_GROUP_ENABLE_MPI ---------------------- **YES**
* VTK_SMP_IMPLEMENTATION_TYPE ---------- **OPENMP**
* VTKm_ENABLE_OPENMP ------------------------- **YES**
* VTKm_ENABLE_MPI --------------------------------- **YES**
* VTKm_ENABLE_TBB -------------------------------- **YES**

While compiling TTK:

* TTK_ENABLE_MPI ----------------------- **ON**

## Compiling the tda-segmentor code:

**cmake** should be used for generating the makefile. Move to the respository directory tda-segmentor/
and type **ccmake .** and one might have to hit configure *'c'* a few times and finally *'g'* to 
generate the makefile. While generating the makefile, the correct paths for vtk, paraview and ttk need
to be provided. If the directories are created as recommended by the TTK website, then these paths are: 

**For Mac Machine**
* **For TTKVTK**: ~/ttk/ttk_install/lib/cmake/ttkVTK 
* **For Paraview**: ~/ttk/ParaView-v5.9.1/build/
* **For VTK**: ~/ttk/ParaView-v5.9.1/build/lib/cmake/paraview-5.9/vtk
* **For VTKm**: ~/ttk/ParaView-v5.9.1/build/lib/cmake/paraview-5.9/vtk/vtkm

**For Linux Machine**
* **For TTKVTK**: /usr/local/lib64/cmake/ttkVTK
* **For Paraview**: ~/ttk/ParaView-v5.11.1/build/
* **For VTK**: ~/ttk/ParaView-v5.11.1/build/lib/cmake/paraview-5.11/vtk
* **For VTKm**: ~/ttk/ParaView-v5.11.1/build/lib/cmake/paraview-5.11/vtk/vtkm

Once the makefile is generated, the code is simply compiled using
**make -j**  and the executable should be located in the same directory
for a linux machine and in tda-segmentor.app/Contents/MacOS/ on a Mac machine.
 

# tda-segmentor documentation

We provide here some examples on the usage of the tda-segmentor software package. 


## Generation of Input Grids

The tda-segmentor accepts two types of input grids - **distance grids** and **energy grids**. Both these grids are scalar functions defined on a discrete grid, from which topological information can be extracted. 

### Distance grid generation

For nanoporous materials, the distance grids can be generated from the software package [Zeo++](http://www.zeoplusplus.org/). The grids are saved in the Gaussian cube format with an extension *.cube* and can be read directly into the tda-segmentor code. 

### Energy grid generation

The energy grids are generated from the interaction of a guest molecule with the nanoporous material. These can be generated from the Julia software package [PorousMaterials.jl](https://github.com/SimonEnsemble/PorousMaterials.jl). This generates too a grid file in the Gaussian cube format (*.cube*) and can be read directly into the tda-segmentor code. Note that the energy grid generated can have points with **Inf** values at certain grid points; such strings are identified and are fixed in  the tda-segmentor code. 

## Usage

The code has many functionalities that are implemented as various modules such as the **persistence curve**, **accessible void graph**, **solid segmentation**, **accessible solid graph**, etc. Once an analysis is run, the code generates a log file on which lots of information such as the number of grid points, grid resolution, and a commentary on the analysis are stored. All the results of an analysis are stored in the folder **segmentor-baseFileName.results/**, where if the input file name is *FAU.cube*, the baseFileName is *FAU*. 

Given an input grid with limited information, the first step in topological analysis would be to plot the persistence curve of the scalar field. 

### Persistence Curve

The persistence curve for the input grid can be obtained through the command : 

    tda-segmentor -pc -useallcores FAU.cube

This generates four .txt files which saves the number of critcal pairs as a function of the persistence; this can be plotted as shown below.

![](https://github.com/AMDatIMDEA/tda-segmentor/blob/main/images/persistence-curve.png?raw=true)

The persistence curve helps in estimating the **persistence threshold** which is used to remove noise in the analysis. All critical pairs below the **persistence threshold** will be ignored, and much cleaner segments will be generated. A good estimate of the persistence threshold is to choose a value just before the plateau, as the plateau separates noise from actual geometric features. The blue vertical line can be chosen for example as a persistence threshold for this example. **NOTE:** If TTK is compiled correctly with MPI, the flag *-useallcores* uses the multicore capabilities of TTK and gives significant speed-up. 

### Morse Smale Complex

This is the main segmentation routine of the code and is run by the command 

    tda-segmentor -msc 0.04 -useallcores FAU.cube 
   
This routine takes an optional **persistence threshold** which is input here as 0.04. If no persistence threshold is given then 1% of the maximum persistence is taken as the persistence threshold. The -msc analysis saves two files : the *FAU_CriticalPoints.vtk* that saves the critical points data, and *FAU_Segmentation.vtk*, that saves the segmentation data. These *.vtk* files can be viewed using *Paraview*. For the above analysis, the segmentation drawn on a contour of distance 1.6 is as shown below. 

![](https://github.com/AMDatIMDEA/tda-segmentor/blob/main/images/FAU-segmentation.png?raw=true)


### Void Segmentation

For a more high-throughput analysis, we might want to output the segmentation data of a certain space of the scalar field as easily readable files. The command
    
    tda-segmentor -vs 0.04 -useallcores FAU.cube

saves segmentation information of the void space written to a .csv file, *FAU_Void_Segments.csv* in addition to the *.vtk* files, *FAU_Segmentation.vtk* and *FAU_CriticalPoints.vtk*. The .csv file has 10 columns viz. regionID (ID of the segment), x y z (coordinates of the grid point), Scalar (value of distance function at the grid point), RegionMaxValue (maximum value of the distance function in the segment), isMaximum, isSaddle, (if the grid point corresponds to a maxima or 2-saddle), numberOfPoints (number of grid points in the segment), Volume (volume of the segment), numberOfConnections (number of segments the current segment is connected to). The column numberOfConnections can be very useful in identifying isolated segments. If numberOfConnections is 0, then they are completely disconnected, and are inaccessible to a guest molecule. **NOTE:** Here the persistence threshold is an optional parameter and *tda-segmentor -vs -useallcores FAU.cube* can also be used; the persistence threshold is taken as 1% of the maximum value, but it is recommended to give the persistence threshold as an input. 

### Solid Segmentation

The solid segmentation routine does the exact same thing as void segmentation, but for the solid space. It is invoked by,

    tda-segmentor -ss 0.04 -useallcores FAU.cube

This saves too a .csv file, *FAU_Solid_Segments.csv*, in addition to the segmenation *.vtk* files. This file too has 10 columns but as the solid space analyzes points with distance function < 0, we have minimas within the segments and we save minimum value in each segment. The 10 colums are : regionID (ID of the segment), x y z (coordinates of the grid point), Scalar (value of distance function at the grid point), RegionMinValue (minimum value of the distance function in the segment), isMinima, isSaddle, (if the grid point corresponds to a minima or 1-saddle), numberOfPoints (number of grid points in the segment), Volume (volume of the segment), numberOfConnections (number of segments the current segment is connected to). **NOTE:** Here the persistence threshold is an optional parameter and *tda-segmentor -ss -useallcores FAU.cube* can also be used; the persistence threshold is taken as 1% of the maximum value, but it is recommended to give the persistence threshold as an input. 

### Accessible Void Space

This routine saves the segment information in a .csv file accessible to a guest molecule to the nanoporous material. This is invoked by:

    tda-segmentor -avs 0.04 1.6 -useallcores FAU.cube
 
Here, the first parameter, the persistence threshold 0.04 is optional, but the second parameter, the radius of the guest molecule is mandatory. Thus, this command also accepts input of the form

    tda-segmentor -avs 1.6 -useallcores FAU.cube

where the persistence threshold is chosen automatically as 1% of the maximum persistence. This routine too saves in a .csv file, *FAU-avs-radiusVal.csv*, the segmentation information. This file also has 10 columns and have exactly the same meaning as Void Segmentation. Note that the difference between this routine and the void segmentation routine, is that void segmentation is for the entire void space (distance function > 0.0), but accessible void space, is only the void space accessible to a guest molecule (distance function > radius of guest molecule). 

### Accessible Void Graph

From the morse-smale complex, this module constructs a graph representation of the accessible void space. This routine is invoked by:

    tda-segmentor -avg 0.04 1.6 -useallcores FAU.cube

Like the accessible void space module, *-avg* too accepts two arguments - an optional persistence threshold and a mandatory radius of the guest molecule. If the persistence threshold is to be chosen automatically as 1% of the max. persistence, then, 

    tda-segmentor -avg 1.6 -useallcores FAU.cube
    
will also work. This command saves the segmentation data as *.vtk* files, i.e. *FAU_Segmentation.vtk* and *FAU_CriticalPoints.vtk* and also another visualization, *FAU_viz_graph.vtk* for the visualization of the graph as shown below.

![](https://github.com/AMDatIMDEA/tda-segmentor/blob/main/images/FAU-avg.png?raw=true)

In the above picture, on the left is the segmentation and on the right is the graph representation. Red spheres are the local maximas, and blue spheres are the 2-saddles. In topological terms, this is the Reeb Graph connecting 3-maximas and 2-saddles restricting to scalar data of distance greater than the radius of the guest atom. The graph representation is periodic and just for visualization the periodic node is plotted in gray spheres in the neighboring box. 

This routine also writes the graph in a *.nt2* file that can be used for post-processing. First all the nodes of the graph are stored with an ID, and its X, Y, Z, coordinates respectively. Then the edges are stored with the first two columns indicating the nodes they connect; let us call first the birth and the second the death. Then the next three columns indicates periodicity in X, Y, Z, directions respectively. If the value is 0, then it is not periodic along this axis. If the value is 1, then the death connects the birth along the positive direction of the axis. If the value is -1, then the death connects the birth along the negative direction of the axis. 

### Accessible Solid Graph

From the morse-smale complex, this constructs the graph representation of the solid space. This is invoked by: 

    tda-segmentor -asg 0.01 random-structure.cube
    
Like the -ss module, this accepts an optional persistence threshold, and if not provided choses an automatic 1% of the maximum persistence. This commands saves the segmentation data as *.vtk* files, i.e., *random-structure_Segmentation.vtk* and *random-structure_CriticalPoints.vtk* and also another visualization *.vtk* file, *random-structure_viz_graph.vtk* for the visualization of the graph as shown below.

![](https://github.com/AMDatIMDEA/tda-segmentor/blob/main/images/random-structure-asg.png?raw=true)

In the above picture, red spheres are the local minimas, and blue spheres are the 1-saddles. The graph representation is periodic and just for visualization the periodic node is plotted in gray spheres in the neighboring box. This routine also writes the graph in a *.nt2* file that can be used for post-processing, and has the exact same format as the accessible void graph module. 

### Other flags

Some other flags can also be given to the tda-segmentor tool. 

* -writefractionalgrid : writes the input grid in the fractional coordinates from 0 to 1. This can be used for testing directly on Paraview, where many of the TTK modules can be applied with the GUI interface. 
* -savelogfile : If an analysis is already executed in a folder, then this flag, preserves the older log file before running a new analysis. 
* -usesupercell : For very small lattices, it might be necessary to use the super cell for analysis. By using this flag, a supercell is generated that is used for segmentation. The invocation looks as follows: 

        tda-segmentor -avg 0.03 0.5 -usesupercell -useallcores NPO.cube

  For the zeolite NPO, the supercell segmentation and its graph representation for a guest molecule of radius 0.5 Angstrom is 
  as shown below
  ![](https://github.com/AMDatIMDEA/tda-segmentor/blob/main/images/NPO-supercell.png?raw=true)
  
* tda-segmentor can also be executed without any options, i.e.

        tda-segmentor FAU.cube
  
  will work. This by default runs *-msc* module, with an automatic choice of persistence threshold. 

* One can use a combination of modules such as :

        tda-segmentor -msc 0.01 -avg 0.01 1.6 -avs 0.01 1.6 -useallcores -usesupercell FAU.cube
        
  Such commands will also work.

* Two macro variables **DEBUG** and **NITERATIONS** are defined in the files *grid.h* and *headers.h* respectively that can be changed if necessary. Setting **DEBUG** to 1, dumps additional information to the log file such as segment connectivity, while **NITERATIONS** by default takes a value 2, which actually smoothens the input data by **NITERATIONS** iterations. 

# Copyright notice

Copyright (c) 2023, Fundacion IMDEA Materiales 
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of tda-segmentor nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# Acknowledgements

This work was supported by the Plan de Recuperación, Transformación y Resiliencia of the European Union, through the project (MAD2D) - Two-dimensional (2D) disruptive materials for the new technological transformation. 

# References

If you have used this package, please consider citing the following in your technical documents or scientific publications:

* TDA Segmentor - Tools to extract and analyze local structure and porosity features in porous materials, A. Vasudevan, J.  Z. Prieto, S. Zorkaltsev, M. Haranczyk. 

# Contact Information

We are always excited to receive feedback on the tool, and if you may have any additional questions, please find the contact information of the developers below: 

* Aditya Vasudevan (adityavv.iitkgp@gmail.com)
* Maciej Haranczyk (maciej.haranczyk@imdea.org)
