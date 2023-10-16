# tda-segmentor documentation

We provide here some examples on the usage of the tda-segmentor tool. 


## Generation of Input Grids

The tda-segmentor accepts two types of input grids - **distance grids** and **energy grids**. Both these grids are scalar functions defined on a discrete grid, from which topological information can be extracted. 

### Distance grid generation

For nanoporous materials, the distance grids can be generated from the software package [Zeo++](http://www.zeoplusplus.org/). These files are saved in the Gaussian cube format with an extension *.cube* and can be read directly into the tda-segmentor code. 

### Energy grid generation

The energy grids are generated from the interaction of a guest molecule with the nanoporous material. These can be generated from the Julia software package [PorousMaterials.jl](https://github.com/SimonEnsemble/PorousMaterials.jl). This generates too a grid file in the Gaussian cube format (*.cube*) and can be read directly into the tda-segmentor code. Note that the energy grid generated can have points with **Inf** values at some certain grid points; such strings are identified and are alread fixed in the reader function of the tda-segmentor code. 

## Usage

The code has many functionalities that are implemented as various modules such as the **persistence curve**, **accessible void graph**, **solid segmentation**, **accessible solid graph**, etc. Once an analysis is run, the code generates a log file on which lot of information such as the number of grid points, grid resolution, and a commentary on the analysis are stored. All the results of an analysis are stored in the folder **segmentor-baseFileName.results/**

Given an input grid with limited information, the first step can be to plot the persistence curve of the scalar field. 

### Persistence Curve

The persistence curve for the input grid can be obtained through the command : 

    tda-segmentor -pc -useallcores FAU.cube

This generates four .txt files which saves the number of critcal pairs as a function of the persistence; this can be plotted as shown below.



The persistence curve helps in estimating the **persistence threshold** which is used to removing noise in the analysis. All critical pairs below the **persistence threshold** will be ignored, and much cleaner segments will be generated. A good estimate of the persistence threshold is to choose a value just before the plateau, as the plateau separates noise from actual geometric features. 