# tda-segmentor
A segmentation tool for the analysis of porous structures 
using the TTK toolkit (https://topology-tool-kit.github.io/)

### Installation Instructions
Compiling the tda-segmentor code should be done in two steps :

* Install the ttk library from source.
* Compile the tda-segmentor code linking the ttk library.

##### Installing the Topology Toolkit (TTK)
* TTK must be installed from the source code and detailed instructions
  can be found here (https://topology-tool-kit.github.io/installation.html)
* To install TTK from source, first Paraview must be compiled from sources, 
  and on a linux machine it is recommended to install ParaView-v5.9.1
  while on a Mac machine it is recommended to install ParaView-v5.10.1
* Once Paraview is installed, to install TTK, we must use the developer version 
  of TTK directly from the GitHub repo (https://github.com/topology-tool-kit/ttk) 
  using the following commands :
  * **git clone https://github.com/topology-tool-kit/ttk.git**
  * Next we checkout to a certain commit : **git checkout 6ff3c19c2** (this is necessary for the periodic box module)
  * The above code will be the source code for installation and the installation should proceed as instructed in the website. 

#### Compiling the tda-segmentor code:

**cmake** should be used for generating the makefile. Move to the respository directory tda-segmentor/
and type **ccmake .** and one might have to hit configure *'c'* a few times and finally *'g'* to 
generate the makefile. While generating the makefile, the correct paths for vtk, paraview and ttk need
to be provided. If the directories are created as recommended by the TTK website, then these paths are: 

* **For TTKVTK**: ~/ttk/ttk_install/lib/cmake/ttkVTK 
* **For Paraview**: ~/ttk/ParaView-v5.9.1/build/
* **For VTK**: ~/ttk/ParaView-v5.9.1/build/lib/cmake/paraview-5.9/vtk
* **For VTKm**: ~/ttk/ParaView-v5.9.1/build/lib/cmake/paraview-5.9/vtk/vtkm

Once the makefile is generated, the code is simply compiled using
**make -j**  and the executable should be located in the same directory
for a linux machine and in tda-segmentor.app/Contents/MacOS/ on a Mac machine.
 

### Usage
* Invocation syntax
Once compiled correctly, the code is called directly from the command-line. 
Different functionalities of the code are implemented as modules and are called as: 

**poretda -module *moduleName* INPUT-FILE**

Input parameters are given by its keyword followed by its value. For example, 
the module **accessiblevoidspace** needs the radius of the probe atom and 
is input as keyword followed by its value. (**-proberadius rad_val**)

Some optional flags are provided such as **-usesupercell**, **-usetbb**, **-savelogfile**
that can serve additional functionality as described below and in Docs/ 

* List of all modules: 
    * **-module segmentation** : This is the segmentation module that segments the 
      structure and voids into many segments. It requires an input parameter persistencethreshold
      and is given as **-persistencethreshold value**. If this module is called without an input
      parameter, then a default value 0.1 is chosen for value. 
    * **-module accessiblevoidspace**: Generates the accessible void space on the segmented voids.
      This module requires the radius of a probe atom and is input as (**-proberadius rad_val**). This
      module implicitly also runs the segmentation module and hence **-persistencethreshold** must also be given.  
    * **-module voidsegmentation** : Segments the void space and dumps results in a .csv file. Requires persistencethreshold parameter. **(explain info on the .csv file in Output section)** 
    * **-module solidsegmentaiton** : Segments the solid space and dumps results in a .csv file. Requires persistencethreshold parameter. **(explain info on the .csv file in Output section)**
    * **-module persistencecurve** : Save the persistence curve without segmenting the structure. **(to be implemented)**

* List of all optional flags: 

    * **-usesupercell** : For some structures, the neighboring 26 unit cells are also used for segmentation
                          (Beware this consumes a lot of memory!)
    * **-usetbb**: If ttk is installed correctly with intel threading blocks, then multithreading flag can be enabled for faster computation. (**needs to be tested**)
    * **-savelogfile** : saves the log file of the previous run, if runs are being made from the same folder.  (**needs to be implemented**)

* Some invocation examples: 

    * **poretda -module segmentation -module accessiblevoidspace -persistencethreshold 0.12 -proberadius 1.6 FAU.cube**
    * **poretda -module segmentation -persistencethreshold 0.12 -usesupercell FAU.cube**
    * **poretda -module persistencecurve FAU.cube**

* Output
    * Lots of information is written in the basefilename.log and all the results are stored in tda-segmentor.results/ folder. 
### Copyright notice

### Licensing

### Acknowledgements
Acknowledge funding sources

### References
Cite preprint

### Contact Information

For additional questions on the code, contact: 

* Aditya Vasudevan (aditya.vasudevan@imdea.org)
* Maciej Haranczyk (maciej.haranczyk@imdea.org)


