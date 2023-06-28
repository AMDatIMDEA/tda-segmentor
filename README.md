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
  * Next we checkout to a certain commit : **git checkout 6ff3c19c2**
  * The above code will be the source code for installation. 

#### Compiling the tda-segmentor code:

**cmake** should be used for generating the makefile. Move to the respository directory tda-segmentor/
and type **ccmake .** and one might have to hit configure *'c'* a few times and finally *'g'* to 
generate the makefile. Once the makefile is generated, the code is simply compiled using
*make -j*
 
* Specific instructions for MAC machine

* Specific instructions for Linux machine

### Usage
* Invocation syntax

* List of modules

* Optional flags

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


