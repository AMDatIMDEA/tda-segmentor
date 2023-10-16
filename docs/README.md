# tda-segmentor documentation

We provide here some examples on the usage of the tda-segmentor tool. 


## Generation of Input Grids

The tda-segmentor at the moment accepts two types of input grids - **distance grids** and **energy grids**. Both these grids are two discrete functions defined on a grid, from which topological information can be extracted. 

### Distance grid generation

For nanoporous materials, the distance grids can be generated from the software package [Zeo++](http://www.zeoplusplus.org/). These files are saved in the Gaussian cube format with an extension *.cube* and can be read directly into the tda-segmentor code. 

### Energy grid generation

The energy grids are generated from the interaction of a guest molecule with the nanoporous material. These can be generated in a high-throughput way from the Julia software package [PorousMaterials.jl](https://github.com/SimonEnsemble/PorousMaterials.jl). This generates too a grid file in the Gaussian cube format (*.cube*) and can be read directly into the tda-segmentor code. Note that the energy grid generated can have points with **Inf** values at a certain grid point; such strings are identified and are alread fixed in the reader function of the tda-segmentor code. 

## Usage


