# NEXD 3D

NEXD is a software package for high order simulation of seismic waves using the nodal discontinuous Galerkin method.
NEXD 3D is the 3D solver of this software package.

To deal with highly complex heterogeneous models, here the Nodal Discontinuous Galerkin Method (NDG) is used to calculate synthetic seismograms. The advantage of this method is that complex mesh geometries can be computed by using triangular or tetrahedral elements for domain discretization together with a high order spatial approximation of the wave field. The simulation tool NEXD is developed which has the capability of simulating elastic and anelastic wave fields for seismic experiments for one-, two- and three- dimensional settings. The implementation of poroelasticity and simulation of slip interfaces are currently in progress and are working for the one dimensional part. External models provided by e.g. Trelis/Cubit can be used for parallelized computations on triangular or tetrahedral meshes. For absorbing boundary conditions either a fluxes based approach or a Nearly Perfectly Matched Layer (NPML) can be used. 

## Authors and License

NEXD 3D and some of its components, as well as documentation and some examples
are available under terms of the [GNU General Public License](LICENSE) (version 3 or higher)
on [github](https://github.com/seismology-RUB/NEXD-3D).
Please find contact addresses [there](https://github.com/seismology-RUB), or visit 
http://www.rub.de/nexd in case you want to get in touch with the authors. If you 
encounter any problems installing or using the software, it will be helpful to 
open (or add to) an "issues" topic at the [github repository](https://github.com/seismology-RUB/NEXD-3D).

The main authors are Lasse Lambrecht, Andre Lamert, and Wolfgang Friederich (Ruhr-Universität Bochum, Germany).


## Documentation

A complete documentation will be available in the future.


## Examples

NEXD 3D comes with two examples: [simulations/brick](simple example) and a [simulations/2layerbox](two layer problem). Please use these examples to familiarize yourself with the software.


## Requirements

* GNU Make
* Fortran compiler (sufficient standard)
* LAPACK libraries
* MPI libraries (e.g., OpenMPI) and METIS (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview; tested with version 4.0.3) for parallel applications


## Installation

0. If not yet done, you should download the source code of the NEXD 3D main package by cloning the master branch of the ASKI repository on gitHub.com:
     ```
     git clone --depth 1 --branch master https://github.com/seismology-RUB/NEXD-3D
     ```

1. Install all software dependencies.

2. Adjust the software to your system and personal requirements by changing the [Makefile](Makefile) appropriately (e.g., change the path to your METIS installation and set your compiler).

3. Run command
     ```
     make all
     ```
   from installation path `NEXD-3D/` to compile NEXD 3D.
