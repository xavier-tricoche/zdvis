# zdvis
**Lagrangian Visualization for Vector, Tensor, and Multifield Data.**

Visualization library offering a unified Lagrangian approach to the structural analysis and visual representation of vector, tensor, and multi-field datasets.

This library contains the code written by Zi'ang Ding as part of his PhD thesis titled **"LAGRANGIAN ANALYSIS OF VECTOR AND TENSOR FIELDS: ALGORITHMIC FOUNDATIONS AND APPLICATIONS IN MEDICAL IMAGING AND COMPUTATIONAL FLUID DYNAMICS"** in the Computer Science department at Purdue University from 2010 to 2016. Examples are included that document the proposed approach in fluid dynamics and medical imaging applications.

This project requires following third-party open source libraries:
* Teem: http://teem.sourceforge.net/
* Eigen: http://eigen.tuxfamily.org/
* VTK: https://www.vtk.org/

In addition, the compilation of some examples requires NetCDF: https://www.unidata.ucar.edu/software/netcdf/

The test files included in the library use some reference datasets that have been used in publications. These datasets include 2D and 3D time-dependent CFD datasets along with a diffusion MRI dataset. They can be downloaded from ftp://ftp.cs.purdue.edu/pub/zdvis/datasets/. 

This work was supported in part by NSF OCI CAREER award 1150000 *"Efficient Structural Analysis of Multivariate Fields for Scalable Visualization"* (Xavier Tricoche, PI)
