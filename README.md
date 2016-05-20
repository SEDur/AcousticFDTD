# AcousticFDTD
An acoustic finite-difference time domain engine based in Python.
-----------------------------------------------------------------
The core FDTD engine will be written in C++, and will feature a Python wrapper.

Future implementation will see a number of features (roughly in order):

     2D solving for homogeneous geometries using predescribed nodes
      3D solving for homogeneous geometries using predescribed nodes
-----Convergeance and solution analysis tools
  |  Flexible reciever design
  |  Multipoint source design
 GUI Directional source design
  |  Heterogeneous solutions within problem space
  |  Viscous losses
  |  Rectilinear meshing of CAD files
  |  PML bounded problem spaces
-----Multiple CUDA GPU implementation

------------------------------------------------------------------
The elements of the program will be prototyped in Matlab
Successful prototype will lead to implementation in a seperate branch
