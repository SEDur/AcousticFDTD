AcousticFDTD
An acoustic finite-difference time domain engine based in Python.
-----------------------------------------------------------------
The core FDTD engine will be written in C++, and will feature a Python wrapper.
<br>
Future implementation will see a number of features (roughly in order):<br>
<ul>
  <li>GUI</li>
  <li>2D solving for homogeneous geometries using predescribed nodes</li>
  <li>3D solving for homogeneous geometries using predescribed nodes</li>
  <li>Convergeance and solution analysis tools</li>
  <li>Flexible reciever design</li>
  <li>Multipoint source design</li>
  <li>Directional source design</li>
  <li>Heterogeneous solutions within problem space</li>
  <li>Viscous losses</li>
  <li>Rectilinear meshing of CAD files</li>
  <li>PML bounded problem spaces</li>
  <li>Parallel CUDA GPU implementation</li>
</ul>
------------------------------------------------------------------
The elements of the program will be prototyped in Matlab
Successful prototype will lead to implementation in a seperate branch
