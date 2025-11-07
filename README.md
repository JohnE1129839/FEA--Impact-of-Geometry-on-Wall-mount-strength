# Description
This repository is a part of an assignment for the ME4291 Module in NUS regarding Finite Element Analysis. 
The objective of this object is to use Finite Element Analysis to evaluate the effect of changing a fixture's geometry on a drilled wall-mount.
This is done by using MATLAB code to conduct a 2D triangular element analysis on different meshes representing different possible geometries of the mount.

# Objective
The objective of this project is to understand how different shapes and sizes of the support can affect the strength of the mount. The relationship may be complex and 
non-linear but it is usefull to have a general sense on how stronger mounts can be achieved using the same type and volume of material.

# Problem Statement
The mount is modeled to be a planar square made of alloy steel. It can be assumed to have a uniform thickness of 5mm and has the size of a (10 x 10) cm hole.
A hole of certain shape and size is drilled onto the mount and a rigid support is passed through the hole. It can be assumed that the support is far stronger and 
more rigid than the mount. The task is to determine the stress distribution within the mount given a load of 1 Newton is applied to the mount.

# Methods used
The main method used in this project is 2D planar stress finite element analysis, more specifically using linear triangular elements. The project uses the 'pde' add-on in MATLAB
to generate meshes for the suggested geometries made in SolidWorks. These meshes are then converted into nodal vectors and stiffness matrices which shall be used to 
determine the nodal displacement and the Von Mises stress distribution throughout the structure.

# Software required
This project already contains some mesh models which can be used for analyses. These meshes can be used in MATLAB without any additional downloads required. 
However, to generate a new mesh for a certain geometry from a STL file, it is necessary to download the 'pde' add-on on MATLAB.

# Instructions on usage
There are 2 main files useful for this project, **simulation.m** and **generateMesh.m**. Simulation.m is used as the main file to run the analysis in the project. 
To run the simulation, simply specify the mesh being used in the simulation via the *fileName* variable in the file.  
Some notable meshes already included within this project are:
- c1, c2, c3 -> Circular holes within the centre of the mount with the diameters 5cm, 4cm, and 3cm
- s1_rounded, s2_rounded, s3_rounded -> Square holes in the centre of the mount with lengths 5cm, 4cm, and 3cm
- s1_rounded_rotated, s2_rounded_rotated, s3_rounded_rotated -> Diamond shape holes with diagonals 5cm, 4cm, and 3cm
- c3_off1, c3_off2, c3_off3, c3_off4 -> Variations of c3 with different offsets ranging from +2cm to -2cm  

To generate a Mesh using the generateMesh.m file, a STL file of the specified geometry must be added to the STLs folder of the project first. The user would then have to specify the
STL file name through the *filename* variable and verify the edge selected. As of now, the only way of verifying the specified internal edges of the hole is visually. The user will
have to specify the internal edges of the mount through the *edgeIDs* variable. After the code has run, three text files containing the information of the mesh should be generated in the
Meshes folder, it is now possible to use simulation.m to test the mesh.
