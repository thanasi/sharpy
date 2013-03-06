sharpy 0.2
==================

Expand a 3D surface as spherical harmonics Ylm using Monte Carlo integration

####Requirements:
1. python (dev using python3)
2. vtk with python bindings

####Input File Formats Supported:
- STL

####ToDo:

1. ~~Import STL~~ (testscripts/testSTLPoints)
2. ~~Check if a point is contained within STL~~ (testscripts/testSTLPoints)
3. ~~Visualize points + color based on if in or out~~ (testscripts/testSTLPoints)
4. ~~Implement as library + class~~ (sharpy/shape3d.py, testscripts/testSTLPoints2)
5. ~~Calculate moment of inertia tensor + diagonalize~~
	- speed up calculation for large N (parallelize or use numpy tricks. or make a c module)
6. Visualize eigenvectors of moment of inertia tensor
7. STL Registration along principal axes of moment of inertia tensor:
	- center of mass at (0,0,0)
	- x = major axis
	- y = middle
	- z = minor axis
8. Monte Carlo integral against different Ylm to get eigenvalues
9. Calculate SphericalHarmonicTransform
10. Calculate reverse transform
11. Visualize decomposition (eigenvector YLM sized based on eigenvalue)
12. Visualize reverse transform and compare to STL