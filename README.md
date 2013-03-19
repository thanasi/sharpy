sharpy 0.5
==================

[https://github.com/thanasi/sharpy](https://github.com/thanasi/sharpy)

Expand a 3D surface as spherical harmonics Ylm using Monte Carlo integration

####Requirements:
1. python (dev using python2.7)
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
6. ~~Visualize eigenvectors of moment of inertia tensor~~
7. ~~STL Registration along principal axes of moment of inertia tensor:~~
	- center of mass at (0,0,0)
	- x = major axis
	- y = middle
	- z = minor axis
	- **this is still slow/fails for meshes with many points
8. ~~Write surface to output STL~~
9. Monte Carlo STintegral against different Ylm to get eigenvalues
10. Calculate SphericalHarmonicTransform
11. Calculate reverse transform
12. Visualize decomposition (eigenvector YLM sized based on eigenvalue)
13. Visualize reverse transform and compare to STL