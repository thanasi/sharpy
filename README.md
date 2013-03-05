sharpy 0.1
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
4. Implement as library
6. STL Registration along principal axes : x = minor axis, y = middle, z = major axis
7. Monte Carlo integral against different Ylm to get eigenvalues
8. Calculate SphericalHarmonicTransform
9. Calculate reverse transform
10. Visualize decomposition
11. Visualize reverse transform and compare to STL