##############################
## 
## sharpy.Shape3D
## 
## Author: A. Athanassiadis
## Date: March, 2013
## 
##############################
import os
from itertools import izip
import numpy as np
from numpy import linalg as LA
import vtk

## tolerance for point merging
TOL = 1e-4

class ImportError(Exception):
    """ Base class for surface import errors. """
    pass


class Shape3D(object):
    """
    Shape3D
        object that holds a 3D surface and can check
        if a point is enclosed by the surface or not.
        
    Shape3D(filename)
        load surface from filename
        removes duplicate points within a tolerance of 1e-4
        Supported file types:
            --STL
    
    Shape3D.is_inside(points)
        check if points in list are inside surface
        will accept single points or a list
        returns single binary value or an array of them
        
    Shape3D.visualize()
        visualize surface and points checked
        
    """
    
    def __init__(self, fn, clean=True):
        self._filename = fn
        
        self._loaded = False
        self._aligned = False
        self._cleaned = False
        self._pointch = False
        self._pointcloudgen = False
        
        self.Surf = None
        self.InertiaTensor = None
        self.PointCloud = None
        self.CoM = None
        
        ext = np.char.lower(os.path.splitext(fn)[-1])
        
        if ext == ".stl":
            self._load_stl()
            
            if clean:
                self._clean_dup_points()
            
        else:
            raise ImportError("Shape3D cannot currently load %s files" % ext)

    def _load_stl(self):
        """" Load an STL file into the Shape3D object """
        
        ## set up reader
        self._reader = vtk.vtkSTLReader()
        self._reader.SetFileName(self._filename)
        self._reader.Update()
        
        ## set up surface
        self.Surf = self._reader.GetOutput()
        self._bounds = self.Surf.GetBounds()
        
        self._loaded = True
        
    def _clean_dup_points(self):
        """ Remove duplicate points within TOL in the loaded surface"""     
        
        scrubber = vtk.vtkCleanPolyData()
        scrubber.SetTolerance(TOL)
        scrubber.SetInput(self.Surf)
        scrubber.Update()
        
        N1 = self.Surf.GetNumberOfPoints()
        N2 = scrubber.GetOutput().GetNumberOfPoints()
        
        if N2<N1:
            print "Removed %d duplicate points" % (N1-N2)
            self.Surf = scrubber.GetOutput()
        else:
            print "No duplicate points within tolerance"

        self._cleaned = True   
        
    def _align_axes(self):
        """ 
        Align the imported surface according to the following rules:
            - move center of mass to (0,0,0)
            - x axis: minor axis of the moment of inertia tensor
            - y axis: middle axis
            - z axis: major  axis
        
        """
        
        # 
        
        # self._aligned = True
        pass
    
    
    def generate_pointcloud(self, N=1000000):
        """ Build an N-point pointcloud representation of the surface """

        x = np.random.uniform(low=self._bounds[0], high=self._bounds[1], size=N)
        y = np.random.uniform(low=self._bounds[2], high=self._bounds[3], size=N)
        z = np.random.uniform(low=self._bounds[4], high=self._bounds[5], size=N)
        
        pc = np.array([x,y,z]).T
        
        inbool = self.is_inside(pc,save=False)
        
        self.PointCloud = pc[inbool==1]
                
    def calc_InertiaTensor(self):
        """ Calculate moment of Inertia Tensor """
        
        self.CoM = self.PointCloud.mean(0)
        
        pc = self.PointCloud.copy() - self.CoM
        
        ## moment of inertia tensor
        I = np.ones((3,3))
        
        gr = np.mgrid[:3,:3]
        delta = np.eye(3)
        
        for p in pc:
            crossmat = np.outer(p,p)
            I += delta * np.sum(p**2) - crossmat

        self.InertiaTensor = I
                    
        # calculate and sort eigenvectors
        eigval, eigvec = LA.eigh(I)
        
        eigsys = zip(eigval, eigvec.T)
        eigsys.sort(key = lambda x: x[0])
        
        eigval, eigvec = zip(*eigsys)
        
        self._eigval = eigval[::-1]
        self._eigvec = eigsys[::-1]
        
        print I
        print eigval
        print eigvec       
        
    
    def is_inside(self,points,save=True):
        """ Check if given points lie inside the surface """
        
        points = np.array(points)
        if points.ndim==1:
            points = np.array([points]) 
        
        assert points.shape[1] == 3, "input point must be x,y,z"

        n = points.shape[0]
        
        ## set up points to check
        vpoints = vtk.vtkPoints()
        for p in points:
            vpoints.InsertNextPoint(p)
        
        checkPoints = vtk.vtkPolyData()
        checkPoints.SetPoints(vpoints)
        
        ## set up point checking object
        pointChecker = vtk.vtkSelectEnclosedPoints()
        pointChecker.SetInput(checkPoints)
        pointChecker.SetSurface(self.Surf)
        pointChecker.Update()
        
        ## check the status for each point
        inout = []
        inPoints0 = vtk.vtkPoints()
        outPoints0 = vtk.vtkPoints()
        
        for i in range(checkPoints.GetNumberOfPoints()):
            inout.append(pointChecker.IsInside(i))
            # print i, inout[-1]
            if inout[-1]:
                inPoints0.InsertNextPoint(vpoints.GetPoint(i))
            else:
                outPoints0.InsertNextPoint(vpoints.GetPoint(i))
        
        if save:
            self._vpoints = vpoints
            self._inPoints = vtk.vtkPolyData()
            self._inPoints.SetPoints(inPoints0)
            self._outPoints = vtk.vtkPolyData()
            self._outPoints.SetPoints(outPoints0)
        
            self._pointch = True
        
        if n==1:
            return inout[0]
        else:
            return np.array(inout)
    
    def visualize(self):
        """ 
            Use VTK to visualize the geometry and 
            points if points have been checked 
        """ 
        
        ## shape mapper, actor
        shapeMapper = vtk.vtkPolyDataMapper()
        shapeMapper.SetInputConnection(self._reader.GetOutputPort())

        shapeActor = vtk.vtkActor()
        shapeActor.SetMapper(shapeMapper)
        shapeActor.GetProperty().SetOpacity(0.5)
    
        ## point mappers, actors
        useinpoints = 0
        useoutpoints = 0
        
        if self._pointch and self._inPoints.GetNumberOfPoints() > 0:
            useinpoints=1
            
            inverteces = vtk.vtkVertexGlyphFilter()
            inverteces.AddInput(self._inPoints)
            inverteces.Update()
    
            inpointsMapper = vtk.vtkPolyDataMapper()
            inpointsMapper.SetInputConnection(inverteces.GetOutputPort())
    
            inpointsActor = vtk.vtkActor()
            inpointsActor.SetMapper(inpointsMapper)
            inpointsActor.GetProperty().SetPointSize(5)
            inpointsActor.GetProperty().SetColor(0, 0, 1)
            inpointsActor.GetProperty().SetOpacity(0.75)    

        if self._pointch and self._outPoints.GetNumberOfPoints() > 0:
            useoutpoints=1
            
            outverteces = vtk.vtkVertexGlyphFilter()
            outverteces.AddInput(self._outPoints)
            outverteces.Update()
    
            outpointsMapper = vtk.vtkPolyDataMapper()
            outpointsMapper.SetInputConnection(outverteces.GetOutputPort())
    
            outpointsActor = vtk.vtkActor()
            outpointsActor.SetMapper(outpointsMapper)
            outpointsActor.GetProperty().SetPointSize(5)
            outpointsActor.GetProperty().SetColor(1, 0, 0)
            outpointsActor.GetProperty().SetOpacity(0.1)    
    
        ## renderer, render window, interactor
        renderer = vtk.vtkRenderer()
        renderWindow = vtk.vtkRenderWindow()
        interactor = vtk.vtkRenderWindowInteractor()
    
        renderWindow.AddRenderer(renderer)
        interactor.SetRenderWindow(renderWindow)
    
        ## add actors
        renderer.AddActor(shapeActor)
        if useinpoints:
            renderer.AddActor(inpointsActor)
        if useoutpoints:
            renderer.AddActor(outpointsActor)
            
        renderer.SetBackground(0,0,0)
    
        renderWindow.SetWindowName("Shape3DPoints")
        renderWindow.Render()
        interactor.Start();
       