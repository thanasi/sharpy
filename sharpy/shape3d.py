##############################
## 
## sharpy.Shape3D
## 
## Author: A. Athanassiadis
## Date: March, 2013
## 
##############################
import os
import numpy as np
import vtk

class ImportError(Exception):
    """ Base class for surface import errors. """
    pass


class Shape3D(object):
    """
    Shape3D
        object that holds a 3D surface and can check
        if a point is enclosed by the surface or not.
        
    Shape3D(filename) : load surface from filename
        
        Supported file types:
            --STL
    
    Shape3D.is_inside(points) : check if points in list are inside surface
        
        will accept single points or a list
        returns single binary value or an array of them
        
    Shape3D.visualize() : visualize surface and points checked
        
    """
    
    def __init__(self, fn):
        self._filename = fn
        
        self._loaded = False
        self._aligned = False
        
        ext = np.char.lower(os.path.splitext(fn)[-1])
        
        if ext == ".stl":
            self.load_stl()
            
        else:
            raise ImportError("Shape3D cannot currently load %s files" % ext)

    def load_stl(self):
        """" Load an STL file into the Shape3D object """
        
        ## set up reader
        self._reader = vtk.vtkSTLReader()
        self._reader.SetFileName(self._filename)
        self._reader.Update()
        
        ## set up surface
        self.surf = self._reader.GetOutput()
        self._bounds = self.surf.GetBounds()
                
        
        self._loaded = True
        
        
    def is_inside(self,points):
        """ Check if given points lie inside the surface """
        
        points = np.array(points)
        if points.ndim==1:
            points = np.array([points]) 
        
        assert points.shape[1] == 3, "input point must be x,y,z"

        n = points.shape[0]
        
        ## set up points to check
        self.vpoints = vtk.vtkPoints()
        for i in range(n):
            self.vpoints.InsertNextPoint(points[i])
        
        checkPoints = vtk.vtkPolyData()
        checkPoints.SetPoints(self.vpoints)
        
        ## set up point checking object
        pointChecker = vtk.vtkSelectEnclosedPoints()
        pointChecker.SetInput(checkPoints)
        pointChecker.SetSurface(self.surf)
        pointChecker.Update()
        
        ## check the status for each point
        inout = []
        inPoints0 = vtk.vtkPoints()
        outPoints0 = vtk.vtkPoints()
        
        for i in range(checkPoints.GetNumberOfPoints()):
            inout.append(pointChecker.IsInside(i))
            # print i, inout[-1]
            if inout[-1]:
                inPoints0.InsertNextPoint(self.vpoints.GetPoint(i))
            else:
                outPoints0.InsertNextPoint(self.vpoints.GetPoint(i))
        
        self._inPoints = vtk.vtkPolyData()
        self._inPoints.SetPoints(inPoints0)
        self._outPoints = vtk.vtkPolyData()
        self._outPoints.SetPoints(outPoints0)
        
        if n==1:
            return inout[0]
        else:
            return np.array(inout)
    
    def align_axes(self):
        """ 
        Align the imported surface according to the following rules:
            - x axis is minor principal axis
            - y axis is middle axis
            - z axis is major principal axis
        
        """
        
        # self._aligned = True
        pass
        
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
        
        if self._inPoints.GetNumberOfPoints() > 0:
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

        if self._outPoints.GetNumberOfPoints() > 0:
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
       