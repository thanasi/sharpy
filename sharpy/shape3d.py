##############################
## 
## sharpy.shape3d
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

np.set_printoptions(precision=3)

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
    
    def __init__(self, fn, clean=True, inertia=True, align=True):
        self._filename = fn
        
        self._loaded = False
        self._aligned = False
        self._cleaned = False
        self._pointchecked = False
        self._pointcloudgen = False
        
        self.Surf = None
        self.InertiaTensor = None
        self.PointCloud = None
        self.CoM = None
        
        self._eigval = None
        self._eigvec = None
        
        ext = np.char.lower(os.path.splitext(fn)[-1])
        
        if ext == ".stl":
            self._load_stl()
            
            if clean:
                self._clean_dup_points()
            
            if inertia:
                print 'generating point cloud '
                self._generate_pointcloud()
                print 'calculating moment of inertia and eigenvecs'
                self._calc_inertia_tensor()
                print 'calculated eigenvecs'
                print self._eigvec
                
            if align:
                print 'aligning principal moments to x,y,z axes'
                self._align_axes()
            
        else:
            raise ImportError("Shape3D cannot currently load %s files" % ext)

    def _load_stl(self):
        """" Load an STL file into the Shape3D object """
        
        ## set up reader
        reader = vtk.vtkSTLReader()
        reader.SetFileName(self._filename)
        reader.Update()
        
        ## set up surface
        self.Surf = reader.GetOutput()
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
            - x axis: major axis of the moment of inertia tensor
            - y axis: middle axis
            - z axis: minor  axis
        
        """
        
        ## calculate rotation axes and angles
        
        ## align major principal axis
        
        rax1 = np.cross(self._eigvec[0],(1,0,0))
        rang1 = np.abs(np.arccos(self._eigvec[0][0]) * 180 / np.pi)
        
        txf = vtk.vtkTransform()
        txf.PostMultiply()
        txf.Translate(-1*self.CoM)
        txf.RotateWXYZ(rang1,rax1)
        
        txfPoly = vtk.vtkTransformPolyDataFilter()
        txfPoly.SetInput(self.Surf)
        txfPoly.SetTransform(txf)
        txfPoly.Update()
        
        self.Surf = txfPoly.GetOutput()
        self._bounds = self.Surf.GetBounds()
    
        print 'aligned surface. recalculating eigenvectors of MoI'
        self._generate_pointcloud()
        self._calc_inertia_tensor()
        
        print "updated eigenvectors"
        print self._eigvec
        
        self._aligned = True
        
    def _generate_pointcloud(self, N=1e4):
        """ Build an N-point pointcloud representation of the surface """

        x = np.random.uniform(low=self._bounds[0], high=self._bounds[1], size=N)
        y = np.random.uniform(low=self._bounds[2], high=self._bounds[3], size=N)
        z = np.random.uniform(low=self._bounds[4], high=self._bounds[5], size=N)
        
        pc = np.array([x,y,z]).T
        
        inbool = self.is_inside(pc,save=False)
        
        self.PointCloud = pc[inbool==1]
        print sum(inbool)
        
        self._pointcloudgen = True
                
    def _calc_inertia_tensor(self):
        """ Calculate moment of Inertia Tensor """
        
        self.CoM = self.PointCloud.mean(0)
        
        pc = self.PointCloud.copy() - self.CoM
        
        ## moment of inertia tensor
        gr = np.mgrid[:3,:3]
        delta = np.eye(3)
        
        f = lambda x: delta * np.sum(x**2) - np.outer(x,x)
        I = np.array(map(f,pc)).sum(0)
        
        self.InertiaTensor = I
                    
        # calculate and sort eigenvectors
        eigval, eigvec = LA.eigh(I)
        
        eigsys = zip(eigval, eigvec.T)
        eigsys.sort(key = lambda x: x[0])
        
        eigval, eigvec = zip(*eigsys)
        
        self._eigval = eigval[::-1]
        self._eigvec = np.array(eigvec)[::-1]
        
    
    def is_inside(self,points,save=True):
        """ Check if given points lie inside the surface """
        
        if points.ndim==1:
            points = np.array([points]) 
        
        assert points.shape[1] == 3, "input point must be x,y,z"

        n = points.shape[0]
        
        ## set up points to check
        vpoints = vtk.vtkPoints()
        map(vpoints.InsertNextPoint,points)
        
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
        
            self._pointchecked = True
        
        if n==1:
            return inout[0]
        else:
            return np.array(inout)
    
    def visualize(self, eig=False, axwidg=False):
        """ 
            Use VTK to visualize the geometry and 
            points if points have been checked 
        """ 
        
        ## shape mapper, actor
        shapeMapper = vtk.vtkPolyDataMapper()
        shapeMapper.SetInput(self.Surf)

        shapeActor = vtk.vtkActor()
        shapeActor.SetMapper(shapeMapper)
        shapeActor.GetProperty().SetOpacity(0.5)
    
        ## point mappers, actors
        useinpoints = 0
        useoutpoints = 0
        
        if self._pointchecked and self._inPoints.GetNumberOfPoints() > 0:
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

        if self._pointchecked and self._outPoints.GetNumberOfPoints() > 0:
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
            
        ## eigenvector mappers, actors
        if eig:
            ## major axis
            cylSource1 = vtk.vtkCylinderSource()
            cylSource1.SetCenter(0,0,0)
            h1 = self._eigval[0] / max(self._eigval)
            cylSource1.SetHeight(h1)
            cylSource1.SetRadius(.05)
            cylMapper1 = vtk.vtkPolyDataMapper()
            cylMapper1.SetInput(cylSource1.GetOutput())
            
            eig1Actor = vtk.vtkActor()
            eig1Actor.SetMapper(cylMapper1)
            eig1Actor.GetProperty().SetColor(1,0,0)
            eig1Actor.GetProperty().SetOpacity(.5)
            
            ## find rotation axis and angle
            rax1 = np.cross((0,1,0), self._eigvec[0])
            rang1 = np.arccos(self._eigvec[0][1]) * 180 / np.pi
            
            tfx1 = vtk.vtkTransform()
            tfx1.PostMultiply()
            tfx1.RotateWXYZ(rang1, rax1)
            # tfx1.Translate(self._eigvec[0] * h1/2)
            tfx1.Translate(self.CoM)
            
            eig1Actor.SetUserTransform(tfx1)
            
            ## middle axis
            cylSource2 = vtk.vtkCylinderSource()
            cylSource2.SetCenter(0,0,0)
            h2 = self._eigval[1] / max(self._eigval)
            cylSource2.SetHeight(h2)
            cylSource2.SetRadius(.05)
            cylMapper2 = vtk.vtkPolyDataMapper()
            cylMapper2.SetInput(cylSource2.GetOutput())
            
            eig2Actor = vtk.vtkActor()
            eig2Actor.SetMapper(cylMapper2)
            eig2Actor.GetProperty().SetColor(0,1,0)
            eig2Actor.GetProperty().SetOpacity(.5)

            
            ## find rotation axis and angle
            rax2 = np.cross((0,1,0), self._eigvec[1])
            rang2 = np.arccos(self._eigvec[1][1]) * 180 / np.pi
            
            tfx2 = vtk.vtkTransform()
            tfx2.PostMultiply()
            tfx2.RotateWXYZ(rang2, rax2)
            # tfx2.Translate(self._eigvec[1] * h2/2)
            tfx2.Translate(self.CoM)

            
            eig2Actor.SetUserTransform(tfx2)
            
            ## minor axis
            cylSource3 = vtk.vtkCylinderSource()
            cylSource3.SetCenter(0,0,0)
            h3 = self._eigval[2] / max(self._eigval)
            cylSource3.SetHeight(h3)
            cylSource3.SetRadius(.05)
            cylMapper3 = vtk.vtkPolyDataMapper()
            cylMapper3.SetInput(cylSource3.GetOutput())
            
            eig3Actor = vtk.vtkActor()
            eig3Actor.SetMapper(cylMapper3)
            eig3Actor.GetProperty().SetColor(0,0,1)
            eig3Actor.GetProperty().SetOpacity(.5)
            
            rax3 = np.cross((0,1,0),self._eigvec[2])
            rang3 = np.arccos(self._eigvec[2][1]) * 180 / np.pi
            
            tfx3 = vtk.vtkTransform()
            tfx3.PostMultiply()
            tfx3.RotateWXYZ(rang3, rax3)
            # tfx3.Translate(self._eigvec[2] * h3/2)
            tfx3.Translate(self.CoM)
            
            eig3Actor.SetUserTransform(tfx3)
            
            
        ## add axes widget
        if axwidg:
            axes = vtk.vtkAxesActor()
            widget = vtk.vtkOrientationMarkerWidget()
            widget.SetOutlineColor(0.9300, 0.5700, 0.1300)
            widget.SetOrientationMarker(axes)


        ## renderer, render window, interactor
        renderer = vtk.vtkRenderer()
        renderWindow = vtk.vtkRenderWindow()
        interactor = vtk.vtkRenderWindowInteractor()
    
        renderWindow.AddRenderer(renderer)
        interactor.SetRenderWindow(renderWindow)

        if axwidg:
            widget.SetInteractor(interactor)
            widget.SetViewport(0.0, 0.0, 0.4, 0.4)
            widget.SetEnabled(1)
            widget.InteractiveOn()
    
        ## add actors
        renderer.AddActor(shapeActor)
        if useinpoints:
            renderer.AddActor(inpointsActor)
        if useoutpoints:
            renderer.AddActor(outpointsActor)
        if eig:
            renderer.AddActor(eig1Actor)
            renderer.AddActor(eig2Actor)
            renderer.AddActor(eig3Actor)
            
        renderer.SetBackground(0,0,0)
    
        renderWindow.SetWindowName("Shape3DPoints")
        renderWindow.Render()
        interactor.Start();
       