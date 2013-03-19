##############################
## 
## sharpy.shape3d
## 
## Author: A. Athanassiadis
## Date: March, 2013
## 
##############################
from __future__ import division
import os
import vtk
import numpy as np
from numpy import linalg as LA
from itertools import izip, imap

np.set_printoptions(precision=3)

## tolerance for point merging
TOL = 1e-4

## shorthand for readably accessing arrays
X = 0
Y = 1
Z = 2

## subexpressions for integrals
def subexpr(w):

    f = np.zeros(3)
    g = np.zeros(3)

    tmp0 = w[0] + w[1]
    tmp1 = w[0]**2
    tmp2 = tmp1 + w[1] * tmp0

    f[0] = w.sum()    
    f[1] = tmp2 + w[2] * f[0]
    f[2] = w[0] * tmp1 + w[1] * tmp2 + w[2] * f[1]
    
    g = f[1] + w * (f[0] + w)
    
    return f,g
    

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
    
    def __init__(self, fn, clean=True, mass=True, align=True):
        
        if not os.path.exists(fn):
            raise IOError("input file does not exists")
        
        self._filename = fn
        
        self._loaded = False
        self._centered = False
        self._aligned = False
        self._cleaned = False
        self._pointchecked = False
        self._pointcloudgen = False
        
        self.Surf = None
        self.nFaces = None
        self.PointCloud = None

        self.mass = 0
        self.CoM = np.zeros(3)
        self.inertia = np.zeros((3,3))
        
        self._eigval = None
        self._eigvec = None
        
        ext = np.char.lower(os.path.splitext(fn)[-1])
        
        if ext == ".stl":
            self._load_stl()
            
            if clean:
                self._clean_dup_points()
            
            if mass:
                print 'calculating mass properties'
                self._calc_mass_prop()
                # print self.inertia
                # print self._eigvec
                
            if align:
                ## works perfectly if you run it twice. what can I say?
                print 'aligning principal moments to x,y,z axes'
                self._align_axes()
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
        self.nFaces = self.Surf.GetNumberOfCells()
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
        
    def _calc_mass_prop(self):
        """
        calculate mass, center of mass, and moment of inertia of mesh
            
        """
        
        mult = np.array([1/6,1/24,1/24,1/24,1/60,1/60,1/60,1/120,1/120,1/120])
        intg = np.zeros(10) ## order: 1, x, y, z, x^2, y^2, z^2, xy, yz, zx
        
        for f in range(self.nFaces):
            ## get vertices
            p = map(self.Surf.GetCell(f).GetPoints().GetPoint, range(3))
            p = np.array(p)
            x,y,z = p.T
            fx = np.zeros(3)
            fy = np.zeros(3)
            fz = np.zeros(3)
            gx = np.zeros(3)
            gy = np.zeros(3)
            gz = np.zeros(3)
            
            ## get edges
            e1 = p[1] - p[0]
            e2 = p[2] - p[0]
            d = np.cross(e1,e2)
            
            ## compute integral terms
            fx,gx = subexpr(x)
            fy,gy = subexpr(y)
            fz,gz = subexpr(z)
            
            ## update integrals
            intg[0] += d[X] * fx[0]
            
            intg[1] += d[X] * fx[1]
            intg[2] += d[Y] * fy[1]
            intg[3] += d[Z] * fz[1]
            
            intg[4] += d[X] * fx[2]
            intg[5] += d[Y] * fy[2]
            intg[6] += d[Z] * fz[2]
            
            intg[7] += d[X] * np.dot(y,gx)
            intg[8] += d[Y] * np.dot(z,gy)
            intg[9] += d[Z] * np.dot(x,gz)      
            
        ## apply weights
        intg *= mult
        
        ## calc mass
        mass = intg[0]
        
        ## calc CoM
        cm = intg[1:4] / mass
        
        inertia = np.zeros((3,3))
        
        ## calc inertia wrt CoM
        inertia[X,X] = intg[5] + intg[6] - mass * (cm[Y]**2 + cm[Z]**2)
        inertia[Y,Y] = intg[4] + intg[6] - mass * (cm[Z]**2 + cm[X]**2)
        inertia[Z,Z] = intg[4] + intg[5] - mass * (cm[X]**2 + cm[Y]**2)
        inertia[X,Y] = -(intg[7] - mass * cm[X] * cm[Y])
        inertia[Y,Z] = -(intg[8] - mass * cm[Y] * cm[Z])
        inertia[X,Z] = -(intg[8] - mass * cm[Z] * cm[X])
        inertia[Y,X] = inertia[X,Y]
        inertia[Z,Y] = inertia[Y,Z]
        inertia[Z,X] = inertia[X,Z]
        
        self.mass = mass
        self.CoM = cm
        self.inertia = inertia   
        
        # calculate and sort eigenvectors of inertia tensor
        eigval, eigvec = LA.eigh(self.inertia)

        eigsys = zip(eigval, eigvec.T)
        eigsys.sort(key = lambda x: x[0])

        eigval, eigvec = zip(*eigsys)

        self._eigval = eigval[::-1]
        self._eigvec = np.array(eigvec)[::-1]
        
    def _center(self):
        """
        center mesh at (0,0,0)
        
        """
        
        txf = vtk.vtkTransform()
        txf.PostMultiply()
        txf.Translate(-1*self.CoM)
        
        txfPoly = vtk.vtkTransformPolyDataFilter()
        txfPoly.SetInput(self.Surf)
        txfPoly.SetTransform(txf)
        txfPoly.Update()
        
        self.Surf = txfPoly.GetOutput()
        self._bounds = self.Surf.GetBounds()
        
        self._calc_mass_prop()
        
        self._centered = True
                      
        
    def _align_axes(self):
        """ 
        Align the imported surface according to the following rules:
            - x axis: major axis of the moment of inertia tensor
            - y axis: middle axis
            - z axis: minor  axis
        
        """
        
        # make sure we're centered
        if not self._centered:
            self._center()
            
        ## align major principal axis
        rax1 = np.cross(self._eigvec[0],(1,0,0))
        rax1 /= LA.norm(rax1)
        rang1 = np.abs(np.arccos(self._eigvec[0][0]) * 180 / np.pi)
        
        c = self._eigvec[0][0]
        s = np.sin(np.arccos(c))
        
        R = np.eye(3) * c + \
            (1-c) * np.outer(rax1,rax1) + \
            s * np.array([ [0, -rax1[2], rax1[1]],
                           [rax1[2], 0, -rax1[0]],
                           [-rax1[1], rax1[0], 0]])
        
        e1rot = np.dot(R,self._eigvec[0])
        e2 = self._eigvec[1]
        e2rot = np.dot(R,e2)
        
        rax2 = np.cross(e2rot,(0,1,0))
        rax2 /= LA.norm(rax2)
        rang2 = np.abs(np.arccos(e2rot[1]) * 180/np.pi)
        
        txf = vtk.vtkTransform()
        txf.PostMultiply()
        txf.RotateWXYZ(rang1,rax1)
        txf.RotateWXYZ(rang2,rax2)
        
        txfPoly = vtk.vtkTransformPolyDataFilter()
        txfPoly.SetInput(self.Surf)
        txfPoly.SetTransform(txf)
        txfPoly.Update()
        
        self.Surf = txfPoly.GetOutput()
        self._bounds = self.Surf.GetBounds()
    
        print 'aligned surface. recalculating eigenvectors of inertia tensor'
        self._calc_mass_prop()
        
        print "updated eigenvectors"
        print self._eigvec
        ## print self.inertia
        
        self._aligned = True
        
    def _generate_pointcloud(self, N=1e5, save=False):
        """ Build an N-point pointcloud representation of the surface """

        x = np.random.uniform(low=self._bounds[0], high=self._bounds[1], size=N)
        y = np.random.uniform(low=self._bounds[2], high=self._bounds[3], size=N)
        z = np.random.uniform(low=self._bounds[4], high=self._bounds[5], size=N)
        
        pc = np.array([x,y,z]).T
        
        inbool = self.is_inside(pc,save=save)
        
        self.PointCloud = pc[inbool==1]
        print sum(inbool)
        
        self._pointcloudgen = True
                
    
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
            outpointsActor.GetProperty().SetOpacity(0.01)    
            
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
       