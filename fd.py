# -*- coding: utf-8 -*-
"""
fd.py

(c) Timothy Drysdale <timothy.d.drysdale@gmail.com>


Finite difference demonstration
Written for readability, not efficiency


Calculates the potentials in a four-sided box, on a unitary aspect-ratio
rectangular grid.

The potentials on each wall are set to be constant, but need not be constant.

Aside from the challenge of numerical performance, one of the biggest 
time-soaks in making a modelling tool is handling the geometry input and mesh 
generation. It's no surprise that competing commercial tools use the same 
geometry library.

For simplicity - and to encourage you to get stuck in yourself - there are no 
such fancy graphical input tools with this code, but some convenience routines
have been provided to get you started and show the way to automating other
shaped fixed potentials.

Potentials can be specified within the domain as well, by adding an optional 
mask that is OR'd with the default update mask - although that particular 
implentation detail is hidden from the user via the fixV() method. Here's an
example doUpdate array for a 3x5 region with potentials defined on all 
boundaries.
FFFFF
FTTTF
FFFFF

Gotchas
Trying to fix the potential in the corners throws an error because it 
won't affect anything (due to the way the rectange grid and cross-shaped stamp
work - the stamp being the name to the update rule that we apply - the name
coming from the idea of using a rubber stamp to highlight on a printed array
which elements are used to update the element at the centre of the stamp) - 
here, X is updated using, N,S,E,W but not NW/SW/NE/SE.
  N
W X E
  S
  
We don't define the size of the grid because it is not incuded in this static
solution. This is an example of scale-invariance, although it only holds under 
the assumption that we are zoomed out enough we can't see the bumps of the
individual atoms in our materials (i.e. macroscopic).


"""
import numpy as np
import matplotlib.pyplot as plt

class Grid:
    
    def update(self):
        """ this is the update stamp ... not much to it! """
        for x in np.arange(self.West, self.East, self.goEast):
            if x not in [self.West, self.East]:
                for y in np.arange(self.South, self.North, self.goNorth):
                    if y not in [self.South, self.North]:
                        if not self.doUpdate[x,y]:
                            continue #skip fixed potentials
                        oldV = self.V[x,y]
                        newV = np.mean([
                                self.V[x + self.goWest, y],
                                self.V[x + self.goEast, y],
                                self.V[x, y + self.goNorth],
                                self.V[x, y + self.goSouth]
                                ])
                        deltaV = np.abs(newV - oldV)
                        self.maxDeltaV = np.max([self.maxDeltaV, deltaV])
                        self.V[x,y] = newV
                        
    def solve(self, accuracyV = 0.1, maxIterations = 9):
        """ here we decide how many times to run the update """
        self.maxDeltaV = 0.
        self.ErrNotConverged = False
        self.update()
        iterationCount = 1
        while (self.maxDeltaV > accuracyV) and (iterationCount < maxIterations):
            self.maxDeltaV = 0
            self.update()
            iterationCount = iterationCount + 1
        self.ErrNotConverged = self.maxDeltaV > accuracyV    
        self.iterationCount = iterationCount
       
  
    def __init__(self, Nx, Ny, southIsZero = True):
        """
        The grid needs setting up for each solution. For ease of writing and
        reading the rest of the code, we define indices and directions in
        terms of (n)orth, (s)outh, (e)ast (w)est.
        """
        if np.min([Nx,Ny]) < 3:
            raise RuntimeError("Grid must be 3x3 or larger")
        
        self.Nx = Nx
        self.Ny = Ny
                
        if southIsZero:
            self.goNorth = 1  #add this to index to go North by one grid
            self.North = self.Ny - 1
            self.South = 0
        else:
            self.goNorth = -1
            self.North = 0
            self.South = self.Ny - 1 
            
        self.goEast = 1 #add this to an index to go East by one grid
        self.West = 0
        self.East = self.Ny - 1
        self.goSouth = self.goNorth * -1
        self.goWest = self.goEast * -1
        
        self.NextToNorth = self.North + (1 * self.goSouth)
        self.NextToSouth = self.South + (1 * self.goNorth)
        self.NextToEast  = self.East  + (1 * self.goWest)
        self.NextToWest  = self.West  + (1 * self.goEast)
        
        """
        these arrays are indexed in the order of Easting, Northing to be 
        consistent with our conventional expection for the x, y dimensions;
        note that this violates our spoken convention of using the northing
        first - impossible to satisfy both conventions at once, unless you 
        violate another (and rotate the field so that north points right...!)
        """
        
        self.doUpdate = np.ones((Nx,Ny),dtype=bool)
        self.V = np.zeros((Nx,Ny),dtype=float)
        self.maxDeltaV = 0.
       
        """
        fix zero potentials on walls by default -
        we need to define boundary conditions on all walls no matter what,
        and a fixed potential is the only boundary condition implemented so far
        """
        self.fixWall('n',0)
        self.fixWall('s',0)
        self.fixWall('e',0)
        self.fixWall('w',0)
        
   
    def isCorner(self,x,y):
        """ True if x,y lie in one of the four grid corners"""
        return (x in [self.West, self.East]) and (y in [self.South, self.North])
    
       
    def fixV(self, x, y, V):
        """
        Fix a potential - except in a corner!
        """
        if self.isCorner(x,y):
            raise RuntimeError("Can't fix potential in the corner")
            
        self.doUpdate[x,y] = False
        self.V[x,y] = V
        
    def fixPoly(self, poly, V):
        for x in np.arange(self.West, self.East+1, self.goEast):
            for y in np.arange(self.South, self.North+1, self.goNorth):
                if poly.inside(x,y):
                    self.V[x,y] = V
                    self.doUpdate[x,y] = False
                    
    def floatPoly(self, poly):
        for x in np.arange(self.West, self.East, self.goEast):
            for y in np.arange(self.South, self.North, self.goNorth):
                if poly.inside(x,y):
                    self.doUpdate[x,y] = True
                    
    def fixWall(self,wall,V):
        """ fix a whole wall's potential at once 
        
        walls are (n/N)orth,(s/S)outh,(e/E)ast,(w/W)est
        
        """
        wallFound = False
        
        if wall[0].lower() == "n":
            wallFound = True
            for x in np.arange(self.West,self.East):
                self.V[x, self.North] = V
                self.doUpdate[x, self.North] = False
        
        if wall[0].lower() == "s":
            wallFound = True
            for x in np.arange(self.West,self.East):
                self.V[x, self.South] = V
                self.doUpdate[x, self.South] = False
                
        
        if wall[0].lower() == "w":
            wallFound = True
            for y in np.arange(self.South, self.North, self.goNorth):
                self.V[self.West, y] = V            
                self.doUpdate[self.West, y] = False

        if wall[0].lower() == "e":
            wallFound = True
            for y in np.arange(self.South, self.North, self.goNorth):
                self.V[self.East, y] = V  
                self.doUpdate[self.East, y] = False

        self.updateCorners()
        
        if not wallFound:
            raise RuntimeError("Wall not found, valid wall names are \
                                   N,S,E,W (not case sensitive)")
        
        
        
    def floatV(self, x, y):
        """ 
            for overwriting a fixed potential so it can float again
            may help in building some complex geometries
        """
        self.doUpdate[x,y] = True
    
    def updateCorners(self):
        """
        Treat corners as floating because
        a voltage fixed in a corner will not be noticed by the 
        update routines; Probably just need to do this when setting 
        the walls - if user is doing anything complicated with boundary 
        potentials, then they can call it themselves
        """
        self.V[self.West, self.North]  = np.mean(
                [self.V[self.NextToWest, self.North],
                 self.V[self.West, self.NextToNorth]])
        self.V[self.West, self.South]  = np.mean(
                [self.V[self.NextToWest, self.South],
                 self.V[self.West, self.NextToSouth]])    
        self.V[self.East, self.North]  = np.mean(
                [self.V[self.NextToEast, self.North],
                 self.V[self.East, self.NextToNorth]])
        self.V[self.East, self.South]  = np.mean(
                [self.V[self.NextToEast, self.South],
                 self.V[self.East, self.NextToSouth]])

"""

local oddNodes = false
local j = #polygon
for i = 1, #polygon do
    if (polygon[i].y < point.y and polygon[j].y >= point.y or polygon[j].y < point.y and polygon[i].y >= point.y) then
        if (polygon[i].x + ( point.y - polygon[i].y ) / (polygon[j].y - polygon[i].y) * (polygon[j].x - polygon[i].x) < point.x) then
            oddNodes = not oddNodes;
        end
    end
    j = i;
end
return oddNodes end

"""
class Poly:
    def __init__(self,points):
        self.points = points
        
    def inside(self,x,y):
        oddNodes= False
        j = len(self.points)-1

        for i in range(0,j+1):
            pix = self.points[i][0]
            piy = self.points[i][1]
            pjx = self.points[j][0]
            pjy = self.points[j][1]
            if ((piy < y and pjy >= y) or
                (pjy < y and piy >= y)):
                if ( pix + (y - piy)/(pjy-piy) * (pjx-pix)) < x:
                    oddNodes = not oddNodes
                    
            j=i
            
        return oddNodes    
            
            
def testPoly():
    
    p = Poly([[0,0],[2,0],[1,2]])
    if p.inside(1,1) == False:
        raise RuntimeError("error: point labelled outside poly when inside")
    if p.inside(0,3) == True:
        raise RuntimeError("error: point labelled inside poly when outside")
    
  
    
    
    
""" Some non-class functions to test  / demo the Grid class """
        
def testGrid():
    """
    Test routine to make sure that grid init and potential
    fixing/floating are working ok
    """
    #check that a new grid has its boundaries fixed at 0V 
    g = Grid(5,5)
    
    expectedV = np.array([
 [0., 0., 0., 0., 0.],
 [0., 0., 0., 0., 0.],
 [0., 0., 0., 0., 0.],
 [0., 0., 0., 0., 0.],
 [0., 0., 0., 0., 0.]])
   
    expectedDoUpdate = np.array([
 [False, False, False, False, False],
 [False,  True,  True, True, False],
 [False,  True, True,  True, False],
 [False,  True,  True,  True, False],
 [False, False, False, False,  True]])
    
    if not g.V.all() == expectedV.all():
        raise RuntimeError("Failed grid init test (non-zero potential found)")
    if not g.doUpdate.all() == expectedDoUpdate.all():
        raise RuntimeError("Failed grid init test (doUpdate incorrect)") 
  
    #fix a potential in the centre, and set non-zero boundary potentials
    g.fixV(2,2,10)
    g.fixWall('n',5)
    g.fixWall('s',10)
    g.fixWall('e',15)
    g.fixWall('w',20)
    
    expectedV = np.array(
 [[15.,  20., 20.,  20.,  12.5],
 [10.,   0.,   0.,   0.,   5. ],
 [10.,   0.,  10.,   0.,   5., ],
 [10.,   0.,   0.,   0.,   5., ],
 [12.5, 15.,  15.,  15.,  10. ]])
   
    expectedDoUpdate = np.array(
 [[False, False, False, False, False],
 [False,  True,  True, True, False],
 [False,  True, False,  True, False],
 [False,  True,  True,  True, False],
 [False, False, False, False,  True]])
    
    if not g.V.all() == expectedV.all():
        raise RuntimeError("Failed potential fixing test")
    if not g.doUpdate.all() == expectedDoUpdate.all():
        raise RuntimeError("Failed potential fixing test")        
 
    
    # test that a fixed potential can be floated
    g.floatV(2,2)   
   
    expectedDoUpdate = np.array(
 [[False, False, False, False, False],
 [False,  True,  True, True, False],
 [False,  True, True,  True, False],
 [False,  True,  True,  True, False],
 [False, False, False, False,  True]])
    
    if not g.doUpdate.all() == expectedDoUpdate.all():
        raise RuntimeError("Failed potential floating test")  
        
    """
    check that the solver is converging 
    """

    g = Grid(25,25)
    for x in np.arange(2,12):
        g.fixV(x,16,10)
        g.fixV(x,8,-10)

    g.solve(accuracyV = 0.1, maxIterations = 20)
    if g.ErrNotConverged:
        raise RuntimeError("Solver did not converge soon enough") 

def demoGrid():
    plt.figure()
    g = Grid(25,25)
    for x in np.arange(2,23):
        g.fixV(x,16,10)
        g.fixV(x,8,-10)
    
    g.solve(accuracyV=0.01, maxIterations=99)
    plt.contourf(g.V,10)
    plt.colorbar(label='Potential [Volts]')
    plt.xlabel('x [unit lengths]')
    plt.ylabel('y [unit lengths]')
    converged = " (Converged)"
    if g.ErrNotConverged:
        converged = " (Not Converged)"
    errV = "[error=%0.3fV]"%g.maxDeltaV
    plt.title('Solution after %d iterations%s%s'%(g.iterationCount,converged,errV))
    plt.savefig('demo.png',dpi=300)
    return g


def showEvolution():
    """ 
    produce some images for my lecture notes, showing the field in a finite
    capacitor evolving with each time step.
    Note - inefficient way to do a step-by-step export!
    """
    for maxIt in np.arange(1,52): 
        g = Grid(25,25)
        for x in np.arange(2,23):
            g.fixV(x,16,10)
            g.fixV(x,8,-10)
        g.solve(accuracyV=0.01, maxIterations=maxIt)
        plt.figure()
        plt.contourf(g.V,10)
        plt.colorbar(label='Potential [Volts]')
        plt.xlabel('x [unit lengths]')
        plt.ylabel('y [unit lengths]')
        converged = " (Converged)"
        if g.ErrNotConverged:
            converged = " (Not Converged)"
        errV = "[error=%0.3fV]"%g.maxDeltaV
        plt.title('Solution after %d iterations%s%s'%(g.iterationCount,converged,errV))
        plt.savefig('evo_step_%d.png'%g.iterationCount,dpi=300)
        print("%d/%d"%(g.iterationCount,maxIt))
        plt.close()

def demoPoly():
    
    plt.figure()
    g = Grid(50,50)
    
    p = Poly([[10,10],[30,10],[20,30]])

    g.fixPoly(p,50)    
    
    q = Poly([[35,35],[40,35],[40,40],[35,40]])
    
    g.fixPoly(q,-50)
    
    g.solve(accuracyV=0.01, maxIterations=500)
    plt.contourf(g.V,10)
    plt.colorbar(label='Potential [Volts]')
    plt.xlabel('x [unit lengths]')
    plt.ylabel('y [unit lengths]')
    converged = " (Converged)"
    if g.ErrNotConverged:
        converged = " (Not Converged)"
    errV = "[error=%0.3fV]"%g.maxDeltaV
    plt.title('Solution after %d iterations%s%s'%(g.iterationCount,converged,errV))
    plt.savefig('demoPoly.png',dpi=300)
    return g


def demoEM3():
    
    plt.figure()
    g = Grid(100,100)
    
    x0 = 10
    x1 = 15
    x2 = 25
    x3 = 35
    y0 = 10
    y1 = 15
    y2 = 25
    y3 = 30
    y4 = 40
    y5 = 45
    p = Poly([
            [x0,y0],
            [x3,y0],
            [x3,y1],
            [x1,y1],
            [x1,y2],
            [x2,y2],
            [x2,y3],
            [x1,y3],
            [x1,y4],
            [x2,y4],
            [x3,y4],
            [x3,y5],
            [x0,y5],
            [x0,y0]])
    
    g.fixPoly(p,-50)
    
    q = Poly([[x1 + 2, y1 + 4],[x3, y1 + 4],[x3, y1+6], [x1 + 2, y1 + 6], [x1 + 2, y1 + 4]])
    g.fixPoly(q,0)
    
    r = Poly([[x1 + 2, y3 + 4],[x3, y3 + 4],[x3, y3+6], [x1 + 2, y3 + 6]])
    g.fixPoly(r,0)  
 
    """            4
    *       *  3
    * *   * *  2
    *   *   *  1
    *       * 
    *       * 0
  x0 x1 x2 x3 x4
    """    
    x0 = 45
    x1 = 50
    x2 = 62
    x3 = 75
    x4 = 80

    y0 = 10
    y1 = 30
    y2 = 25
    y3 = 35
    y4 = 45
    
    ml = Poly([
            [x0,y0],
            [x1,y0],
            [x1,y4],
            [x0,y4],
            [x0,y0]])
    
    g.fixPoly(ml,50)
    
    mr = Poly([
            [x3,y0],
            [x4,y0],
            [x4,y4],
            [x3,y4],
            [x3,y0]])
    
    g.fixPoly(mr,50)    
    
    mv = Poly([
            [x1,y3],
            [x2,y2],
            [x3,y3],
            [x3,y4],
            [x2,y3],
            [x1,y4],
            [x1,y3]
            ])
    g.fixPoly(mv,50)
    
    mlm = Poly([
            [x0+7,y0],
            [x0+8,y0],
            [x0+8,y3-2],
            [x0+7,y3-2],
            [x0+7,y0]])
    
    g.fixPoly(mlm,0)

    mrm = Poly([
            [x4-7,y0],
            [x4-8,y0],
            [x4-8,y3-2],
            [x4-7,y3-2],
            [x4-7,y0]])
    
    g.fixPoly(mrm,0)
    
    g.solve(accuracyV=0.1, maxIterations=200)
    plt.contourf(np.flipud(np.rot90(g.V)),10)
    plt.colorbar(label='Potential [Volts]')
    plt.xlabel('x [unit lengths]')
    plt.ylabel('y [unit lengths]')
    converged = " (Converged)"
    if g.ErrNotConverged:
        converged = " (Not Converged)"
    errV = "[error=%0.3fV]"%g.maxDeltaV
    plt.title('Solution after %d iterations%s%s'%(g.iterationCount,converged,errV))
    plt.savefig('demoEM3.png',dpi=300)
    

    
    return g

if __name__ == "__main__":
    
    doTestPoly = False
    doTest = False
    
    doDemoPoly = False
    doDemoEM3 = False

    doDemo = True
    doShowEvolution = False #Slow-ish!
 
   
    if doDemoEM3:
        demoEM3()
        
    if doTestPoly:
        testPoly()

    if doDemoPoly:
        demoPoly()
        
    if doTest:
        testGrid()
        
    if doDemo:
        g = demoGrid() #returns g so you can explore the results yourself
        
    if doShowEvolution:
        showEvolution()