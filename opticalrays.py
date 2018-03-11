# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 17:34:01 2018

@author: fsk16
"""

"""
Module to create and keep track of rays 
"""
import numpy as _np
import pylab as pl

def normalise(vector):
    norm = _np.linalg.norm(vector)
    if norm == 0: 
       return vector
    return vector / norm

norm = _np.linalg.norm

def dot_prod(v1,v2):
    return sum( v1 * v2 )

def index(array,value):  #to find the index of a certain element in an array
    x = 0
    for i in array:
#        if round(i,2) == value:
#            return x
        if round(i,1) == value:
            return x
        x += 1
        
def SnellsLaw(incident , normal , n_1 , n_2):
    normal = normalise(normal)
    incident = normalise(incident)
    n = n_1 / n_2
    c = -dot_prod(incident,normal)
    if 1 - n**2 *(1 - c**2) < 0:    # reflection etc.
        return None
    return n * incident  +  ( n*c - _np.sqrt(1 - n**2 *(1 - c**2))  )  * normal

class Ray:
    """
    Ray object stores the path and direction of a ray. Rays move at 0.1 u/s and can be changed
    by Ray.dr = new_velocity.
    
    """
    dr = 0.1
    def __init__(self, point = [0,0,0], direction = [0,0,0]):
        self.__position = _np.array([point[0],point[1],point[2]])
        self.__direction = normalise(_np.array([direction[0],direction[1],direction[2]]))
        self.__P = [self.__position]
        self.__D = [self.__direction]
        
    def pos(self):
        return self.__position
    def direction(self):
        return self.__direction

    def append(self, pos, direction):
        """
        Add new position and direction to ray, ray continues to move from this point
        """
        self.__position = pos # _np.array([pos[0],pos[1],pos[2]])
        self.__direction = direction # _np.array([direction[0],direction[1],direction[2]])
        self.__P.append(pos)
        self.__D.append(direction)

    def vertices(self):
        x = []
        y = []
        z = []
        for i in self.__P:
            x.append(i[0])
            y.append(i[1])
            z.append(i[2])
        A =_np.array([x,y,z])
        return A

    def  __repr__(self):
        return "Ray object at ( %g,%g,%g) " %(self.pos()[0],self.pos()[1],self.pos()[2])


class OpticalElement:
    """
    Represents the optical surfaces that refract/reflect the rays.
    """
    def propagate_ray(self, ray):
        """ 
        propagate a ray through the optical element assuming the 
        ray travels at +z direction with arbitrary directions in x and y 
        """
        dr = 0.1
        intercept = self.intercept(ray)
        if intercept != None:   # intercept exists
            while ray.pos()[2] < intercept[2]:
                ray.append(ray.pos() + ray.direction() * dr , ray.direction())      # this also updates position and direction of ray
            if self.refract == "positivecurvature":
                normal = intercept - self.focalpoint 
                new_direction = SnellsLaw(ray.direction(), normal ,self.n1 ,self.n2)
                ray.append(ray.pos() , new_direction)
            elif self.refract == "negativecurvature":
                normal = self.focalpoint - intercept
                new_direction = SnellsLaw(ray.direction(), normal ,self.n1 ,self.n2)
                ray.append(ray.pos() , new_direction)
            elif self.refract == "plane":
                normal = [0,0,-1]
                new_direction = SnellsLaw(ray.direction(), normal ,self.n1 ,self.n2)
                ray.append(ray.pos() , new_direction)
            return ray             #  "No interseption with surface, cannot propogate"   #Terminate the ray? IDK see worksheet
    
        
class SphericalRefraction(OpticalElement):
    """
    """
    def __init__(self, z_0 , curvature , n_1 , n_2 , aperture ):
        if curvature != 0:
            self.__radius = _np.abs(1/curvature)
            self.focalpoint = _np.array([0,0, z_0 + 1/curvature])
            self.refract = "curvature"
        else:
            self.__radius = 0
            self.focalpoint = _np.array([0,0,z_0])
            self.refract = "plane"
        self.__curvature = curvature
        self.__aperture = aperture
        self.n1 = n_1
        self.n2 = n_2

    def intercept(self,ray): 
        """
        From a given position and propogation direction of ray, finds the intersection with the surface of optical element.
        """
        k = ray.direction()  
        R = self.__radius
        r = ray.pos() - self.focalpoint
        
        if self.__curvature > 0:
            self.refract = "positivecurvature"
            if dot_prod(r , k )**2 - (norm(r)**2 - R**2) >= 0:
                length = -dot_prod(r , k) - _np.sqrt( dot_prod(r , k )**2 - (norm(r)**2 - R**2) )
            else:
                return None
        elif self.__curvature < 0:
            self.refract = "negativecurvature"
            if dot_prod(r , k )**2 - (norm(r)**2 - R**2) >= 0:
                length = -dot_prod(r , k) + _np.sqrt( dot_prod(r , k )**2 - (norm(r)**2 - R**2) )
            else:
                return None
        else:       #curv = 0 , plane surface
            length = norm(r[2]/k[2])
            
            
        intercept = ray.pos() + length * k
        
        if _np.abs(intercept[0]) > self.__aperture/2 or _np.abs(intercept[1]) > self.__aperture/2:  #If ray intersects but out of bounds due to aperture
            return None
        else:
            return intercept
        
class OutputPlane(OpticalElement):
    def __init__(self , z_f):
        self.__position = z_f
        self.pos = z_f
        self.refract = False
    def intercept(self,ray):
        k = ray.direction()
        length = norm((self.__position - ray.pos()[2] )  /  k[2])
        return length * k + ray.pos()
    
class RayBundle:
    """
    """
    def __init__(self, radius, origin, direction):
        self.__origin = origin
        self.__direction = direction
        self.rays = [Ray( self.__origin , self.__direction)]
        n=1
        while n <= 5:
            r = radius * n / 5
            for angle in  _np.arange( 1 , 6*n+1 , 1 ):
                x = r * _np.cos( 2*_np.pi/(6*n)  *angle)
                y = r * _np.sin( 2*_np.pi/(6*n)  *angle)
                z = self.__origin[2]
                ray_i = Ray([x , y , z] , self.__direction)
                self.rays.append(ray_i)
            n += 1

    def spotgraph(self, z):
        for ray in self.rays:
            ind = index(ray.vertices()[2],z)
            x = ray.vertices()[0]
            y = ray.vertices()[1]
            pl.scatter(x[ind],y[ind],s=10,c='b')
    def RMS(self, z):
        rsqr_sum = 0
        ind = index(self.rays[0].vertices()[2],z)
        xorigin_at_z = (self.rays[0].vertices()[0])[ind]
        y_at_z = (self.rays[0].vertices()[1])[ind]
        for ray in self.rays:
            ind = index(ray.vertices()[2],z)
            x = (ray.vertices()[0])[ind] - xorigin_at_z
            y = (ray.vertices()[1])[ind] - y_at_z
            rsqr_sum += x*x + y*y
        rms = _np.sqrt(rsqr_sum/91)
        return rms
            
def optimize(curvatures, focus , radius , n_1 , n_2):
    curvature2 = curvatures[1]
    curvature1 = curvatures[0]
    Obj1 = SphericalRefraction(100 , curvature1, n_1 , n_2 , 30)
    Obj2 = SphericalRefraction(105 , curvature2, n_2 , n_1 , 30)
    Output = OutputPlane(focus)
    bundle = RayBundle( radius, [0,0,0] , [0,0,1])
    for i in bundle.rays:
        Obj1.propagate_ray(i)
        Obj2.propagate_ray(i)
        Output.propagate_ray(i)
  #  print (bundle.RMS(focus))
    return bundle.RMS(focus)

    #Parameters of the model
if __name__ == "__main__":
    z_0 = 100
    z_1 = 105
    curv = 0.03
    n_1 = 1
    n_2 = 1.5
    aperture = 10
    z_output = 200
    radius_bundle = 1
    bundle_origin = [0,0,0]
    bundle_direction = [0,0,1]
    lambd = 400e-9
    #
    diffraction_scale = lambd * (z_output-z_0) / (radius_bundle*2)
    #
    ##    #Initialise objects
    convex = SphericalRefraction(z_0 , curv , n_1 , n_2 , aperture)
    plane = SphericalRefraction(z_1, 0 , n_2 , n_1 , aperture)
    M = SphericalRefraction(5 ,0, n_1 , n_2 , aperture)
    P = OutputPlane(z_output)
    bundle = RayBundle(radius_bundle,bundle_origin,bundle_direction)
    
#    print(optimise(0.03,0.02))
    #    #Propagate and plot the rays
    pl.figure(1)
    for i in bundle.rays:
        convex.propagate_ray(i,500)
        plane.propagate_ray(i,100)
        P.propagate_ray(i,500)
#        pl.scatter(i.vertices()[2],i.vertices()[0],s=0.0005,c='b')
#    #pl.axis([-1,220,-5,5])
    pl.title("x-z z0%.3f radius %.3f " %(z_0,radius_bundle))
    pl.xlabel("z(mm)")
    pl.ylabel("x(mm)")
    #pl.savefig("x-z z0%.3f radius%.3f .png" %(z_0,radius_bundle))
    print(bundle.RMS(z_output), diffraction_scale)
