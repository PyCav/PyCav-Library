import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from mpl_toolkits.mplot3d import Axes3D

import pycav.display as display

def cross(a,b):
    return np.array((a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]))

class Ray():

    def __init__(self,x,y,f,dfdx,dfdy,h,n_arr):
        self.x = x
        self.y = y
        self.vec = np.array((0.,0.,-1.))
        self.f = f
        self.dfdx = dfdx
        self.dfdy = dfdy
        self.h = h
        self.n_arr = n_arr

    def findnormal(self,t):   
        self.norm  = np.array((-self.dfdx(self.x,self.y,t),-self.dfdy(self.x,self.y,t),1))
        self.norm /= np.linalg.norm(self.norm)

    def snellsangle(self):
        self.theta_1 = np.arccos(-self.vec.dot(self.norm))
        self.theta_2 = np.arcsin(self.n_arr[0]*np.sin(self.theta_1)/self.n_arr[1])
        
    def transformation_matrix(self):
        
        t_12 = self.theta_1-self.theta_2
        c = np.cos(t_12)
        s = np.sin(t_12)
        R_z = np.array(((c ,s ,0.),
                        (-s,c ,0.),
                        (0.,0.,1.)))
        
        axis  = cross(self.vec,self.norm)
    
        if any(axis != np.array((0.,0.,0.))):
    
            axis /= np.linalg.norm(axis)
    
            x_vec = cross(self.vec,axis)

            sin   = np.linalg.norm(x_vec)
            cos   = axis.dot(self.vec)
        
            if sin == 0.:
                a = 0.
            else:
                a = (1-cos)/sin**2

            skew_nz = np.array(((0.       ,-x_vec[2],x_vec[1] ),
                                (x_vec[2] ,0.       ,-x_vec[0]),
                                (-x_vec[1],x_vec[0] ,0.       )))

            R_nz = np.identity(3)+skew_nz+a*skew_nz.dot(skew_nz)

            self.R = np.linalg.inv(R_nz).dot(R_z.dot(R_nz))
        
        else:
            self.R = np.identity(3)
                          
    def refract(self,t):
        self.findnormal(t)
        self.snellsangle()
        self.transformation_matrix()
        
        self.r_x = self.x+(self.h+self.f(self.x,self.y,t))*self.R[2,0]/self.R[2,2]
        self.r_y = self.y+(self.h+self.f(self.x,self.y,t))*self.R[2,1]/self.R[2,2]

def ray_grid(N_x,N_y,n_arr,h,f,dfdx,dfdy,x_lims = [0.,1.], y_lims = [0.,1.]):

    x = np.linspace(x_lims[0],x_lims[1],N_x)
    y = np.linspace(y_lims[0],y_lims[1],N_y)

    X,Y = np.meshgrid(x,y)
    x_coord = X.ravel()
    y_coord = Y.ravel()

    return x,y,[Ray(x_coord[i],y_coord[i],f,dfdx,dfdy,h,n_arr) for i in range(N_x*N_y)]

def rays_refract(rays,t):
    for ray in rays:
        ray.refract(t)
    rays_x = np.asarray([ray.r_x for ray in rays])
    rays_y = np.asarray([ray.r_y for ray in rays])
    return rays_x,rays_y

def ray_count(rays_x,rays_y,boxsize_x,boxsize_y):    
    x_max = rays_x.max()
    x_min = rays_x.min()
    y_max = rays_y.max()
    y_min = rays_y.min()
    
    x_steps = int((x_max-x_min)/boxsize_x)
    y_steps = int((y_max-y_min)/boxsize_y) 
    
    I, x_edges, y_edges = np.histogram2d(rays_x, rays_y, bins=[x_steps, y_steps], 
                                         range=[[x_min, x_max], [y_min, y_max]])
    
    x = 0.5*(x_edges[:-1] + x_edges[1:])
    y = 0.5*(y_edges[:-1] + y_edges[1:])
    
    XX,YY = np.meshgrid(x,y)
                
    return XX,YY,I.T

def evolve(rays,t,boxsize_x,boxsize_y):
    len_t = len(t)
    
    XX_t = []
    YY_t = []
    I_t = []    

    for i in range(len_t):
        rays_x,rays_y = rays_refract(rays,t[i])
        XX,YY,I = ray_count(rays_x,rays_y,boxsize_x,boxsize_y)
        XX_t.append(XX)
        YY_t.append(YY)
        I_t.append(I)
        
    return XX_t,YY_t,I_t

def single_time_image(rays,boxsize_x,boxsize_y):
    rays_x,rays_y = rays_refract(rays,0.)
    X,Y,I = ray_count(rays_x,rays_y,boxsize_x,boxsize_y)
    
    return rays_x,rays_y,X,Y,I

######################
# Plotting functions #
######################

def surface_calc(surface,x,y,t,N,h,f,disturbance_height,plot_height):
    N_x = N[0]
    N_y = N[1]
    for i in range(N_x):
        for j in range(N_y):
            surface[j,i] = plot_height*h*f(x[i],y[j],t)/disturbance_height
    return surface

def caustic_image(x,y,N,XX,YY,II,h,f,disturbance_height,plot_height,c_map = 'Blues_r'):
    fig = plt.figure(figsize = (9,9))
    ax = fig.gca(projection='3d')
    
    X,Y = np.meshgrid(x,y)
    
    surface = np.zeros_like(X)
    max_I = II.max()
    
    surface = surface_calc(surface,x,y,0.,N,h,f,disturbance_height,plot_height)
            
    ax.plot_surface(X,Y,surface, color = 'w', alpha = 0.5)
    ax.contourf(XX, YY, II, 100, zdir='z', offset=-h, cmap=c_map,vmin=0., vmax=max_I)
    ax.set_zlim(-h,2*plot_height*h);
    
def caustic_anim(x,y,t,N,XX_t,YY_t,II_t,h,f,disturbance_height,plot_height,c_map='Blues_r',interval = 100,fname = None):
    fig = plt.figure(figsize = (9,9))
    ax = fig.gca(projection='3d')
    
    X,Y = np.meshgrid(x,y)    
    
    surface = np.zeros_like(X)

    max_I = max([I.max() for I in II_t])
    
    surface = surface_calc(surface,x,y,t[0],N,h,f,disturbance_height,plot_height)    
    ax.plot_surface(X,Y,surface, color = 'w', alpha = 0.5)
    ax.contourf(XX_t[0], YY_t[0], II_t[0], 100,
                zdir='z', offset=-h, cmap=c_map,vmin=0., vmax=max_I)
    ax.set_zlim(-h,0.5*h);
    
    def nextframe(arg):
        surface[:,:] = surface_calc(surface,x,y,t[arg],N,h,f,disturbance_height,plot_height)
        ax.clear()
        ax.plot_surface(X,Y,surface, color = 'w', alpha = 0.5)
        ax.contourf(XX_t[arg], YY_t[arg], II_t[arg], 100,
                zdir='z', offset=-h, cmap=c_map,vmin=0., vmax=max_I)
        ax.set_zlim(-h,0.5*h);
    
    num_frames = len(t)
    animate = anim.FuncAnimation(fig,nextframe,interval = interval, frames = num_frames, repeat = False)
    if fname == None:
        animate = display.create_animation(animate, temp = True)
    else:
        animate = display.create_animation(animate, fname = fname)

    return animate