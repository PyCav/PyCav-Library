Optics
====================================

Introduction to the Optics module
--------------------------

This module contains functions for use in geometric optics problems. The example usage is for caustics, such those observed at the bottom of a swimming pool. As we are working in the geometric optics limit, this simply requires manipulation of vectors. Refraction is a rotation about a axis perpendicular to the surface normal and incoming wave vector by an angle determined by Snell's law.

In-depth Documentation
-------

.. toctree::
   :maxdepth: 1
   
   ray_object
   caustics
   dispersion


The list below summarises the functions, their input arguments and their outputs for quick reference for the informed user:

Functions
---------

ray_grid(N_x,N_y,n_arr,h,f,dfdx,dfdy,x_lims = [0.,1.], y_lims = [0.,1.])
^^^^^^^^^^^^^^^^^^^^^^^^^^

   **Parameters:**

   *N_x: integer*
   
   **Returns:**

rays_refract(rays,t)
^^^^^^^^^^^^^^

ray_count(rays_x,rays_y,boxsize_x,boxsize_y)
^^^^^^^^^^^^^^^^^^^^^

evolve(rays,t,boxsize_x,boxsize_y)
^^^^^^^^^^^^^^^^^^^^^

single_time_image(rays,boxsize_x,boxsize_y)
^^^^^^^^^^^^^^^^^^^^^^

surface_calc(surface,x,y,t,N,h,f,disturbance_height,plot_height)
^^^^^^^^^^^^^^^^^^^^^^

caustic_image(x,y,N,XX,YY,II,h,f,disturbance_height,plot_height,c_map = 'Blues_r')
^^^^^^^^^^^^^^^^^^^^^^^^^

caustic_anim(x,y,t,N,XX_t,YY_t,II_t,h,f,disturbance_height,plot_height,c_map='Blues_r',interval = 100,fname = None)
^^^^^^^^^^^^^^^^^^^^^^^^

Ray Object
---------

Ray(x,y,f,dfdx,dfdy,h,n_arr) (__init__ function)
^^^^^^^^^^^^

findnormal(self,t)
^^^^^^^^^^^^^

snellsangle(self)
^^^^^^^^^^^^^

transformation_matrix(self)
^^^^^^^^^^^^^

refract(self,t)
^^^^^^^^^^^^