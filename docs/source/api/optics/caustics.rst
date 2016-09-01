Caustics
============================

.. image:: ../images/refraction.png

Argument list
^^^^^^^^^^^^^^^^

ray_grid(N_x,N_y,n_arr,h,f,dfdx,dfdy,x_lims = [0.,1.], y_lims = [0.,1.])

   **Parameters:**

   *N_x: integer*
   
   **Returns:**

rays_refract(rays,t)


ray_count(rays_x,rays_y,boxsize_x,boxsize_y)


evolve(rays,t,boxsize_x,boxsize_y)


single_time_image(rays,boxsize_x,boxsize_y)


surface_calc(surface,x,y,t,N,h,f,disturbance_height,plot_height)


caustic_image(x,y,N,XX,YY,II,h,f,disturbance_height,plot_height,c_map = 'Blues_r')


caustic_anim(x,y,t,N,XX_t,YY_t,II_t,h,f,disturbance_height,plot_height,c_map='Blues_r',interval = 100,fname = None)
