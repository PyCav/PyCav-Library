Particle
========

Particle is a class that represents a particle. It is not tied to any visualization method by design. By default, it implements a gravitational field, which can be changed by subclassing.

Functions
---------
__init__(pos=None, v=None, radius=1., inv_mass=0., color=None, alpha=1.,   fixed=False, applied_force=no_force, q = 0., make_trail = False)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  Initialises the Particle object.

  **Parameters:**

  *pos: numpy array*

  Initial position of particle, default 0, 0, 0

  *v: numpy array*

  Initial velocity of particle, default 0, 0, 0

  *radius: float*

  Radius of particle

  *inv_mass: float*

  Inverse mass of particle

  *color: array*

  Color of particle, given in form [R G B], default 1, 0, 0

  *alpha: float*

  Alpha of particle, 1 is completely opaque, 0 is completely transparent, used in visualisation


  *fixed: boolean*

  Whether particle can move or not

  *applied_force: function taking arguments of: particle(Particle) and time(float), that returns a numpy array*

  Gives the applied force based on the particle's properties and time. By default, no force applied on particle.

  *q: float*

  Charge on particle

  *make_trail: boolean*

  Whether the particle will make a trail or not

update(dt)
^^^^^^^^^^^^^^^^
  Updates the position of the particle using the velocity Verlet method.

  **Parameters:**

  *dt: float*

  Size of time step to take

force_on(other, if_at=None)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  Gives force which another particle will feel from this particle (if at means that this can also give the force that the other particle would feel if it were at some other position). Default implementation gives gravitational force. Subclass particle and override this function to implement custom forces.

  **Parameters:**

  *other: Particle* 

  The particle which feels the force

  *if_at: numpy array*

  If this parameter is used, the function gives the force the 'other' particle   would feel if it were at this position

  **Returns:**

  A *numpy array* with 3 elements giving the vector force on the other particle from this particle.

Properties
----------
pos
^^^
  *numpy array*

  3 element array giving position of particle in 3D space.

v
^
  *numpy array*

  3 element array giving the velocity of the particle.

fixed
^^^^^
  *boolean*

  If True, particle will be fixed and not move around however much force is applied to it. If False, particle will move due to interactions.

q
^^
  *float*

  Charge on particle.

color
^^^^^
  *array*

  Gives the color of the particle as an array, given in the form [R G B]. Each element of the array should be between 0 and 1.

radius
^^^^^^
  *float*

  Gives the radius of the particle.

alpha
^^^^^
  *float*

  Float between 0 and 1, giving the opacity of the particle.

make_trail
^^^^^^^^^^
  *boolean*

  Decides whether the particle will make a trail or not when visualized.

amplitude
^^^^^^^^^
  *float, read only*

  Gives the amplitude of oscillations. Depends on the system class the particle is in to update.

prev_pos
^^^^^^^^
  *numpy array, read only*

  3 element array giving the previous position of the particle.




