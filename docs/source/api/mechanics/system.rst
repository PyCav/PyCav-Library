System
======

System is a class used for the physical simulation of a collection of Particles, Springs, and Containers. It is also responsible for the visualization of such a collection. It does this by translating the properties of Particles, Springs, and Containers, to properties that vpython can use.

To simulate any new types of physical objects, subclass System and extend the simulate function. To make these new objects visualize in vpython, extend the create_vis and update_vis functions too. 

To visualize the system using some different system than vpython, subclass System and override the create_vis and update_vis functions.

Functions
-----------

__init__(collides, interacts, visualize, particles=None, springs=None, container=None, visualizer_type="vpython", canvas=None, stop_on_cycle=False, record_amplitudes=False, display_forces=False, record_pressure=False)
^^^^^^^^^^^^^^^^^
	
	Initialises a System class

	**Parameters**

	*collides: boolean*

	Whether particles in this system collide with each other or not

	*interacts: boolean*

	Whether particles in this system interact with each other via fields

	*visualize: boolean*

	Whether the things in this system are visualized

	*particles: array of Particles*

	An array of all the particles in this system

	*springs: array of springs*

	An array of all the springs in this system

	*container: Container*

	Container within which the simulation takes place

	*visualizer_type: string*

	What type of visualizer is used. Used to make sure that some vpython-specific non-essential things don't run if different visualization method is used in a subclass. Set to anything you like if using another visualization method

	*canvas: vpython cavas*

	The canvas within which the system is visualized. Can be any type of view if System is subclassed

	*stop_on_cycle: boolean*

	Whether the simulation stops runnign when a full cycle is done. Only tested on relatively simple 1D systems, and only works when using the function run_for instead of simulate

	*record_amplitudes: boolean*

	Whether the system records the amplitudes of any oscillations. Only tested with simple 1D systems

	*display_forces: boolean*

	Whether the forces applied on the particle are displayed. By default, these vectors are displaced from the particles by 2*particle radius along the 2-axis. To change, need to subclass and edit simulate method.

	*record_pressure: boolean*

	Whether to record the pressure on the walls of the container or not

create_vis(canvas=None)
^^^^^^^^^^^

	Creates a visualization. Override this method to change the visualization method.

	**Parameters:**

	*canvas: vpython canvas*

	The canvas into which the visualization is drawn. For this implementation, must be a vpython canvas

update_vis()
^^^^^^^^^^^
	Updates the visualization. Override this method to change the visualization method.

run_for(time, dt=0.01, on_step=None)
^^^^^^^^^^^
	Run simulation for a certain amount of time(as measured in the simulated system's time). Recommended to use this instead of simulate(dt) for most situations, unless need some mechanism to stop simulation on some external condition. 

	**Parameters:**

	*time: float*

	Time for which the simulation will run for in the system's time

	*dt: float*

	Size of each step taken in time

	*on_step: function taking one unnamed argument of System*

	This system is passed to the function, and the defined function will be performed at the end of every step


simulate(dt = 0.01)
^^^^^^^^^^^
	Simulates a time-step with a step size of dt. Collision detection, etc. happen here, so when adding new classes to simulate, extend this to add logic to simulate them.

	**Parameters:**

	*dt: float*

	Size of time step taken

create_particles_in_container(number=0, speed=0, radius=0, inv_mass=1.)
^^^^^^^^^^^^^^^^^^
	Creates the given number of particles, with the given parameters, in random locations within the container. If the system has no container, this method will raise a RuntimeError.

	**Parameters:**

	*number: integer*

	The number of particles to create

	*speed: float*

	The speed of these particles

	*radius: float*

	The radius of these particles

	*inv_mass: float*

	The inverse mass of these particles

Properties
-----------

*particles: array of Particles*

An array of all the particles in this system

*springs: array of springs*

An array of all the springs in this system

*container: Container*

Container within which the simulation takes place

*collides: boolean*

Whether particles in this system collide with each other or not

*interacts: boolean*

Whether particles in this system interact with each other via fields

*visualize: boolean*

Whether the things in this system are visualized

*visualizer_type: string*

What type of visualizer is used. Used to make sure that some vpython-specific non-essential things don't run if different visualization method is used in a subclass. Set to anything you like if using another visualization method

*canvas: vpython cavas*

The canvas within which the system is visualized. Can be any type of view if System is subclassed

*stop_on_cycle: boolean*

Whether the simulation stops runnign when a full cycle is done. Only tested on relatively simple 1D systems, and only works when using the function run_for instead of simulate

*record_amplitudes: boolean*

Whether the system records the amplitudes of any oscillations. Only tested with simple 1D systems

*display_forces: boolean*

Whether the forces applied on the particle are displayed. By default, these vectors are displaced from the particles by 2*particle radius along the 2-axis. To change, need to subclass and edit simulate method.

*record_pressure: boolean*

Whether to record the pressure on the walls of the container or not

*speeds: Array of floats, read only*

3D speed distribution of system as an unsorted array

*one_d_velocities: Array of floats, read only*

1D velocity distribution of system as an unsorted array




