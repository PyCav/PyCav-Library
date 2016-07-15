"""
An object oriented library to allow quick prototyping of simulations.
"""
from __future__ import division, print_function
import numpy as np
from random import *
# This is here so that vpython can be imported later on if not yet imported but needed,
# and this method is faster than importing every time it's needed.
vpython = None


def no_force(pos, time):
    """
    We can apply some time/position-dependent force to each particle.
    This is the default function, where there is no force applied on the particle.

    Parameters
    ----------
    pos : numpy array
        Position at which the force is felt
    time: float
        Time at which force is felt
    """
    return np.array([0., 0., 0.])


def element_mult(vec_1, vec_2):
    """
    Multiplies each element of vector vec_1 with corresponding element of vec_2

    Parameters
    ----------
    vec_1: numpy array
        First vector for which elements are multiplied
    vec_2: numpy array
        Second vector for which elements are multiplied
    """
    vec = []
    vec.append(vec_1[0] * vec_2[0])
    vec.append(vec_1[1] * vec_2[1])
    vec.append(vec_1[2] * vec_2[2])
    return np.array(vec)


def vector_from(arr=np.array([0, 0, 0])):
    """
    Creates a vpython vector from a numpy array

    Parameters
    ----------
    arr: numpy array
        Array that is converted into a vpython vector
    """
    global vpython
    if vpython is None:
        import vpython
    return vpython.vector(arr[0], arr[1], arr[2])


def normalized(arr):
    """
    Returns normalised numpy array.

    Parameters
    ----------
    arr: numpy array
        Array to be normalised
    """
    return arr / np.sqrt(np.inner(arr, arr))


def _duplicate_vector(arr):
    """
    Duplicates a vector and returns it. Faster than copy.deepcopy
    """
    return np.array([arr[0], arr[1], arr[2]])


def _perpendicular_vector(vec):
    """
    Gets an arbitary vector that is perpendicular to the vector given
    Parameters
    ----------
    vec : numpy array
        Vector whose perpendicular you want to find
    """
    if vec[0] == 0 and vec[1] == 0:
        return _cross(vec, np.array([0, 1, 0]))
    return _cross(vec, np.array([0, 0, 1]))


def _cross(vec1, vec2):
    """
    Method to cross two vectors, represented as arrays.
    Use instead of numpy cross method for slight performance gain.
    Parameters
    ----------
    vec1: numpy array
        First vector in cross product
    vec2: numpy array
        Second vector in cross product
    """
    x = ((vec1[1] * vec2[2]) - (vec1[2] * vec2[1]))
    y = ((vec1[2] * vec2[0]) - (vec1[0] * vec2[2]))
    z = ((vec1[0] * vec2[1]) - (vec1[1] * vec2[0]))
    return (np.array([x, y, z]))


class Particle(object):
    """
    Class representing a particle, which can be attached to a spring.
    Not tied to any visualisation method.
    """


    
    # Getters and setters for certain properties so that visualization is only updated if it has changed.

    @property
    def color(self):
        return self._color
    @color.setter
    def color(self, color):
        self._color = color
        self._vis_values_changed = True

    @property
    def radius(self):
        return self._radius
    @radius.setter
    def radius(self, radius):
        self._radius = radius
        self._vis_values_changed = True

    @property
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha
        self._vis_values_changed = True

    @property
    def make_trail(self):
        return self._make_trail
    @make_trail.setter
    def make_trail(self, make_trail):
        self._make_trail = make_trail
        self._vis_values_changed = True
    

    def __init__(self, pos=None, v=None, radius=1.,
        inv_mass=0., color=None, alpha=1., fixed=False, applied_force=no_force,
        q = 0., make_trail = False):
        """
        Parameters
        ----------
        pos: numpy array
            Initial position of particle, default 0,0,0
        v: numpy array
            Initial velocity of particle, default 0,0,0
        radius: float
            Radius of particle
        inv_mass: float
            Inverse mass of particle
        color: array
            Color of particle, given in form [R G B], default 1, 0, 0
        alpha: float
            Alpha of particle, 1 is completely opaque, 0 is completely transparent, used in visualisation
        fixed: boolean
            Whether particle can move or not
        applied_force: function taking arguments of: pos(numpy array) and time(float). Returns a numpy array
            Gives the applied force when the particle is at a certain location and time. By default, no force applied on particle.
        q: Float
            Charge on particle
        make_trail: boolean
            Whether the particle will make a trail or not
        """
        if pos is not None:
            self.pos = pos
        else:
            self.pos = np.array([0, 0, 0])
        if v is not None:
            self.v = v
        else:
            self.v = np.array([0, 0, 0])
        self.inv_mass = inv_mass
        self._radius = radius
        if color is not None:
            self._color = color
        else:
            self._color = [1, 0, 0]
        self._make_trail = make_trail
        self._alpha = alpha
        self._vis_values_changed = False
        self.fixed = fixed
        self.applied_force = applied_force
        self.total_force = np.array([0, 0, 0])
        self.max_point = np.array([None])
        self.min_point = np.array([None])
        self._visualized = False  # Used in visualisation
        self._pointer_assigned = False  # Used when looking at forces using pointers
        self.q = q
        self.initial_v = _duplicate_vector(self.v)
        self.initial_pos = _duplicate_vector(self.pos)
        self.prev_pos = _duplicate_vector(self.pos)
        self.prev_v = _duplicate_vector(self.v)

    def update(self, dt):
        """
        Updates the position of the particle.
        Parameters
        ----------
        dt: float
            Time step to take
        """
        self.v += (dt * self.total_force) * self.inv_mass
        self.pos += (self.v * dt) + ((0.5 * self.total_force * self.inv_mass) * (dt**2))

    def force_on(self, other, if_at=None):
        """
        Gives force which another particle will feel from this particle.
        Default implementation gives gravitational force.
        Subclass and change this class to implement custom forces.
        Parameters
        ----------
        other: Particle
            The particle which feels the force
        if_at: numpy array
            If this parameter is used, the function gives the force the 'other' particle would feel if it were at this position
        """
        if not if_at:
            if_at = other.pos
        position_difference = if_at - self.pos
        determinant = np.sqrt(np.inner(position_difference, position_difference))
        if determinant == 0:
            return np.array([0., 0., 0.])
        g_force_scalar = (-1) / (self.inv_mass * other.inv_mass * (determinant**3))
        g_force_vector = g_force_scalar * position_difference
        return g_force_vector

    @property
    def amplitude(self):
        """
        Property which gives the amplitude of oscillation. Depends on external thing to set max_point and min_point.
        """
        if self.max_point.any() and self.min_point.any():
            _amplitude = np.linalg.norm(self.max_point - self.min_point)
            return _amplitude
        return None


class Container(object):
    """
    Class describing a cubic box which particles can be in.
    """
    @property
    def dimension(self):
        """Property which stores the dimensions of the cube"""
        return self._dimension

    @dimension.setter
    def dimension(self, dimension):
        """Setter so that we know when dimensions have changed"""
        self._dimension = dimension
        self.size_changed = True


    # Getters and setters for certain properties so that visualization is only updated if it has changed.

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, color):
        self._color = color
        self._vis_values_changed = True

    @property
    def _origin_pos(self):
        """
        Property which gives position of the 'origin' of the cube,
        i.e. the position of the bottom left corner
        """
        l2 = self.dimension / 2
        return (self.pos - np.array([l2, l2, l2]))

    @property
    def surface_area(self):
        """Property which describes the total surface area of the container"""
        return 6 * (self.dimension)**2

    def __init__(self, dimension, pos=None, color=None, alpha=0.3):
        """
        Parameters
        ----------
        dimension: float
            The dimension of the cube, i.e. the length, width, height of the cube
        pos: numpy array
            Position of centre of cube, default 0,0,0
        color: array
            Color of particle, given in form [R G B], default 1,1,1
        alpha: float
            Alpha of particle, 1 is completely opaque, 0 is completely transparent, used in visualisation
        """
        if pos is not None:
            self.pos = pos
        else:
            self.pos = np.array([0, 0, 0])
        self._dimension = dimension
        self.size_changed = False
        if color is not None:
            self._color = color
        else:
            self._color = [1, 1, 1]
        self._alpha = alpha
        self._visualized = False
        self._vis_values_changed = False

    def contains(self, particle):
        """
        Returns True if is inside container.
        Returns the index of the axis along which it is outside the container if it is outside the container
        Parameters
        ----------
        particle: Particle
            Particle which is being checked to see if inside container or not
        """
        relative_pos = particle.pos - self._origin_pos
        for i in range(0, 3):
            if relative_pos[i] > self.dimension - particle.radius or relative_pos[i] < particle.radius:
                return i
        return True


class _Domain(Container):
    def __init__(self, pos, dimension):
        """
        Class used for collision detection, by splitting up containers into small domains.
        Parameters
        ----------
        pos: numpy array
            Position of 'origin' of domain, i.e. lower left corner of domain.
        dimension : float
            The length of each side of a domain.
        """
        l2 = dimension / 2
        Container.__init__(self, dimension, pos=pos + np.array([l2, l2, l2]), alpha=0.)
        self.particles = []

    def contains(self, particle):
        """
        Returns True if is inside container.
        Returns the index of the axis along which it is outside the container if it is outside the container
        Parameters
        ----------
        particle: Particle
            Particle which is being checked to see if inside container or not
        """
        relative_pos = particle.pos - self._origin_pos
        for i in range(0, 3):
            if relative_pos[i] > self.dimension or relative_pos[i] < 0:
                return i
        return True


class Spring(object):
    """
    Class representing a spring. Not tied to any visualisation method.
    """

    @property
    def axis(self):
        """
        Property giving the axis, i.e. the vector showing the orientation and length of the spring.
        """
        return (self.particle_2.pos - self.particle_1.pos)

    @property
    def pos(self):
        """
        Property giving the position of one of the ends of the spring.
        """
        return self.particle_1.pos

    # Getters and setters for certain properties so that visualization is only updated if it has changed.
    @property
    def radius(self):
        return self._radius
    @radius.setter
    def radius(self, radius):
        self._radius = radius
        self._vis_values_changed = True

    @property
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha
        self._vis_values_changed = True

    @property
    def color(self):
        return self._color
    @color.setter
    def color(self, color):
        self._color = color
        self._vis_values_changed = True


    def __init__(self, particle_1, particle_2, k, l0=None, radius=0.5, color=None, alpha=1.):
        """
        Parameters
        ----------
        particle_1: Particle object
            Particle on one end of spring
        particle_2: Particle object
            Particle on other end of spring
        k: float
            The spring constant of the spring (F = kx)
        l0: float
            Original length of the spring
        radius: float
            Radius of spring.
        color: array
            Color of particle, given in form [R G B]
        alpha: float
            Alpha of particle, 1 is completely opaque, 0 is completely transparent, used in visualisation
        """
        self.particle_1 = particle_1
        self.particle_2 = particle_2
        self.k = k
        self._alpha = alpha
        if color is not None:
            self._color = color
        else:
            self._color = [1, 1, 1]
        self._visualized = False
        self._radius = radius
        if l0:
            self.l0 = l0
        else:
            self.l0 = np.linalg.norm(self.particle_1.pos - self.particle_2.pos)
        self._vis_values_changed = False

    def force_on(self, particle, if_at=np.array([None])):
        """
        Given an arbitary particle, gives the force on that particle. No force if the spring isn't connected to that particle.
        Parameters
        ----------
        particle: Particle object
            Particle which feels the force
        if_at: NumPy Array
            If this parameter is used, this gives the force felt if the particle were at this position.
        """
        if not if_at.any():
            if_at = particle.pos
        if particle == self.particle_1:
            x = np.linalg.norm(if_at - self.particle_2.pos) - self.l0
            return normalized(self.particle_2.pos - if_at) * self.k * x
        elif particle == self.particle_2:
            x = np.linalg.norm(if_at - self.particle_1.pos) - self.l0
            return normalized(self.particle_1.pos - if_at) * self.k * x
        return np.array([0, 0, 0])


class PointerArrow(object):
    """
    Class representing an arrow/pointer, that isn't tied to any specific visualization method/library.
    """


    # Getters and setters for certain properties so that visualization is only updated if it has changed.

    @property
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha
        self._vis_values_changed = True

    @property
    def color(self):
        return self._color
    @color.setter
    def color(self, color):
        self._color = color
        self._vis_values_changed = True

    @property
    def shaftwidth(self):
        return self._shaftwidth
    @shaftwidth.setter
    def shaftwidth(self, shaftwidth):
        self._shaftwidth = shaftwidth
        self._vis_values_changed = True

    def __init__(self, pos, axis, shaftwidth=1, color=None, alpha=1):
        """
        Parameters
        ----------
        pos: numpy array
            Location of tail end of pointer
        axis: numpy array
            A vector giving the length and orientation of the pointer
        shaftwidth: float
            The width of the pointer's shaft
        color: array
            Color of particle, given in form [R G B], default 1, 1, 1
        alpha: float
            Alpha of particle, 1 is completely opaque, 0 is completely transparent, used in visualisation
        """
        self.pos = pos
        self.axis = axis
        self._shaftwidth = shaftwidth
        self._alpha = alpha
        if color is not None:
            self._color = color
        else:
            self._color = [1, 1, 1]
        self._vis_values_changed = False
        self._visualized = False


class System(object):
    """
    Class representing a collection of particles, springs, pointers, and a container(not yet implemented).
    Also can handle visualisation.
    If want to use anything other than vpython, subclass and change create_vis, update_vis, (necessary),
    and maybe run_for (If want a convenience function, but can just use simulate).
    Alternatively, could just put visualize = True and extract information from particles, etc. and visualize.
    """
    @property
    def one_d_velocities(self):
        """
        Property giving 1d velocity distribution of system as unsorted array.
        """
        self._one_d_velocities = np.zeros(len(self.particles))
        for index, particle in enumerate(self.particles):
            self._one_d_velocities[index] = particle.v[0]
        return self._one_d_velocities

    @property
    def speeds(self):
        """
        Property giving 3d speed distribution of system as unsorted array.
        """
        self._speeds = np.zeros(len(self.particles))
        for index, particle in enumerate(self.particles):
            self._speeds[index] = np.sqrt(np.inner(particle.v, particle.v))
        return self._speeds

    def __init__(self, collides, interacts, visualize,
        particles=None, springs=None, container=None,
        visualizer_type="vpython", canvas=None,
        stop_on_cycle=False, record_amplitudes=False, display_forces=False,
        record_pressure=False):
        """
        Parameters
        ----------
        particles: array of Particles
            An array of all the particles in this system.
        interacts: boolean
            Whether particles interact with each other via fields they emit.
        springs: array of Springs
            An array of all the springs in this system.
        container: Container
            The container containing the particles in this system.
        collides: boolean
            Determines whether particles will collide with each other or not.
        visualize: boolean
            Whether system visualizes itself.
        visualizer_type: string
            What type of visualizer is used. Used to make sure that some vpython-specific non-essential things don't run.
            Set to anything you like if using another visualization method.
        canvas: some view
            Some view to draw everything into. By default, a vpython canvas. Depends upon visualization method used.
        stop_on_cycle: boolean
            Whether the system stops running when a full cycle is done. Only tested on relatively simple 1-D systems, but seems to work.
        record_amplitudes: boolean
            Whether the system records the amplitudes of oscillations. Only works on simple 1-D systems going along the x-axis.
        display_forces: boolean
            Whether the forces applied onto a particle are displayed.
            By default, these vectors are displaced from the particles by 2*particle radius in y-direction.
            Look in simulate method to change this.
        record_pressure: boolean
            Whether to record pressure on walls or not
        """
        self.visualize = visualize
        self.interacts = interacts
        if particles is not None:
            self.particles = particles
        else:
            self.particles = []
        if springs is not None:
            self.springs = springs
        else:
            self.springs = []
        self.collides = collides
        self.pointerarrows = []  # Only used for displaying forces
        # visual representations of the more abstract classes particles, springs, and pointerarrows.
        # By default, vpython is used for these.
        self.container = container
        self.spheres = []
        self.helices = []
        self.arrows = []
        self.box = None
        self.scene = None
        self.stop_on_cycle = stop_on_cycle
        self.domains = []
        self.record_amplitudes = record_amplitudes
        self.time = 0.
        self.display_forces = display_forces
        self.record_pressure = record_pressure
        self.pressure = 0           # Pressure is set to 0 at the start
        self.steps = 0              # Number of steps taken, is set to 0 at the start
        self.pressure_history = []  # History of instantaneous values of pressure
        if display_forces:
            self._assign_pointers()
        if visualize:
            self.visualizer_type = visualizer_type
            self.create_vis(canvas=canvas)
        else:
            self.visualizer_type = None

    def create_vis(self, canvas=None):
        """
        Creates a visualisation. By default, uses a vpython canvas, so imports vpython.
        Subclass class and change this to change the way visualisations are made.
        Parameters
        ----------
        canvas: vpython canvas
            display into which the visualization is drawn
        """
        global vpython
        if vpython is None:
            import vpython
        # self.scene stores the display into which the system draws
        if not canvas:
            if not self.scene:
                self.scene = vpython.canvas()
        if canvas:
            self.scene = canvas

        # Draw particles if they aren't drawn yet
        for particle in self.particles:
            if not particle._visualized:
                self.spheres.append(vpython.sphere(pos=vector_from(particle.pos),
                    radius=particle.radius,
                    color=vector_from(particle.color),
                    opacity=particle.alpha,
                    display=self.scene,
                    make_trail=particle.make_trail))
                particle._visualized = True
        # Draw springs if they aren't drawn yet
        for spring in self.springs:
            if not spring._visualized:
                self.helices.append(vpython.helix(pos=vector_from(spring.particle_1.pos),
                    axis=vector_from(spring.axis),
                    radius=spring.radius,
                    opacity=spring.alpha,
                    color=vector_from(spring.color),
                    display=self.scene))
                spring._visualized = True
        # Draw pointers if they aren't drawn yet
        for pointer in self.pointerarrows:
            if not pointer._visualized:
                self.arrows.append(vpython.arrow(pos=vector_from(pointer.pos),
                    axis=vector_from(pointer.axis),
                    shaftwidth=pointer.shaftwidth,
                    opacity=pointer.alpha,
                    color=vector_from(pointer.color),
                    display=self.scene))
                pointer._visualized = True
        # Draw container if exists
        if not self.box and self.container:
            self.box = vpython.box(pos=vector_from(self.container.pos),
                length=self.container.dimension,
                width=self.container.dimension,
                height=self.container.dimension,
                color=vector_from(self.container.color),
                opacity=self.container.alpha)

    def update_vis(self):
        """
        Function which updates the visualization. In this case, updates the vpython visualization
        """
        global vpython
        if vpython is None:
            import vpython
        # For each type of object, only updates certain properties if they have changed.
        # Update display of particles(rendered as spheres)
        for (index, __) in enumerate(self.spheres):
            self.spheres[index].pos = vector_from(self.particles[index].pos)
            if self.particles[index]._vis_values_changed:
                self.particles[index]._vis_values_changed = False
                self.spheres[index].radius = self.particles[index].radius
                self.spheres[index].color = vector_from(self.particles[index].color)
                self.spheres[index].opacity = self.particles[index].alpha
                self.spheres[index].make_trail = self.particles[index].make_trail
        # Update display of springs(rendered as helices)
        for (index, __) in enumerate(self.helices):
            self.helices[index].pos = vector_from(self.springs[index].pos)
            self.helices[index].axis = vector_from(self.springs[index].axis)
            if self.springs[index]._vis_values_changed:
                self.springs[index]._vis_values_changed = False
                self.helices[index].radius = self.springs[index].radius
                self.helices[index].color = vector_from(self.springs[index].color)
                self.helices[index].opacity = self.springs[index].alpha
        # Update display of pointers(rendered as arrows)
        for (index, __) in enumerate(self.arrows):
            self.arrows[index].pos = vector_from(self.pointerarrows[index].pos)
            self.arrows[index].axis = vector_from(self.pointerarrows[index].axis)
            if self.pointerarrows[index]._vis_values_changed:
                self.arrows[index].shaftwidth = self.pointerarrows[index].shaftwidth
                self.arrows[index].color = vector_from(self.pointerarrows[index].color)
                self.arrows[index].opacity = self.pointerarrows[index].alpha
        if self.box:
            self.box.pos = vector_from(self.container.pos)
            if self.container._vis_values_changed:
                self.container._vis_values_changed = False
                self.box.color = vector_from(self.container.color)
                self.box.opacity = self.container.alpha
            if self.container.size_changed:
                self.container.size_changed = False
                self.box.width = self.container.dimension
                self.box.length = self.container.dimension
                self.box.height = self.container.dimension

    def run_for(self, time, dt=0.01):
        """
        Run simulation for a certain amount of time(as measured in the simulated system's time).
        Recommended to use this instead of simulate(dt) for most situations.
        Parameters
        ----------
        time: float
            Time for which the simulation will go on for in the system's time
        dt: float
            Time steps taken.
        """
        # Make pointer objects if not already created
        if self.display_forces:
            self._assign_pointers()
        # Create visualization if necessary
        if self.visualize:
            self.create_vis()
        # Import vpython if necessary
        if self.visualize and self.visualizer_type == "vpython":
            global vpython
            if vpython is None:
                import vpython
        # Simulate for given time
        while self.time < time:
            if self.visualize and self.visualizer_type == "vpython":
                vpython.rate(150)
            self.simulate(dt)
            if self.stop_on_cycle:
                if self._cycle_completed():
                    break
                for particle in self.particles:
                    particle.prev_pos = _duplicate_vector(particle.pos)

    def create_particles_in_container(self, number=0, speed=0, radius=0, inv_mass=1.):
        """
        Creates the given number of particles, with the given parameters, in random locations within the container.
        If the system has no container, this method will raise an error.
        Parameters
        ----------
        number: integer
            The number of particles to create.
        speed: float
            The speed of these particles.
        radius: float
            The radius of these particles
        inv_mass: float
            The inverse mass of these particles.
        """
        if not self.container:
            raise RuntimeError("No container in system")
        else:
            for i in range(0, number):
                l = self.container.dimension - radius * 2
                position = np.array([l * random() - l / 2, l * random() - l / 2, l * random() - l / 2]) + self.container.pos
                there_is = self._has_a_particle_at(position, radius)
                while there_is:
                    position = np.array([l * random() - l / 2, l * random() - l / 2, l * random() - l / 2]) + self.container.pos
                    there_is = self._has_a_particle_at(position, radius)
                velocity = speed * normalized(np.array([1 * random() - 0.5, 1 * random() - 0.5, 1 * random() - 0.5]))
                particle = Particle(pos=position, v=velocity, inv_mass=inv_mass, radius=radius)
                self.particles.append(particle)

    def simulate(self, dt=0.01):
        """
        Simulates a time-step with a step size of dt.
        Parameters
        ----------
        dt: float
            Size of time step taken
        """
        # Set up system if only just starting, in case things have changed since system was created.
        if self.time == 0 or self.steps == 0:
            if self.display_forces:
                self._assign_pointers()
            if self.visualize:
                self.create_vis()
            if self.collides and self.container:
                self._setup_domains()
                self._assign_particles_to_domains()

        # Make pointers appropriate sizes according to forces on particles.
        if self.display_forces:
            for index, pointer in enumerate(self.pointerarrows):
                pointer.pos = self.particles[index].pos + np.array([0., self.particles[index].radius + self.particles[index].radius, 0.])
                # Add radius to itself as slightly faster performance that way than doing 2*
                pointer.axis = self.particles[index].applied_force(self.particles[index].pos, self.time)

        # Forces on particles depending on their applied force
        for particle in self.particles:
            if not particle.fixed:
                particle.total_force = particle.applied_force(particle.pos, self.time)
                for spring in self.springs:
                    particle.total_force += spring.force_on(particle)
            else:
                particle.total_force = np.array([0., 0., 0.])
        # If particles interact with each other, add forces from fields.
        if self.interacts:
            for particle in self.particles:
                if not particle.fixed:
                    for other_particle in self.particles:
                        if particle != other_particle:
                            particle.total_force += other_particle.force_on(particle)

        if self.container:
            # Collision detection between particles using domains if has container.
            if self.collides:
                self._collision_detection_with_domains()

            # Collision detection with walls of container if has one.
            momenta_change = 0.
            # Go through all the particles
            for (index, particle) in enumerate(self.particles):
                # Check for collisions with walls
                wall_collision_index = self.container.contains(particle)
                if wall_collision_index is not True:
                    particle.v[wall_collision_index] = -particle.v[wall_collision_index]
                    if particle.pos[wall_collision_index] > 0:
                        particle.pos[wall_collision_index] = particle.pos[wall_collision_index] - (particle.radius / 10)
                    else:
                        particle.pos[wall_collision_index] = particle.pos[wall_collision_index] + particle.radius / 10
                        momenta_change += abs(particle.v[wall_collision_index] / particle.inv_mass)
            # Record the pressure, but only after a certain number of steps have been taken, when the system will be in equilibrium
            if self.record_pressure and self.steps > 200:
                instantaneous_pressure = (momenta_change / dt) / self.container.surface_area
                self._update_pressure(instantaneous_pressure)

        else:
            # If the system does not have a container, still does collision detection, but a lot less efficient.
            # (O(n^2) as opposed to O(nlog(n)))
            if self.collides:
                self._collision_detection()
        # Update particle positions according to the forces.
        for particle in self.particles:
            particle.update(dt)
        # record amplitudes/visualize if required.
        if self.record_amplitudes:
            self._get_amplitudes()
        # Update visualisation
        if self.visualize:
            self.update_vis()
        self.time += dt
        self.steps += 1

    def _setup_domains(self):
        """
        Sets up domains in the container.
        So that the more efficient collision detection algorithm functions:
        don't make the domain_size any less than 2r,
        and the step should be smaller than domain_size by at least 1 particle_radius.
        The actual values chosen here are just what seem to work well, do tinker with them
        """
        self.domain_size = self.particles[0].radius * 17
        step = self.particles[0].radius * 15.5
        x1 = np.arange(- self.container.dimension/2 + self.container.pos[0], (1.05 * self.container.dimension + self.container.pos[0])/2, step)
        x2 = np.arange(- self.container.dimension/2 + self.container.pos[1], (1.05 * self.container.dimension + self.container.pos[1])/2, step)
        x3 = np.arange(- self.container.dimension/2 + self.container.pos[2], (1.05 * self.container.dimension + self.container.pos[2])/2, step)
        x, y, z = np.meshgrid(x1, x2, x3)
        self._setup_domains_helper(self,x,y,z)
        self._assign_particles_to_domains()

    @np.vectorize
    def _setup_domains_helper(self, x, y, z):
        """
        Helper function for setup_domains(), made so that @np.vectorize can be used for faster execution
        """
        pos = np.array([x, y, z])
        self.domains.append(_Domain(pos, self.domain_size))

    def _assign_particles_to_domains(self):
        """
        Assigns a list of particles that are in each domain to the list of domains.
        """
        for domain in self.domains:
            domain.particles = []
            for particle in self.particles:
                if domain.contains(particle) is True:
                    domain.particles.append(particle)

    def _collided(self, particle_1, particle_2):
        """
        Checks if two particles have collided, returns true if they have, false if they haven't.
        Parameters
        ----------
        particle_1: Particle
            First particle that is being checked for collision
        particle_2: Particle
            The particle that particle_1 is being checked for a collision with
        """
        diff = particle_2.pos - particle_1.pos
        distance = np.sqrt(np.inner(diff, diff))  # more efficient than np.linalg.norm()
        # If the two particles overlap, then a collision is detected:
        if distance <= particle_2.radius + particle_1.radius:
            return True
        return False

    def _collision(self, particle_1, particle_2):
        """
        Defines how to change velocities when two particles, particle_1 and particle_2 collide.
        Parameters
        ----------
        particle_1: Particle
            One of the particles that has collided
        particle_2: Particle
            The other particle that has collided
        """
        diff = particle_1.pos - particle_2.pos
        axis_1 = normalized(diff)
        particle_1.pos += axis_1 * np.sqrt(np.inner(diff,diff)) / 20    # more efficient than np.linalg.norm()
        axis_2 = normalized(_perpendicular_vector(axis_1))
        axis_3 = normalized(_cross(axis_1, axis_2))
        axes = [axis_1, axis_2, axis_3]
        particle_1_v_after = np.array([0., 0., 0.])
        particle_2_v_after = np.array([0., 0., 0.])
        for i in range(1, 3):
            particle_1_v_after += np.dot(particle_1.v, axes[i]) * axes[i]
            particle_2_v_after += np.dot(particle_2.v, axes[i]) * axes[i]
        u1_parallel_axis_1 = np.dot(particle_1.v, axes[0])
        u2_parallel_axis_1 = np.dot(particle_2.v, axes[0])
        Z = (((u1_parallel_axis_1 * particle_2.inv_mass) + (u2_parallel_axis_1 * particle_1.inv_mass))
            / (particle_1.inv_mass + particle_2.inv_mass))
        particle_1_v_after += axes[0] * (-1* u1_parallel_axis_1 + Z*2)
        particle_1.v = particle_1_v_after
        axis_1 = normalized(particle_1.pos - particle_2.pos)
        particle_2_v_after += axes[0] * (-1* u2_parallel_axis_1 + Z*2)
        particle_2.v = particle_2_v_after

    def _has_a_particle_at(self, pos, radius):
        """
        Checks whether a particle alread exists at a certain position.
        Parameters
        ----------
        pos: vector
            Position that is being checked for the existence of a particle
        """
        diam = radius * 2
        for particle in self.particles:
            diff = particle.pos - pos
            distance = np.sqrt(np.inner(diff, diff))    # more efficient than np.linalg.norm()
            if distance <= diam:
                return True
        return False

    def _cycle_completed(self):
        """
        Function to see if an oscillation cycle has been completed. Only should work for normal modes.
        """
        if self.time == 0:
            for particle in self.particles:
                particle.initial_v = _duplicate_vector(particle.v)
                particle.initial_pos = _duplicate_vector(particle.pos)
                particle.prev_pos = _duplicate_vector(particle.pos)
        if self.time <= 0.01:
            # Don't expect any oscillations to finish faster than 0.01 seconds. Change if necessary.
            return False
        # If any of the particles have gone past their original position/are on their original positions,
        # and the velocity is in the same direction as originally, then cycle is completed.
        for particle in self.particles:
            if not particle.fixed:
                diff = element_mult(particle.pos - particle.initial_pos,
                    particle.prev_pos - particle.initial_pos)
                if (diff[0] <= 0. and diff[1] <= 0. and diff[2] <= 0.):
                    vel_diff = element_mult(particle.initial_v, particle.v)
                    if vel_diff[0] < 0. or vel_diff[1] < 0. or vel_diff[2] < 0.:
                        return False
                    return True
                return False
        return False

    def _get_amplitudes(self):
        """
        Gets the amplitudes for all the particles system. Only verified to work for 1D oscillations.
        """
        for particle in self.particles:
            vel_diff = element_mult(particle.v, particle.prev_v)
            if vel_diff[0] < 0. or vel_diff[1] < 0. or vel_diff[2] < 0.:
                if not particle.max_point.any():
                    if not particle.min_point.any():
                        particle.max_point = _duplicate_vector(particle.pos)
                        particle.min_point = _duplicate_vector(particle.pos)
                else:
                    if particle.pos[0] > particle.max_point[0]:
                        particle.max_point = _duplicate_vector(particle.pos)
                    else:
                        particle.min_point = _duplicate_vector(particle.pos)
            particle.prev_v = _duplicate_vector(particle.v)

    def _assign_pointers(self):
        """
        Assign Pointers to the forces on each particle.
        """
        for particle in self.particles:
            if not particle._pointer_assigned:
                self.pointerarrows.append(PointerArrow(pos=particle.pos,
                    axis=particle.applied_force(particle.pos, 0)))
                particle._pointer_assigned = True

    def _collision_detection(self):
        """
        Inefficient collision detection algorithm that works by checking each particle with all other particles after it in particles.
        """
        for index1, particle in enumerate(self.particles):
            for index2, other_particle in enumerate(self.particles):
                if index2 > index1:
                    if self._collided(particle, other_particle):
                        self._collision(particle, other_particle)

    def _collision_detection_with_domains(self):
        # Remake domains if the container's size has changed
        if self.container.size_changed:
            self._setup_domains()
            self._assign_particles_to_domains()
            self.container.size_changed = False
        # Assign particles to domains every 3 steps.
        # (Don't need to do every step, as particles are unlikely to move to different domains in 3 steps time.
        # Should be fine as if particles are moving so quickly, simulation is inacurrate anyway.)
        if self.steps % 3 == 0:
            self._assign_particles_to_domains()

        # The actual collision detection
        for domain in self.domains:
            for index_1, particle_1 in enumerate(domain.particles):
                for index_2, particle_2 in enumerate(domain.particles):
                    if (index_2 > index_1):
                        collided = self._collided(particle_1, particle_2)
                        if collided:
                            self._collision(particle_1, particle_2)

    def _update_pressure(self, instantaneous_pressure):
        """
        Updates pressure of system.
        Parameters
        ----------
        instantaneous_pressure: float
            Pressure at a certain time
        """
        self.pressure_history.append(instantaneous_pressure)
        self.pressure = np.mean(self.pressure_history)
