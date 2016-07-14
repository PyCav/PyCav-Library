import numpy as np

def explicit_wave_equation(psi_0, dx, N, a = 1., c = 1., bound_cond = 'periodic',init_grad = None, init_vel = None):
	# CFL condition
	dt = a*dx/c

	t = np.linspace(0.,(N-1)*dt,N)

	dim = len(psi_0.shape)
	len_x = psi_0.shape[0]
	if dim > 1:
		len_y = psi_0.shape[1]

	if dim > 2 or dim < 1:
		print('1D or 2D only, adjust initial wave array')
		return 0

	if dim == 1:
		psi = np.zeros((len_x,N))
		psi[:,0] = psi_0

	if dim == 2:
		psi = np.zeros((len_x,len_y,N))
		psi[:,:,0] = psi_0

	if init_vel == None:
		def _vel_0():
			def _V(psi_0):
				return np.zeros_like(psi_0)
			return _V
		V = _vel_0()
	else:
		V = init_vel

	if init_grad == None:
		def _grad_0():
			def _grad(psi_0):
				if dim == 1:
					central_grad = np.zeros_like(psi_0)
					psi_right = psi_0[2:]
					psi_left = psi_0[:-2]
					central_grad[1:-1] = (psi_right-psi_left)/(2.*dx)

					if bound_cond == 'reflective':
						central_grad[0] = 0.0
						central_grad[-1] = 0.0

					else:
						central_grad[0] = (psi_0[0]-psi_0[1])/dx
						central_grad[-1] = (psi_0[-2]-psi_0[-1])/dx

					return central_grad

				if dim == 2:
					central_grad_x = np.zeros_like(psi_0)
					central_grad_y = np.zeros_like(psi_0)
					psi_x_right = psi_0[2:,:]
					psi_x_left = psi_0[:-2,:]
					central_grad_x[1:-1,:] = (psi_x_right-psi_x_left)/(2.*dx)

					psi_y_right = psi_0[:,2:]
					psi_y_left = psi_0[:,:-2]
					central_grad_y[:,1:-1] = (psi_y_right-psi_y_left)/(2.*dx)

					if bound_cond == 'reflective':
						central_grad_x[0,:] = 0.0
						central_grad_x[-1,:] = 0.0
						central_grad_y[:,0] = 0.0
						central_grad_y[:,-1] = 0.0

					else:
						central_grad[0] = (psi_0[0]-psi_0[1])/dx
						central_grad[-1] = (psi_0[-2]-psi_0[-1])/dx

					return central_grad_x,central_grad_y

			return _grad
		D = _grad_0()
	else:
		D = init_grad

	def _initialise(r,s,l = None):
		if dim == 1:
			s[:,0] = V(psi_0)
			r[:,0] = D(psi_0)
			return r,s
		if dim == 2:
			s[:,:,0] = V(psi_0)
			r[:,:,0],l[:,:,0] = D(psi_0)
			return r,s,l			

	def _derivatives(psi):
		if dim == 1:
			r = np.zeros((len_x,2))
			s = np.zeros((len_x,2))
		if dim == 2:
			r = np.zeros((len_x,len_y,2))
			s = np.zeros((len_x,len_y,2))
			l = np.zeros((len_x,len_y,2))			

		for i in range(0,N-1):
			if dim == 1:
				if i == 0:
					r,s = _initialise(r,s)

				r_right = r[2:,0]
				r_left  = r[:-2,0]
				s_right = s[2:,0]
				s_left  = s[:-2,0]

				r[1:-1,1] = r[1:-1,0]+0.5*a*(s_right-s_left+a*(r_right+r_left-2*r[1:-1,0]))
				s[1:-1,1] = s[1:-1,0]+0.5*a*(r_right-r_left+a*(s_right+s_left-2*s[1:-1,0]))

				if bound_cond == 'reflective':
					r[0,:] = 0.0
					r[-1,:] = 0.0
					s[0,1] = s[0,0]+a*r[1,0]
					s[-1,1] = s[-1,0]-a*r[-2,0]

				if bound_cond == 'periodic':
					r[0,:] = r[-2,:]
					s[0,:] = s[-2,:]
					r[-1,:] = r[1,:]
					s[-1,:] = s[1,:]					

				elif bound_cond == 'fixed':
					s[0,:] = 0.0
					s[-1,:] = 0.0
					r[0,1] = r[0,0]+a*s[1,0]
					r[-1,1] = r[-1,0]-a*s[-2,0]

				psi[:,i+1] = psi[:,i]+0.5*dt*(s[:,1]+s[:,0])

				r[:,0] = r[:,1]
				s[:,0] = s[:,1]

			if dim == 2:
				if i == 0:
					r,s,l = _initialise(r,s,l = l)

				r_right = r[2:,:,0]
				r_left  = r[:-2,:,0]

				l_right = l[:,2:,0]
				l_left  = l[:,:-2,0]

				s_x_right = s[2:,:,0]
				s_x_left  = s[:-2,:,0]

				s_y_right = s[:,2:,0]
				s_y_left  = s[:,:-2,0]

				r[1:-1,:,1] = r[1:-1,:,0]+0.5*a*(s_x_right-s_x_left+a*(r_right+r_left-2*r[1:-1,:,0]))
				l[:,1:-1,1] = l[:,1:-1,0]+0.5*a*(s_y_right-s_y_left+a*(l_right+l_left-2*l[:,1:-1,0]))

				s[1:-1,1:-1,1] = s[1:-1,1:-1,0]+0.5*a*(r_right[:,1:-1]-r_left[:,1:-1]+l_right[1:-1,:]-l_left[1:-1,:]+
								 a*(s_x_right[:,1:-1]+s_x_left[:,1:-1]+s_y_right[1:-1,:]+s_y_left[1:-1,:]-4*s[1:-1,1:-1,0]))

				if bound_cond == 'reflective':
					r[0,:,:] = 0.0
					r[-1,:,:] = 0.0

					l[:,0,:] = 0.0
					l[:,-1,:] = 0.0

					# x boundaries
					s[0,1:-1,1] = s[0,1:-1,0]+a*r[1,1:-1,0]+0.5*a*(l_right[0,:]-l_left[0,:])
					s[-1,1:-1,1] = s[-1,1:-1,0]-a*r[-2,1:-1,0]+0.5*a*(l_right[-1,:]-l_left[-1,:])

					# y boundaries
					s[1:-1,0,1] = s[1:-1,0,0]+a*l[1:-1,1,0]+0.5*a*(r_right[:,0]-r_left[:,0])
					s[1:-1,-1,1] = s[1:-1,-1,0]-a*l[1:-1,-2,0]+0.5*a*(r_right[:,-1]-r_left[:,-1])

					# corners
					s[0,0,1] = s[0,0,1]+a*r[1,0,0]+a*l[0,1,0]
					s[0,-1,1] = s[0,-1,1]+a*r[1,-1,0]-a*l[0,-2,0]
					s[-1,0,1] = s[-1,0,1]-a*r[-2,0,0]+a*l[-1,1,0]
					s[-1,-1,1] = s[-1,-1,1]-a*r[-2,-1,0]-a*l[-1,-2,0]

				if bound_cond == 'periodic':
					psi[0,:,i+1] = psi[-1,:,i+1]
					psi[:,0,i+1] = psi[:,-1,i+1]

				elif bound_cond == 'fixed':
					psi[0,:,i+1] = 0.0
					psi[-1,:,i+1] = 0.0
					psi[:,0,i+1] = 0.0
					psi[:,-1,i+1] = 0.0

				psi[:,:,i+1] = psi[:,:,i]+0.5*dt*(s[:,:,1]+s[:,:,0])

				r[:,:,0] = r[:,:,1]
				l[:,:,0] = l[:,:,1]
				s[:,:,0] = s[:,:,1]


		return psi

	def _main_step(psi):
		psi = _derivatives(psi)

	_main_step(psi)

	return psi, t