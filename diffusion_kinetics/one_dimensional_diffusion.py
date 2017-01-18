"""
solve the one dimensional diffusion equation for some species
dy/dt = 1/rho d(K rho dy/dz)/dz
where
y: mass fractions, y(i,z), where i is the chemical species 
rho: density
K: diffusion coefficient
z: uniform grid of location
all quantities can be a function of z
"""

import numpy as np

class ExplicitDiffusion:
	def __init__(self,K,grid):
		self.N = grid.N
		self.rho = grid.rho
		self.rho_K_boundary = (self.rho[0:self.N-1]*K[0:self.N-1] + self.rho[1:self.N]*K[1:self.N])/2.
		self.dz = grid.z[1] - grid.z[0]

	def advance(self,y,n_species,dt):
		"""
		advance a small timestep dt
		"""
		y_new = np.zeros_like(y)
		for k in range(n_species):
			for i in range(self.N):
				if i==0:
					phi_l = 0.
					phi_u = self.rho_K_boundary[i]*(y[k,i+1]-y[k,i])/self.dz
				elif i==self.N-1:
					phi_u = 0.
					phi_l = self.rho_K_boundary[i-1]*(y[k,i]-y[k,i-1])/self.dz;
				else:
					phi_u = self.rho_K_boundary[i]*(y[k,i+1]-y[k,i])/self.dz
					phi_l = self.rho_K_boundary[i-1]*(y[k,i]-y[k,i-1])/self.dz;
				y_new[k,i]= y[k,i] + (phi_u-phi_l)*dt/(self.dz*self.rho[i])
		return y_new
		
