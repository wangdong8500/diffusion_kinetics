"""
This module is responsible for interpolating a one dimensional uniform grid of atmosphere
"""
import numpy as np 
from scipy import interpolate

class UniformGrid:
	def __init__(self,N,PG_min,PG_max,PG_ref, T_ref, z_ref, rho_ref):
		self.N = N
		self.T_ref = T_ref
		self.PG_ref = PG_ref
		z_fit = interpolate.interp1d(np.log10(PG_ref),z_ref,kind = 'cubic')
		z_min = z_fit(np.log10(PG_min))
		z_max = z_fit(np.log10(PG_max))
		# construct uniform z grid
		self.z = np.linspace(z_min,z_max,N)
		T_fit = interpolate.interp1d(z_ref,T_ref,kind='cubic')
		self.T = T_fit(self.z)
		logP_fit = interpolate.interp1d(z_ref,np.log10(PG_ref),kind='cubic')
		self.PG = 10**logP_fit(self.z)
		rho_fit = interpolate.interp1d(z_ref,np.log10(rho_ref),'cubic')
		self.rho = 10**rho_fit(self.z)

	def Plot(self):
		import matplotlib.pyplot as plt 
		fig = plt.figure()
		ax1 = fig.add_subplot(2,2,1)
		ax1.plot(self.T,self.PG,lw = 2,marker = '*')
		ax1.set_xlabel('T')
		ax1.set_ylabel('P')
		ax1.set_yscale('log')
		ax2 = fig.add_subplot(2,2,2)
		ax2.plot(self.z,self.PG,lw = 2, marker = '*')
		ax2.set_xlabel('z')
		ax2.set_ylabel('P')
		ax2.set_yscale('log')
		ax3 = fig.add_subplot(2,2,3)
		ax3.plot(self.rho,self.PG,lw = 2, marker = '*')
		ax3.set_xlabel(r'$\rho$')
		ax3.set_ylabel('P')
		ax3.set_yscale('log')
		ax3.set_xscale('log')
		plt.show()
