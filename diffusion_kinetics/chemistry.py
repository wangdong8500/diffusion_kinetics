"""
This module constructs chemistry objects for time evolution
"""
import numpy as np 
import cantera as ct 

class Chemistry:
	def __init__(self,mechanism_path,composition,grid,atol=None,rtol=None):
		self.gas = ct.Solution(mechanism_path)
		self.gas.X = composition
		self.T = grid.T
		self.PG = grid.PG
		self.N = grid.N
		self.reactor = ct.IdealGasReactor(contents = self.gas, energy = 'off')
		self.network=ct.ReactorNet([self.reactor])
		if atol != None:
			self.network.atol = atol
		if rtol != None:
			self.network.rtol = rtol

	def equilibrate(self):
		n_species = self.gas.n_species
		y = np.zeros((n_species,self.N))
		# set up the initial state using equilibrium states
		for i in range(self.N):
			self.gas.TP = self.T[i],self.PG[i]
			self.gas.equilibrate('TP')
			y[:,i]=self.gas.Y
		return y


	def advance(self,y,dt):
		# chemical evolution for dt
		for i in range(self.N):
			self.gas.TPY = self.T[i],self.PG[i],y[:,i]
			self.reactor.syncState()
			self.network.reinitialize()
			self.network.set_initial_time(0.0)
			self.network.advance(dt)
			y[:,i]=self.gas.Y

		return y

