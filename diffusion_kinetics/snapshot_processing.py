"""
This is a module for plotting and processing the mass fractions data
"""
import numpy as np 
import matplotlib.pyplot as plt 

class SnapshotProcessing:
	def __init__(self,y,chem_obj,grid):
		self.y = y
		self.y_eq = chem_obj.equilibrate()
		self.gas = chem_obj.gas  
		self.grid = grid
		self.x = np.zeros((self.gas.n_species,self.grid.N))
		self.x_eq = self.x.copy()
		for i in range(self.grid.N):
			self.x[:,i] = self.y[:,i]/self.gas.molecular_weights*self.gas.mean_molecular_weight
			self.x_eq[:,i] = self.y_eq[:,i]/self.gas.molecular_weights*self.gas.mean_molecular_weight

	def Plot(self,species_names,output_fig_path):
		n = len(species_names)
		species_indices = []
		for species_name in species_names:
			species_indices.append(self.gas.species_index(species_name))
		
		fig = plt.figure(figsize=(10,5))
		ax1 = fig.add_subplot(1,2,1) 
		ax1.plot(self.grid.T,self.grid.PG,'-',self.grid.T_ref,self.grid.PG_ref,'--')
		ax1.set_yscale('log')
		ax1.invert_yaxis()
		ax1.set_xlabel('T(K)')
		ax1.set_ylabel('P(Pa)')
		
		colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
		color_index = 0
		ax2 = fig.add_subplot(1,2,2)
		for species_index, species_name in zip(species_indices,species_names):
			ax2.plot(self.x[species_index,:],self.grid.PG,linestyle='-',color = colors[color_index],label = species_name)
			ax2.plot(self.x_eq[species_index,:],self.grid.PG,linestyle='--',color = colors[color_index])
			color_index = color_index + 1
		ax2.set_xscale('log')
		ax2.set_yscale('log')
		ax2.invert_yaxis()
		ax2.set_xlabel('X')
		ax2.set_ylabel('P(Pa)')
		ax2.legend(loc=3)
		plt.savefig(output_fig_path)
		print "snapshot saved successfully"

	def SaveData(self,species_names,output_file_path,header_variable_names):
		# change to mole fractions
	    x_selected_species = np.zeros((len(species_names),self.grid.N))
	    i=0
	    for species_name in species_names: 
	        x_selected_species[i,:] = self.x[self.gas.species_index(species_name),:]
	        i=i+1
	    data = np.array([self.grid.z,self.grid.PG,self.grid.T])
	    data = np.transpose(np.concatenate((data,x_selected_species),axis=0))
	    np.savetxt(output_file_path,data,fmt='%.8e',header = header_variable_names)
	    print "File written successfully"





