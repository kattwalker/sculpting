import hashlib
import os
import os.path
import time
import random
import subprocess as sub
from file_write_up import write_voxelyze_file
from file_write_up import read_voxelyze_file
from file_write_up import read_voxelyze_file_un
from file_write_hoz import write_voxelyze_file_hoz
from file_write_hoz import read_voxelyze_file_hoz
from update_k import update_k
import numpy as np



def test_population_slope(population):
	pop_size=5
	t=0
	fitness=[]
	pop_morph=[] 
	pop_stiff=[]
	for t in range(pop_size):
		new_distance=[]
		morphology=[]
		stiffness=[]
		for x in range(9):
			l1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
			stiff=[5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,5000000,]
			morphology.append(l1)
			stiffness.append(stiff)
			
		
		pop_morph.append(morphology)
		pop_stiff.append(stiffness)

	old_distance=np.zeros(pop_size)
        #old_distance=old_distance+1.5
	for episode in range(15):
		print(old_distance)
		after_episode=[]
		new_distance=[]
		for t in range(pop_size):
			write_voxelyze_file(t,pop_morph[t],pop_stiff[t])
			ta=str(t)
			my_file=sub.Popen("./voxelyze  -f  katt"+ta+".vxa", shell=True)
			
		new_dist=np.zeros(pop_size)
		time.sleep(100)		
		count=0
		while count!=pop_size:
			for t in range(pop_size):
				ta=str(t)				
				if os.path.exists("pressuresmy_fitness"+ta+".xml.csv") and os.path.exists("my_fitness"+ta+".xml"):
						print('i am here katt!')
						results=read_voxelyze_file_un(ta)
						new_dist[t]=results[0]- results[1]*3
                                                if new_dist[t]==0:
                                                   time.sleep(50)
                                                   results=read_voxelyze_file_un(ta)
                                                   print(results)
						   new_dist[t]=results[0]- results[1]*3
						results=update_k(population[t],new_dist[t],old_distance[t],pop_morph[t],pop_stiff[t],ta,episode)
						pop_morph[t]=results[0]
						pop_stiff[t]=results[1]
						os.remove("pressuresmy_fitness"+ta+".xml.csv")
						os.remove("my_fitness"+ta+".xml")
						count=count+1
						
		old_distance=new_dist
              
   
	for t in range(pop_size):       
		#fitness.append(new_dist[t])
		fitness=new_dist             
	return fitness
  
