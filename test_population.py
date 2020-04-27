import hashlib
import os
import os.path
import time
import random
import subprocess as sub
from file_write_up import write_voxelyze_file
from file_write_up import read_voxelyze_file
from file_write_up import read_voxelyze_file_un

from update_k import update_k
import numpy as np



def test_population(population,pop_size):

	t=0 #start at episode zero
	fitness=[]
	pop_morph=[] 
	pop_stiff=[]
	#This for loop initializes starting robot - aka it creates the morphology and stiffness of lump that is then read by
	#the .vxa file
	for t in range(pop_size):

		morphology=[]
		stiffness=[]
		for x in range(7):
			l1=[3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,]
			stiff=[500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,500000,]
			morphology.append(l1)
			stiffness.append(stiff)
			
		
		pop_morph.append(morphology)
		pop_stiff.append(stiffness)

	old_distance=np.zeros(pop_size) #we say that the starting distance is zero
	#start the loop - robot is allowed to sculpt after each episode and 'try again' to better locomote
	for episode in range(15):


		for t in range(pop_size):
			write_voxelyze_file(t,pop_morph[t],pop_stiff[t])
			ta=str(t)
			my_file=sub.Popen("./voxelyze  -f  katt"+ta+".vxa", shell=True) #here we actually do the testing. 
			#voxelyze will output a fitness file we then read 
			
		new_dist=np.zeros(pop_size)
		time.sleep(10)		
		count=0 #the count allows us to read the fitness files in any order rather than waiting for
		#a chronological file to be produced
		while count!=pop_size:
			for t in range(pop_size):
				ta=str(t)
				#only attempt to read if the files are there!				
				if os.path.exists("pressuresmy_fitness"+ta+".xml.csv") and os.path.exists("my_fitness"+ta+".xml"):
						
					results=read_voxelyze_file_un(ta)
						#make sure that it is going in a straight line - penalize if not
					new_dist[t]=results[1]-results[0]*3
					#sometimes we get an error saying results ==0. Recheck if thats the case
					if new_dist[t]==0:
						time.sleep(50)
						results=read_voxelyze_file_un(ta)
						new_dist[t]=results[1]- results[0]*3
                    
					#this function updates the morphology and does the sculpting
					results=update_k(population[t],new_dist[t],old_distance[t],pop_morph[t],pop_stiff[t],ta)
					pop_morph[t]=results[0]
					pop_stiff[t]=results[1]
					os.remove("pressuresmy_fitness"+ta+".xml.csv")
					os.remove("my_fitness"+ta+".xml")
					count=count+1
						
		old_distance=new_dist
    #our fitness is thedistance at the end of the final episode            
	for t in range(pop_size):       
		fitness.append(new_dist[t])
                print(new_dist[t])
             
	return fitness
  
