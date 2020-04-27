import hashlib
import os
import time
import random
import subprocess as sub
import csv
import math
import numpy
from shutil import copy2
from neural_net import neural_net
import numpy as np

def update_k(population,new_dist,old_distance,morphology,stiffness,ta):

       da=1
       g1=population

       #calculate gradient of the landscape - this will be used to calculate how much we should change our stiffness
       d=(new_dist-old_distance)/10
       #make sure we have negative value if we are getting worse
       if old_distance>new_dist and d>0:
          d=d*-1

       delta_d=d


       time_steps=0
       average_force=[]
       for z in range(216):
              average_force.append(0)
       #read the outputted kinetic energy and get it into a format we can use
       with open("kemy_fitness"+ta+".xml.csv") as csvfile:
              reader = csv.reader(csvfile)
              for row in reader:
                     time_steps=time_steps+1
                     length=len(row)-1
                     for x in range(length):
                            number=float(row[x])
                            average_force[x]=average_force[x]+number*10 #N.B average KE is very tiny 
                    
       length=len(average_force)-1
      
       ultimate_average=0
       #clculate the 'ultimate average' aka the average K.E in ALL the voxels 
       for i in range(length):
          average_force[i]=average_force[i]/(length)
	   ultimate_average=ultimate_average+average_force[i]/length       
 
      
       new_morphology=[]
       new_stiffness=[]
       size=0
       #here we do the sculpting, we loop through every voxel
       for x in range(6):
              
              line=morphology[x]
              line_s=stiffness[x]
              
              for y in range(36):
                     offset=(x*36)+y
                     #compare the k.e in each voxel with the average fore across the robot
                     pressure=10*(ultimate_average-average_force[offset])   
                     #it says pressure but its actully k.e, sorry.
                     # now we calculate the stiffness change using the neural network. remember that the genome is
                     # the weights of the neural network                  
                     stiffness_ind=neural_net(pressure,delta_d,g1)
                     #print(stiffness_ind)
                     line_s[y]=line_s[y]+(stiffness_ind*100000)
                     #if the stiffness is below a certain value, remove it
                     if line_s[y]<100000:
                            line_s[y]=100000
                            line[y]=0
                            size=size+1
                     elif line_s[y]>50000000:
                            line_s[y]=50000000
                            line[y]=2
                     else:
                            line[y]=3
	      #a little hacky bug fix
              if x==0:
                     line_s[35]=line_s[29]
                     line[35]=line[29]
              new_morphology.append(line)
              new_stiffness.append(line_s)
       size=216-size

       morphology=[]
       stiffness=[] 

       #update the morphology and stiffness      
       for x in range(6):
              line=new_morphology[x]
              line_s=new_stiffness[x]

              morphology.append(line)
              stiffness.append(line_s)
       results=[]
       results.append(morphology)
       results.append(stiffness)
       
       return results
        
