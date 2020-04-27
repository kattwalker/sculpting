import hashlib
import os
import time
import random
import subprocess as sub
import numpy as np

def assemble_population(population,genome):
    population.append(genome[:])
    return population

def new_population(population,fitness,pop_size,cpop):

    for i in range(pop_size):

        if i<1:
            new_genome=[]
            new_genome=cpop[i]
            new_genome=new_genome[0]
            print('this is my new genome!')
            print(new_genome)
            
        
        elif i<3:
            new_genome=[]
            which_genome=np.random.randint(0,3)
            selected_genome=cpop[which_genome]
            selected_genome=selected_genome[0]

            for t in range(10):
                my_num=float(np.random.randint(-10,10))
                my_location=(np.random.randint(0,28))              
                selected_genome[my_location]=my_num/10
                
            new_genome=selected_genome
            print('this is my new genome!')
            print(new_genome)

            
            
        elif i>3:
            new_genome=[]
            for t in range(28):
                my_num=float(np.random.randint(-10,10))
                my_num=my_num/10
                new_genome.append(my_num)
            print('this is my new genome!')
            print(new_genome)

    
            
        population=assemble_population(population,new_genome)

    return population
