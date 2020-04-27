import hashlib
import os
import time
import random
import subprocess as sub
import numpy as np
from test_population import test_population

from sort_population import sort_population
from new_population import new_population




if __name__ == "__main__":
    
    #how big is our population?
    pop_size=15
    population=[]
  
    #create a random population of genomes for testing
    for i in range(pop_size):     
            my_genome=[]
            for i in range(9):
                my_num=float(np.random.randint(-20,20))
                my_num=my_num/10
                my_genome.append(my_num)
            population.append(my_genome)  
 
    print(population)


    #now we are ready to start our GA
    num_gen=50 #how long will we run it for?
    for time in range(num_gen):
        
        print('GENERATION NUMBER:')
        print(time)

        new_distance2=test_population(population,pop_size) #this function does the testing of the population
        print(new_distance2)
        
        cpop=[] #allocate a place for our combined populaton and fitness - required for sorting
        fitness=new_distance2
        
        cpop=sort_population(population,fitness,pop_size,cpop) #now we sort the population to see how good it is
        population=[]
        population=new_population(population,fitness,pop_size,cpop) # create updated population for next round of GA
        print(population)
