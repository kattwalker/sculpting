import hashlib
import os
import time
import random
import subprocess as sub

def take_second(elm):
    return elm[1]
def sort_population(population,fitness,pop_size,cpop):

    for i in range(pop_size):
        genome=population[i]
        fit=fitness[i]
        cpop_genome=[]
        cpop_genome.append(genome)
        cpop_genome.append(fit)
        cpop.append(cpop_genome)
    cpop.sort(key=take_second,reverse=True)
    print("showing the sorted population!")
    for i in range(pop_size):
        best=cpop[i]
        
        print(best)

    best_genome=cpop[0]
    print("this is my best genome")
    print(best_genome)

    return cpop
