import hashlib
import os
import time
import random
import subprocess as sub
import csv
import math
import numpy
from shutil import copy2

def neural_net(pressure,delta_d,g1):
   g1=[1,1.2,1.4,-0.1,-1.4,0.9,-1.4]
   node1=pressure*g1[0]+delta_d*g1[1]+g1[2]
   node2=pressure*g1[3]+delta_d*g1[4]+g1[5]
   if node1>100:
      node1=100
   if node1<-100:
      node1=-100
   if node2>100:
      node2=100
   if node2<-100:
      node2=-100
   node1=1/(1+math.exp(-node1))-0.5
   node2=1/(1+math.exp(-node2))-0.5
   node3=g1[6]*(node1+node2)
   return node3
