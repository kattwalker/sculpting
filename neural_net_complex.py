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
   g1=[-0.6, 0.6, 0.0, 0.3, 0.0, 0.5, 0.9, -0.2, -0.2, -0.5, -0.3, 0.7, 0.7, -0.1, -1.0, -0.9, -0.4, 0.9, -0.3, 0.1, 0.2, 0.3, -0.4, -0.1, 0.0, 0.9, -0.6, -0.8]
   bias1=1
   bias2=1

   node1=bias1*g1[0]+pressure*g1[1]+delta_d*g1[2]
   node2=bias1*g1[3]+pressure*g1[4]+delta_d*g1[5]
   node3=bias2*g1[6]+pressure*g1[7]+delta_d*g1[8]
   node4=bias2*g1[9]+pressure*g1[10]+delta_d*g1[11]
   if node1>100:
      node1=100
   if node1<-100:
      node1=-100
   if node2>100:
      node2=100
   if node2<-100:
      node2=-100
   if node3>100:
      node3=100
   if node3<-100:
      node3=-100
   if node4>100:
      node4=100
   if node4<-100:
      node4=-100   
   node1=1/(1+math.exp(-node1))-0.5
   node2=1/(1+math.exp(-node2))-0.5
   node3=1/(1+math.exp(-node3))-0.5
   node4=1/(1+math.exp(-node4))-0.5



   node5=node1*g1[12]+node2*g1[13]+node3*g1[14]+node4*g1[15]
   node6=node1*g1[16]+node2*g1[17]+node3*g1[18]+node4*g1[19]
   node7=node1*g1[20]+node2*g1[21]+node3*g1[22]+node4*g1[23]

   if node5>100:
      node5=100
   if node5<-100:
      node5=-100
   if node6>100:
      node6=100
   if node6<-100:
      node6=-100
   if node7>100:
      node7=100
   if node7<-100:
      node7=-100

   node6=1/(1+math.exp(-node6))-0.5
   node7=1/(1+math.exp(-node7))-0.5
   node5=1/(1+math.exp(-node5))-0.5

   node8=node5*g1[24]+node6*g1[25]+node7*g1[26]
   if node8<-100:
      node8=-100
   node8=1/(1+math.exp(-node8))-0.5

   node8=node8*g1[27]

   return node8
