#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TDA-Segmentor plotting persistence curve: 
    
Developers:       Aditya Vasudevan
                  Maciej Haranzcyk
                  
                  
                  IMDEA Materiales Institue,
                  Getafe, Spain



If tda-segmentor is used with the -module persistencecurve,
   then a 4 files are created:
       basefilename-minSaddlePairs.txt
       basefilename-saddleSaddlePairs.txt
       basefilename-SaddleMaxPairs.txt
       basefilename-allPairs.txt
This script, plots the number of critical pairs 
       as function of the persistence

One must then choose a persistence value that corresponds
   to the plateau as the persistence threshold for the
   morse smale complex routine.    
   
   Provide as input the basefilename - FAU, MFI, etc.
    
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


plt.ioff();

if (len(sys.argv) <= 1) : 
    print("Please provide the basefile name as input! (MFI, FAU etc.)");
    sys.exit();

basefilename = sys.argv[1];



#######################################################
##### Plotting the minSaddlePairs vs persistence
#######################################################


inputfilename = basefilename + "-minSaddlePairs.txt"
persistence, minSaddlePairs = np.loadtxt(inputfilename, unpack=True, skiprows=1)


fntsize = 9.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(2.0, 2.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel("persistence",fontsize=fntsize)
plt.ylabel("min-saddle-pairs",fontsize=fntsize)
plt.plot(persistence,minSaddlePairs,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')


outputfilename = "minSaddlePairs.png"

plt.savefig(outputfilename, format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);



#######################################################
##### Plotting the SaddleSaddlePairs vs persistence
#######################################################

inputfilename = basefilename + "-SaddleSaddlePairs.txt"
persistence, SaddleSaddlePairs = np.loadtxt(inputfilename, unpack=True, skiprows=1)


fntsize = 9.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(2.0, 2.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel("persistence",fontsize=fntsize)
plt.ylabel("saddle-saddle-pairs",fontsize=fntsize)
plt.plot(persistence,SaddleSaddlePairs,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')


outputfilename = "SaddleSaddlePairs.png"

plt.savefig(outputfilename, format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);


#######################################################
##### Plotting the maxSaddlePairs vs persistence
#######################################################


inputfilename = basefilename + "-maxSaddlePairs.txt"
persistence, maxSaddlePairs = np.loadtxt(inputfilename, unpack=True, skiprows=1)


fntsize = 9.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(2.0, 2.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel("persistence",fontsize=fntsize)
plt.ylabel("max-saddle-pairs",fontsize=fntsize)
plt.plot(persistence,maxSaddlePairs,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')


outputfilename = "maxSaddlePairs.png"

plt.savefig(outputfilename, format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);

#######################################################
##### Plotting allPairs vs persistence
#######################################################

inputfilename = basefilename + "-allPairs.txt"
persistence, allPairs = np.loadtxt(inputfilename, unpack=True, skiprows=1)


fntsize = 9.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(2.0, 2.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel("persistence",fontsize=fntsize)
plt.ylabel("all-pairs",fontsize=fntsize)
plt.plot(persistence,allPairs,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')


outputfilename = "allPairs.png"

plt.savefig(outputfilename, format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);
