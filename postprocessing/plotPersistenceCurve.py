#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TDA-Segmentor postprocessing scripts: 
    
Developers:       Aditya Vasudevan
                  Maciej Haranzcyk
                  
                  
                  IMDEA Materiales Institue,
                  Getafe, Spain



If tda-segmentor is used with the -module persistencecurve,
   then a file basefilename_persistencecurve.vtk is created, 
   that stores as a vtk table, the persistence and the 
   number of minimum saddle pairs.
This script, takes this file as an input and plots the 
   number of minimum saddle pairs vs the persistence value

One must then choose a persistence value that corresponds
   to the plateau as the persistence threshold for the
   morse smale complex routine.    
    
"""

import os
import vtk
import math
import numpy as np
import matplotlib.pyplot as plt
import sys


plt.ioff();

if (len(sys.argv) <= 1) : 
    print("Please provide the basefilename_persistencecurve.vtk input file!");
    sys.exit();

inputfilename = sys.argv[1];

tableReader = vtk.vtkTableReader();
tableReader.SetFileName(inputfilename);
tableReader.Update();

table = vtk.vtkTable();
table.DeepCopy(tableReader.GetOutput(0));
DataColumnName1 = table.GetColumnName(0);
DataColumnName2 = table.GetColumnName(1);

persistenceData = table.GetColumnByName(DataColumnName1);
MinSaddlePairsData = table.GetColumnByName(DataColumnName2);

persistenceRange = persistenceData.GetRange();


fntsize = 13.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(5.0, 5.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel(r'persistence',fontsize=fntsize)
plt.ylabel(r'number of minimum-saddle pairs',fontsize=fntsize)
plt.plot(persistenceData,MinSaddlePairsData,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')

plt.xlim(1e-5,persistenceRange[1])

plt.savefig("persistenceCurve.png", format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);
