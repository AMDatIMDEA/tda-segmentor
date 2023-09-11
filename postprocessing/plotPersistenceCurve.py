#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TDA-Segmentor postprocessing scripts: 
    
Developers:       Aditya Vasudevan
                  Maciej Haranzcyk
                  
                  
                  IMDEA Materiales Institue,
                  Getafe, Spain



If tda-segmentor is used with the -module persistencecurve,
   then a 4 files are created:
       basefilename_minSaddlePairs.vtk
       basefilename_saddleSaddlePairs.vtk
       basefilename_SaddleMaxPairs.vtk
       basefilename_allPairs.vtk
   These stores as a vtk table, the persistence and the 
   number of respective critical pairs.
This script, takes this file as an input and plots the 
   number of critical pairs as function of the persistence

One must then choose a persistence value that corresponds
   to the plateau as the persistence threshold for the
   morse smale complex routine.    
   
   Provide as input the basefilename - FAU, MFI, etc.
    
"""

import os
import vtk
import math
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


inputfilename = basefilename + "_minSaddlePairs.vtk"

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

xlimMax = persistenceRange[1];
if (xlimMax > 100.0) : 
    xlimMax = 100.0


fntsize = 9.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(5.0, 5.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel(DataColumnName1,fontsize=fntsize)
plt.ylabel(DataColumnName2,fontsize=fntsize)
plt.plot(persistenceData,MinSaddlePairsData,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')

plt.xlim(1e-5,xlimMax)

outputfilename = DataColumnName2 + ".png"
outputfilename.replace(" ", "")

plt.savefig(outputfilename, format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);



#######################################################
##### Plotting the SaddleSaddlePairs vs persistence
#######################################################


inputfilename = basefilename + "_saddleSaddlePairs.vtk"

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

xlimMax = persistenceRange[1];
if (xlimMax > 100.0) : 
    xlimMax = 100.0

fntsize = 9.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(5.0, 5.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel(DataColumnName1,fontsize=fntsize)
plt.ylabel(DataColumnName2,fontsize=fntsize)
plt.plot(persistenceData,MinSaddlePairsData,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')

plt.xlim(1e-5,xlimMax)

outputfilename = DataColumnName2 + ".png"
outputfilename.replace(" ", "")

plt.savefig(outputfilename, format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);


#######################################################
##### Plotting the maxSaddlePairs vs persistence
#######################################################

inputfilename = basefilename + "_SaddleMaxPairs.vtk"

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

xlimMax = persistenceRange[1];
if (xlimMax > 100.0) : 
    xlimMax = 100.0

fntsize = 9.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(5.0, 5.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel(DataColumnName1,fontsize=fntsize)
plt.ylabel(DataColumnName2,fontsize=fntsize)
plt.plot(persistenceData,MinSaddlePairsData,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')

plt.xlim(1e-5,xlimMax)

outputfilename = DataColumnName2 + ".png"
outputfilename.replace(" ", "")

plt.savefig(outputfilename, format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);


#######################################################
##### Plotting allPairs vs persistence
#######################################################

inputfilename = basefilename + "_allPairs.vtk"

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

xlimMax = persistenceRange[1];
if (xlimMax > 100.0) : 
    xlimMax = 100.0

fntsize = 9.0;lwdth = 2.0;
fntsizelgd = 8.0;
fig = plt.figure(figsize=(5.0, 5.0)); # Figure handle

ax = plt.gca();

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.xticks(fontsize=fntsize)
plt.yticks(fontsize=fntsize)
plt.xlabel(DataColumnName1,fontsize=fntsize)
plt.ylabel(DataColumnName2,fontsize=fntsize)
plt.plot(persistenceData,MinSaddlePairsData,'-', color='BLACK', linewidth=lwdth)
ax.set_yscale('log')
ax.set_xscale('log')

plt.xlim(1e-5,xlimMax)

outputfilename = DataColumnName2 + ".png"
outputfilename.replace(" ", "")

plt.savefig(outputfilename, format="png", dpi=1000,bbox_inches="tight")
plt.close(fig);
