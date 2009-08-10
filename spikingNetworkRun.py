#!/usr/bin/env python
#PROGRAM spikingNetworkRun.py
#10 July 2009, Nicholas Cain
#UW Department of Applied Mathematics under Dr. Eric Shea-Brown

import sys
import os
import time
import random

# Declare settings:
numberOfTrials=1				# Number of calls of the spikingNetwork.c file
jobNameBase='Temp'				# Unique base-name for saving results and plots
dt = .02;						# Time-step, in ms.
simDuration = 200;				# Simulation duration, in ms.
stimOnset = 100;				# Time for stimulus onset, in ms.

# Declare strings that vary between OS:
if sys.platform == 'darwin':
	compileString='make'
	matlabSettings = '-nosplash -nojvm -r'
	runPrefix = './'
elif sys.platform == 'win32':
	compileString='cl spikingNetwork.c'
	matlabSettings = '-wait -automation -r'
	runPrefix = ''

# Compile current build of spikingNetwork program:
print 'Beginning Run:'
print '  Compiling spikingNetwork.c ...'
os.system(compileString)

# Run Simulations, each time creating a data output file:
print '  Running Simulations ...'
tBegin=time.mktime(time.localtime())
for i in range(1,numberOfTrials + 1):
	jobName=jobNameBase + '_' + str(i)
	print '    Simulation # ' + str(i) + ' now running...'
	os.system(runPrefix + 'spikingNetwork ' + str(simDuration/dt) + ' ' + str(stimOnset/dt) + ' ' + str(dt))
	print '    Compiling results ...'
	os.system('matlab ' + matlabSettings + ' "cd(\'' + os.getcwd() + '\'); getFR_simple(\'' + jobName + '\'' + ',0,' + str(simDuration) + ',' + str(dt) + ');exit"')

# Find AIC for this trial:
print '  Performing model comparison ...'
os.system('matlab ' + matlabSettings + ' "cd(\'' + os.getcwd() + '\'); FRFit(\'' + jobNameBase + '\'' + ',' + str(numberOfTrials) + ',1,' + str(stimOnset) + ');exit"')

# Finalize
tEnd = time.mktime(time.localtime())
secondsToCompute=tEnd-tBegin
print ('  Total Computation Time: ', time.strftime("H:%H M:%M S:%S",time.gmtime(secondsToCompute)))


