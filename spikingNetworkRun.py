#!/usr/bin/env python
#PROGRAM spikingNetworkRun.py
#10 July 2009, Nicholas Cain
#UW Department of Applied Mathematics under Dr. Eric Shea-Brown

import sys
import os
import time
import random

# Simulation settings:
trialRangeBegin=101				# trialRangeBegin - trialRangeEnd = Number of calls of the spikingNetwork.c file
trialRangeEnd=250	
jobNameBase='2sim_1on'				# Unique base-name for saving results and plots
dt = .02;						# Time-step, in ms.
simDuration = 2000;				# Simulation duration, in ms.
stimOnset = 1000;				# Time for stimulus onset, in ms.

# Analysis settings:
analyze = 0;
x_label = 'S1';
y_label = 'S1';
numberOfTrials = 12;

# Declare strings that vary between OS:
if sys.platform == 'darwin':
	compileString='gcc -o spikingNetwork spikingNetwork.c'
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
for i in range(trialRangeBegin,trialRangeEnd + 1):
	jobName=jobNameBase + '_' + str(i)
	print '    Simulation # ' + str(i) + ' now running...'
	os.system(runPrefix + 'spikingNetwork ' + str(simDuration/dt) + ' ' + str(stimOnset/dt) + ' ' + str(dt))
	print '    Compiling results ...'
	os.system('matlab ' + matlabSettings + ' "cd(\'' + os.getcwd() + '\'); getFR(\'' + jobName + '\'' + ',0,' + str(simDuration) + ',' + str(dt) + ');exit"')

# Find AIC for this trial:
if analyze == 1:
	print '  Performing model comparison ...'
	os.system('matlab ' + matlabSettings + ' "cd(\'' + os.getcwd() + '\'); spikeRateAnalysis(\'' + jobNameBase + '\'' + ',' + str(numberOfTrials) + ',' + str(stimOnset) + ',\'' + x_label + '\',\'' + y_label + '\',1);exit"')
else:
	print '  Skipping model comparison'
# Finalize
tEnd = time.mktime(time.localtime())
secondsToCompute=tEnd-tBegin
print ('  Total Computation Time: ', time.strftime("H:%H M:%M S:%S",time.gmtime(secondsToCompute)))


