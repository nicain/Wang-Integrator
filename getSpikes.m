function getSpikes(outputFileName,Time,dt)

% Declare settings:
Nall = 2000;
Eall = 1600;
Iall = 400;
f = 0.15;

% Initialize relevant variables
spikeMatrix = sparse(Time/dt,Nall);

% Upload data set:
load SN_spikes.txt

% Create population spike data matrix:
spikeMatrix(sub2ind([Time/dt,Nall],SN_spikes(:,1),SN_spikes(:,2)+1))=true;
clear('SN_spikes');

% Create single-cell spike data matrix:
E1cell=1;
E2cell=1;
spikeMatrixSingleE1=spikeMatrix(:,E1cell);
spikeMatrixSingleE2=spikeMatrix(:,E1cell+E2cell);

% Save outpus
save([pwd,'/savedResults/',outputFileName],'spikeMatrix','spikeMatrixSingleE1','spikeMatrixSingleE2');

return