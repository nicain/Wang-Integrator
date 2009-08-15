function getFR_simple(outputFileName,plotsOn,Time,dt)
% I added a comment.
% Declare settings:
Nall = 2000;
Eall = 1600;
Iall = 400;
%dt = 0.02; 
%Time = 12000; % Set the simulation duration (= maxStep/50)
f = 0.15;

% Initialize relevant variables
spikeMatrix = sparse(Time/dt,Nall);
E1 = Eall * f;
E2 = Eall * f;

% Upload data set:
load SN_spikes.txt
spikeMatrix(sub2ind([Time/dt,Nall],SN_spikes(:,1),SN_spikes(:,2)+1))=true;
clear('SN_spikes');

w=20; %ms
increment=5; %ms

% Calculating the firing rate, sp/sec: 
totalLength=floor((Time/dt-w/dt-1)/(increment/dt));
FRE1=zeros(totalLength,1);
FRE2=zeros(totalLength,1);
plotT=zeros(totalLength,1);
for i=1:totalLength
    FRE1(i) = sum(sum(spikeMatrix((1+(increment/dt)*(i-1)):((increment/dt)*(i-1)+w/dt),1:E1)))./(E1*w/1000);
    FRE2(i) = sum(sum(spikeMatrix((1+(increment/dt)*(i-1)):((increment/dt)*(i-1)+w/dt),(E1+1):(E1+E2))))./(E2*w/1000);
    plotT(i) = dt+i*increment;
end

clear('spikeMatrix');
if plotsOn==1
    % Plot firing rate of each population
    figure(2)
    clf
    hold on
        plot(plotT,FRE1,'-r','LineWidth',1); 
        plot(plotT,FRE2,'-r','LineWidth',1); 
    hold off
    legend('Selective1','Selective2','Location','NorthWest');
    xlabel('Time (ms)');
    ylabel('Firing Rate (Spikes/sec)');
    ylim([0 50])
    
    % Save figures
    saveas(gcf,['./savedFigures/',outputFileName,'_','FR','.eps'],'eps')
    saveas(gcf,['./savedFigures/',outputFileName,'_','FR','.fig'],'fig')
    unix(['open ./savedFigures/',outputFileName,'_','FR','.eps']);
end

save([pwd,'/savedResults/',outputFileName],'FRE1','FRE2','plotT','increment')

return
