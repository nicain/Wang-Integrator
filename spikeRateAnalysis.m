function spikeRateAnalysis(jobNameBase,numberOfJobs,switchOver,x_label,y_label,plotsOn)

% Set up directory structure:
workingDirectory=pwd;
saveDirectory=[workingDirectory,'/savedFigures'];
loadDirectory=[workingDirectory,'/savedResults'];

% Declare settings:
maxEval=5000;
maxIter=5000;

% Upload data:
for n=1:numberOfJobs
    currFullJobName=[loadDirectory,'/',jobNameBase,'_',num2str(n),'.mat'];
    load(currFullJobName);
    FRES1(n,:)=FRE1;
    FRES2(n,:)=FRE2;
    plotTMat(n,:)=plotT;
end

% Create a design matrix (dm)
A = [0 0; 1 0; 0 1; 1 1];
B = [ 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];

dm = [];
for i = 1:size(B,1)
    B_temp = ones(size(A,1),1)*B(i,:);
    dm = [dm;A,B_temp];
end

numberOfModels = size(dm,1);

% Models to compare:
for n=1:numberOfJobs
    initial_condition = [0.1,0.01,rand(1,3)];
    while initial_condition(3) < 0.1;
        initial_condition(3) = rand(1);
    end
    for i=1:numberOfModels;
        mu{i} = @(X,Y,theta,dv) theta(1).*dv(1) + theta(2).*dv(2).*X;
        g{i} = @(X,Y,theta,dv) theta(3).*dv(3) + theta(4).*dv(4).*abs(Y) + theta(5).*dv(5).*sqrt(abs(Y));
        IC{n,i} = initial_condition;
        theta{n,i} = logical(dm(i,:));
        description{i} = modelDescription(dm(i,:));
        likelihoodPre(n,i) = NaN;
        likelihoodPost(n,i) = NaN;
        AICPre(n,i) = NaN;
        AICPost(n,i) = NaN;
    end
end

% Initialize relevant variables
options = optimset('MaxFunEvals', maxEval,'MaxIter',maxIter);

dt=plotT(2)-plotT(1);
onIndex=floor(switchOver/increment);

pre_candidates  = 1:size(dm,1); % [1 2 3 7 8 9];
post_candidates = 1:size(dm,1);

%% Pre-stimulus epoch

% Minimize the negative log-likelyhood for each model, selective1-selective2, before input:
warning off all
for n=1:numberOfJobs
    % Choose xCurr
    if strcmp(x_label,'S1-S2')
        xCurr = FRES1(n,1:(onIndex))-FRES2(n,1:(onIndex));
    elseif strcmp(x_label,'S1+S2')
        xCurr = FRES1(n,1:(onIndex))+FRES2(n,1:(onIndex));
    elseif strcmp(x_label,'S1')
        xCurr = FRES1(n,1:(onIndex));
    end

    % Choose yCurr
    if strcmp(y_label,'S1-S2')
        yCurr = FRES1(n,1:(onIndex))-FRES2(n,1:(onIndex));
    elseif strcmp(y_label,'S1+S2')
        yCurr = FRES1(n,1:(onIndex))+FRES2(n,1:(onIndex));
    elseif strcmp(y_label,'S1')
        yCurr = FRES1(n,1:(onIndex));
    end
    
    for i=pre_candidates
        dv = dm(i,:);
        % Develop f:
        f{i}=@(xNext,xCurr,yCurr,theta,dv) normpdf(xNext,xCurr+mu{i}(xCurr,yCurr,theta,dv)*dt,abs(g{i}(xCurr,yCurr,theta,dv))*sqrt(dt));

        % Develop Log-Likelyhood function, using the looping method to help:
        L{i}=@(theta) likelyhoodFunction(theta,f{i},xCurr,yCurr,dv);
        
        % Do the search:
        [temp1,temp2,temp3]=fminsearch(L{i},IC{n,i},options);
        thetaEstPre{n,i}=temp1;
        likelihoodPre(n,i)=temp2;
        exitflag=temp3;
        
        % Determine the Akaiki Information for each model, provided convergence:
        thetaVector = nan(1,5);
        if exitflag==1
            AICPre(n,i)=2*likelihoodPre(n,i)+2*length(thetaEstPre{n,i});
            thetaVector(logical(dm(i,:))) = thetaEstPre{n,i}(logical(dm(i,:)));
        else
            [thetaEstPre{n,i},likelihoodPre(n,i),exitflag]=fminunc(L{i},IC{n,i},options);
            if exitflag==1
                AICPre(n,i)=2*likelihoodPre(n,i)+2*length(thetaEstPre{n,i});
                thetaVector(logical(dm(i,:))) = thetaEstPre{n,i}(logical(dm(i,:)));
            else
                AICPre(n,i)=Inf;
                likelihoodPre(n,i)=-Inf;
            end
        end
        thetaPreMatrix(i,:,n) = thetaVector;
    end
end
optimalModelCountPre=zeros(1,length(description));
for n=1:numberOfJobs
    [val,ind]=min(AICPre(n,:));
    optimalModelCountPre(ind)=optimalModelCountPre(ind)+1;
end

%% Post-stimulus epoch

% Minimize the negative log-likelyhood for each model, selective1-selective2, after input:
for n=1:numberOfJobs
    % Choose xCurr
    if strcmp(x_label,'S1-S2')
        xCurr = FRES1(n,(onIndex+1):end)-FRES2(n,(onIndex+1):end);
    elseif strcmp(x_label,'S1+S2')
        xCurr = FRES1(n,(onIndex+1):end)+FRES2(n,(onIndex+1):end);
    elseif strcmp(x_label,'S1')
        xCurr = FRES1(n,(onIndex+1):end);
    end

    % Choose yCurr
    if strcmp(y_label,'S1-S2')
        yCurr = FRES1(n,(onIndex+1):end)-FRES2(n,(onIndex+1):end);
    elseif strcmp(y_label,'S1+S2')
        yCurr = FRES1(n,(onIndex+1):end)+FRES2(n,(onIndex+1):end);
    elseif strcmp(y_label,'S1')
        yCurr = FRES1(n,(onIndex+1):end);
    end
    for i=post_candidates
        dv = dm(i,:);
        % Develop f:
        f{i}=@(xNext,xCurr,yCurr,theta,dv) normpdf(xNext,xCurr+mu{i}(xCurr,yCurr,theta,dv)*dt,abs(g{i}(xCurr,yCurr,theta,dv)*sqrt(dt)));

        % Develop Log-Likelyhood function, using the looping method to help:
        L{i}=@(theta) likelyhoodFunction(theta,f{i},xCurr,yCurr,dv);
        
        % Do the search:
        [temp1,temp2,temp3]=fminsearch(L{i},IC{n,i},options);
        thetaEstPost{n,i}=temp1;
        likelihoodPost(n,i)=temp2;
        exitflag=temp3;
        
        thetaVector = nan(1,5);
        
        % Determine the Akaiki Information for each model, provided convergence:
        if exitflag==1
            AICPost(n,i)=2*likelihoodPost(n,i)+2*length(thetaEstPost{n,i});
            thetaVector(logical(dm(i,:))) = thetaEstPost{n,i}(logical(dm(i,:)));
        else
            [thetaEstPost{n,i},likelihoodPost(n,i),exitflag]=fminunc(L{i},IC{n,i},options);
            if exitflag==1
                AICPost(n,i)=2*likelihoodPost(n,i)+2*length(thetaEstPost{n,i});
                thetaVector(logical(dm(i,:))) = thetaEstPost{n,i}(logical(dm(i,:)));
            else
                AICPost(n,i)=Inf;
                likelihoodPost(n,i)=-Inf;
            end
        end
        thetaPostMatrix(i,:,n) = thetaVector;
    end
end
optimalModelCountPost=zeros(1,length(description));
for n=1:numberOfJobs
    [val,ind]=min(AICPost(n,:));
    optimalModelCountPost(ind)=optimalModelCountPost(ind)+1;
end
warning on all

thetaPreMatrix=gPosNegFilter(thetaPreMatrix);
thetaPostMatrix=gPosNegFilter(thetaPostMatrix);

%% Plotting figures

if plotsOn==1
    % Plot firing rate of selective populations
    figure(1)
    clf
    plot(plotT,FRES1-FRES2,'-r','LineWidth',1);
    hold on
        plot([switchOver,switchOver],[min(min(FRES1-FRES2)) max(max(FRES1-FRES2))],'--')
    hold off
    xlabel('Time (ms)');
    ylabel('Firing Rate (Spikes/sec)');
    ylim([min(min(FRES1-FRES2)) max(max(FRES1-FRES2))])
    title('Difference in FR')
    
    % Bar chart to compare Akaiki Information for each model;
    figure(2)
    clf
    ind=find(optimalModelCountPre~=0);
    barh(optimalModelCountPre(ind));
    set(gca,'YTick',[(1/7):(1/7):length(ind)]+(.5-1/7));
    set(gca,'YTickLabel',getDescriptionWithTheta(description,thetaPreMatrix,ind))
    title(['AIC Minimization Histogram, Pre-Stimulus: X=',x_label,', Y=',y_label]);
    
    figure(3)
    clf
    ind=find(optimalModelCountPost~=0);
    barh(optimalModelCountPost(ind));
    set(gca,'YTick',[(1/7):(1/7):length(ind)]+(.5-1/7));
    set(gca,'YTickLabel',getDescriptionWithTheta(description,thetaPreMatrix,ind))
    title(['AIC Minimization Histogram, Post-Stimulus: X=',x_label,', Y=',y_label]);

    % Save figures
    saveFigureBase=[saveDirectory,'/',jobNameBase,'_',x_label,'_',y_label];
    saveas(1,[saveFigureBase,'_MultiFR.eps'],'eps')
    saveas(1,[saveFigureBase,'_MultiFR.fig'],'fig')
    
    saveas(2,[saveFigureBase,'_AICPre.eps'],'eps')
    saveas(2,[saveFigureBase,'_AICPre.fig'],'fig')
    
    saveas(3,[saveFigureBase,'_AICPost.eps'],'eps')
    saveas(3,[saveFigureBase,'_AICPost.fig'],'fig')
    
    save([workingDirectory,'/savedResults/',jobNameBase,'_',x_label,'_',y_label,'_results.mat'],'AICPre','AICPost','thetaPreMatrix','thetaPostMatrix');
end 

return

%% getDescriptionWithTheta

function descriptionNew=getDescriptionWithTheta(description,thetaMatrix,ind)

c=1;
for i=1:length(ind)
    descriptionNew{c}=' ';
    c=c+1;
    nOfTheta=0;
    for j=5:-1:1
        tmp=squeeze(thetaMatrix(ind(i),j,:));
        ind2=find(~isnan(tmp));
        if length(ind2)~=0
            currentMean=mean(tmp(ind2));
            currentStdDev=std(tmp(find(~isnan(tmp))));
            currentString=[num2str(currentMean,'%10.2f'),' (',num2str(currentStdDev,'%10.2f'),')  '];
            switch j
                case 1
                    descriptionNew{c}=['a = ',currentString];
                    nOfTheta=nOfTheta+1;
                    c=c+1;
                case 2
                    descriptionNew{c}=['b = ',currentString];
                    nOfTheta=nOfTheta+1;
                    c=c+1;
                case 3
                    descriptionNew{c}=['c = ',currentString];
                    nOfTheta=nOfTheta+1;
                    c=c+1;
                case 4
                    descriptionNew{c}=['d = ',currentString];
                    nOfTheta=nOfTheta+1;
                    c=c+1;
                case 5
                    descriptionNew{c}=['e = ',currentString];
                    nOfTheta=nOfTheta+1;
                    c=c+1;
            end
        end
    end
    descriptionNew{c}=description{ind(i)};
    c=c+1;
    for k=1:5-nOfTheta
        descriptionNew={descriptionNew{1:(end-nOfTheta-1)},' ',descriptionNew{(end-nOfTheta):end}};
        c=c+1;
    end
end

return

%% gPosNegFilter:

function A=gPosNegFilter(A)

for n=1:size(A,3)
    for i=1:size(A,1)
        if sum(isnan(A(i,:,n)))~=5
            if ~isnan(A(i,3,n))
                if A(i,3,n)<0
                    A(i,3:5,n)=-1*A(i,3:5,n);
                end
            elseif ~isnan(A(i,4,n))
                if A(i,4,n)<0
                    A(i,4:5,n)=-1*A(i,4:5,n);
                end
            elseif ~isnan(A(i,5,n))
                if A(i,5,n)<0
                    A(i,5,n)=-1*A(i,5,n);
                end
            end
        end
    end
end

return
%% likelyhoodFunction:

function returnMe=likelyhoodFunction(theta,f,x,y,dv)
    for i=1:size(theta,1)
        returnMe(1,i)=-sum(log(f(x(2:end),x(1:(end-1)),y(1:(end-1)),theta(i,:),dv)));
    end
return

%% modelDescription:

function title = modelDescription(designVector)
    title = 'dX = ';
    switch mat2str(designVector(1:2))
        case '[1 0]'
            title = [title,'a*dt'];
        case '[0 1]'
            title = [title,'b*X*dt'];
        case '[1 1]'
            title = [title,'(a + bX)*dt + '];
    end

    switch mat2str(designVector(3:5))
        case '[0 0 1]'
            title = [title,'e*sqrt(X)*dW'];
        case '[0 1 0]'
            title = [title,'d*X*dW'];
        case '[0 1 1]'
            title = [title,'(d*X + e*sqrt(X))*dW'];
        case '[1 0 0]'
            title = [title,'(c*dW'];
        case '[1 0 1]'
            title = [title,'(c + e*sqrt(X))*dW'];
        case '[1 1 0]'
            title = [title,'(c + d*X)*dW'];
        case '[1 1 1]'
            title = [title,'(c + d*X + e*sqrt(X))*dW'];
    end

return