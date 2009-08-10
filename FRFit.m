function FRFit(jobNameBase,numberOfJobs,plotsOn,switchOver,x_label,y_label)
global model n i

% jobNameBase = '2sim_1on';
% numberOfJobs = 10;
% plotsOn = 1;
% switchOver = 1000;

% Declare settings:
maxEval=5000;
maxIter=5000;
%switchOver=10000; %MS when input turns on

% Upload data:
for n=1:numberOfJobs
    load([jobNameBase,'_',num2str(n),'.mat']);
    FRES1(n,:)=FRE1;
    FRES2(n,:)=FRE2;
    plotTMat(n,:)=plotT;
end

% Create a design matrix (dm)
A = [0 0; 1 0; 1 1];
B = [1 0 0; 1 1 0; 1 0 1; 1 1 1];

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
        model(n).mu{i} = @(X,Y,theta,dv) theta(1).*dv(1) + theta(2).*dv(2).*X;
        model(n).g{i} = @(X,Y,theta,dv) abs(theta(3).*dv(3) + theta(4).*dv(4).*Y + theta(5).*dv(5).*sqrt(abs(Y)));
        model(n).IC{i} = initial_condition;
        model(n).theta{i} = logical(dm(i,:));
        model(n).description{i} = modelDescription(dm(i,:));
        model(n).likelihoodPre{i} = NaN;
        model(n).likelihoodPost{i} = NaN;
        model(n).AICPre(i,1) = NaN;
        model(n).AICPost(i,1) = NaN;
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
    n
    model(n).IC{i}
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
        i
        dv = dm(i,:);
        % Develop f:
        model(n).f{i}=@(xNext,xCurr,yCurr,theta,dv) normpdf(xNext,xCurr+model(n).mu{i}(xCurr,yCurr,theta,dv)*dt,abs(model(n).g{i}(xCurr,yCurr,theta,dv)*sqrt(dt)));

        % Develop Log-Likelyhood function, using the looping method to help:
        model(n).L{i}=@(theta) likelyhoodFunction(theta,model(n).f{i},xCurr,yCurr,dv);
        
        % Do the search:
        [temp1,temp2,temp3]=fminsearch(model(n).L{i},model(n).IC{i},options);
        model(n).thetaEstPre{i}=temp1;
        model(n).likelihoodPre{i}=temp2;
        exitflag=temp3;
        
        % Define the approximate log-likelyhood functions for computing AIC:
        % u=@(t,x,y,theta) -(1/2)*log(2*pi*t)-log(model(n).g{i}(y,theta))-S(x,y,theta).^2./(2*t)+H(x,y,theta)+t*gTilde(x,y,theta);
        % un=@(theta) sum(u(dt,FRES1(n,1:((onIndex)-1))-FRES2(n,1:((onIndex)-1)),FRES1(n,2:(onIndex))-FRES2(n,2:(onIndex)),theta));
        
        thetaVector = nan(1,5);
        
        % Determine the Akaiki Information for each model, provided convergence:
        if exitflag==1
            % model(n).likelihoodPre{i}=un(model(n).thetaEstPre{i});
            model(n).AICPre(i,1)=2*model(n).likelihoodPre{i}+2*length(model(n).thetaEstPre{i});
            thetaVector(logical(dm(i,:))) = model(n).thetaEstPre{i}(logical(dm(i,:)));
        else
            [model(n).thetaEstPre{i},model(n).likelihoodPre{i},exitflag]=fminunc(model(n).L{i},model(n).IC{i},options);
            if exitflag==1
                % model(n).likelihoodPre{i}=un(model(n).thetaEstPre{i});
                model(n).AICPre(i,1)=2*model(n).likelihoodPre{i}+2*length(model(n).thetaEstPre{i});
                thetaVector(logical(dm(i,:))) = model(n).thetaEstPre{i}(logical(dm(i,:)));
            else
                model(n).AICPre(i,1)=Inf;
                model(n).likelihoodPre{i}=-Inf;
            end
        end
        thetaPreMatrix(i,:,n) = thetaVector;
        % thetaPreMatrix(i,:,n) = model(n).thetaEstPre{i};
    end
end
optimalModelCountPre=zeros(1,length(model(n).description));
for n=1:numberOfJobs
    AICPreMatrix(n,:) = model(n).AICPre';
    [val,ind]=min(model(n).AICPre);
    optimalModelCountPre(ind)=optimalModelCountPre(ind)+1;
end

%% Post-stimulus epoch

% Minimize the negative log-likelyhood for each model, selective1-selective2, after input:
for n=1:numberOfJobs
    n
    model(n).IC{i}
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
        i
        dv = dm(i,:);
        % Develop f:
        model(n).f{i}=@(xNext,xCurr,yCurr,theta,dv) normpdf(xNext,xCurr+model(n).mu{i}(xCurr,yCurr,theta,dv)*dt,abs(model(n).g{i}(xCurr,yCurr,theta,dv)*sqrt(dt)));

        % Develop Log-Likelyhood function, using the looping method to help:
        model(n).L{i}=@(theta) likelyhoodFunction(theta,model(n).f{i},xCurr,yCurr,dv);
        
        % Do the search:
        [temp1,temp2,temp3]=fminsearch(model(n).L{i},model(n).IC{i},options);
        model(n).thetaEstPost{i}=temp1;
        model(n).likelihoodPost{i}=temp2;
        exitflag=temp3;
        
        % Define the approximate log-likelyhood functions for computing AIC:
        % u=@(t,x,y,theta) -(1/2)*log(2*pi*t)-log(model(n).g{i}(y,theta))-S(x,y,theta).^2./(2*t)+H(x,y,theta)+t*gTilde(x,y,theta);
        % un=@(theta) sum(u(dt,FRES1(n,1:((onIndex)-1))-FRES2(n,1:((onIndex)-1)),FRES1(n,2:(onIndex))-FRES2(n,2:(onIndex)),theta));
        
        thetaVector = nan(1,5);
        
        % Determine the Akaiki Information for each model, provided convergence:
        if exitflag==1
            % model(n).likelihoodPost{i}=un(model(n).thetaEstPost{i});
            model(n).AICPost(i,1)=2*model(n).likelihoodPost{i}+2*length(model(n).thetaEstPost{i});
            thetaVector(logical(dm(i,:))) = model(n).thetaEstPost{i}(logical(dm(i,:)));
        else
            [model(n).thetaEstPost{i},model(n).likelihoodPost{i},exitflag]=fminunc(model(n).L{i},model(n).IC{i},options);
            if exitflag==1
                % model(n).likelihoodPost{i}=un(model(n).thetaEstPost{i});
                model(n).AICPost(i,1)=2*model(n).likelihoodPost{i}+2*length(model(n).thetaEstPost{i});
                thetaVector(logical(dm(i,:))) = model(n).thetaEstPost{i}(logical(dm(i,:)));
            else
                model(n).AICPost(i,1)=Inf;
                model(n).likelihoodPost{i}=-Inf;
            end
        end
        thetaPostMatrix(i,:,n) = thetaVector;
        % thetaPostMatrix(i,:,n) = model(n).thetaEstPost{i};
    end
end
optimalModelCountPost=zeros(1,length(model(n).description));
for n=1:numberOfJobs
    AICPostMatrix(n,:) = model(n).AICPost';
    [val,ind]=min(model(n).AICPost);
    optimalModelCountPost(ind)=optimalModelCountPost(ind)+1;
end

%% Plotting figures

warning on all

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
    barh(optimalModelCountPre)
    set(gca,'YTickLabel',model(1).description)
    title('AIC Minimization Histogram, Pre-Stimulus')
    
    figure(3)
    clf
    barh(optimalModelCountPost)
    set(gca,'YTickLabel',model(1).description)
    title('AIC Minimization Histogram, Post-Stimulus')

    % Save figures
    saveas(1,['./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_MultiFR.eps'],'eps')
    saveas(1,['./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_MultiFR.fig'],'fig')
    unix(['open ./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_MultiFR.eps']);
    
    saveas(2,['./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_AICPre.eps'],'eps')
    saveas(2,['./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_AICPre.fig'],'fig')
    unix(['open ./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_AICPre.eps']);
    
    saveas(3,['./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_AICPost.eps'],'eps')
    saveas(3,['./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_AICPost.fig'],'fig')
    unix(['open ./savedFigures/',jobNameBase,'_',x_label,'_',y_label,'_AICPost.eps']);
    
    save([jobNameBase,'_',x_label,'_',y_label,'_results.mat'],'AICPreMatrix','AICPostMatrix','thetaPreMatrix','thetaPostMatrix');
end 

return

%%
function SVal=S(x,y,theta)
global model n i
    SFun=@(arg) 1./model(n).g{i}(arg,theta);
    for k=1:length(x)
        if x(k)==y(k)
            SVal(k)=0;
        else
            SVal(k)=quad(SFun,x(k),y(k));
        end
    end
return

%%
function HVal=H(x,y,theta)
global model n i
    HFun=@(arg) B(arg,theta)./model(n).g{i}(arg,theta);
    for k=1:length(x)
        if x(k)==y(k)
            HVal(k)=0;
        else
            HVal(k)=quad(HFun,x(k),y(k));
        end
    end
return

%%
function gTildeVal=gTilde(x,y,theta)
    gTildeVal=-(1/2)*(C(x,theta)+C(y,theta)+(1/3)*B(x,theta).*B(y,theta));
return

%%
function CVal=C(x,theta)
global model n i
    CVal=(1/2)*B(x,theta).^2 + (1/2)*Bx(x,theta).*model(n).g{i}(x,theta);
return

%%
function BVal=B(x,theta)
global model n i
    BVal=model(n).mu{i}(x,theta)./model(n).g{i}(x,theta)-(1/2)*model(n).gx{i}(x,theta);
return

%%
function BxVal=Bx(x,theta)
global model n i
    BxVal=model(n).mux{i}(x,theta)./model(n).g{i}(x,theta)-model(n).mu{i}(x,theta).*model(n).gx{i}(x,theta)./model(n).g{i}(x,theta).^2-(1/2)*model(n).gxx{i}(x,theta);
return
    