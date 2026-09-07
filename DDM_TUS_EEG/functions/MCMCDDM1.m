function []=MCMCDDM1(uVP)
%%
% uVP = 7;
path('/home/afischer/Documents/MATLAB/myfunctions/',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DDM/',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DDM/gwmcmc',path);

load('data/All Behavior.mat')
SetupGlobal(uVP,AD,0);
global VP

%%
VP.fitfunction = @DDM_llf;
VP.fitmode     = 3; %1 = fit to individual quantiles, 2 = fit to all VP equal sized bins, 3 = fit to all VP quantiles
VP.CutLog=1e-5;
VP.FitMethod = 1; %1 = chi2, 2 = log L

%get function handle
DDM = @fDDM;

%     DriftRate = 0.5;            %1 (v)    - speed of evidence accumulation for the indicated response (moves the mean of the deviation to derive noise from)
%     DRvariance = 0.2;           %2 (sv)   - variance of the drift rate distribution from trial to trial
%     Boundary = 0.15;            %3 (±a)   - in DDM terminology
%     StartVar = 0.2;             %4 (sz)   - starting point variability, reflects the range around start value as borders of a uniform distribution. 0.4 would mean that start values can be between ±0.2.
%     NondecisionTime = 120;      %5 (Ter)  - time in ms that visual information need to arrive at our decision maker
%     st = 10;                    %6 (st)   - variance of non-decision time (0 = no variance here)
%     CollapseDegree = inf;%2     %7 Boun

VP.HyperPriors  = [0 5;         0 1;           0.01 1;     0 1;      0.05 0.4;  0 0.3     ]; %;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
VP.ParamNames   = {'DriftRate' 'DRvariance'   'Boundary'  'StartVar' 'Ter'      'st' }; %'Collapse'
VP.Params       = [1           2               3           4           5        6      ]; %this maps input to model order of factors 
VP.OnlyCongruent = 0;

PriorSting = '';
for c = 1 : length(VP.Params) %formulate hard priors (hyperpriors)
    PriorSting = [PriorSting '(DDM(' num2str(c) ')>' num2str(VP.HyperPriors(c,1)) ')&&(DDM(' num2str(c) ')<' num2str(VP.HyperPriors(c,2)) ') && '];
end
PriorSting(end-3:end)='';
DDMprior =@(DDM) eval(PriorSting);
%%
%generate other priors
VP.prior_present = []; %indicate here the factors of the model for which a prior exists (these map onto the input parameters only!)
VP.prior_distributions = [makedist('normal', 0.17, 0.17)];
%normalized functions
VP.prior_functions = {@(x) (pdf(VP.prior_distributions(2),x)/pdf(VP.prior_distributions(2),VP.prior_distributions(2).mu))};
% %actual PDF
% VP.prior_functions = {@(x) (pdf(VP.prior_distributions(1),x))
%      @(x) (pdf(VP.prior_distributions(2),x))};

DriftRate = 0.8;            %1 (v)    - speed of evidence accumulation for the indicated response (moves the mean of the deviation to derive noise from)
DRvariance = 0.2;           %2 (sv)   - variance of the drift rate distribution from trial to trial
Boundary = 0.15;            %3 (±a)   - in DDM terminology
StartVar = 0.2;             %4 (sz)   - starting point variability, reflects the range around start value as borders of a uniform distribution. 0.4 would mean that start values can be between ±0.2.
NondecisionTime = 0.15;     %5 (Ter)  - time in ms that visual information need to arrive at our decision maker
st = 0.1;                   %6 (st)   - variance of non-decision time (0 = no variance here)

init = [DriftRate DRvariance Boundary StartVar NondecisionTime st]';
%when run first time, clear local state of persistent variables in fDDM and the likelihood function
clear fDDM VP.fitfunction
% VP.weight_congruent_errors=10
%get ML fit
VP.Maximize = -1; %minimze with fmincon
init=fmincon(DDM,init,[],[],[],[],VP.HyperPriors(:,1),VP.HyperPriors(:,2));

[StartLL, D] = DDM(init)

%%
VP.Maximize = 1;
%MCMC walkers should be initialized somewhere not too far away from each other
n_walker = length(init)*10; %specifiy the number of walkers (rule of thumb: n x parameters)
n_params = length(init);
initScale = init./4;

%ensure that all parameters are within the limits set by the hyperpriors
lastNchar = 0;
%initialize walkers in gaussian blob around start values
init_walker=bsxfun(@plus,init,randn(n_params,n_walker).*repmat(initScale,1,n_walker)); %start with an initial set of walkers
for c = 1 : n_params %confine to hyperpriors
    init_walker(c,init_walker(c,:)<=VP.HyperPriors(c,1)) = VP.HyperPriors(c,1)+1e-3; %substitute too low values
    init_walker(c,init_walker(c,:)>=VP.HyperPriors(c,2)) = VP.HyperPriors(c,2)-1e-3; %substitute too high values
end
for w = 1 : n_walker
    progressmsg=sprintf('\nOptimizing walker no: %3.0f%',w);
    fprintf('%s%s',repmat(char(8),1,lastNchar),progressmsg); drawnow;
    lastNchar=length(progressmsg);
    ThisWalkerCheck = 0; n_check = 0;
    while ~ThisWalkerCheck
        n_check = n_check +1;
        init_walker(:,w);
        try
            [TestLL] = DDM(init_walker(:,w));
        catch
            TestLL = -inf;
        end
        if TestLL==-inf || ~DDMprior(init_walker(:,w))
            use_scale = initScale / (1 + n_check / 100); %slowly reduce scaling parameter over time to ensure that a finite likelihood is found for each walker
            init_walker(:,w) = bsxfun(@plus,init,randn(n_params,1).*repmat(use_scale,1,1)); %draw a new start point for this walker
        else
            ThisWalkerCheck = 1; %accept this walker
        end
    end
end
disp(' ')
for c = 1 : n_params %display: range of starting walkers
    disp(['Parameter no ' num2str(c) ' ' VP.ParamNames{c} ' min = ' num2str(min(init_walker(c,:))) ' max = ' num2str(max(init_walker(c,:)))])
end

%%
%try out the intital fit
% tic
% VP.Maximize = 1;
% L1 = DDM(init);
% L2 = DDM(ALL.modeParams);
% toc
% L1
% L2
% %%
% for c = 1 : n_walker
%     disp(['LL = ' num2str(DDM(init_walker(:,c))) ' and params: ' num2str(init_walker(:,c)')])
% end

%%
%drop the MCMC Hammer
VP.Maximize = 1;
logP = [];
while isempty(logP)
    try
        tic
        [m,logP]=gwmcmc(init_walker,{DDMprior DDM},75000,'ThinChain',5,'burnin',.1,'StepSize',2,'ProgressBar',true);
        ALL.runtime_minutes = toc/60;
    catch
        disp('Start positions note finite, retry...')
    end
end

m2=m(:,:)'; %reshape matrix to collapse the ensemble member dimension
logP2=logP(:,:)';
logP2=logP2(:,2);

%estimate density functions for parameters and their 2D projections
[~,~,p.ess]=eacorr(m);
p.ess=mean(p.ess);
p.support=VP.HyperPriors; 
M=size(m2,2);
Np=size(m2,1);
%%
close all; figure
for r=1:M
    for c=1:max(r,M)
        if c==r
            [F,X,bw]=ksdensity(m2(:,r),linspace(p.support(r,1),p.support(r,2),500),'support',p.support(r,:)); %TODO: use ESS 
            if p.ess<Np
                [F,X,bw]=ksdensity(m2(:,r),linspace(p.support(r,1),p.support(r,2),500),'width',bw*(Np/p.ess)^.2,'support',p.support(r,:)); %(the power 1/5 comes from examining the bandwidth calculation in ksdensity)
            end
            ALL.ProbabilityMass(:,c) = F;
            ALLsame.ProbabilityMassX(:,c) = X;
            [~,y(r)]=max(F);
        end
    end
    subplot(1, M, r)
    plot(X,F)
end
 
ALL.modeParams = diag(ALLsame.ProbabilityMassX(y,:))';
ALL.meanParams = mean(m2(:,:));
ALL.medianParams = median(m2(:,:));

[~,y]=max(logP2);
for c = 1 : 100
    ll(1,c)=DDM(ALL.modeParams);
    ll(2,c)=DDM(ALL.meanParams);
    ll(3,c)=DDM(ALL.medianParams);
    ll(4,c)=DDM(init);
end
[~,D]=DDM(ALL.modeParams);
D.logl = mean(ll);
DDM_plot(D)
figure
ecornerplot(m,'ks',true,'color',[.6 .35 .3],'names',VP.ParamNames)

ALL.StartParams = D.params;
ALL.m = m;
ALL.logP = logP;
ALL.MCMC_Mode = mean(ll(1,:));
ALL.MCMC_Mean = mean(ll(2,:));
ALL.MCMC_Median = mean(ll(3,:));
ALL.ML = mean(ll(4,:));

ALL.ParamNames = VP.ParamNames;
CPU = cpuinfo;
ALL.CPU = CPU.Name;
ALL.Clock = CPU.Clock;

% save(['/home/data/scratch/adrian/roger/flanker/DDM/Model01/VP' num2str(uVP) '.mat'], 'ALL' )
ALL.Countour2D_N=[];
VP.congr=[]; VP.flanker=[]; VP.target=[]; VP.error=[]; VP.RT=[]; VP.resp=[];
VP.weibull=[]; VP.prior_distributions=[]; VP.prior_functions=[];VP.gompertz=[];
ALL.VP = VP;
save(['/home/data/scratch/adrian/roger/flanker/DDM/Model01_slim/VP' num2str(uVP) '.mat'], 'ALL' )
if uVP == 1
    save('/home/data/scratch/adrian/roger/flanker/DDM/Model01_slim/all same.mat', 'ALLsame' )
end
return;








