function []=emFitAllDDM(uVP)
%Fit the DDM for the Modul Flanker Task
%%
path('/home/afischer/Documents/MATLAB/myfunctions/',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DDM/',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DDM/gwmcmc',path);

LicPresent=0;
while ~LicPresent
    LP = license('checkout','GADS_Toolbox');
    if LP 
        LicPresent = 1;
        disp('License has been checked out...')
    else
        disp('Waiting for license...')
        pause(7)
    end;
end;

% uVP = 1;
load('data/All Behavior.mat')
AD.rt = DDM_prune_RT(AD.rt, 0.02, [80 1000]);
SetupGlobal(uVP,AD,0);


global VP

%%
%define the likelihood function
VP.fitfunction = @DDM_llf;
VP.fitmode     = 3; %1 = fit to individual quantiles, 2 = fit to all VP equal sized bins, 3 = fit to all VP quantiles
VP.CutLog=1e-6;
VP.FitMethod = 2; %1 = chi2, 2 = log L
VP.FitDirection = 2; 
VP.Simulation=0;
%get function handle
DDM = @fDDM;
sRT = sort(VP.RT);

%     DriftRate = 0.5;            %1 (v)    - speed of evidence accumulation for the indicated response (moves the mean of the deviation to derive noise from)
%     DRvariance = 0.2;           %2 (sv)   - variance of the drift rate distribution from trial to trial
%     Boundary = 0.15;            %3 (±a)   - in DDM terminology
%     StartVar = 0.2;             %4 (sz)   - starting point variability, reflects the range around start value as borders of a uniform distribution. 0.4 would mean that start values can be between ±0.2.
%     NondecisionTime = 120;      %5 (Ter)  - time in ms that visual information need to arrive at our decision maker
%     st = 10;                    %6 (st)   - variance of non-decision time (0 = no variance here)
%     CollapseDegree = inf;%2     %7 Boun

LowRT = round(max(sRT(1:round(VP.nTrials*0.05))))/1000;
VP.HyperPriors  = [0 10;          0 2;     0 1;      0.05 LowRT]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
VP.ParamNames   = {'DriftRate'  'sv'     'StartVar' 'Ter'       }; %   'st''Collapse'
VP.Params       = [1             2           4           5       ]; % 6  this maps input to model order of factors 
VP.OnlyCongruent = 0;
VP.GenConstrains = @(v,sv,a,sa,ter,ster,col,favl) (v>sv & a>sa/2); %we dont want more drift rate variance than the actual drift rate and prevent that boundaries are crossed before trial starts

PriorString = '';
for c = 1 : length(VP.Params) %formulate hard priors (boundaries for the model)
    PriorString = [PriorString '(DDM(' num2str(c) ')>' num2str(VP.HyperPriors(c,1)) ')&&(DDM(' num2str(c) ')<' num2str(VP.HyperPriors(c,2)) ') && '];
end
PriorString(end-3:end)='';
DDMprior =@(DDM) eval(PriorString);
%%
%generate other priors
VP.prior_present = []; %indicate here the factors of the model for which a prior exists (these map onto the input parameters only!)
VP.prior_distributions = [makedist('normal', 0.3, 0.01)];
%normalized functions
% VP.prior_functions = {@(x) (pdf(VP.prior_distributions(2),x)/pdf(VP.prior_distributions(2),VP.prior_distributions(2).mu))};
% %actual PDF
VP.prior_functions = {@(x) (pdf(VP.prior_distributions(1),x))};
%      @(x) (pdf(VP.prior_distributions(2),x))};

DriftRate = 0.8;            %1 (v)    - speed of evidence accumulation for the indicated response (moves the mean of the deviation to derive noise from)
DRvariance = 0.45;          %2 (sv)   - variance of the drift rate distribution from trial to trial
Boundary = 0.3;             %3 (±a)   - in DDM terminology
StartVar = 0.2;             %4 (sz)   - starting point variability, reflects the range around start value as borders of a uniform distribution. 0.4 would mean that start values can be between ±0.2.
NondecisionTime = 0.1;      %5 (Ter)  - time in ms that visual information need to arrive at our decision maker
st = 0.1;                   %6 (st)   - variance of non-decision time (0 = no variance here)
collapse = 0.8; 
init = [DriftRate DRvariance StartVar NondecisionTime ]';
%when run first time, clear local state of persistent variables in fDDM and the likelihood function
clear fDDM VP.fitfunction
% VP.weight_congruent_errors=10
%get ML fit
VP.Maximize = -1; %minimze with fmincon
init2=fmincon(DDM,init,[],[],[],[],VP.HyperPriors(:,1),VP.HyperPriors(:,2));

%%
%particle swarm method
% Mname = 'MLE_Model01';
% VP.PD = ParticleDistributions(Mname);
% options = optimoptions('particleswarm', 'Display', 'iter','SwarmSize',100,'CreationFcn',@swarm_priors);
VP.Maximize = -1; %minimze with fmincon
tic        
options = optimoptions('particleswarm', 'Display', 'iter','SwarmSize',75);
% clear fDDM VP.fitfunction
[x,fval] = particleswarm(DDM,length(VP.Params),VP.HyperPriors(:,1),VP.HyperPriors(:,2),options)
ALL.runtime_minutes = toc/60;

%%
%evaluate both fits
for c = 1 : 100
    LL_simplex(c) = DDM(init2);
    LL_particle(c) = DDM(x);
end
ALL.BIC_simplex = -2*(mean(LL_simplex)/VP.Inflation) + length(x)*log(VP.nTrials/VP.Inflation);
ALL.BIC_particle = -2*(mean(LL_particle)/VP.Inflation) + length(x)*log(VP.nTrials/VP.Inflation);
ALL.MLE_simplexParams = init2;
ALL.MLE_particleParams = x;

ALL.ParamNames = VP.ParamNames;
CPU = cpuinfo;
ALL.CPU = CPU.Name;
ALL.Clock = CPU.Clock;

VP.congr=[]; VP.flanker=[]; VP.target=[]; VP.error=[]; VP.RT=[]; VP.resp=[];
VP.weibull=[]; VP.prior_distributions=[]; VP.prior_functions=[];VP.gompertz=[];
VP.fitfunction=[]; VP.GenConstrains=[];
ALL.VP = VP;
save(['/home/data/scratch/adrian/roger/flanker/DDM/MLE_Model03/VP' num2str(uVP) '.mat'], 'ALL' )
return;
