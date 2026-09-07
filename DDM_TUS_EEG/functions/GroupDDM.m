function []=GroupDDM(uVP,FitMethod)
%Fit the DDM for the Modul Flanker Task
%%
if FitMethod==1 %test if GAD license is present when using particle swarm
    NeedLicense = 'statistics_toolbox'; FitString='DiffEvol';
else
    NeedLicense = 'GADS_Toolbox'; FitString='PartSwarm';
end
    
LicPresent=0;
while ~LicPresent
    LP = license('checkout',NeedLicense);
    if LP 
        LicPresent = 1;
        disp([NeedLicense ' has been checked out...'])
    else
        disp(['Waiting for ' NeedLicense ' license...'])
        pause(7)
    end;
end;

path('/home/afischer/Documents/MATLAB/myfunctions/',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DDM/',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DDM/gwmcmc',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DiffEvol',path);

pause(0.001*uVP); %avoid identical start times for random number generators
%%
if ~exist('uVP')
    uVP=1;
end
load('data/All Behavior.mat')
AD.rt = AD.rt + 83;
AD.rt = DDM_prune_RT(AD.rt, 0.02, [130 1000]);
global VP;
VP.Inflation=5;
SetupGlobal(1,AD,0);

%determine which model to run
n_per_mod = 50;
ModelSelect = ceil(uVP/50);
n_model = uVP-n_per_mod*(ModelSelect-1);

disp(['Running Model no ' num2str(ModelSelect) ' in iteration no ' num2str(n_model)])

%define the likelihood function
VP.fitfunction = @DDM_llf;
VP.fitmode     = 2; %1 = fit to all VP equal sized bins, 2 = fit to all VP quantiles, 3 = fit to individual quantiles, 
VP.FitGroup    = 1; %1 = fit group, 0 = fit VP
VP.CutLog=1e-6;
VP.FitMethod = 2; %1 = chi2, 2 = approximate log L
VP.FitDirection = 2; 
VP.Simulation=0;
VP.prior_present = []; %indicate here the factors of the model for which a prior exists (these map onto the input parameters only!)
VP.Maximize = -1; %minimze with fmincon
VP.eta = 0.1;

%get function handle
DDM = @fDDM;
sRT = sort(VP.RT);

LowRT = round(min(sRT(1:round(VP.nTrials*0.05))))/1000;
VP.OnlyCongruent = 0;
VP.GenConstrains = @(v,sv,a,sa,ter,ster,col,favl) (a>sa*0.8 & v>sv & ter>ster/20); %& a>sa/2 we dont want more drift rate variance than the actual drift rate and prevent that boundaries are crossed before trial starts

PR.v = [0 2.5];
PR.sv = [0 1];
PR.a = [0.01 0.45];
PR.sa = [0.05 0.3];
PR.ter = [0.1 0.4]; %100 to 400 ms
PR.st = [0 2]; %0 to 100 ms
PR.k = [10e-4 2]; 
PR.f = [0.5 1.6];
StandVal = [1 0.5 0.22 0.11 0.25 0.8 1.1 1];
StandQua = [0 0 0 0 0.001 0.01 0 0]; %quantization for time (to ms, the unit of DDM here) speeds up convergence
%%

if ModelSelect == 1
    VP.HyperPriors  = [PR.v;        PR.sv; PR.a;  PR.sa;    PR.ter; PR.st]; %fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'  'sv'  'a'     'StartVar' 'Ter'   'st'}; %
    VP.Params       = [1             2     3      4           5       6]; %this maps input to model order of factors 
elseif ModelSelect == 2
    %Model with variable collapsing boundaries
    VP.HyperPriors  = [PR.v;        PR.sv; PR.a;  PR.sa;    PR.ter; PR.st; PR.k];
    VP.ParamNames   = {'DriftRate'  'sv'  'a'     'StartVar' 'Ter'   'st'  'k'};
    VP.Params       = [1             2     3      4           5       6     7];
elseif ModelSelect == 3
    %Model with variable flanker value
    VP.HyperPriors  = [PR.v;        PR.sv; PR.a;  PR.sa;    PR.ter; PR.st; PR.f];
    VP.ParamNames   = {'DriftRate'  'sv'  'a'     'StartVar' 'Ter'   'st'   'f'};
    VP.Params       = [1             2     3      4           5       6      8];
elseif ModelSelect == 4
    %Model with variable flanker value and collapsing boundaries
    VP.HyperPriors  = [PR.v;        PR.sv; PR.a;  PR.sa;    PR.ter; PR.st; PR.k; PR.f];
    VP.ParamNames   = {'DriftRate'  'sv'  'a'     'StartVar' 'Ter'   'st'  'k'   'f'};
    VP.Params       = [1             2     3      4           5       6     7     8];
elseif  ModelSelect == 5
    %Model with variable flanker value and urgency signal
    VP.UrgencyOn = 1;
    VP.HyperPriors  = [PR.v;        PR.sv; PR.a;  PR.sa;    PR.ter; PR.st; PR.k; PR.f];
    VP.ParamNames   = {'DriftRate'  'sv'  'a'     'StartVar' 'Ter'   'st'  'k'   'f'};
    VP.Params       = [1             2     3      4           5       6     7     8];
end

clear fDDM VP.fitfunction
tic
if FitMethod == 1
    % define parameter names, ranges, quantizations and initial values :
    paramDefCell = {'', VP.HyperPriors, StandQua(VP.Params)', StandVal(VP.Params)'};
    objFctSettings={};
    % get default DE parameters
    DEParams = getdefaultparams;
    DEParams.NP = length(VP.Params)*10;

    % use a subfunction to check parameter vectors for validity
    DEParams.validChkHandle = [];%@demo3_constraint;
    DEParams.saveHistory = 0;
    
    % set times
    DEParams.maxiter  = 4000;
    DEParams.maxtime  = 5*3600;  % in seconds
    DEParams.maxclock = [];

    % set display options
    DEParams.infoIterations = 1;
    DEParams.infoPeriod     = 60;  % in seconds

    optimInfo.title = ['Model ' num2str(ModelSelect) ' with ' num2str(length(VP.Params)) ' free parameters'];
    [x, fval] = differentialevolution(DEParams, paramDefCell, DDM, [], [], [], optimInfo)

%%
elseif FitMethod == 2
    %particle swarm method
    rng('shuffle')
    options = optimoptions('particleswarm', 'Display', 'iter','SwarmSize',length(VP.Params)*10);
    [x,fval] = particleswarm(DDM,length(VP.Params),VP.HyperPriors(:,1),VP.HyperPriors(:,2),options)
end

ALL.runtime_minutes = toc/60;
ALL.BestSwarm = fval;
ALL.BestSwarmBIC = -2*(mean(fval)) + length(VP.Params)*log(VP.nTrials)
ALL.HyperPrior = VP.HyperPriors;
ALL.ParamNames = VP.ParamNames;
ALL.Params = VP.Params;

%evulate function at best value 100x
for c = 1 : 100
    LL_particle(c) = DDM(x);
end
ALL.Mean_LL = mean(LL_particle);
ALL.BIC_particle = -2*(mean(LL_particle)) + length(x)*log(VP.nTrials);
ALL.MLE_particleParams = x;
disp(['Saving to : Model fit method ' FitString ' model no ' num2str(ModelSelect) '-' num2str(n_model) '.mat'])
save(['/home/data/scratch/adrian/roger/flanker/DDM/GroupModels/Larger Model fit method ' FitString ' model no ' num2str(ModelSelect) '-' num2str(n_model) '.mat'], 'ALL' )
return
