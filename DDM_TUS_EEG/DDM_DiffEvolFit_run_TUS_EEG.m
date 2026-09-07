function []=DDM_DiffEvolFit_run_TUS_EEG(uVP, ModelSelect,TrialSelect)
%Fit the DDM for the Modul Flanker Task
% uVP           = which participant 
% ModelSelect   = which model to run 
% TrialSelect   = on which trials will the modeling be run

%% to-do
% fix variance parameters 
% collaps models with no variance in drift rate

%%

path('functions/',path);
path('functions/helpers/',path);
path('functions/DDM/',path);
path('functions/DDM/gwmcmc',path);
path('functions/DiffEvol',path);

%%

datapath = 'data/TUS_EEG/';

F = dir ([datapath '*.mat']);

%%
TrialString = {'all' 'all_ce_rem' 'post_cor' 'post_err' 'all_ce_rem_1st_third'};

n_models_to_run = 1;

if ~exist('ModelSelect');ModelSelect=3;end
if ~exist('TrialSelect'); TrialSelect =1; end
if ~exist('uVP'); uVP =1; end

SaveName = ['results/DDM' num2str(ModelSelect) '_TUS-EEG_' TrialString{TrialSelect+1} '_'  F(uVP).name(1:7) '.mat'];
if exist(SaveName)
    disp('Model exists...')
    return;
end

% load data set

% load('data/All Behavior4.mat')
load([datapath F(uVP).name]);

% this is how the mapping should look like
% AD.error —> 1 = correct; 2 = error
% AD.congr —> 1 = congruent; 2 = incongruent 
% AD.resp —> 1 = left; 2 = right 
% AD.Target —> 1 = left; 2 = right 
% AD.Flanker —> 1 = left; 2 = right 

% logg info Camilia
% 1 = left incongruent (expected left button press)
% 2 = left congruent
% 3 = right incongruent (expected right button press)
% 4 = right congruent
% 
% congruence: congruent = 1 and incongruent = -1,
% Resp: 1 = right and -1 = left

AD.congr = [logg.Congruence]';
AD.congr(AD.congr==-1)=2;
AD.congr(AD.congr==1)=1;

AD.error = [logg.ACC]'; 
AD.error(AD.error==0)=2;

AD.rt = [logg.RT]'*1000; 

AD.resp = [logg.Resp]'; 
AD.resp(AD.resp==1)=2;
AD.resp(AD.resp==-1)=1;

StimInfo = [logg.Trial_code]'; 

AD.Target  = ones(length(StimInfo),1); 
AD.Target(StimInfo>2)=2;


AD.Flanker =  ones(length(StimInfo),1);
AD.Flanker(StimInfo==1 | StimInfo==4)=2;

AD.rt = DDM_prune_RT(AD.rt, 0.02, [80 1100]);

posterror = nan(length(StimInfo),1);
posterror(2:end) = AD.error(1:end-1)==2;


global VP;

VP.Inflation=5;
VP.quant_space = 0.2;

n_PostErrorVP = sum(posterror==1 & ~(AD.congr==1 & AD.error==2)); 
n_PostCorrVP = sum(posterror==0 & ~(AD.congr==1 & AD.error==2)); 

MinErrorTrial = 30;
%TrialSelect = 3;
if TrialSelect == 1 %Remove all congruent errors
    AD.rt(AD.congr==1 & AD.error==2) = nan;
    AD.rt(isnan(posterror)) = nan;
elseif TrialSelect == 2 %remove congruent errors and use post correct trials only
    AD.rt(AD.congr==1 & AD.error==2) = nan;
    AD.rt(posterror~=0) = nan;
    VP.Inflation=round(5000/n_PostCorrVP);
    VP.quant_space = 0.2;
    if n_PostErrorVP < MinErrorTrial; disp(['Not enough total trials: ' num2str(n_PostErrorVP)]);return; end
elseif TrialSelect == 3 %remove congruent errors and use post error trials only
    AD.rt(AD.congr==1 & AD.error==2) = nan;
    AD.rt(posterror~=1) = nan; 
    VP.Inflation=round(5000/n_PostErrorVP);
    VP.quant_space = 0.2;
    if n_PostErrorVP < MinErrorTrial; disp(['Not enough total trials: ' num2str(n_PostErrorVP)]);return; end
elseif TrialSelect == 4 %remove congruent errors and use post error trials only
    AD.rt(AD.congr==1 & AD.error==2) = nan;
    halfTime = round(length(AD.rt)/3)+1;
    AD.rt(halfTime:end) = nan;
else
    AD.rt(isnan(posterror)) = nan;
end

disp(['Selecting trials ' TrialString{TrialSelect+1}])

SetupGlobal(1,AD,0);

disp(['Running Model no ' num2str(ModelSelect) ' TUS-EEG for VP ' num2str(uVP) ' with ' num2str(VP.nTrials) ' simulations.'])

%define the likelihood function
VP.fitfunction = @DDM_llf;
VP.fitmode     = 3;     %1 = fit to all VP equal sized bins, 2 = fit to all VP quantiles, 3 = fit to individual quantiles, 
VP.FitGroup    = 0;     %1 = fit group, 0 = fit VP
VP.CutLog=1e-6;         %minimum L per observation (irrelevant if mixture model is used)
VP.FitMethod = 2;       %1 = chi2, 2 = approximate log L
VP.FitDirection = 2; 
VP.Simulation=0;
VP.prior_present = []; %indicate here the factors of the model for which a prior exists (these map onto the input parameters only!)
VP.Maximize = -1; %minimze with fmincon
VP.eta = 0.1;     %noise parameter (here used as scaling factor)
VP.MCMC=0;
VP.MixturePercent = 2; %proportion of trials that can be drawn from the mixture model
VP.MixtureBounds = [min(VP.RT) max(VP.RT)]; %range for the mixture model
VP.UrgencyOn = 0; % if on k is used as a multiplicative gain (urgency signal) instead of boundary collapse



%get function handle
DDM = @fDDM;
sRT = sort(VP.RT);

VP.OnlyCongruent = 0;
VP.GenConstrains = @(v,sv,a,sa,ter,ster,col,favl,sf) (a>sa*0.55); %& a>sa/2 we dont want more drift rate variance than the actual drift rate and prevent that boundaries are crossed before trial starts

tic   
% define hyperparameter bounds
PR.v = [0.1 8.5];       %1 = driftrate
PR.sv = [0 1.5];        %2 = trial-by-trail variance of driftrate
PR.a = [0.1 0.45];      %3 = boundary
PR.sa = [0.05 0.3];     %4 = Start variance
PR.ter = [0.1 0.4];     %5 = non-decision-time (100 to 400 ms) 
PR.st = [0.01 2];       %6 = trial-by-trail variance of non decision time (0 to 100 ms)
PR.k = [10e-4 2.5];     %7 = bounds collapse according to a cumulative weibull distribution (fit the lambda in the weibull distribution)
PR.f = [0.1 1.1];       %8 = flanker surpressor weighting 
PR.sf = [0 1];          %9 = trial-by-trail variance of surpressor weighting

VP.StandVal = [1 0.5 0.22 0.11 0.25 0.8 inf 1 0.45];
VP.FixToGroupMean = 0; VP.GroupTrial = nan;
% StandQua = [0.01 0.001 0.001 0.001 0.001 0.01 0.01 0.001 0.001]; %
StandQua = [0 0 0 0 0.001 0.01 0 0 0]; %quantization for time speeds up convergence
if ModelSelect == 3 %8 parameter DDM
    %Model with variable flanker value and variance
    VP.HyperPriors  = [PR.v;       PR.sv;   PR.a;     PR.sa;       PR.ter;   PR.st;       PR.f;     PR.sf  ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate' 'sv'    'a'        'StartVar'   'Ter'     'st'         'f'       'sf'   }; %   'st''Collapse'
    VP.Params       = [1           2        3         4            5         6            8         9      ];
    GroupBestFit    = [7.3191      0.5461   0.3396     0.1358       0.2336    1.0391       0.3356    0.2759   ]';    

elseif ModelSelect == 4 %weibull model with no variance in sv and sf
    VP.UrgencyOn = 0;
    PR.k = [10e-4 3];
    VP.HyperPriors  = [PR.v;            PR.a;      PR.sa;       PR.ter;   PR.st;  PR.k;   PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'       'a'      'StartVar'     'Ter'     'st'    'k'      'f'    }; %   'st''Collapse'
    VP.Params       = [1                  3          4             5         6      7       8      ];
    VP.StandVal     = [1   0             0.22      0.11           0.25     0.8       inf    1   0  ];
    GroupBestFit    = [4.5542         0.36815     0.18205       0.30014   1.1825  2.2922  0.4637   ]';

elseif ModelSelect == 5 % 8 parameter DDM variance parameters were fixed to the values of the group mean maximumlikelihood
                        % parameters for DDM 3
    VP.HyperPriors  = [PR.v;                PR.a;                  PR.ter;                      PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'         'a'                    'Ter'                         'f'      }; %   'st''Collapse'
    VP.Params       = [1                    3                        5                           8        ];
    VP.FixToGroupMean = 3; VP.GroupTrial = 1;

    PreMLP = readtable('data/DDMParameterReadout.csv');
    PreMLP = table2array(PreMLP);
    disp(['Setting fixed parameters to mean of model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}])
    VP.StandVal =  mean(PreMLP(PreMLP(:,2)== str2double(F(uVP).name(7)),3:10))'; GroupBestFit = VP.StandVal(VP.Params);
    VP.StandVal = [ VP.StandVal(1:6)' inf VP.StandVal(7:8)']'; 
 
elseif ModelSelect == 6 % 8 parameter DDM variance parameters were fixed to the values of the across group mean maximumlikelihood
                        % parameters for DDM 3
    VP.HyperPriors  = [PR.v;                PR.a;                  PR.ter;                      PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'         'a'                    'Ter'                         'f'      }; %   'st''Collapse'
    VP.Params       = [1                    3                        5                           8        ];
    VP.FixToGroupMean = 3; VP.GroupTrial = 1;

    PreMLP = readtable('data/DDMParameterReadout.csv');
    PreMLP = table2array(PreMLP);
    disp(['Setting fixed parameters to mean of model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}])
    VP.StandVal =  mean(PreMLP(:,3:10))'; GroupBestFit = VP.StandVal(VP.Params);
    VP.StandVal = [ VP.StandVal(1:6)' inf VP.StandVal(7:8)']'; 
 
elseif ModelSelect == 7 % 8 parameter DDM variance parameters were fixed to the values of the across group mean maximumlikelihood
                        % parameters for DDM 3
    VP.HyperPriors  = [PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'f'      }; %   'st''Collapse'
    VP.Params       = [8        ];
    VP.FixToGroupMean = 3; VP.GroupTrial = 1;

    PreMLP = readtable('data/DDMParameterReadout.csv');
    PreMLP = table2array(PreMLP);
    disp(['Setting fixed parameters to mean of model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}])
    VP.StandVal =  mean(PreMLP(:,3:10))'; GroupBestFit = VP.StandVal(VP.Params);
    VP.StandVal = [ VP.StandVal(1:6)' inf VP.StandVal(7:8)']'; 
 
end
VP
%%
clear fDDM VP.fitfunction
tic;

% define parameter names, ranges, quantizations and initial values :
paramDefCell = {'', VP.HyperPriors, StandQua(VP.Params)', GroupBestFit};
objFctSettings={};VP.fitfunction

% get default DE parameters
DEParams = getdefaultparams;
DEParams.NP = length(VP.Params)*10;

% use a subfunction to check parameter vectors for validity
DEParams.validChkHandle = [];%@demo3_constraint;
DEParams.saveHistory = 0;

% set times
DEParams.maxiter  = 5000;
DEParams.maxtime  = 7200;  % in seconds - should be at least 7200
DEParams.maxclock = [];

% set display options
DEParams.infoIterations = 1;
DEParams.infoPeriod     = 60;  % in seconds

VP.ExtractSTE       = false; %extract single trial estimates?

optimInfo.title = ['Model ' num2str(ModelSelect) ' with ' num2str(length(VP.Params)) ' free parameters'];
fval = inf;
for c = 1 : n_models_to_run
    [xn, fvaln, ~, nIter] = differentialevolution(DEParams, paramDefCell, DDM, [], [], [], optimInfo)
    if fvaln < fval
        fval = fvaln; x = xn;
    end
    ALL.RunFval(c) = fvaln;
    ALL.RunParams(c,:) = xn;
    disp(['Finished diff evolution fitting run ' num2str(c)])
end



%%
ALL.iterations = nIter*c;
ALL.runtime_minutes = toc/60;
ALL.BestFit = fval;
ALL.DE_Params = x;
ALL.BIC = -2*(mean(fval)) + length(VP.Params)*log(VP.nTrials/VP.Inflation);
ALL.HyperPrior = VP.HyperPriors;
ALL.ParamNames = VP.ParamNames;
ALL.Params = VP.Params;

%evaluate function at best value
for c = 1 : 10
    LL_particle(c) = DDM(x);
end
ALL.BIC_evaluated   = -2*(mean(LL_particle)) + length(x)*log(VP.nTrials/VP.Inflation);
ALL.BIC_evaluatedV2 = -2*(mean(LL_particle)*-1) + length(x)*log(VP.nTrials/VP.Inflation);

% Reference
% Stephan, K. E., Penny, W. D., Daunizeau, J., Moran, R. J., and Friston, K. J. (2009). Bayesian model
% selection for group studies. NeuroImage, 46(4):1004–1017.
ALL.BIC_evaluatedStephan = -mean(LL_particle) - ((length(x)/2)*log(VP.nTrials/VP.Inflation));

ALL.LL_evaluated = mean(LL_particle);

%keyboard

%Get model predictions
VP.Debug            = 0;
VP.ExtractSTE       = false; %extract single trial estimates?
[~, D]              = DDM(x);
ALL.Simrt           = D.rt;
ALL.Simerror        = D.error;
ALL.Congr           = VP.congr;
ALL.HumRT           = VP.RT;
ALL.HumErr          = VP.error;
ALL.D               = D;        

%Get fixed parameters values
ALL.UsedStandVal = VP.StandVal;
ALL.FixedParams = setdiff(1:9, ALL.Params);
ALL.FixedParVal = VP.StandVal(ALL.FixedParams);

rmfield(VP,{'RT' 'congr' 'flanker' 'target' 'error' 'resp' 'fitfunction'...
    'gompertz' 'weibull' 'weibull2' 'urgency' 'GenConstrains'});
ALL.VP = VP;

disp(['Saving to : ' SaveName])
save(SaveName, 'ALL' )

%DDM_plot(D)

return

