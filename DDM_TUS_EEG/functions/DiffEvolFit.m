function []=DiffEvolFit(uVP, ModelSelect,TrialSelect)
%Fit the DDM for the Modul Flanker Task
%%
% LicPresent=0;
% while ~LicPresent
%     LP = license('checkout','GADS_Toolbox');
%     if LP 
%         LicPresent = 1;
%         disp('License has been checked out...')
%     else
%         disp('Waiting for license...')
%         pause(7)
%     end;
% end;

path('/home/afischer/Documents/MATLAB/myfunctions/',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DDM/',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DDM/gwmcmc',path);
path('/home/afischer/Documents/MATLAB/myfunctions/DiffEvol',path);
%%
%%%%%%%%%%%%CHANGE ON SERVER%%%%%%%%%%
% ModelSelect=8
% TrialSelect=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

TrialString = {'all' 'all ce_rem' 'post cor' 'post err'};

n_models_to_run = 1;
if ~exist('uVP')
    uVP=1;
end
if uVP < 10
    AddS = '00';
elseif uVP < 100
    AddS = '0';
else
   AddS = ''; 
end

if ~exist('ModelSelect');ModelSelect=1;end
if ~exist('TrialSelect'); TrialSelect =1; end

SaveName = ['/home/data/scratch/adrian/roger/flanker/DDM/DiffEvolInidividual7/model no ' num2str(ModelSelect) ' ' TrialString{TrialSelect+1} '-' AddS num2str(uVP) '.mat'];
if exist(SaveName)
    disp('Model exists...')
    return;
end
load('data/All Behavior4.mat')
% AD.rt = DDM_prune_RT(AD.rt, 0.02, [130 1000]);
posterror = AD.postE_NEW; 

global VP;

VP.Inflation=5;
VP.quant_space = 0.1;

n_PostErrorVP = sum(posterror(:,uVP)==1 & ~(AD.congr(:,uVP)==1 & AD.error(:,uVP)==2)); 
n_PostCorrVP = sum(posterror(:,uVP)==0 & ~(AD.congr(:,uVP)==1 & AD.error(:,uVP)==2)); 

MinErrorTrial = 30;
% TrialSelect = 3;
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
else
    AD.rt(isnan(posterror)) = nan;
end
disp(['Selecting trials ' TrialString{TrialSelect+1}])

SetupGlobal(uVP,AD,0);

% clear TestRT TestACC 
% for c = 1 : 863
%     uVP = c;
%     SetupGlobal(uVP,AD,0);
%     TestRT(c) = mean(VP.RT);
%     TestACC(c) = 1-(mean(VP.error-1));
%     NTrialsVP(c) = VP.nTrials;
%     disp([c mean(TestRT) mean(TestACC) NTrialsVP(c) n_post_err_trials(c)])
% end
% 
% [n_post_err_trials' NTrialsVP']
% CORRECT = [863.0000  443.3341    0.8745]
% ERROR   = [863.0000  482.9861    0.9031]
% DIFF = ERROR-CORRECT
% %True PIA over all trials when congruent errors are removed = 2.68% (mean ACC post cor = 87.41%, post error 90.09%)
% %True PES over all trials when congruent errors are removed = 32.67 ms (mean RT post cor = 442.69, post error 480.2)


disp(['Running Model no ' num2str(ModelSelect) ' for VP ' num2str(uVP) ' with ' num2str(VP.nTrials) ' simulations.'])

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
VP.UrgencyOn = 0;

%get function handle
DDM = @fDDM;
sRT = sort(VP.RT);

VP.OnlyCongruent = 0;
VP.GenConstrains = @(v,sv,a,sa,ter,ster,col,favl,sf) (a>sa*0.55); %& a>sa/2 we dont want more drift rate variance than the actual drift rate and prevent that boundaries are crossed before trial starts

tic   
PR.v = [0.1 8.5];       %1
PR.sv = [0 1.5];        %2
PR.a = [0.1 0.45];       %3 
PR.sa = [0.05 0.3];     %4
PR.ter = [0.1 0.4];     %5   100 to 400 ms
PR.st = [0.01 2];       %6   0 to 100 ms
PR.k = [10e-4 2.5];     %7
PR.f = [0.1 1.1];       %8
PR.sf = [0 1];          %9

VP.StandVal = [1 0.5 0.22 0.11 0.25 0.8 inf 1 0.45];
VP.FixToGroupMean = 0; VP.GroupTrial = nan;
% StandQua = [0.01 0.001 0.001 0.001 0.001 0.01 0.01 0.001 0.001]; %
StandQua = [0 0 0 0 0.001 0.01 0 0 0]; %quantization for time speeds up convergence
if ModelSelect == 1  %6 Parameter DDM
    VP.StandVal(8) = 0; 
    VP.HyperPriors  = [PR.v;        PR.sv;   PR.a;    PR.sa;    PR.ter;   PR.st                   ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'  'sv'     'a'     'StartVar' 'Ter'     'st'                    }; %   'st''Collapse'
    VP.Params       = [1             2        3       4           5        6                      ]; % 6  this maps input to model order of factors 
    GroupBestFit = [2.3032         0.06      0.2197   0.2139     0.3040    1.5900                 ]';
    VP.StandVal = [ 1               0.5      0.22     0.11       0.25      0.8      inf  1       0];
elseif ModelSelect == 2 %7 Parameter DDM
    %Model with variable flanker value, no urgency or collapse
    VP.HyperPriors  = [PR.v;        PR.sv;   PR.a;    PR.sa;      PR.ter;  PR.st;       PR.f      ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'  'sv'    'a'      'StartVar' 'Ter'      'st'         'f'        }; %   'st''Collapse'
    VP.Params       = [1             2       3        4           5        6             8        ];
    GroupBestFit = [2.1543         1.0000    0.1590   0.0572      0.2910   1.1800        0.5920   ]';
    VP.StandVal = [ 1               0.5      0.22     0.11        0.25     0.8      inf  0.5920  0];
elseif ModelSelect == 3 %8 parameter DDM
    %Model with variable flanker value and variance
    VP.HyperPriors  = [PR.v;       PR.sv;   PR.a;     PR.sa;        PR.ter;   PR.st;       PR.f;     PR.sf]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'  'sv'    'a'       'StartVar'    'Ter'     'st'         'f'       'sf'   }; %   'st''Collapse'
    VP.Params       = [1           2        3          4            5         6            8         9];
    GroupBestFit = [1.4            0.06     0.1590     0.2139       0.2910    1             0.5920   0.45]';
    VP.StandVal =   [2.1543       0.7917    0.1590     0.0572       0.2910    1.1800  inf   0.5920   0.45]';
elseif ModelSelect == 4
    %Everything including urgency
    VP.UrgencyOn = 1;
    PR.k = [10e-4 3];
    VP.HyperPriors  = [PR.v;       PR.sv;   PR.a;     PR.sa;        PR.ter;   PR.st;  PR.k;  PR.f;     PR.sf]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'  'sv'    'a'       'StartVar'    'Ter'     'st'    'k'    'f'       'sf'   }; %   'st''Collapse'
    VP.Params       = [1           2        3          4            5         6        7      8         9];
    GroupBestFit = [4.5542     0.94905      0.36815    0.18205      0.30014  1.1825  2.2922  0.4637     0.64585]';
elseif  ModelSelect == 5
    %Fixed Variance, free urgency
    VP.UrgencyOn = 1;
    PR.k = [10e-4 3];
    VP.HyperPriors  = [PR.v;                PR.a;                  PR.ter;            PR.k;      PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'         'a'                    'Ter'                'k'       'f'      }; %   'st''Collapse'
    VP.Params       = [1                    3                        5                 7         8        ];
    GroupBestFit = [4.5542                0.36815                   0.30014           2.2922    0.4637    ]';
    VP.StandVal =  [4.5542     0.94905     0.36815     0.18205     0.30014   1.1825   2.2922    0.4637     0.64585]';
elseif ModelSelect == 6
    %Fixed Variance, no urgency
    VP.HyperPriors  = [PR.v;                PR.a;                  PR.ter;                      PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'         'a'                    'Ter'                         'f'      }; %   'st''Collapse'
    VP.Params       = [1                    3                        5                           8        ];
    VP.FixToGroupMean = 3; VP.GroupTrial = 1;
    if VP.FixToGroupMean
        PreMLP = load(['data/DE7 - MLE Params model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}]);
        disp(['Setting fixed parameters to mean of model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}])
        VP.StandVal =  mean(PreMLP.DE_Params,2)'; GroupBestFit = VP.StandVal(VP.Params)';
    else
        GroupBestFit = [5.1946                0.37415                 0.30111                       0.4637    ]';
        VP.StandVal =  [5.1946    0.99298     0.37415     0.18267     0.30111   1.1755    inf       0.6207    0.64585]';
    end
    
elseif ModelSelect == 7
    %Model with variable flanker value and variance --> larger space
    PR.sv = [0 2];  PR.a = [0.01 1]; PR.sf = [0 2]; PR.k = [10e-4 5];  
    %try quantiles
    StandQua = [0.01 0.001 0.001 0.001 0.001 0.01 0.01 0.001 0.001]; %quantization for time speeds up convergence
    VP.HyperPriors  = [PR.v;       PR.sv;   PR.a;     PR.sa;        PR.ter;   PR.st;       PR.f;     PR.sf]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'  'sv'    'a'       'StartVar'    'Ter'     'st'         'f'       'sf'   }; %   'st''Collapse'
    VP.Params       = [1           2        3          4            5         6            8         9];
    GroupBestFit = [1.4            0.06     0.1590     0.2139       0.2910    1             0.5920   0.45]';
    VP.StandVal =   [2.1543       0.7917    0.1590     0.0572       0.2910    1.1800  inf   0.5920   0.45]';
elseif ModelSelect == 8
    %Fixed Variance, no urgency
    VP.HyperPriors  = [PR.a;                  PR.ter;                      PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'a'                    'Ter'                         'f'      }; %   'st''Collapse'
    VP.Params       = [3                        5                           8        ];
    
    VP.FixToGroupMean = 6; VP.GroupTrial = 2; VP.PPF = 'DE6 - ';
    if VP.FixToGroupMean
        PreMLP = load(['data/' VP.PPF 'MLE Params model no ' num2str(VP.FixToGroupMean) ' '  TrialString{VP.GroupTrial+1} '.mat']);
        disp(['Setting fixed parameters to mean of model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}])
        VP.StandVal =  [mean(PreMLP.DE_Params,2)' 0.63]; GroupBestFit = VP.StandVal(VP.Params)';
    else
        GroupBestFit = [5.1946                0.37415                 0.30111                       0.4637    ]';
        VP.StandVal =  [5.1946    0.99298     0.37415     0.18267     0.30111   1.1755    inf       0.6207    0.64585]';
    end
elseif ModelSelect == 9
    %Fixed Variance, no urgency
    VP.HyperPriors  = [PR.v;                  PR.ter;                      PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'            'Ter'                        'f'      }; %   'st''Collapse'
    VP.Params       = [1                        5                           8        ];
    
    VP.FixToGroupMean = 6; VP.GroupTrial = 2; VP.PPF = 'DE6 - ';
    if VP.FixToGroupMean
        PreMLP = load(['data/' VP.PPF 'MLE Params model no ' num2str(VP.FixToGroupMean) ' '  TrialString{VP.GroupTrial+1} '.mat']);
        disp(['Setting fixed parameters to mean of model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}])
        VP.StandVal =  [mean(PreMLP.DE_Params,2)' 0.63]; GroupBestFit = VP.StandVal(VP.Params)';
        VP.StandVal(3) =  0.3;
    else
        GroupBestFit = [2.7                   0.37415                 0.30111                       0.4637    ]';
        VP.StandVal =  [5.1946    0.99298     0.37415     0.18267     0.30111   1.1755    inf       0.6207    0.64585]';
    end
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
DEParams.maxiter  = 4000;
DEParams.maxtime  = 900;  % in seconds
DEParams.maxclock = [];

% set display options
DEParams.infoIterations = 1;
DEParams.infoPeriod     = 60;  % in seconds

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
ALL.BIC_evaluated = -2*(mean(LL_particle)) + length(x)*log(VP.nTrials/VP.Inflation);
ALL.LL_evaluated = mean(LL_particle);

%Get model predictions
VP.Debug=0;
[~, D] = DDM(x);
ALL.Simrt = D.rt;
ALL.Simerror = D.error;
ALL.Congr = VP.congr;
ALL.HumRT = VP.RT;
ALL.HumErr = VP.error;

%Get fixed parameters values
ALL.UsedStandVal = VP.StandVal;
ALL.FixedParams = setdiff(1:9, ALL.Params);
ALL.FixedParVal = VP.StandVal(ALL.FixedParams);

rmfield(VP,{'RT' 'congr' 'flanker' 'target' 'error' 'resp' 'fitfunction'...
    'gompertz' 'weibull' 'weibull2' 'urgency' 'GenConstrains'});
ALL.VP = VP;

disp(['Saving to : ' SaveName])
save(SaveName, 'ALL' )

return



%     %particle swarm method
%     rng('shuffle')
%     options = optimoptions('particleswarm', 'Display', 'iter','SwarmSize',100,'MaxStallIterations', 20);
%     [x,fval] = particleswarm(DDM,length(VP.Params),VP.HyperPriors(:,1),VP.HyperPriors(:,2),options)
