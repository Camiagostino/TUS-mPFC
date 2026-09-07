function []=ParameterRecovery(uVP, ModelSelect)
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
TrialString = {'all' 'all ce_rem' 'post cor' 'post err'};
%%%%%%%%%CHANGE ON SERVER%%%%%%%%%%%%%%%
ModelSelect = 6; 
%%%%%%%%%CHANGE ON SERVER%%%%%%%%%%%%%%%
pause(uVP/1000) %prevent identical random bumbers for iterations

n_models_to_run = 1;
if ~exist('uVP')
    uVP=1;
end
if uVP < 10
    AddS = '000';
elseif uVP < 100
    AddS = '00';
elseif uVP < 1000
    AddS = '0';
else
   AddS = ''; 
end

if ~exist('ModelSelect');ModelSelect=1;end
if ~exist('TrialSelect'); TrialSelect =1; end

SaveName = ['/home/data/scratch/adrian/roger/flanker/DDM/ParameterRecovery01/model no ' num2str(ModelSelect) '-' AddS num2str(uVP) '.mat'];
load('data/All Behavior4.mat')
% AD.rt = DDM_prune_RT(AD.rt, 0.02, [130 1000]);
posterror = AD.postE_NEW; 

global VP;
VP.Inflation=5;
VP.quant_space = 0.1;
SetupGlobal(1,AD,0);
disp(['Running Model no ' num2str(ModelSelect) ' for Iteration ' num2str(uVP)])

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
VP.GenConstrains = @(v,sv,a,sa,ter,ster,col,favl,sf) (a>sa*0.8); %& a>sa/2 we dont want more drift rate variance than the actual drift rate and prevent that boundaries are crossed before trial starts

tic   
PR.v = [0.1 8.5];       %1
PR.sv = [0 1.5];        %2
PR.a = [0.01 0.45];     %3 
PR.sa = [0.05 0.3];     %4
PR.ter = [0.1 0.4];     %5   100 to 400 ms
PR.st = [0.01 2];       %6   0 to 100 ms
PR.k = [10e-4 2.5];     %7
PR.f = [0.1 1.1];       %8
PR.sf = [0 1];          %9

VP.StandVal = [1 0.5 0.22 0.11 0.25 0.8 inf 1 0.45];
StandQua = [0 0 0 0 0.001 0.01 0 0 0]; %quantization for time speeds up convergence
VP.FixToGroupMean = 0; VP.GroupTrial = nan;



if ModelSelect == 1  %6 Parameter DDM
elseif ModelSelect == 2 %7 Parameter DDM
elseif ModelSelect == 3 %8 parameter DDM
    %Model with variable flanker value and variance
    VP.HyperPriors  = [PR.v;       PR.sv;   PR.a;     PR.sa;        PR.ter;   PR.st;       PR.f;     PR.sf]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'  'sv'    'a'       'StartVar'    'Ter'     'st'         'f'       'sf'   }; %   'st''Collapse'
    VP.Params       = [1           2        3          4            5         6            8         9];
    GroupBestFit = [1.4            0.06     0.1590     0.0572       0.2910    1             0.5920   0.45]';
    VP.StandVal =   [2.1543       0.7917    0.1590     0.0572       0.2910    1.1800  inf   0.5920   0.45]';
elseif ModelSelect == 4
elseif ModelSelect == 5
elseif ModelSelect == 6
        %Fixed Variance, no urgency
    VP.HyperPriors  = [PR.v;                PR.a;                  PR.ter;                      PR.f;    ]; %;  0 2     ;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate'         'a'                    'Ter'                         'f'      }; %   'st''Collapse'
    VP.Params       = [1                    3                        5                           8        ];
    VP.FixToGroupMean = 3; VP.GroupTrial = 1;
    if VP.FixToGroupMean
        PreMLP = load(['data/DE6 - MLE Params model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}]);
        disp(['Setting fixed parameters to mean of model no ' num2str(VP.FixToGroupMean) ' ' TrialString{VP.GroupTrial+1}])
        VP.StandVal =  mean(PreMLP.DE_Params,2)'; GroupBestFit = VP.StandVal(VP.Params)';
    else
        GroupBestFit = [5.1946                0.37415                 0.30111                       0.4637    ]';
        VP.StandVal =  [5.1946    0.99298     0.37415     0.18267     0.30111   1.1755    inf       0.6207    0.64585]';
    end
elseif ModelSelect == 7
elseif ModelSelect == 8
elseif ModelSelect == 9
end
%% Simulate model data
VP.Simulation = 0;
TrialString = {'all' 'all ce_rem' 'post cor' 'post err'};
clear fDDM VP.fitfunction
PreMLP = load(['data/DE6 - MLE Params model no ' num2str(3) ' ' TrialString{2}]);
disp(['Setting fixed parameters to mean of model no ' num2str(3) ' ' TrialString{2}])
VP.StandVal =  mean(PreMLP.DE_Params,2)'; GroupBestFit = VP.StandVal(VP.Params)';

StdDevParams = std(PreMLP.DE_Params,0,2);
StdDevParams = StdDevParams(VP.Params);
MeanParams = VP.StandVal(VP.Params);
rng('shuffle');

for c = 1 : length(VP.Params)
   pd(c) = makedist('Normal','mu',MeanParams(c),'sigma',StdDevParams(c));
end

FieldStr = {'v' 'sv' 'a' 'sa' 'ter' 'st' 'k' 'f' 'sf'};
FieldStr = FieldStr(VP.Params);

ResultPlausible = false; n_draw_count=1;
while ~ResultPlausible
    %initialize
    v = -inf;sv = -inf;a = -inf;sa = -inf;ter = -inf;st = -inf;k = -inf;f = -inf;sf = -inf;
    ap = [v sv a sa ter st k f sf];
    ap = ap(VP.Params);
    
    for c = 1 : length(VP.Params)
    %draw samples
        while ~all(ap(c) >= PR.(FieldStr{c})(1) & ap(c) <= PR.(FieldStr{c})(2))
            ap(c) = random(pd(c));
        end
    end
    
    disp('All parameters fullfil range restrictions...')
    
    [LL, D] = DDM(ap');
    
    if LL~=inf && sum(D.error==2)>30
        ResultPlausible=1; 
        disp('Simulation without failure and enough errors...')
    else
        disp('Simulation produced error of error number not sufficient, repeating...')
    end
end

VP.error    = D.error;
VP.RT       = D.rt;
VP.resp     = D.resp;
VP.MixtureBounds = [min(VP.RT) max(VP.RT)]; %range for the mixture model

ALL.GT_Params = ap;

DDM_quantilizer(D.rt,D.error);
% clf
% plot(VP.Quant_Edge_VP_cor);hold on;
% plot(VP.Quant_Edge_VP_err)

GroupBestFit = MeanParams';

disp(['GROUND TRUTZ PARAMETERS: ' num2str(ap)])

LL = DDM(ap');
ALL.LL_GT = LL;

disp(['GROUND TRUTZ LL: ' num2str(LL)])


%%
clear fDDM VP.fitfunction
tic;
% define parameter names, ranges, quantizations and initial values :
paramDefCell = {'', VP.HyperPriors, StandQua(VP.Params)', GroupBestFit};
objFctSettings={};
% get default DE parameters
DEParams = getdefaultparams;
DEParams.NP = length(VP.Params)*10;

% use a subfunction to check parameter vectors for validity
DEParams.validChkHandle = [];%@demo3_constraint;
DEParams.saveHistory = 0;

% set times
DEParams.maxiter  = 4000;
DEParams.maxtime  = 3600;  % in seconds
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
ALL.DE_Params = x';
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
ALL.FixedParams = setdiff(1:8, ALL.Params);
ALL.FixedParVal = VP.StandVal(ALL.FixedParams);

rmfield(VP,{'RT' 'congr' 'flanker' 'target' 'error' 'resp' 'fitfunction'...
    'gompertz' 'weibull' 'weibull2' 'urgency' 'GenConstrains'});
% ALL.VP = VP;

disp(['Saving to : ' SaveName])
save(SaveName, 'ALL' )

return





ResultPlausible = false; n_draw_count=1;
while ~ResultPlausible
    %initialize
    DriftRate = -inf;
    DRvariance = -inf;
    Boundary = -inf;
    StartVar = -inf;
    NondecisionTime = -inf;
    st = -inf;
    
    %We draw our test values randomly, yet constrained to plausible values
    pd1 = makedist('Normal','mu',0.8,'sigma',0.6);
    pd2= makedist('HalfNormal','mu',0,'sigma',1);
    pd3= makedist('Normal','mu',0.15,'sigma',0.2);
    pd4= makedist('Normal','mu',0.2,'sigma',0.15);
    pd5= makedist('Normal','mu',0.2,'sigma',0.3);
    pd6= makedist('Normal','mu',0.3,'sigma',0.5);
    
    %draw samples
    while DriftRate<0
        DriftRate = random(pd1);
    end
    DRvariance = random(pd2);
    while Boundary<0
        Boundary = random(pd3);
    end
    while StartVar<0
        StartVar = random(pd4);
    end
    while NondecisionTime<0
        NondecisionTime = random(pd5);
    end
    while st<0
        st = random(pd6);
    end
    

%     DriftRate = 0.5;            %1 (v)    - speed of evidence accumulation for the indicated response (moves the mean of the deviation to derive noise from)
%     DRvariance = 0.4;           %2 (sv)   - variance of the drift rate distribution from trial to trial
%     Boundary = 0.2;            %3 (±a)   - in DDM terminology
%     StartVar = 0.4;             %4 (sz)   - starting point variability, reflects the range around start value as borders of a uniform distribution. 0.4 would mean that start values can be between ±0.2.
%     NondecisionTime = 0.12;      %5 (Ter)  - time in ms that visual information need to arrive at our decision maker
%     st = 0.1;                    %6 (st)   - variance of non-decision time (0 = no variance here)
%     CollapseDegree = inf;%2     %7 Boun

    VP.HyperPriors  = [0 5;         0 1;           0.01 1;     0 1;      0.05 0.4;  0 2     ]; %;    0 5 fixed borders for all parameters (these are evaluated before the model call, which speeds up MCMC considerably)   
    VP.ParamNames   = {'DriftRate' 'DRvariance'   'Boundary'  'StartVar' 'Ter'      'st' }; %'Collapse'
    VP.Params       = [1           2               3           4           5        6      ]; %this maps input to model order of factors 
    VP.OnlyCongruent = 0;
    VP.Maximize = 1;
    VP.CountError = 1;
    VP.nTrials = nTrials;
    VP.rng              = []; %rng
    VP.UseLog           = 0; %log transform pseudo chi2 statistic (1) or not (0)
    VP.fitfunction=@no_action;
    VP.CutLog=1e-5;
    VP.FitMethod = 2; 

    PriorSting = '';
    for c = 1 : length(VP.Params) %formulate hard priors (hyperpriors)
        PriorSting = [PriorSting '(DDM(' num2str(c) ')>' num2str(VP.HyperPriors(c,1)) ')&&(DDM(' num2str(c) ')<' num2str(VP.HyperPriors(c,2)) ') && '];
    end
    PriorSting(end-3:end)='';
    DDMprior =@(DDM) eval(PriorSting);

    %derive RT distribution from these parameters
    VP.Simulation=1;
    ParamVector = [DriftRate DRvariance Boundary StartVar NondecisionTime st]; 
    [~, D] = DDM(ParamVector);
    if ~isempty(D) %model did not aboart (i.e., due to negative rt or something else)
        VP.Simulation=0;
        VP.error    = D.error;
        VP.RT       = D.rt;
        VP.resp     = D.resp;
        nanmean(D.rt)
        nanmean(VP.error-1)
        ParamVector
        if nanmean(D.rt)>lower_rt && nanmean(D.rt)<higher_rt && nanmean(VP.error-1)<lower_er && nanmean(VP.error-1)>higher_er
            ResultPlausible = true;
        else
           disp('Results violate participant restrictions, repeating the draw of model parameters...') 
           n_draw_count=n_draw_count+1
        end
    end
end