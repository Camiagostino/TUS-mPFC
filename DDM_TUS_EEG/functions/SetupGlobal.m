function [ ] = SetupGlobal(uVP, AD, QR)
%%
%This function sets the global variable VP to default fields.
global VP

%setup global VP if it does not exist in the workspace
if ~isfield(VP,'UseLog')
    VP.congr= [];
    VP.flanker= [];
    VP.target= [];
    VP.error= [];
    VP.RT= [];
    VP.resp= [];
    if ~isfield(VP,'Inflation')
        VP.Inflation=3; %inflate a VP's data n times (reduced noise due to random number generation)
    end
    if VP.Inflation > 1
        disp(['Increasing trial number to ' num2str(VP.Inflation) 'x original value.'])
    end
    VP.Debug=0; %activate to plot the decision process (drift rate and RT) per trial
    
    VP.Params= [];
    VP.ParamNames= {};
    VP.OnlyCongruent = 0; %boundary collapse (if set) confined  to congruent trials
    VP.UrgencyOn = 0;              %0 = no urgency, 1 = urgency signal (multiplicative gain on full diffusion which scales both the mean and variance of the distribution up)
                                   %note: the model assumes that increased gain makes only sense after target presentation, thus this is only added to the target diffusion process
                                   %if 1, speed parameter of collapse (lamda) becomes gain of urgency signal    
    VP.Maximize= 1;
    VP.Simulation= 0; %1 = only produce data (no fit calculated
    VP.FitMethod= 1; % 1 = chi2, 2 = loglikelihood
    VP.fitmode     = 2; %1 = fit to individual quantiles, 2 = fit to all VP equal sized bins, 3 = fit to all VP quantiles
    VP.FitGroup    = 0; %0 = fit to data from VP, 1 = fit to data from group (VP.fitmode then has to be either 2 or 3) 
    VP.CountError= 1;
    VP.nTrials= 1088;
    VP.rng= [];
    VP.UseLog= 0;
    VP.CutLog= 1e-5;   %lower limit for calculation accuracy
    VP.CutPrior = 1e-10;
    VP.FitDirection = 1; %1 = return P(M|D), 2 = return P(D|M) (Bayesian)
    VP.Disp = 1; %display stuff
    VP.DispSetup = 0; %display setup of RT
    VP.eta = 0.2; %drift noise within;
    if ~isfield(VP,'quant_space')
        VP.quant_space = 0.1; %space between quantiles
    end
    
    %define mixture model
    VP.MixturePercent = 2; %proportion of trials that can be drawn from the mixture model
    VP.MixtureBounds = [130 1000]; %range for the mixture model
    
    VP.fitfunction= @DDM_llf;
    
    VP.fit_ratio_congruent_error= 0; %set to 1 to include congruent errors in chi2 stats, 2 = incongruent errors as well
    VP.weight_congruent_errors= 1; %this multiplies SSE to increase importance for model fit
    VP.useExp= 1;
    VP.symbmath= 0;
    VP.UseGroupErrQuant = 0;
    
    VP.EqBin_Edge_Group_cor= [0 110 220 330 440 550 660 770 880 990 Inf];
    VP.EqBin_Edge_Group_err= [0 110 220 330 440 550 660 770 880 990 Inf];
    VP.Quant_Edge_Group_cor= [0 267 297 323 346 368 390 413 444 496 Inf];
    VP.Quant_Edge_Group_err= [0 197 229 245 257 269 281 296 316 351 Inf];
    VP.Quant_Edge_VP_cor= [ ];
    VP.Quant_Edge_VP_err= [ ];
    VP.n_Group_Quantile_cor= [ ];
    VP.n_Group_Quantile_err= [ ];
    VP.x2df= [ ];

    VP.n_VP_EqBin_Edge_Group_cor= [ ];
    VP.n_VP_EqBin_Edge_Group_err= [ ];
    VP.n_VP_Quant_Edge_Group_cor= [ ];
    VP.n_VP_Quant_Edge_Group_err= [ ];
    VP.n_VP_Quant_Edge_VP_cor= [ ];
    VP.n_VP_Quant_Edge_VP_err= [ ];

    VP.prior_present= [];
    VP.prior_distributions= [];
    VP.prior_functions= [];

    %define functions for boundary collapse
    VP.gompertz = @(x,b,y) (1-exp(-exp(-y*(x-b)))); %gompertz in the range of 0:1 (a is implicitly set to 1)
    VP.weibull = @(x,a,l) (a - (1-exp(-(x./l).^3))*a);
    VP.weibull2 = @(x,a,l) (a - (1-exp(-(x./l).^3))*0.5*a);%a can be fixed to 1 (scale of collapse if total collapse is assumed), l reflects the quickness of the collapse
    VP.urgency = @(t,d,sx,sy) (sy*exp(sx*(t-d))/1+exp(sx*(t-d)))+(1+(1-sy)*exp(-sx*d)/1+exp(-sx*d));
end
% keyboard
%%
%specific setup for this VP
% if isempty(AD)
%     load('data/All Behavior3.mat')
% end

NanInd = ~isnan(AD.rt(:,uVP)); %remove all nan values of RTs

VP.congr    = AD.congr(NanInd,uVP);
VP.error    = AD.error(NanInd,uVP);
VP.RT       = AD.rt(NanInd,uVP);
VP.resp     = AD.resp(NanInd,uVP);
VP.flanker  = 2*(AD.Flanker(NanInd,uVP)-1)-1; %needs to be -1 (left) and +1 (right)
VP.target   = 2*(AD.Target(NanInd,uVP)-1)-1;
if VP.Inflation>1
    VP.congr = repmat(VP.congr,VP.Inflation,1);
    VP.error = repmat(VP.error,VP.Inflation,1);
    VP.RT = repmat(VP.RT,VP.Inflation,1);
    VP.resp = repmat(VP.resp,VP.Inflation,1);
    VP.flanker = repmat(VP.flanker,VP.Inflation,1);
    VP.target = repmat(VP.target,VP.Inflation,1);
end
VP.nTrials  = length(VP.congr);

if QR==1 %quick setup: do not calculate this again, just load the subject data
    return;
end

% %Do some pruning: RT must not be negative (participant pressed button before target onset)
rt = [AD.rt]; error = [AD.error];
% rt(rt<0) = nan;
% rt(rt>1100) = nan; %we also discard extreme RT values over 1100 ms (0.3% of data)

%What is the cut off for fast guesses?? (around 210 ms)
% Quantile_Spacing = [0.0:0.02:1];
% Quantile_Number = length(Quantile_Spacing)-1;
% [Quant_Edge_Group_cor]=lequantile(nanrem(reshape(rt(error==1),[],1)),'qq', Quantile_Spacing, 'eQuant', [0 inf]); %as we do not know the highest RT in the DDM, the border of the last quantile should be inf
% [Quant_Edge_Group_err]=lequantile(nanrem(reshape(rt(error==2),[],1)),'qq', Quantile_Spacing, 'eQuant', [0 inf]);
% close all; figure
% subplot(231)
% His1=histogram(nanrem(reshape(rt(error==1),[],1)),Quantile_Number); hold %,'Normalization','Probability'
% His2=histogram(nanrem(reshape(rt(error==2),[],1)),Quantile_Number);
% keyboard
%%
VP.n_VP_EqBin_Edge_Group_cor = [];
VP.n_VP_EqBin_Edge_Group_err = [];
VP.n_VP_Quant_Edge_Group_cor = [];
VP.n_VP_Quant_Edge_Group_err = [];
VP.n_VP_Quant_Edge_VP_cor    = [];
VP.n_VP_Quant_Edge_VP_err    = [];

%split the global rt data into n quantiles (to calculate X2 statistics)
Quantile_Spacing = [0:VP.quant_space:1];
Quantile_Number = length(Quantile_Spacing)-1;

QuantileSet.space = Quantile_Spacing;
QuantileSet.extremes = [0 inf];

Group_cor_RT = nanrem(reshape(rt(error==1),[],1));
Group_err_RT = nanrem(reshape(rt(error==2),[],1));

VP_cor_RT = nanrem(reshape(VP.RT(VP.error==1),[],1));
VP_err_RT = nanrem(reshape(VP.RT(VP.error==2),[],1));

VP.NoErr = 0;
if isempty(VP_err_RT); VP.NoErr = 1; end

[VP.Quant_Edge_Group_cor,VP.n_Group_Quantile_cor] = agf_quantile(Group_cor_RT,  QuantileSet);
[VP.Quant_Edge_Group_err,VP.n_Group_Quantile_err] = agf_quantile(Group_err_RT, QuantileSet);


[~,VP.EqBin_Edge_Group_cor]     = histcounts(Group_cor_RT,Quantile_Number); VP.EqBin_Edge_Group_cor(end)    = inf;

%Individual VP quantiles and N
VP.Quant_Edge_VP_cor = agf_quantile(VP_cor_RT, QuantileSet); 

%All errror trials
if ~VP.NoErr && numel(VP_err_RT)/VP.Inflation > Quantile_Number
    VP.Quant_Edge_VP_err = agf_quantile(VP_err_RT, QuantileSet); 
    %very particular case: sometimes quantiles can not be evenly space, we have to use correct quantiles then
    if sum(hist(VP.Quant_Edge_VP_err(1:end-1),unique(VP.Quant_Edge_VP_err(1:end-1)))>1)>0 %one value is more than once a quantile edge
        VP.Quant_Edge_VP_err        = VP.Quant_Edge_Group_err; VP.UseGroupErrQuant = 1;
    end
    [~,VP.EqBin_Edge_Group_err]     = histcounts(Group_err_RT,Quantile_Number); VP.EqBin_Edge_Group_cor(end)    = inf;
else %VP has no errors, quantiles are the same as group error quantiles but count will be zero
    VP.UseGroupErrQuant = 1;
    VP.Quant_Edge_VP_err        = VP.Quant_Edge_Group_err;
    VP.EqBin_Edge_Group_err     = VP.EqBin_Edge_Group_err;
end

VP.rng                      = []; %rng
VP.UseLog                   = 0; %log transform pseudo chi2 statistic (1) or not (0)

%VP quantiles based on group edges must be calcualted
for c =  1 : Quantile_Number
    VP.n_VP_EqBin_Edge_Group_cor(c)      = numel(VP.RT(VP.error==1 & VP.RT>=VP.EqBin_Edge_Group_cor(c) & VP.RT<VP.EqBin_Edge_Group_cor(c+1)));
    VP.n_VP_EqBin_Edge_Group_err(c)      = numel(VP.RT(VP.error==2 & VP.RT>=VP.EqBin_Edge_Group_err(c) & VP.RT<VP.EqBin_Edge_Group_err(c+1)));
    VP.n_VP_Quant_Edge_Group_cor(c)      = numel(VP.RT(VP.error==1 & VP.RT>=VP.Quant_Edge_Group_cor(c) & VP.RT<VP.Quant_Edge_Group_cor(c+1)));
    VP.n_VP_Quant_Edge_Group_err(c)      = numel(VP.RT(VP.error==2 & VP.RT>=VP.Quant_Edge_Group_err(c) & VP.RT<VP.Quant_Edge_Group_err(c+1)));
    VP.n_VP_Quant_Edge_VP_cor(c)         = numel(VP.RT(VP.error==1 & VP.RT>=VP.Quant_Edge_VP_cor(c) & VP.RT<VP.Quant_Edge_VP_cor(c+1)));
    VP.n_VP_Quant_Edge_VP_err(c)         = numel(VP.RT(VP.error==2 & VP.RT>=VP.Quant_Edge_VP_err(c) & VP.RT<VP.Quant_Edge_VP_err(c+1)));
end

return;






%%
% [VP.Quant_Edge_Group_cor,VP.n_Group_Quantile_cor] = lequantile(nanrem(reshape(rt(error==1),'qq', Quantile_Spacing, 'eQuant', [0 inf]); %as we do not know the highest RT in the DDM, the border of the last quantile should be inf
% [Quant_Edge_Group_err,~,~,VP.n_Group_Quantile_err]=lequantile(nanrem(reshape(rt(error==2),[],1)),'qq', Quantile_Spacing, 'eQuant', [0 inf]);
% VP.n_Group_Quantile_cor(1)=[];VP.n_Group_Quantile_err(1)=[];
% 
% Quant_Edge_Group_cor
% VP.n_Group_Quantile_cor
% ert
%get quantiles / bins
%All correct trials









[Quant_Edge_VP_cor]         = lequantile(nanrem(reshape(VP.RT(VP.error==1),[],1)),'qq', Quantile_Spacing, 'eQuant', [0 inf]); 
[~,edges1]                  = histcounts(nanrem(reshape(rt(error==1),[],1)),Quantile_Number);
VP.EqBin_Edge_Group_cor     = edges1; VP.EqBin_Edge_Group_cor(end)   = inf;
VP.Quant_Edge_Group_cor     = Quant_Edge_Group_cor; %quantile bin borders of globa RT distribution of whole dataset

%All errror trials
ErrRTVP = nanrem(reshape(VP.RT(VP.error==2),[],1));
if ~isempty(ErrRTVP) && numel(ErrRTVP)/VP.Inflation > Quantile_Number
    [Quant_Edge_VP_err]         = lequantile(ErrRTVP,'qq', Quantile_Spacing, 'eQuant', [0 inf]); 
    [~,edges2]                  = histcounts(ErrRTVP,Quantile_Number);
else %VP has no errors, quantiles are the same as correct but count will be zero
    [Quant_Edge_VP_err]         = Quant_Edge_VP_cor;
    edges2                      = edges1;
end
VP.EqBin_Edge_Group_err     = edges2; VP.EqBin_Edge_Group_err(end)   = inf;
VP.Quant_Edge_Group_err     = Quant_Edge_Group_err;

VP.Quant_Edge_VP_cor        = Quant_Edge_VP_cor; %quantile bin borders of RT based on this VP
VP.Quant_Edge_VP_err        = Quant_Edge_VP_err;
VP.rng                      = []; %rng
VP.UseLog                   = 0; %log transform pseudo chi2 statistic (1) or not (0)

VP.n_VP_EqBin_Edge_Group_cor = [];
VP.n_VP_EqBin_Edge_Group_err = [];
VP.n_VP_Quant_Edge_Group_cor = [];
VP.n_VP_Quant_Edge_Group_err = [];
VP.n_VP_Quant_Edge_VP_cor    = [];
VP.n_VP_Quant_Edge_VP_err    = [];

for c =  1 : Quantile_Number
    VP.n_VP_EqBin_Edge_Group_cor(c)      = numel(VP.RT(VP.error==1 & VP.RT>=VP.EqBin_Edge_Group_cor(c) & VP.RT<VP.EqBin_Edge_Group_cor(c+1)));
    VP.n_VP_EqBin_Edge_Group_err(c)      = numel(VP.RT(VP.error==2 & VP.RT>=VP.EqBin_Edge_Group_err(c) & VP.RT<VP.EqBin_Edge_Group_err(c+1)));
    VP.n_VP_Quant_Edge_Group_cor(c)      = numel(VP.RT(VP.error==1 & VP.RT>=VP.Quant_Edge_Group_cor(c) & VP.RT<VP.Quant_Edge_Group_cor(c+1)));
    VP.n_VP_Quant_Edge_Group_err(c)      = numel(VP.RT(VP.error==2 & VP.RT>=VP.Quant_Edge_Group_err(c) & VP.RT<VP.Quant_Edge_Group_err(c+1)));
    VP.n_VP_Quant_Edge_VP_cor(c)         = numel(VP.RT(VP.error==1 & VP.RT>=VP.Quant_Edge_VP_cor(c) & VP.RT<VP.Quant_Edge_VP_cor(c+1)));
    VP.n_VP_Quant_Edge_VP_err(c)         = numel(VP.RT(VP.error==2 & VP.RT>=VP.Quant_Edge_VP_err(c) & VP.RT<VP.Quant_Edge_VP_err(c+1)));
end

if VP.DispSetup
    close all; figure
    subplot(231)
    His1=histogram(nanrem(reshape(rt(error==1),[],1)),Quantile_Number); hold %,'Normalization','Probability'
    His2=histogram(nanrem(reshape(rt(error==2),[],1)),Quantile_Number);
    title(['Median RT: ' num2str(mean(nanmedian(rt(error==1)))) ' ' num2str(mean(nanmedian(rt(error==2))))])

    subplot(232)
    histogram(nanrem(reshape(VP.RT(VP.error==1),[],1)),20,'Normalization','Probability'); hold
    histogram(nanrem(reshape(VP.RT(VP.error==2),[],1)),20,'Normalization','Probability');
    title(['Median RT: ' num2str(nanmedian(VP.RT(VP.error==1))) ' ' num2str(nanmedian(VP.RT(VP.error==2)))])

    subplot(234)
    bar([VP.n_VP_EqBin_Edge_Group_cor; VP.n_VP_EqBin_Edge_Group_err]);
    title('VP binned by all data bins')

    subplot(235)
    bar([VP.n_VP_Quant_Edge_Group_cor; VP.n_VP_Quant_Edge_Group_err]);
    title('VP binned by all data quantiles')

    subplot(236)
    bar([VP.n_VP_Quant_Edge_VP_cor; VP.n_VP_Quant_Edge_VP_err]);
    title('VP binned within')
end






end
