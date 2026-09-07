function [varargout]=fDDM(TheInput, varargin)
%%
global VP
%%
%free
%always fixed:
if isempty(TheInput) 
    TheInput=varargin{:};
end

persistent DriftRate DriftNoiseWithin DRvariance NondecisionTime StartVar Boundary  D Prior Priors st StartVal SOA k StepNumber FlankerValue RTcutoff LTargetD logl InputSort TerD MaxGain sf
if isempty(DriftRate) %set standard values for fittable parameters
    DriftRate       = VP.StandVal(1);           %1 (v)    - speed of evidence accumulation for the indicated response (moves the mean of the deviation to derive noise from)
    DRvariance      = VP.StandVal(2)^2;         %2 (sv)   - variance of the drift rate distribution from trial to trial
    Boundary        = VP.StandVal(3);           %3 (±a)   - in DDM terminology
    StartVar        = VP.StandVal(4);           %4 (sz)   - starting point variability, reflects the range around start value as borders of a uniform distribution. 0.4 would mean that start values can be between ±0.2.
    NondecisionTime = VP.StandVal(5)*1000;      %5 (Ter)  - time in ms that visual information need to arrive at our decision maker
    st              = VP.StandVal(6)*10;        %6 (st)   - variance of non-decision time (0 = no variance here)
    k               = VP.StandVal(7);           %7 Bound collapse = lambda of weibull distribution, inf = no collapse
    FlankerValue    = VP.StandVal(8);           % degree to which flanker relative to target influence accumulation (1 = equal influence, most simple assumption)
    sf              = VP.StandVal(9)^2;         % variance in flanker drift rate
    
    %non-free settings
    StartVal = 0;               %(z)    - (or bias) Usually starts each trial at 0 + noise (can reflect a bias of a trial, e.g., beta rebound)
    DriftNoiseWithin = VP.eta;  %(s/eta)    - amount of noise per computation step of the diffusion (usually fixed, often to 0.01)
    SOA = 80;                   %Offset between flanker and target
    StepNumber = 0.001;         %reflects delta parameter in normal DDM nomenclature, this simply scales parameters to keep them in a nicer range
    MaxGain = 5;                %this is the maximal gain added to the signal in the urgency on condition
    
    %how long does the target drift maximally have to be? 
    RTcutoff = 1100; %everything longer than this will be quantified in the last bin and is no longer explicitly simulated (sign of last datapoint determines response)
    LTargetD = RTcutoff -SOA-1;
    logl = nan;
    Priors = nan(length(VP.prior_present),1);
    Prior=0;

    %retrieve parameters
    InputSort = eval(get_sorthandle);
    [DriftRate DRvariance  Boundary  StartVar NondecisionTime st k  FlankerValue sf] = InputSort(TheInput);
    if VP.Disp
        %%%%%%print a short summary of model settings for easier control and bug fixing
        FixVsFree = ismember(1:9, VP.Params);
        p = [[DriftRate DRvariance  Boundary  StartVar NondecisionTime st k  FlankerValue sf]; FixVsFree];
        printtext=sprintf('v:\t%.2f\t\tfree:\t%d\n sv:\t%.2f\t\tfree:\t%d\n a:\t%.2f\t\tfree:\t%d\n sa:\t%.2f\t\tfree:\t%d\n ter:\t%.2f\t\tfree:\t%d\n st:\t%.2f\t\tfree:\t%d\n coll:\t%.2f\t\tfree:\t%d\n fval:\t%.2f\t\tfree:\t%d\n  sf:\t%.2f\t\tfree:\t%d\n',p(:));
        fprintf('\n\n\nStarting fDDM with the following settings:\n\n %s\n\nWithin trial noise is fixed to: %.2f\n',printtext, DriftNoiseWithin);
        if VP.UrgencyOn
            disp('k is used as a multiplicative gain (urgency signal) instead of boundary collapse.')
        end
        drawnow;
    end
end
%%
%reset values on novel iteration
[DriftRate DRvariance  Boundary  StartVar NondecisionTime st k  FlankerValue sf] = InputSort(TheInput);
% keyboard
%%
% %calculate prior probabilities
if ~isempty(VP.prior_present)
    for c = 1 : length(VP.prior_present)
%         Prior = VP.prior_functions{c}(TheInput(VP.prior_present(c))) * Prior;
        Priors(c,1) = VP.prior_functions{c}(TheInput(VP.prior_present(c)));
    end
    Priors(Priors<VP.CutPrior) = VP.CutPrior;
    Prior = sum(log(Priors))*VP.Inflation;
end
% if Prior==0 %stop here if prior is 0
%     varargout{1}=-inf;varargout{2}=[];return;
% end
%%
%test general constraints (if these are violated, script terminates to save time)
if ~VP.GenConstrains(DriftRate, DRvariance, Boundary, StartVar, NondecisionTime, st, k, FlankerValue, sf)
    varargout{1} = -inf * VP.Maximize; varargout{2} = []; return;
end
%%
%save time for variables that need to be calculated and do this only on the first function call
persistent CollapseTime Bound uBound TrialNoise TrialMu rt resp  RandomDriftRate RandomNonDecisionTime RandomStartPoint Error Temp FlankerDrift nTrials Gain RandomFlanker
if isempty(rt)
    if ~isfield(VP,'nTrials')
        nTrials = length(VP.RT);    %number of trials to simulate
    else
        nTrials = VP.nTrials;
    end
    
    %Simulate random trial sequences in a persistent variable
    Error = ones(1,nTrials);
%     RandomFlanker = ones(1,nTrials);
    rt = nan(1,nTrials);
    resp=rt;
    FlankerDrift = nan(nTrials, RTcutoff+1);
    TrialNoise = nan(nTrials, RTcutoff-1);
    TrialMu = nan(nTrials, RTcutoff-1);
    
    Temp = FlankerDrift';
    CollapseTime = linspace(0,2,RTcutoff+1); %bit faster to use additional var
    %calculate the decision boundaries that can collapse exponentially
%     Bound = [repmat(TreshA, 1, NonCollapse) TreshA.*exp(-[1 : RTcutoff+1 - NonCollapse].*StepNumber.*k)];
    Bound = VP.weibull(CollapseTime,1,inf).*Boundary;
    uBound = repmat(Bound,[nTrials 1])';
    Gain = ones(1,length(CollapseTime)); %Gain of urgency signal (1 means no change)
end
%%
%these are fittable random sequences that need to change upon model iteration
if ~isempty(VP.rng); rng(VP.rng); end
RandomDriftRate  = normrnd(DriftRate,DRvariance,nTrials,1);if ~isempty(VP.rng); rng(VP.rng); end
RandomNonDecisionTime = round(unifrnd(NondecisionTime-st/2,NondecisionTime+st/2,nTrials,1));if ~isempty(VP.rng); rng(VP.rng); end %needs to be integer
RandomStartPoint=unifrnd(StartVal-StartVar/2,StartVal+StartVar/2,nTrials,1); % starting point for all trials
%draw random flanker modulation
RandomFlanker  = normrnd(FlankerValue,sf,nTrials,1); 

if st ~= 0
    TerD    = zeros(1, max(RandomNonDecisionTime));
    for c = 1 : nTrials %if Ter is variable, could not find fast vectorized solution that beats simple loop
        TrialMu(c,:) = [TerD(1:RandomNonDecisionTime(c)) repmat(RandomDriftRate(c).*VP.flanker(c).*RandomFlanker(c).*StepNumber,1,SOA)  repmat(RandomDriftRate(c).*VP.target(c).*StepNumber,1,LTargetD-RandomNonDecisionTime(c))];
    end
else %~3x faster if Ter is always the same    
    TrialMu = [zeros(nTrials,RandomNonDecisionTime(1)) repmat(RandomDriftRate.*VP.flanker.*RandomFlanker.*StepNumber,1,SOA) repmat(RandomDriftRate.*VP.target.*StepNumber,1,LTargetD-RandomNonDecisionTime(1))];
end

%calculate boundary collapse
if any(VP.Params==7) %boundary collappse as a free parameter requires that boundaries are newly calculated on every iteration
    if k>VP.HyperPriors(VP.Params==7,2)-1e-5 %close to upper boundary, do not use collapse at all
        k = inf;
    end
    if ~VP.UrgencyOn %boundary collapse
        Bound = VP.weibull(CollapseTime,1,k).*Boundary;
        uBound = repmat(Bound,[nTrials 1])';
        if VP.OnlyCongruent %do not collapse boundaries on incongruent trials
            uBound(:,VP.congr==2)=Boundary; %congr ==2 = incongruent trials, replace with normal threshold
        end
        TrialNoise = repmat(DriftNoiseWithin*sqrt(StepNumber),nTrials,RTcutoff-1);
    else %urgency signal
        uBound(:) = Boundary;
        Gain = 1+(1./VP.weibull2(CollapseTime,1,k))*MaxGain-MaxGain;
        TrialNoise = repmat(repmat(DriftNoiseWithin*sqrt(StepNumber),1,RTcutoff-1).*Gain(1:RTcutoff-1),[nTrials,1]);%scale noise with gain
        TrialMu = TrialMu.*repmat(Gain(1:RTcutoff-1),[nTrials,1]);%scale diffusion with gain
    end
else %just use the regular boundaries
    uBound(:) = Boundary;
    TrialNoise = repmat(DriftNoiseWithin*sqrt(StepNumber),nTrials,RTcutoff-1);
end
% RandomDriftRate(RandomDriftRate<0)=0; % This can happen if drift rate variance is extremely high, and drift rate low. This leads to negativ drift --> somewhat unintuitive, but should not lead to problems
% keyboard
%%
% try
    FlankerDrift(:,1:RTcutoff) = cumsum([RandomStartPoint normrnd(TrialMu,TrialNoise)],2);
% catch ME
%     [DriftRate DRvariance  Boundary  StartVar NondecisionTime st k  FlankerValue sf]
%     disp('Error in size of TrialMu and TrialNoise'); varargout{1} = -inf * VP.Maximize; varargout{2} = []; ert; return;
% end

%a problem for a vectorized solution is if the diffusion did not terminate
%we add a last row that will always cross any boundary (representing an instant collapse, which punishes the model fit for too slow drifts which should be very rare in the Flanker case anyways)
FlankerDrift(:,end) = FlankerDrift(:,end-1).*inf;

%find when borders are crossed
Temp = abs(FlankerDrift')>=uBound; %not sure if there is a way around a temp variable...
[~,rt] = max(Temp(:,any(Temp)));

resp = sign(FlankerDrift(sub2ind(size(FlankerDrift),1:nTrials, rt)))';
if any(rt<=0) %this can happen when Ter is very low and variance extremely large
    varargout{1} = -inf * VP.Maximize; varargout{2} = []; return;
end
Error(:) = 1;
Error(resp ~= VP.target) = 2;

%%
%keyboard

if VP.Debug
    h=figure;
    set(h, 'Position', [99 265 1231 677])%
    set(h,'PaperPositionMode','Auto')
    %%
    for c = 1 : nTrials
        clf
        col = 'g';
        if Error(c) == 2
            col = 'r';
        end
        subplot(211)
        yyaxis left
        plot(FlankerDrift(c,:),col); hold on;
        

        plot(uBound(:,1)','b')
        plot(-uBound(:,1)','m')
%         gridxy([], [-Boundary ])
%         gridxy([], Boundary)
        title(['RT: ' num2str(rt(c)) ' DR: ' num2str(RandomDriftRate(c))])
        set(gca, 'YLim', [-uBound(1,1) uBound(1,1)].*1.5) 
        gridxy(rt(c)); 
        gridxy(RandomNonDecisionTime(c));
        if VP.UrgencyOn
             yyaxis right
           plot((Gain),'r')
        end
        
        subplot(212)
        plot(Temp(:,c))
        pause
    end
end

%%
if ~VP.Simulation %simulation just generates data, no fitting
    %determine the likelihood based on this distribution
    [logl] = VP.fitfunction(rt', Error');
    D.prior = Prior;
%     D.observation = logl.chi_2;
%     D.posterior = (logl.chi_2 + Prior*VP.Maximize);
    D.posterior = VP.LL_tot + Prior;
    varargout{1} = D.posterior;
    if isnan(D.posterior)
        keyboard
    end
end
% D.posterior
D.params.DriftRate = DriftRate;
D.params.DRvariance = sqrt(DRvariance);
D.params.Boundary = Boundary;
D.params.StartVar = StartVar;
D.params.Ter = NondecisionTime/1000;
D.params.st = st/100;
D.params.k = k;
D.params.DriftNoiseWithin = DriftNoiseWithin;
D.params.StartPointBias = StartVal;
D.params.FlankerValue = FlankerValue;
D.params.sf = sqrt(sf);

D.rt = rt';
D.error = Error';
D.resp = resp;

if VP.ExtractSTE
    %keyboard
    
    % For display purposes add 700ms baseline and model a consecutive return of the decision variable to baseline similar to an Ornstein-
    % Uhlenbeck process to facilitate comparisons between model and BPL
    % fixed params for OU
    
    Baseline = nan(VP.nTrials,700);

    %Number of simulations
    nSims   = 1;
    
    mu      = 0;% Reverting level
    Vol     = 0.001;% volatility
    theta   = 0.04;% Speed of mean reversion
    
    for i = 1:VP.nTrials
        
        if resp(i) == 1 %flip, so that the later chosen response is always plotted downwards
            FlankerDrift(i,:)= FlankerDrift(i,:)*-1;
        end
   
        if resp(i) == 1
            Baseline(i,:) = normrnd(RandomStartPoint(i)*-1,.001,1,700);
        else
            Baseline(i,:) = normrnd(RandomStartPoint(i),.005,1,700);
        end
        % timepoint to simulate forward
        nTP = length([rt(i)+1:RTcutoff]);
        s0    = FlankerDrift(i,rt(i));% Initial value

        FlankerDrift(i,rt(i)+1:RTcutoff) = OrnsteinUhlenbeck(nTP, nSims, mu, s0, Vol, theta);

    end
    
    D.TimeInfo = -700:1100;

    %extract drift trajectory   
    D.FlankerDrift                  = FlankerDrift;
    D.STE_Bound                     = uBound;
    D.STE_RandomDrift               = RandomDriftRate;
    D.STE_RandomStart               = RandomStartPoint;
    D.STE_RandomNonDecisionTime     = RandomNonDecisionTime;

    D.FullMultistageDrift           = [Baseline FlankerDrift];

    %trial info
    D.congr = VP.congr;
    
end


varargout{2} = D;
return;

D.prior
D.posteriorX2
D.params
varargout{1}
%%
if Prior%if there is a prior, combine model evidence and prior probability
    D.ModelP = 1-chi2cdf(logl.chi_2,VP.x2df);
    D.JointP = D.ModelP*Prior;
    D.posteriorX2 = chi2inv(1-D.JointP,VP.x2df);
else
    D.posteriorX2 = logl.chi_2;
end
if D.posteriorX2 == Inf

end













