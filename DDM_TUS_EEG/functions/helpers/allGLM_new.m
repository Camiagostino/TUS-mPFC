function [r, p, R2] = allGLM_new( data, regressoren, zeitfenster, varargin )
%%
%Function that performs singlepointcorrelations for singletrial EEG data on
%all channels. Returns regression weights or correlation coeffictients, p
%values for the model and an R2 or R2 erstimate for the whole model per
%datapoint. 
%
%INPUT:
%'indices'      Equals an array of selected epochs. Only these selected
%               epochs will be used for calculation.
%'mode'         Can be set to 'GLM', 'regress', 'robust1/2' or 'COR', will use either glmfit, regress, robustfit or
%               corr function of matlab (default is 'GLM'). Note: Function
%               returns beta-values when set to GLM, robust and Regress and
%               correlation-coefficients when set to COR. 
%               Note: Regress only sopports return of p statistics for the full model.
%               If robust1 is specified, error variance is returned instead of p value for later normalization.  
%               If robust2 is specified, error variance is automatically
%               used to normalize b values and normalized values are returned.
%'Disp'         Toggle display of percentage [0/1].


%Set defaults
Fmode       = 'glm';    % Default will use glmfit.
Fdisp       = 0;        % Default: Display is deactivated.

%Get Variables from varargin array
nargs = nargin-3;
if nargs > 1 
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
  end
end;

for i = 1:2:length(varargin)
    Param = varargin{i};
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
        case 'indices'
            Findices=varargin{i+1};
            vorher=size(data,2);
            nachher=length(Findices);
            data=data(Findices);
            disp(['Selecting ' num2str(nachher) ' epochs out of ' num2str(vorher) ' total epochs.'])
        case 'mode'
            Fmode=varargin{i+1};
        case 'disp'
            Fdisp=varargin{i+1};
        otherwise
            disp(['Unknown argument ' varargin{i} '...'])
            pause
    end;
end; 

Fzeitfenster=zeitfenster;
Fdaten=shiftdim(reshape([data.EEG_activity]',length(Fzeitfenster),size(data,2),size([data(1).EEG_activity],1)),2);

%Fzeitfenster=750:1100;
Fregressoren=regressoren;
%Fregressoren(:,1)=random('exp',1,1,size(Fdaten,3));
%Fregressoren(:,2)=random('exp',1,1,size(Fdaten,3));
nReg=size(Fregressoren,2);
r(size(Fdaten,1),length(Fzeitfenster),nReg) = NaN;
p(size(Fdaten,1),length(Fzeitfenster),nReg) = NaN;
R2(size(Fdaten,1),length(Fzeitfenster))     = NaN;

if length(size(Fdaten))<3 %only one channel data
    Fdaten=shiftdim(Fdaten,-1); %add singleton dimension
end;

disp(['Performing singlepoint correlations over ' num2str(size(Fdaten,1)) ' channels and ' num2str(length(Fzeitfenster)) ' datapoints in ' num2str(size(Fdaten,3)) ' epochs with ' num2str(nReg) ' different regressors using mode ' Fmode '.'])
%pause(1)

%Check if regressor fits data
if size(Fregressoren,1)~=size(Fdaten,3)
    disp('Error: Regressor length must equal input data epoch number.');pause;
end;
fprintf(1,['Progress: Percent of regressions done:   ']);

            
for DP=1:length(Fzeitfenster) %loop through all selected datapoints
    do=floor((DP/length(Fzeitfenster))*100);
    if do < 10
        fprintf(2,'\b%d', do);
    else
        fprintf(2,'\b\b%d', do);
    end;
    for channel=1:size(Fdaten,1) %loop through all the channels
        if strcmp(Fmode, 'glm') %Use GLM FIT
            [b,dev,stats]   = glmfit(Fregressoren, squeeze(Fdaten(channel,DP,:)));
            R2(channel,DP)  = stats.sfit;           %save output without constant term for fitting
            r(channel,DP,:) = stats.beta(2:end);    %save output without constant term for beta values
            p(channel,DP,:) = stats.p(2:end);       %save output without constant term for p statistics
            
        elseif strcmp(Fmode, 'regress') %Use 'regress'
            [B,BI,R,RI,STATS]   = regress(squeeze(Fdaten(channel,DP,:)), [ones(size(Fregressoren,1),1) Fregressoren]);
            p(channel,DP,:)     = STATS(3);         %save output for p statistics
            R2(channel,DP)      = STATS(1);         %save output for fitting
            r(channel,DP,:)     = B(2:end);         %save output without constant term for beta values
            
        elseif strcmp(Fmode(1:6), 'robust') %Use 'robustfit1', return error variance for each regressor separately
            [B,STATS] = robustfit([Fregressoren], squeeze(Fdaten(channel,DP,:)));
             if strcmp(Fmode, 'robust1')
                r(channel,DP,:)=B(2:end);
                p(channel,DP,:)=STATS.se(2:end);     %save output without constant term for p statistics
            else
                r(channel,DP,:)=STATS.t(2:end);     %save output without constant term for beta values
                p(channel,DP,:)=STATS.p(2:end);     %save output without constant term for p statistics
            end;
            
            %Calculate R2 estimate
            sse = STATS.dfe * STATS.robust_s^2;
            phat = B(1);
            for pc = 1 : length(B)-2
                phat = phat + B(pc+1)*Fregressoren(:,pc+1);
            end;
            ssr = norm(phat-mean(phat))^2;
            R2(channel,DP) = 1 - sse / (sse + ssr);

            
        else %Use simple correlation coefficients (which leads to the same result, but is not scaled to microvolt)
            [r(channel,DP,:) p(channel,DP,:)]=corr(squeeze(Fdaten(channel,DP,:)), [Fregressoren]);
        end;
    end;
end;fprintf('\n')
return