function [r p ft] = allGLM( data, regressoren, zeitfenster, varargin )
%%
%Function that performs singlepointcorrelations for singletrial EEG data on
%all channels.
%INPUT:
%'indices'      Equals an array of selected epochs. Only these selected
%               epochs will be used for calculation.
%'mode'         Can be set to 'GLM', 'regress', 'robust1/2' or 'COR', will use either glmfit, regress, robustfit or
%               corr function of matlab (default is 'GLM'). Note: Function
%               returns beta-values when set to GLM, robust and Regress and
%               correlation-coefficients when set to COR. 
%               Note: Regress only sopports return of p statistics for the full model.
%               If robust1 is specified, error variance is returned for later normalization.  
%               If robust2 is specified, error variance is automatically
%               used to normalize b values and normalized values are returned.
%'Nreturn'      Number of Regressors to be included for returnin, default
%               only first. 
%'Disp'         Toggle display of percentage [0/1].
%'R2'           Returns also an R2 or R2 estimate for the data [0/1].

%Set defaults
Fmode       = 'glm';    % Default will use glmfit.
Fnreturn    = 'first';  % Default returns only first regressor (should be changed).
Fdisp       = 0;        % Default: Display is deactivated.
Fr2         = 0;        % Default: No r2 is computed.

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
        case 'nreturn'
            Fnreturn=varargin{i+1};
        case 'disp'
            Fdisp=varargin{i+1};
        case 'r2'
            Fr2=varargin{i+1};
            if Fr2 ~= 1 && Fr2 ~= 0
                disp('Error: To activate R2 calculation set value to 1, else set to 0. Will deactivate R2...')
                Fr2 = 0;
            end;
    end;
end; 

Fzeitfenster=zeitfenster;
Fdaten=shiftdim(reshape([data.EEG_activity]',length(Fzeitfenster),size(data,2),size([data(1).EEG_activity],1)),2);

%Fzeitfenster=750:1100;
Fregressoren=regressoren;
%Fregressoren(:,1)=random('exp',1,1,size(Fdaten,3));
%Fregressoren(:,2)=random('exp',1,1,size(Fdaten,3));
nReg=size(Fregressoren,2);
if strcmp(Fnreturn,'first')
    r(size(Fdaten,1),length(Fzeitfenster),1)=NaN;
    p(size(Fdaten,1),length(Fzeitfenster),1)=NaN;
    ft(size(Fdaten,1),length(Fzeitfenster))=NaN;
else
    r(size(Fdaten,1),length(Fzeitfenster),nReg)=NaN;
    p(size(Fdaten,1),length(Fzeitfenster),nReg)=NaN;
    ft(size(Fdaten,1),length(Fzeitfenster))=NaN;
end;
if strcmp(Fmode(1:6), 'robust')
    ft=r;
    t=ft;
end;

if Fr2
    R2 = r;
end;

if length(size(Fdaten))<3 %only one channel data
    Fdaten=shiftdim(Fdaten,-1); %add singleton dimension
end;

% size(Fdaten)
% figure
% plot(Fdaten(22,:,7))
% brzk

disp(['Performing singlepoint correlations over ' num2str(size(Fdaten,1)) ' channels and ' num2str(length(Fzeitfenster)) ' datapoints in ' num2str(size(Fdaten,3)) ' epochs with ' num2str(nReg) ' different regressors using mode ' Fmode '.'])
pause(2)

%Check if regressor fits data
if size(Fregressoren,1)~=size(Fdaten,3)
    disp('Error: Regressor length must equal input data epoch number.');pause;
end;

for DP=1:length(Fzeitfenster) %loop through all selected datapoints
    if Fdisp==1
        disp([num2str(DP/length(Fzeitfenster)*100) '% done.'])
    end;
    for channel=1:size(Fdaten,1) %loop through all the channels
        if strcmp(Fmode, 'glm') %Use GLM FIT
            [b,dev,stats] = glmfit(Fregressoren, squeeze(Fdaten(channel,DP,:)));
            ft(channel,DP)=stats.sfit;      %save output without constant term for fitting
            if strcmp(Fnreturn,'first')
                r(channel,DP,:)=stats.beta(2);  %save output without constant term for beta values
                p(channel,DP,:)=stats.p(2);     %save output without constant term for p statistics
            else
                r(channel,DP,:)=stats.beta(2:end);  %save output without constant term for beta values
                p(channel,DP,:)=stats.p(2:end);     %save output without constant term for p statistics
            end;
            
        elseif strcmp(Fmode, 'regress') %Use 'regress'
            [B,BINT,R,RINT,STATS] = regress(squeeze(Fdaten(channel,DP,:)), [ones(size(Fregressoren,1),1) Fregressoren]);
            p(channel,DP,:)=STATS(3);           %save output for p statistics
            ft(channel,DP)=STATS(1);            %save output for fitting
            if strcmp(Fnreturn,'first')
                r(channel,DP,:)=B(2);           %save output without constant term for beta values
            else
                r(channel,DP,:)=B(2:end);       %save output without constant term for beta values
            end;
            
        elseif strcmp(Fmode(1:6), 'robust') %Use 'robustfit1', return error variance for each regressor separately
            [B,STATS] = robustfit([Fregressoren], squeeze(Fdaten(channel,DP,:)));
            if strcmp(Fnreturn,'first')
                r(channel,DP,:)=B(2);           %save output without constant term for beta values
                p(channel,DP,:)=STATS.p(2);     %save output without constant term for p statistics
                ft(channel,DP,:)=STATS.se(2);   %Standard error of residuals
                t(channel,DP,:)=STATS.t(2);   %Standard error of residuals
            else
                r(channel,DP,:)=B(2:end);           %save output without constant term for beta values
                p(channel,DP,:)=STATS.p(2:end);     %save output without constant term for p statistics
                t(channel,DP,:)=STATS.t(2:end);   %Standard error of residuals
            end;
            
            if Fr2 %R2 should be calculated
                sse = stats_rob.dfe * stats_rob.robust_s^2;
                phat = B(1);
                for pc = 2 : length (B)
                    phat = phat + B(pc)*x;
                end;
                ssr = norm(phat-mean(phat))^2;
                R2(channel,DP,:) = 1 - sse / (sse + ssr);
            end;
            
        else %Use simple correlation coefficients (which leads to the same result, but is not scaled to microvolt)
            [r(channel,DP,:) p(channel,DP,:)]=corr(squeeze(Fdaten(channel,DP,:)), [Fregressoren]);
        end;
    end;
end;
if strcmp(Fmode, 'robust2') %Automatically normalize b with error variance 
    %r=r./ft;
    r=t; %Same effect but saves time
end;
return