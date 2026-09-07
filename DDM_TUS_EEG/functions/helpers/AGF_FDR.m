function [ corrected_data, crit_p ] = AGF_FDR( data, varargin )
%Function that will perform subject wise ttests over a given dataset and correct results for multiple comparisons using
%selected methods. All nonsignificant values will be set to 0 (or NaN if selected), the rest is left unchanged.
%AGF, 2011
%__________________________________________________________________________
%Input:
%'data'     -   Required input has to be a 3d Matrix, test will be
%               performed over the last dimension.
%'method'   -   ['BH' / 'BY' / 'BKY' / 'Bonferroni' / 'Sidak' / n / 't_map']
%               Method for correction and limit, n=any number <1 that is used as threshold. tmap simply plots maps of t values without correction.
%'inputN'   -   Can be set to any numer, will overwrite the number of
%               comparisons used to estimate threshhold for Bonferroni or Sidak correction.
%'out_Value'-   Can be set to 0, NaN or any number. Will be the value all
%               nonsignificant value are set to.
%'alpha'    -   Set alpha level to different value (standard = 0.05), also used as q value if FDR is selected.
%'display'  -   Toggels display of percentage done [0 = off | 1 = on].
%'range'    -   Vector of datapoints to be used for thresholding p-values.
%               Not that all other datapoints will be corrected with the threshold in the
%               given timewindow (example: Epoch = 0ms-600ms, timewindow for correction:
%               200-400ms, but everything will be thresholded. FWER correction is only
%               valid in this timewindow!)
%'uplim'    -   Upper limit for threshold (ex. 0.001). Can be necessary if a
%               lot of tests reject 0 hypothesis which may increase threshold to
%               undesirably high values (only supported for BH and BY).
%
%Output:
%'corrected_data'    is a vector of the same size as input 'data', but all
%                    nonsignificant fields are set to value 0 (or 'out_Value' if specified).
%'thresh'            Threshold for p-values to be determined as significant (H0 rejected).
%__________________________________________________________________________


%Set standard values and adjust custom input parameters.
nargs      = nargin-1;
Fmethod    = 'bh';                         % Standard method is FDR BH.
FinputN    = size(data,1)*size(data,2);    % Number of multiple comparisons equals all input data.
Fout_Value = 0;                            % Replace nonsignificant elemets with 0.
Falpha     = 0.05;                         % Standard alpha is 0.05.
Fdisp      = 1;                            % Default: Does not display progress.
Frange     = [];                           % Default range for inculsion in multiple comparisons is the whole epoch.
Fuplim     = [];                           % No upper limit for critical p value is set.

if nargs > 1
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
  end
  for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
        case 'method'
            Fmethod = lower(Value);
            if ~strcmp(Fmethod, 'bh') && ~strcmp(Fmethod, 'by') && ~strcmp(Fmethod, 'bky') && ~strcmp(Fmethod, 'bonferroni') && ~strcmp(Fmethod, 'sidak') && ~strcmp(Fmethod, 't_map') && ~isnumeric(Fmethod)
                disp(['Error: Correction method: ' Fmethod ' unknown.'])
                return
            end;
        case 'inputn'
            FinputN = lower(Value);
        case 'out_value'
            Fout_Value = lower(Value);  
        case 'display'
            Fdisp = lower(Value);  
        case 'alpha'
            Falpha = lower(Value);
        case 'uplim'
            Fuplim = lower(Value);    
        case 'range'
            Frange = lower(Value);
            if length(Frange)>size(data,2)
                disp(['Error: Size of range must be smaller than epoch length!'])
                return
            end;
    end;
  end;
end;
disp(['Correcting input dataset with ' num2str(size(data,1)) ' channels, ' num2str(size(data,2)) ' datapoints and ' num2str(size(data,3)) ' subjects for: '...
    num2str(FinputN) ' multiple comparisons using method: ' upper(Fmethod) '.'])

%Calculat n-ttests against 0.
for c1=1:size(data,1)
    if Fdisp
        prct=(c1/size(data,1))*100;
        disp([num2str(prct) ' % of ttests done.'])
    end;
    for c2=1:size(data,2)
        [h,p(c1,c2),ci, stats]=ttest(data(c1,c2,:));
        t(c1,c2)=stats.tstat;
    end;
end;

if Frange
    disp('Reducing datapoint range to input range')
    save_p=p;
    p=p(:,Frange);
    save_t=t;
    t=t(:,Frange);
end;
    
%Adjust threshold. 
if isnumeric(Fmethod)
    if Fmethod>1
        disp('Error: Input ''n'' must be smaller 1!'); return;
    end;
    h=p<=Fmethod; crit_p=Fmethod;
    proc_acc=size(find(h==1),1)/numel(h)*100;
elseif strcmp(Fmethod, 't_map')
    corrected_data=t;return;
elseif strcmp(Fmethod, 'bonferroni')
    %Caculate desired threshold for p
    crit_p = Falpha / FinputN;
    h=p<=crit_p;
    proc_acc=size(find(h==1),1)/numel(h)*100;
elseif strcmp(Fmethod, 'sidak')    
    %Caculate desired threshold for p
    crit_p = 1-(1-Falpha)^(1/FinputN);
    h=p<=crit_p;
    proc_acc=size(find(h==1),1)/numel(h)*100;
elseif strcmp(Fmethod, 'bh')
    if Fdisp
        disp('Determine p-value threshold...')
        [h, crit_p]=fdr_bh(p, Falpha, 'yes', 'pdep', Fuplim);
    else
        [h, crit_p]=fdr_bh(p, Falpha, 'no', 'pdep', Fuplim);
    end;
    proc_acc=size(find(h==1),1)/numel(h)*100;
elseif strcmp(Fmethod, 'by')
    if Fdisp
        disp('Determine p-value threshold...')
        [h, crit_p]=fdr_bh(p, Falpha, 'yes', 'dep', Fuplim);
    else
        [h, crit_p]=fdr_bh(p, Falpha, 'no', 'dep', Fuplim);
    end;
    proc_acc=size(find(h==1),1)/numel(h)*100;
elseif strcmp(Fmethod, 'bky')
    if Fdisp
        disp('Determine p-value threshold...')
        [h, crit_p]=fdr_bky(p, Falpha, 'yes');
    else
        [h, crit_p]=fdr_bky(p, Falpha, 'no');
    end;
    proc_acc=size(find(h==1),1)/numel(h)*100;
end;

%if range for FDR determination is differet from epoch length all other values will be corrected using the same critical p value as determined
%(this is for plotting purpose, the non included timerange is not guaranteed to be corrected for FWER of course!)
if Frange
   if Fdisp
      disp(['Correcting all datapoints with p value determined in specified range from ' num2str(Frange(1)) ' to ' num2str(Frange(end)) '. WARNING: FWER correction not guaranteed!'])
   end    
   h=save_p<=crit_p;
end;

%Reject data at nonsignificant points.
for c1=1:size(data,1)
    if Fdisp
        prct=(c1/size(data,1))*100;
        %disp([num2str(prct) ' % of data have been analysed for rejection.'])
    end;
    for c2=1:size(data,2)
        if ~h(c1,c2)
            data(c1,c2,:)=Fout_Value;
        end;
    end;
end;

corrected_data=data;
disp([num2str(100 - proc_acc) '% of the data have been rejected and set to ' num2str(Fout_Value) '.'])
return

    
    
    
    