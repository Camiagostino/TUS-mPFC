function [ cor_data, dif_data, t_map ] = pop_differencewaves( indata1, indata2, varargin)
%function that creates differencewaves of the averaged ERPs of two input
%datasets. Inputs must be before averaging! Function can automatically
%equalize the epochnumbers if wanted (random trials are picked)... Function will always
%subtract dataset2 from dataset1 (i.e. if 1 is more negative, result stays negative).
%
%Required inputs
%'indata1'                   | EEG structure
%'indata2'                   | EEG structure
%
%Optional inputs
%'method' = ['subtract'/'p'] | Defines the method of creating the difference-datasets: 
%                            | 'subtract' just calculates the difference between both datasets.
%                            | 'p' calculates a paired/two sample ttest at significancelevel 'alpha'. Nonsignificant datapoints are set to 0.
%'alpha' = 0 - 1             | Significance level for tests (default: 0.05)
%
%Outputs
%'cor_data'                  | Differencewaves correct using predefined alpha values and method.
%'dif_data'                  | Simple complete differencewave.
%'t_map'                     | Returns the t-values at every point in time and electrode for the difference between the datasets.
%AGF, 2013

%defaults
Falpha=0.05;        %Default: Significance level = 0.05.
Fmethod='subtract'; %Default: Datasets be subtracted.

%Check n-args in
nargs = nargin-2;
if nargs > 1 
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
    return
  end;
end;
%Read arg in
for i = 1:2:length(varargin) %Get input parameters
    Param = varargin{i};
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
     case 'alpha'
         Falpha=varargin{i+1};
     case 'method'
         Fmethod=varargin{i+1};
    end;
end;

%Epoch sizes
EP1=size(indata1.epoch,2);
EP2=size(indata2.epoch,2);

%Subtract datasets
if EP1==EP2
    dif_data=indata1.data-indata2.data;
else
    dif_data=zeros(size(indata1.data,1),size(indata1.data,2),2);
    dif_data(:,:,1)=mean(indata1.data,3)-mean(indata2.data,3);
    dif_data(:,:,2)=mean(indata1.data,3)-mean(indata2.data,3);
end;
cor_data=dif_data;
h=zeros(size(indata1.data,1),size(indata1.data,2));
p=zeros(size(h));
t_map=zeros(size(indata1.data,1),size(indata1.data,2),2);

%Generate Differencewaveform and save new dataset
%Display choosen parameters
if strcmp(Fmethod,'subtract')
    disp('Dataset1 will be subtracted from dataset2...')

elseif strcmp(Fmethod,'p')
    if EP1==EP2
        disp(['Performing ' num2str(numel(p)) ' paired ttests of dataset1 and dataset2 at alpha: ' num2str(Falpha) ' - all other datapoints will be set to zero.'])
        %Calculat n-ttests.
        for c1=1:size(indata1.data,1)
            if mod(c1,round(mod(size(indata1.data,1)/10,10)))==0
                %progr=strcat(progr,'.');
                fprintf('.');
            end;
            for c2=1:size(indata1.data,2)
                [h(c1,c2),p(c1,c2), temp, stats]=ttest(indata1.data(c1,c2,:), indata2.data(c1,c2,:), Falpha);
                t_map(c1,c2,1)=stats.tstat;
                t_map(c1,c2,2)=stats.tstat; %Add second epoch, otherwise plotting from EEGLAB gui does not work
            end;
        end;
        fprintf('\n')
    else
        disp(['Performing ' num2str(numel(p)) ' two sample ttests of dataset1 and dataset2 at alpha: ' num2str(Falpha) ' - all other datapoints will be set to zero.'])
        %Calculat n-ttests.
        for c1=1:size(indata1.data,1)
            if mod(c1,round(mod(size(indata1.data,1)/10,10)))==0
                %progr=strcat(progr,'.');
                fprintf('.');
            end;
            for c2=1:size(indata1.data,2)
                [h(c1,c2),p(c1,c2), temp, stats]=ttest2(indata1.data(c1,c2,:), indata2.data(c1,c2,:), Falpha);
                t_map(c1,c2)=stats.tstat;
                t_map(c1,c2,2)=stats.tstat; %Add second epoch, otherwise plotting from EEGLAB gui does not work
            end;
        end;
    end;
    %Reject data at nonsignificant points.
    for c1=1:size(indata1.data,1)
        for c2=1:size(indata1.data,2)
            if ~h(c1,c2)
                cor_data(c1,c2,:)=0;
            end;
        end;
    end;
    fprintf('\n')
else
    disp(['Input argument `method´ set to invalid value: ' Fmethod '.'])
    return
end;

return;