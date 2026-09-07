function [EEG1, rmeans] = AGF_rmbase(EEG1, Fzeit, varargin)
%Function that performs baselinecorrection and allows to take the baseline
%of a separate Dataset. Note: To verify identity of epochs, EEG.epoch.trnr
%is used. Epochs that are not shared between the datasets are deleted.
%
%Required:
%'EEG1' =   Standard EEG Dataset that will be corrected (standard base_rm if
%               no other dataset is provided).
%'zeit' =   Timewindow in ms (example: [-200 -50], default is whole epoch: []). 
%               Can also be set to vector of datapoints ([100:225]).
%
%Optional: 
%'EEG2' =   Dataset to take baseline from. 'zeit' has to be specified for
%               this dataset.
%
%Output: 
%dataout=   Corrected dataset.
%rmeans =   Values subtracted on each channel and epoch.
%
%AG Fischer, 2011

%SET DEFAULT VALUES
nargs   = nargin-2;
EEG2    = [];        %No second dataset.


if nargs > 1 
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
  end
end;

for i = 1:2:length(varargin)                                            %Get Variables from varargin array
    Param = varargin{i};
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
     case 'eeg2'
         EEG2=varargin{i+1};
    end;
end; 


if isempty(EEG2)
    if numel(Fzeit)==0                                                      %Check if 'zeit' is time in ms, whole epoch or vector.
        disp('Using whole epoch...')
        Fzeit = 1:EEG1.pnts;
    elseif numel(Fzeit)==2
        disp(['Subtracting baseline from ' num2str(Fzeit(1)) ' to ' num2str(Fzeit(2)) 'ms.'])  
        Fzeit = round((Fzeit(1)/1000-EEG1.xmin)*EEG1.srate+1):round((Fzeit(2)/1000-EEG1.xmin)*EEG1.srate);
    elseif numel(Fzeit)>2
        disp(['Subtracting baseline from datapoint ' num2str(Fzeit(1)) ' to datapoint ' num2str(Fzeit(end)) '.'])
    end;
    disp('Standard baselinecorrection for one dataset.')
    [EEG1]=pop_rmbase(EEG1, [], Fzeit); rmeans=[]; return;

else
    if numel(Fzeit)==0                                                      %Check if 'zeit' is time in ms, whole epoch or vector.
        disp('Using whole epoch of dataset2...')
        Fzeit = 1:EEG2.pnts;
    elseif numel(Fzeit)==2
        disp(['Subtracting baseline from ' num2str(Fzeit(1)) ' to ' num2str(Fzeit(2)) 'ms taken from dataset2.'])  
        Fzeit = round((Fzeit(1)/1000-EEG2.xmin)*EEG2.srate+1):round((Fzeit(2)/1000-EEG2.xmin)*EEG2.srate);
    elseif numel(Fzeit)>2
        disp(['Subtracting baseline from datapoint ' num2str(Fzeit(1)) ' to datapoint ' num2str(Fzeit(end)) ' taken from dataset2.'])
    end;

    %Index frames of dataset2
    for cepoch=1:length(EEG1.epoch)
        Aevent0(cepoch)=find([EEG1.epoch(cepoch).eventlatency{:}]==0);  %Get index for event at timepoint 0 in this epoch.
        Atrial(cepoch)=EEG1.epoch(cepoch).eventtrnr{Aevent0(cepoch)};   %Read out original trialnumber for this event.        
    end;
    %Same thing for dataset2
    for cepoch=1:length(EEG2.epoch)
        Bevent0(cepoch)=find([EEG2.epoch(cepoch).eventlatency{:}]==0);  %Get index for event at timepoint 0 in this epoch.
        Btrial(cepoch)=EEG2.epoch(cepoch).eventtrnr{Bevent0(cepoch)};   %Read out original trialnumber for this event.        
    end;
    
    miss = setdiff(Atrial, Btrial);                                     %Find trials that cannot be found in EEG2
    if numel(miss)>0
        for c=1:numel(miss)
            missI(c)=find(Atrial==miss(c));                             %Index these trials.
        end;
        [EEG1]=pop_rejepoch( EEG1, missI ,0);                           %Reject indexed epochs.
        Atrial(missI)=[];                                               %Remove missing trials in Atrial index
    else
        disp('All trials are found in dataset2...')
    end;
    for getind=1:numel(Atrial)
        Bind(getind)=find(Atrial(getind)==Btrial);                      %Index trials in EEG2 (=Btrial) to access them
    end;

    disp('Removing baseline taken from dataset2 from dataset1...')    
    for e=1:size(EEG1.data,3) %count epochs
        for c=1:size(EEG1.data,1) %count channels
            rmeans(c,e)=nan_mean(double(EEG2.data(c,Fzeit,Bind(e))));
            EEG1.data(c,:,e)=[EEG1.data(c,:,e)]-rmeans(c,e);
        end;
    end;
    EEG1.icaact = [];
    return;
end;


    