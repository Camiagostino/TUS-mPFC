function [EEG, indices, com] = AGF_epoch_by_number( EEG, epochs, newname, time ) 
%Function that simply calls pop_epoch eeglab function but selects only
%those epochs, that are indicated in 'epochs' input.
%AGF, 2013
%__________________________________________________________________________
%Input:
%'loadpath'     - Input dataset. Data may already be epoched; in this case,
%                 extract (shorter) subepochs time locked to epoch events.
%'epochs'       - Array of epochs to use. Function reads out the events at
%                 latency 0 and extracts these events.
%__________________________________________________________________________
Ev_array=[];
if isempty(time)
    term1='[EEG.xmin EEG.xmax]';
else
    term1=[ '[' num2str(time(1)) ' ' num2str(time(2)) ']' ];
end;

%Check if dataset has initial epochs, otherwise create them 
if isempty(EEG.epoch)
    EEG = pop_epoch( EEG, {  }, [EEG.xmin EEG.xmax], 'newname', EEG.setname, 'epochinfo', 'yes');
end;

%Get events at latency 0 in epoch structure
for c = 1 : length(EEG.epoch)
    if iscell(EEG.epoch(c).eventlatency(1)) % Check if data is cell or integer array and read that out
        Ev0=find([EEG.epoch(c).eventlatency{:}]==0);
    else
        Ev0=find([EEG.epoch(c).eventlatency(:)]==0);
    end;
    Ev_array=[Ev_array EEG.epoch(c).event(Ev0)];
end;
if isempty(newname)
    [EEG, indices, com] = pop_epoch( EEG, {  }, eval(term1), 'newname', EEG.setname, 'epochinfo', 'yes', 'eventindices', Ev_array(epochs));
else
    [EEG, indices, com] = pop_epoch( EEG, {  }, eval(term1), 'newname', newname, 'epochinfo', 'yes', 'eventindices', Ev_array(epochs));
end;
return