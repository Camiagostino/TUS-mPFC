function [ ] = new_average( inname, outname, inlocation, outlocation, varargin )
%Will average an inputvektor of names

%defaults
Fbaserem=NaN;       %Default: No baseline correction.
L=length(inname);
Fsetname=outname;   %Default: Setname is the same as the outname.

%Check n-args in
nargs = nargin-6;
if nargs > 1 
  if ~(round(nargs/2) == nargs/2)
    error('Odd number of input arguments??')
    return
  end;
end;
%Read arg in
for i = 1:2:length(varargin)
    Param = varargin{i};
    if ~isstr(Param)
      error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
     case 'baserem'
         Fbaserem=varargin{i+1};
     case 'setname'
         clear Fsetname
         Fsetname=varargin{i+1};
    end;
end;

if L~=length(outname)
    disp('Error: Array sizes do not match...')
    return;
end;

%Adjust input pathes
if size(inlocation,1)==1
    disp(['All Datasets are read from folder: ' inlocation '.'])
    for i=1:L
        Finlocation{i,1}=inlocation;
    end;
elseif size(inlocation,1)>1
    if size(inlocation,1)~=L
        disp('Dimensions of inlocations do not match input datasets numbers.')
        return;
    end;
    disp('Reading datasets from different folders.')
    Finlocation=inlocation;
end;

%Adjust output pathes
if size(outlocation,1)==1
    disp(['All Datasets are saved in folder: ' outlocation '.'])
    for i=1:L
        Foutlocation{i,1}=outlocation;
    end;
elseif size(outlocation,1)>1
    if size(outlocation,1)~=L
        disp('Dimensions of outlocations do not match input datasets numbers.')
        return;
    end;
    disp('Saving datasets in different folders.')
    Foutlocation=outlocation;
end;

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for c=1:L
    EEG = pop_loadset( 'filename', [inname{c} '.set'], 'filepath', Finlocation{c});
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    %If baserem is set, perform baserem before creating Differencewaves
    if isempty(Fbaserem)
        disp('De-meaning the whole epoch.')
        ALLEEG(c) = pop_rmbase( ALLEEG(c), [],1:size(ALLEEG(c).data,2));
    elseif ~isnan(Fbaserem)
        disp(['Removing Baseline from ' num2str(Fbaserem(1)) 'ms to ' num2str(Fbaserem(2)) 'ms.'])
        ALLEEG(c) = pop_rmbase( ALLEEG(c), Fbaserem);
    end;
    EEG = pop_grandaverage(ALLEEG(c), 'datasets', 1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, c, 'overwrite', 'on', 'gui', 'off','setname',Fsetname{c}); 
    EEG = pop_saveset( ALLEEG(c),  'filename', [outname{c} '.set'], 'filepath', Foutlocation{c});
end;

