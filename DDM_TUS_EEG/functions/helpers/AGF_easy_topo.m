function [ ] = AGF_easy_topo( loadpath, savepath, steps, varargin )
%Function that plots topos and es very easy to use, nothing more niothing less.
%AGF, 2011
%__________________________________________________________________________
%Input:
%'loadpath'     -   Folder with datasets to plot
%'savepath'     -   Folder where to save plots
%'steps'        -   Stepsize for interval of plotting, can be set to string
%                   values that will force this exact (or lower) number of plots.
%'cmap'         -   Path to a colormap to load and use, default AGF_cmap
%'addtit'       -   String to attach to the title for all files!
%'FixLim'       -   Specify maplimits here (note: these have to be symmetrical, otherwise colormaps are incorrect)
%                   Example: 'FixLim', [-1 1].
%'Mask'         -   Specify a range of values that is set to 0 (example [-0.05 0.05])
%'filetype'     -   Specify fieltype to plot graphic in (epsc, jpg, bmp...)
%Output:        -   Just saves the topos!
%__________________________________________________________________________

%Set standard values and adjust custom input parameters.
nargs      = nargin-1;
Fcmap = load('/Volumes/AGF work/SMAC/1auswertung/0myfunctions/data/AGF_cmap.mat');
Fcmap=Fcmap.AGF_cmap;

%Set default values
Ftitle  = [];       %No default title addition.
Ffixlim = [];       %Estimate maplimits from the data.
Fmask   = [];       %Do not mask out any datapoints.
Fsavetyp= 'jpg';    %Default filetype: jpg.
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
     case 'cmap'
      Fcmap = lower(Value);
     case 'addtit'
      Ftitle = lower(Value);
    case 'fixlim'
      Ffixlim = lower(Value);
    case 'mask'
      Fmask = lower(Value);  
      if length(Fmask) > 2
          disp('Error: Mask has to have two or one entries');return;
      end;
      if Fmask(1)>Fmask(2)
          disp('Error: Entry 1 of mask argument has to be the smaller one');return;
      end;
    case 'filetype'
      Fsavetyp = lower(Value); 
    end;
  end;
end;

current=pwd;
for CC=1:length(loadpath)
    clear liste
    if exist(loadpath{CC})==7 %Savepath is existent
        cd(loadpath{CC})
        liste=dir('*.set');
        cd(current)
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        for c=1:size(liste,1)
            if size(ALLEEG,2)>0
                ALLEEG = pop_delset( ALLEEG, [1 2] );
            end
            if Ftitle
                TIT=[liste(c).name(1:end-4) Ftitle];
            else
                TIT=liste(c).name(1:end-4);
            end;
            EEG = pop_loadset( 'filename', [liste(c).name], 'filepath', loadpath{CC});
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );
            
            if ~isempty(Fmask)
                if length(Fmask) == 1
                    EEG.data([EEG.data]==Fmask(1)) = 0;
                    disp(['Masking out value: ' num2str(Fmask(1))])
                else
                    EEG.data([EEG.data]<Fmask(2) & [EEG.data]>Fmask(1)) = 0;
                    disp(['Masking out values between ' num2str(Fmask(1)) ' and ' num2str(Fmask(2))])
                end;
            end;
            
            if ischar(steps)
                Fsteps = str2num(steps);
                if EEG.pnts <= Fsteps*(1000/EEG.srate)
                    Fsteps = 1000/EEG.srate;
                else
                    thresh = 1000/EEG.srate; n_pic = inf;
                    while n_pic > Fsteps
                        %find stepsize that leads to no more than STEPS plots
                        n_pic = EEG.pnts / thresh;
                        thresh = thresh + 1000/EEG.srate;
                    end;
                    Fsteps=thresh;
                end;
            else
                Fsteps = steps;
            end;
            times=[EEG.xmin*1000:Fsteps:EEG.xmax*1000];
            %EEG = pop_rmbase( EEG, [1 size(EEG.data,2)*2]);%size(EEG.data,2)*2
            close all; 
            if isempty(Ffixlim)
                pop_topoplot(EEG,1, times , TIT , 0, 'electrodes', 'on');
            else
                pop_topoplot(EEG,1, times , TIT , 0, 'electrodes', 'on', 'maplimits', fixlim);
            end;
            colormap(Fcmap)
            sname = EEG.setname;
            sname(sname=='.')=','; %remove dots from savenames
            saveas(gca, [savepath sname], Fsavetyp);
        end;
    else
        disp(['Error: Specified input path (' loadpath{CC} ') is not existent...'])
    end;
end;
