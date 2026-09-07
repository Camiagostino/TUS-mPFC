function [ Odata ] = AGF_TRIALPLOT( indata, Trials, Categ, time_w, Ele, varargin)
%%
%Function that runs matlab integrated SVM on an EEG dataset.
%
%INPUT:
%'indata'       either EEG structure or path to eeg file that will be
%               openend.
%'trial_ind'    indexes the trials of interest.
%'time_w'       time window (ms) that should be classified.
%
%OPTIONAL INPUT:
%'options'      A struct that can contain the fields:
%'limits'       plotlimits [-y y]
%'reverse'      1 = reverse Y axis
%               .

%%
if strcmp(computer, 'MACI64')
    Fcmap = load('/Users/Erdnase/Documents/Arbeit/0myfunctions/data/AGF_cmap.mat');
else
    Fcmap = load('/home/afischer/Documents/MATLAB/myfunctions/data/AGF_cmap.mat');
end;
Fcmap=Fcmap.AGF_cmap;
UseE    = find(strcmpi({indata.chanlocs.labels},Ele));
colors={'m' 'g' 'b' 'c' 'r' 'k'};
transp = 0;
Confidence=0.99;
Flabels = {};
Fvis = 1;
savep = '';
Freverse = 1;
Flimits=[];

if isempty(Trials)
    Trials = 1 : size(indata.data,3);
end;

if Confidence
    ConfFactor = norminv(1+(0.5-(1-Confidence/2)),0,1);
    TitStr = [' (shade = ' num2str(Confidence*100) '% CI)'];
else
    ConfFactor = 1;
    TitStr = [' (shade = sem)'];
end;

for hide = 1:1 %Get Variables from varargin array
    nargs = nargin-5;
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
            case 'labels'
                Flabels=varargin{i+1};
            case 'confidence'
                Confidence=varargin{i+1};
            case 'savep'
                savep=varargin{i+1};
            case 'reverse'
                Freverse=varargin{i+1};
            case 'limits'
                Flimits=varargin{i+1};
            case 'visible'
                Fvis=varargin{i+1};
            otherwise
                disp(['Unknown argument ' varargin{i} '...'])
                pause
        end;
    end; 
end;
u = unique(Categ);
u = u(~isnan(u));
n = length(u);
Time_in_DP = find([indata.times]==time_w(1)) : find([indata.times]==time_w(2));
PlotTime=indata.times(Time_in_DP);

%Prepare data and plot it
This_EEG    = squeeze(nanmean(indata.data(UseE,Time_in_DP,Trials(find(Categ==u(1)))),3));
This_SD     = ConfFactor*(squeeze(nanstd(indata.data(UseE,Time_in_DP,Trials(find(Categ==u(1)))),0,3))./sqrt(length(find(Categ==u(1)))-1));
hFig = figure;
if Fvis == 0
    set(gcf,'Visible', 'off'); 
end;
AGF_shadedErrorBar(PlotTime, This_EEG, This_SD, colors{1}, transp);
title(['\bf' num2str(n) ' categories at electrode ' Ele '\rm' TitStr], 'FontSize', 15)
hold;
if n > 1
    for c = 2 : n
        This_EEG    = squeeze(nanmean(indata.data(UseE,Time_in_DP,Trials(find(Categ==u(c)))),3));
        This_SD     = ConfFactor*(squeeze(nanstd(indata.data(UseE,Time_in_DP,Trials(find(Categ==u(c)))),0,3))./sqrt(length(find(Categ==u(c)))-1));
        AGF_shadedErrorBar(PlotTime, This_EEG, This_SD, colors{c}, transp);
    end;
end;
A=gca; A.XLim=[PlotTime(1) PlotTime(end)]


if length(Flabels)>0
    L = legend(Flabels );
    %Append trial number information to legend
    for c = 1 : length(L.String)
        L.String{c}=[L.String{c} ' (n = ' num2str(length(find(Categ==u(c)))) ')'];
    end;
end;

if Freverse
    set(gca,'YDir','reverse');
end;

if ~isempty(Flimits)
    a = ver;
    if str2num(a(1).Release(5:6)) < 14
        set(gca,'ylim', Flimits);
    else
        A = gca;
        A.YLim = Flimits;
    end;
end;

if ~isempty(savep)
    saveas (hFig,savep,'epsc');
end;

