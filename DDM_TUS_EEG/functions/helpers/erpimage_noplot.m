% erpimage() - Plot a colored image of a collection of single-trial data epochs, optionally
%              sorted on and/or aligned to an input sorting variable and smoothed across
%              trials with a Gaussian weighted moving-average. (To return event-aligned data
%              without plotting, use eegalign()).  Optionally sort trials on value, amplitude
%              or phase within a specified latency window. Optionally plot the ERP mean
%              and std. dev.and moving-window spectral amplitude and inter-trial coherence
%              at aselected or peak frequency. Optionally 'time warp' the single trial
%              time-domain (potential) or power data to align the plotted data to a series
%              of events with varying latencies that occur in each trial. Click on
%              individual figures parts to examine them separately and zoom (using axcopy()).
% Usage:
%            >> figure; erpimage(data,[],times); % image trials as colored lines in input order
%
%            >> figure; [outdata,outvar,outtrials,limits,axhndls, ...
%                        erp,amps,cohers,cohsig,ampsig,outamps,...
%                        phsangls,phsamp,sortidx,erpsig] ...
%                            = erpimage(data,sortvar,times,'title',avewidth,decimate,...
%                                  'key', 'val', ...); % use options
% Required input:
%   data     = [vector or matrix] Single-channel input data to image.
%               Formats (1,frames*trials) or (frames,trials)
%
% Optional ordered inputs {with defaults}:
%
%   sortvar  = [vector | []] Variable to sort epochs on (length(sortvar) = nepochs)
%              Example: sortvar may by subject response time in each epoch (in ms)
%              {default|[]: plot in input order}
%   times    = [vector | []] vector of latencies (in ms) for each epoch time point.
%               Else [startms ntimes srate] = [start latency (ms), time points
%               (=frames) per epoch, sampling rate (Hz)]. Else [] -> 0:nframes-1
%               {default: []}
%  'title'   = ['string'] Plot title {default: none}
%   avewidth = [positive scalar (may be non-integer)]. If avg_type is set to 'boxcar'
%               (the default), this is the number of trials used to smooth
%               (vertically) with a moving-average. If avg_type is set to
%               'Gaussian,' this is the standard deviation (in units of
%               trials) of the Gaussian window used to smooth (vertically)
%               with a moving-average.  Gaussian window extends three
%               standard deviations below and three standard deviations above window
%               center (trials beyond window are not incorporated into average). {default: no
%               smoothing}
%   decimate = Factor to decimate|interpolate ntrials by (may be non-integer)
%               Else, if this is large (> sqrt(ntrials)), output this many epochs.
%               {default|0->1}
%
% Optional unordered 'keyword',argument pairs:
%
% Re-align data epochs:
%   'align'  = [latency] Time-lock data to sortvar. Plot sortvar at given latency
%               (in ms). Else Inf -> plot sortvar at median sortvar latency
%               {default: do not align}
%   'timewarp' = {[events], [warpms], {colors}} Time warp ERP, amplitude and phase
%               time-courses before smoothing. 'events' is a matrix whose columns
%               specify the latencies (in ms) at which a series of successive events occur
%               in each trial. 'warpms' is an optional vector of latencies (in ms) to which
%               the series of events should be time locked. (Note: Epoch start and end
%               should not be declared as events or warpms}. If 'warpms' is absent or [],
%               the median of each 'events' column will be used. {colors} contains a
%               list of Matlab linestyles to use for vertical lines marking the occurrence
%               of the time warped events. If '', no line will be drawn for this event
%               column. If fewer colors than event columns, cycles through the given color
%               labels.  Note: Not compatible with 'vert' (below).
%   'renorm' = ['yes'|'no'| formula] Normalize sorting variable to epoch
%               latency range and plot. 'yes'= autoscale. formula must be a linear
%               transformation in the format 'a*x+b'
%               Example of formula: '3*x+2'. {default: 'no'}
%               If sorting by string values like event type, suggested formulas for:
%                 letter string: '1000*x', number string: '30000*x-1500'
%   'noplot' = ['on'|'off'] Do not plot sortvar {default: Plot sortvar if in times range}
%   'NoShow' = ['on'|'off'] Do not plot erpimage, simply return outputs {default: 'off'}
%
% Sort data epochs:
% 'nosort'       = ['on'|'off'] Do not sort data epochs. {default: Sort data epochs by
%                  sortvar (see sortvar input above)}
% 'replace_ties' = ['yes'|'no'] Replace trials with the same value of
%                  sortvar with the mean of those trials.  Only works if sorting trials
%                  by sortvar. {default: 'no'}
% 'valsort'      = [startms endms direction] Sort data on (mean) value
%                  between startms and (optional) endms. Direction is 1 or -1.
%                  If -1, plot max-value epoch at bottom {default: sort on sortvar}
% 'phasesort'    = [ms_center prct freq maxfreq topphase] Sort epochs by phase in
%                  a 3-cycle window centered at latency ms_center (ms).
%                  Percentile (prct) in range [0,100] gives percent of trials
%                  to reject for (too) low amplitude. Else, if in range [-100,0],
%                  percent of trials to reject for (too) high amplitude;
%                  freq (Hz) is the phase-sorting frequency. With optional
%                  maxfreq, sort by phase at freq of max power in the data in
%                  the range [freq,maxfreq] (Note: 'phasesort' arg freq overrides
%                  the frequency specified in 'coher'). With optional topphase,
%                  sort by phase, putting topphase (degrees, in range [-180,180])
%                  at the top of the image. Note: 'phasesort' now uses circular
%                  smoothing. Use 'cycles' (below) for wavelet length.
%                  {default: [0 25 8 13 180]}
%  'ampsort'     = [center_ms prcnt freq maxfreq]  Sort epochs by amplitude.
%                  (See 'phasesort' above). If ms_center is 'Inf', then sorting
%                  is by mean power across the time window specified by 'sortwin'
%                  below. If third arg, freq, is < 0, sort by mean power in the range
%                  [ abs(freq)   maxfreq ].
%  'sortwin'     = [start_ms end_ms] If center_ms == Inf in 'ampsort' arg (above), sorts
%                  by mean amplitude across window centers shifted from start_ms
%                  to end_ms by 10 ms.
%  'showwin'     = ['on'|'off'] Show sorting window behind ERP trace. {default: 'off'}
%
% Plot time-varying spectral amplitude instead of potential:
% 'plotamps'     = ['on'|'off'] Image amplitudes at each trial and latency instead of 
%                  potential values. Note: Currently requires 'coher' (below) with alpha signif.
%                  Use 'cycles' (see below) > (its default) 3 for better frequency specificity,
%                  {default: plot potential, not amplitudes, with no minimum}. The average power
%                  (in log space) before time 0 is automatically removed. Note that the 
%                  'baseline' parameter has no effect on 'plotamps'. Instead use
%                  change "baselinedb" or "basedB" in the 'limits' parameter. By default
%                  the baseline is removed before time 0.
%
% Specify plot parameters:
%   'limits'         = [lotime hitime minerp maxerp lodB hidB locoher hicoher basedB]
%                      Plot axes limits. Can use NaN (or nan, but not Nan) for missing items
%                      and omit late items. Use last input, basedB, to set the
%                      baseline dB amplitude in 'plotamps' plots {default: from data}
%   'sortvar_limits' = [min max] minimum and maximum sorting variable
%                      values to image. This only affects visualization of
%                      ERPimage and ERPs (not smoothing).  Cannot be used
%                      if sorting by any factor besides sortvar (e.g.,
%                      phase).
%   'signif'         = [lo_dB, hi_dB, coher_signif_level] Use precomputed significance
%                      thresholds (as from outputs ampsig, cohsig) to save time. {default: none}
%   'caxis'          = [lo hi] Set color axis limits. Else [fraction] Set caxis limits at
%                      (+/-)fraction*max(abs(data)) {default: symmetrical in dB, based on data limits}
%
% Add epoch-mean ERP to plot:
%   'erp'      = ['on'|'off'|1|2|3|4] Plot ERP time average of the trials below the
%                image.  If 'on' or 1, a single ERP (the mean of all trials) is shown.  If 2,
%                two ERPs (super and sub median trials) are shown.  If 3, the trials are split into
%                tertiles and their three ERPs are shown.  If 4, the trials are split into quartiles
%                and four ERPs are shown. Note, if you want negative voltage plotted up, change YDIR
%                to -1 in icadefs.m.  If 'erpalpha' option is used, any values of 'erp' greater than
%                1 will be reset to 1. {default no ERP plotted}
%   'erpalpha' = [alpha] Visualizes two-sided significance threshold (i.e., a two-tailed test) for the
%                null hypothesis of a zero mean, symmetric distribution (range: [.001 0.1]). Thresholds
%                are determined via a permutation test. Requires 'erp' to be a value other than 'off'.
%                If 'erp' is set  to a value greater than 1, it is reset to 1 to increase plot readability.
%                {default: no alpha significance thresholds plotted}
%   'erpstd'   = ['on'|'off'] Plot ERP +/- stdev. Requires 'erp' {default: no std. dev. plotted}
%   'erp_grid' = If 'erp_grid' is added as an option voltage axis dashed grid lines will be
%                 added to the ERP plot to facilitate judging ERP amplitude
%   'rmerp'    = ['on'|'off'] Subtract the average ERP from each trial before processing {default: no}
%
% Add time/frequency information:
%   'coher'  = [freq] Plot ERP average plus mean amplitude & coherence at freq (Hz)
%               Else [minfrq maxfrq] = same, but select frequency with max power in
%               given range. (Note: the 'phasesort' freq (above) overwrites these
%               parameters). Else [minfrq maxfrq alpha] = plot coher. signif. level line
%               at probability alpha (range: [0,0.1]) {default: no coher, no alpha level}
%   'srate'  = [freq] Specify the data sampling rate in Hz for amp/coher (if not
%               implicit in third arg, times) {default: as defined in icadefs.m}
%   'cycles' = [float] Number of cycles in the wavelet time/frequency decomposition {default: 3}
%
% Add plot features:
%   'cbar'           = ['on'|'off'] Plot color bar to right of ERP-image {default no}
%   'cbar_title'     = [string] The title for the color bar (e.g., '\muV' for
%                       microvolts).
%   'topo'           = {map,chan_locs,eloc_info} Plot a 2-D scalp map at upper left of image.
%                       map may be a single integer, representing the plotted data channel,
%                       or a vector of scalp map channel values. chan_locs may be a channel locations
%                       file or a chanlocs structure (EEG.chanlocs). See '>> topoplot example'
%                       eloc_info (EEG.chaninfo), if empty ([]) or absent, implies the 'X' direction
%                       points towards the nose and all channels are plotted {default: no scalp map}
%   'spec'           = [loHz,hiHz] Plot the mean data spectrum at upper right of image.
%   'specaxis'       = ['log'|'lin] Use 'lin' for linear spectrum frequency scaling, 
%                       else 'log' for log scaling {default: 'log'}
%   'horz'           = [epochs_vector] Plot horizontal lines at specified epoch numbers.
%   'vert'           = [times_vector] Plot vertical dashed lines at specified latencies
%   'auxvar'         = [size(nvars,ntrials) matrix] Plot auxiliary variable(s) for each trial
%                       as separate traces. Else, 'auxvar',{[matrix],{colorstrings}}
%                       to specify N trace colors.  Ex: colorstrings = {'r','bo-','','k:'}
%                       (see also: 'vert' and 'timewarp' above). {default: none}
%   'sortvarpercent' = [float vector] Plot percentiles for the sorting variable
%                       for instance, [0.1 0.5 0.9] plots the 10th percentile, the median
%                       and the 90th percentile.
% Plot options:
% 'noxlabel'          = ['on'|'off'] Do not plot "Time (ms)" on the bottom x-axis
% 'yerplabel'         = ['string'] ERP ordinate axis label (default is ERP). Print uV with '\muV'
% 'avg_type'          = ['boxcar'|'Gaussian'] The type of moving average used to smooth
%                        the data. 'Boxcar' smoothes the data by simply taking the mean of
%                        a certain number of trials above and below each trial.
%                        'Gaussian' does the same but first weights the trials
%                        according to a Gaussian distribution (e.g., nearby trials
%                        receive greater weight).  The Gaussian is better than the
%                        boxcar in that it rather evenly filters out high frequency
%                        vertical components in the ERPimage. See 'avewidth' argument
%                        description for more information. {default: boxcar}
% 'img_trialax_label' = ['string'] The label of the axis corresponding to trials in the ERPimage
%                        (e.g., 'Reaction Time').  Note, if img_trialax_label is set to something
%                        besides 'Trials' or [], the tick marks on this axis will be set in units
%                        of the sorting variable.  This is a useful alternative to plotting the
%                        sorting variable when the sorting variable is not in milliseconds. This
%                        option is not effective if sorting by amplitude, phase, or EEG value. {default: 'Trials'}
% 'img_trialax_ticks' = Vector of sorting variable values at which tick marks (e.g., [300 350 400 450]
%                        for reaction time in msec) will appear on the trial axis of the erpimage. Tick mark
%                        values should be given in units img_trialax_label (e.g., 'Trials' or msec).
%                        This option is not effective if sorting by amplitude, phase, or EEG value.
%                        {default: automatic}
% 'baseline'          = [low_boundary high_boundary] a time window (in msec) whose mean amplitude in
%                        each trial will be removed from each trial (e.g., [-100 0]) after filtering.
%                        Useful in conjunction with 'filt' option to re-basline trials after they have been
%                        filtered. Not necessary if data have already been baselined and erpimage
%                        processing does not affect baseline amplitude {default: no further baselining
%                        of data}.
% 'baselinedb'        = [low_boundary high_boundary] a time window (in msec) whose mean amplitude in
%                        each trial will be removed from each trial (e.g., [-100 0]). Use basedB in limits
%                        to remove a fixed value. Default is time before 0. If you do not want to use a 
%                        baseline for amplitude plotting, enter a NaN value.
% 'filt'              = [low_boundary high_boundary] a two element vector indicating the frequency
%                        cut-offs for a 3rd order Butterworth filter that will be applied to each
%                        trial of data.  If low_boundary=0, the filter is a low pass filter.  If
%                        high_boundary=srate/2, then the filter is a high pass filter.  If both
%                        boundaries are between 0 and srate/2, then the filter is a bandpass filter.
%                        If both boundaries are between 0 and -srate/2, then the filter is a bandstop
%                        filter (with boundaries equal to the absolute values of low_boundary and
%                        high_boundary).  Note, using this option requires the 'srate' option to be
%                        specified and the signal processing toolbox function butter.m.  You should
%                        probably use the 'baseline' option as well since the mean prestimulus baseline
%                        may no longer be 0 after the filter is applied {default: no filtering}
%
% Optional outputs:
%    outdata  = (times,epochsout) data matrix (after smoothing)
%     outvar  = (1,epochsout) actual values trials are sorted on (after smoothing).
%               if 'sortvarpercent' is used, this variable contains a cell array with
%               { sorted_values { sorted_percent1 ... sorted_percentN } }
%   outtrials = (1,epochsout)  smoothed trial numbers
%     limits  = (1,10) array, 1-9 as in 'limits' above, then analysis frequency (Hz)
%    axhndls  = vector of 1-7 plot axes handles (img,cbar,erp,amp,coh,topo,spec)
%        erp  = plotted ERP average
%       amps  = mean amplitude time course
%      coher  = mean inter-trial phase coherence time course
%     cohsig  = coherence significance level
%     ampsig  = amplitude significance levels [lo high]
%    outamps  = matrix of imaged amplitudes (from option 'plotamps')
%   phsangls  = vector of sorted trial phases at the phase-sorting frequency
%     phsamp  = vector of sorted trial amplitudes at the phase-sorting frequency
%    sortidx  = indices of input data epochs in the sorting order
%     erpsig  = trial average significance levels [2,frames]
%
% Example:  >> figure;
%              erpimage(data,RTs,[-400 256 256],'Test',1,1,...
%                            'erp','cbar','vert',-350);
%    Plots an ERP-image of 1-s data epochs sampled at 256 Hz, sorted by RTs, with
%    title ('Test'), and sorted epochs not smoothed or decimated (1,1). Overplots
%    the (unsmoothed) RT latencies on the colored ERP-image. Also plots the
%    epoch-mean (ERP), a color bar, and a dashed vertical line at -350 ms.
%
% Authors: Scott Makeig, Tzyy-Ping Jung & Arnaud Delorme,
%          CNL/Salk Institute, La Jolla, 3-2-1998 -
%
% See also: phasecoher, rmbase, cbar, movav

% Copyright (C) Scott Makeig & Tzyy-Ping Jung, CNL / Salk Institute, La Jolla 3-2-98
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Uses external toolbox functions: phasecoher(), rmbase(), cbar(), movav()
% Uses included functions:         plot1trace(), phasedet()

% UNIMPLEMENTED - 'allcohers',[data2] -> image the coherences at each latency & epoch.
%                 Requires arg 'coher' with alpha significance.
%                 Shows projection on grand mean coherence vector at each latency
%                 and trial. {default: no}

%% LOG COMMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,outsort,outtrials,limits,axhndls,erp,amps,cohers,cohsig,ampsig,allamps,phaseangles,phsamp,sortidx,erpsig] = erpimage(data,sortvar,times,titl,avewidth,decfactor,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24,arg25,arg26,arg27,arg28,arg29,arg30,arg31,arg32,arg33,arg34,arg35,arg36,arg37,arg38,arg39,arg40,arg41,arg42,arg43,arg44,arg45,arg46,arg47,arg48,arg49,arg50,arg51,arg52,arg53,arg54,arg55,arg56,arg57,arg58,arg59,arg60,arg61,arg62,arg63,arg64,arg65,arg66)

%
%% %%%%%%%%%%%%%%%%% Define defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize optional output variables:
warning off;
erp = []; amps = []; cohers = []; cohsig = []; ampsig = [];
allamps = []; phaseangles = []; phsamp = []; sortidx = [];
auxvar = []; erpsig = []; winloc = [];winlocs = [];
timeStretchColors = {};

YES = 1;  % logical variables
NO  = 0;

DEFAULT_BASELINE_END = 0; % ms
TIMEX = 1;          % 1 -> plot time on x-axis;
% 0 -> plot trials on x-axis

BACKCOLOR = [0.8 0.8 0.8]; % grey background
try, icadefs; catch, end;
% read BACKCOLOR for plot from defs file (edit this)
% read DEFAULT_SRATE for coher,phase,allamps, etc.
% read YDIR for plotting ERP
% Fix plotting text and line style parameters
SORTWIDTH = 2.5;    % Linewidth of plotted sortvar
ZEROWIDTH = 3.0;    % Linewidth of vertical 0 line
VERTWIDTH = 2.5;    % Linewidth of optional vertical lines
HORZWIDTH = 2.1;    % Linewidth of optional vertical lines
SIGNIFWIDTH = 1.9;  % Linewidth of red significance lines for amp, coher
DOTSTYLE   = 'k--'; % line style to use for vetical dotted/dashed lines
LINESTYLE = '-';    % solid line
LABELFONT = 10;     % font sizes for axis labels, tick labels
TICKFONT  = 10;

PLOT_HEIGHT = 0.2;  % fraction of y dim taken up by each time series axes
YGAP = 0.03;        % fraction gap between time axes
YEXPAND = 1.3;      % expansion factor for y-axis about erp, amp data limits

DEFAULT_SDEV  = 1/7; % smooth trials with this window size by default if Gaussian window
DEFAULT_AVEWIDTH  = 1; % smooth trials with this window size by default
DEFAULT_DECFACTOR = 1; % decimate by this factor by default
DEFAULT_CYCLES    = 3; % use this many cycles in amp,coher computation window
cycles = DEFAULT_CYCLES;
DEFAULT_CBAR      = NO;% do not plot color bar by default
DEFAULT_PHARGS = [0 25 8 13]; % Default arguments for phase sorting
DEFAULT_ALPHA     = 0.01;
alpha     = 0;      % default alpha level for coherence significance

MIN_ERPALPHA = 0.001; % significance bounds for ERP
MAX_ERPALPHA = 0.1;

NoShowVar = NO;     % show sortvar by default
Nosort    = NO;     % sort on sortvar by default
Caxflag   = NO;     % use default caxis by default

timestretchflag = NO; % Added -JH

mvavg_type='boxcar'; % use the original rectangular moving average -DG
erp_grid = NO; % add y-tick grids to ERP plot -DG
cbar_title = []; % title to add above ERPimage color bar (e.g., '\muV') -DG
img_ylab = 'Trials'; % make the ERPimage y-axis in units of the sorting variable -DG
img_ytick_lab = []; % the values at which tick marks will appear on the trial axis of the ERPimage (the y-axis by default).
%Note, this is in units of the sorting variable if img_ylab~='Trials', otherwise it is in units of trials -DG
baseline   = []; %time window of each trial whose mean amp will be used to baseline the trial -DG
baselinedb = []; %time window of each trial whose mean power will be used to baseline the trial
flt=[]; %frequency domain filter parameters -DG
sortvar_limits=[]; %plotting limits for sorting variable/trials; limits only affect visualization, not smoothing -DG
replace_ties = NO; %if YES, trials with the exact same value of a sorting variable will be replaced by their average -DG
erp_vltg_ticks=[]; %if non-empty, these are the voltage axis ticks for the ERP plot

Caxis     = [];
caxfraction = [];
Coherflag = NO;     % don't compute or show amp,coher by default
Cohsigflag= NO;     % default: do not compute coherence significance
Allampsflag=NO;     % don't image the amplitudes by default
Allcohersflag=NO;   % don't image the coherence amplitudes by default
Topoflag  = NO;     % don't plot a topoplot in upper left
Specflag  = NO;     % don't plot a spectrum in upper right
SpecAxisflag = NO;  % don't change the default spectrum axis type from default
SpecAxis = 'log';   % default log frequency spectrum axis (if Specflag)
Erpflag   = NO;     % don't show erp average by default
Erpstdflag= NO;
Erpalphaflag= NO;
Alignflag = NO;     % don't align data to sortvar by default
Colorbar  = NO;     % if YES, plot a colorbar to right of erp image
Limitflag = NO;     % plot whole times range by default
Phaseflag = NO;     % don't sort by phase
Ampflag   = NO;     % don't sort by amplitude
Sortwinflag = NO;   % sort by amplitude over a window
Valflag   = NO;     % don't sort by value
Srateflag = NO;     % srate not given
Vertflag  = NO;
Horzflag  = NO;
titleflag = NO;
NoShowflag  = NO;
Renormflag = NO;
Showwin = NO;
yerplabel = 'ERP';
yerplabelflag = NO;
verttimes = [];
horzepochs = [];
NoTimeflag= NO;     % by default DO print "Time (ms)" below bottom axis
Signifflag= NO;     % compute significance instead of receiving it
Auxvarflag= NO;
plotmodeflag= NO;
plotmode = 'normal';
Cycleflag = NO;
signifs   = NaN;
coherfreq = nan;    % amp/coher-calculating frequency
freq = 0;           % phase-sorting frequency
srate = DEFAULT_SRATE; % from icadefs.m
aligntime = nan;
timelimits= nan;
topomap   = [];     % topo map vector
lospecHz  = [];     % spec lo frequency
topphase = 180;     % default top phase for 'phase' option
renorm    = 'no';
NoShow    = 'no';
Rmerp     = 'no';
percentiles = [];
percentileflag = NO;
erp_ptiles = 1;

minerp = NaN; % default limits
maxerp = NaN;
minamp = NaN;
maxamp = NaN;
mincoh = NaN;
maxcoh = NaN;
baseamp =NaN;
allamps = []; % default return

ax1    = NaN; % default axes handles
axcb   = NaN;
ax2    = NaN;
ax3    = NaN;
ax4    = NaN;


timeStretchRef = [];
timeStretchMarks = [];
tsurdata = []; % time-stretched data before smoothing (for timestretched-erp computation)

% If time-stretching is off, this variable remains empty
%
%% %%%%%%%%%%%%%%%%% Test, fill in commandline args %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin < 1
    help erpimage
    return
end

data = squeeze(data);
if nargin < 3 | isempty(times)
    if size(data,1)==1 || size(data,2)==1
        fprintf('\nerpimage(): either input a times vector or make data size = (frames,trials).\n')
        return
    end
    times = 1:size(data,1);
    NoTimesPassed= 1;
end

if nargin < 2 | isempty(sortvar)
    sortvar = 1:size(data,2);
    NoShowVar = 1; % don't plot the dummy sortvar
end

framestot = size(data,1)*size(data,2);
ntrials = length(sortvar);
if ntrials < 2
    help erpimage
    fprintf('\nerpimage(): too few trials.\n');
    return
end

frames = floor(framestot/ntrials);
if frames*ntrials ~= framestot
    help erpimage
    fprintf(...
        '\nerpimage(); length of sortvar doesn''t divide number of data elements??\n')
    return
end

if nargin < 6
    decfactor = 0;
end
if nargin < 5
    avewidth = DEFAULT_AVEWIDTH;
end
if nargin<4
    titl = ''; % default no title
end
if nargin<3
    times = NO;
end
if (length(times) == 1) | (times == NO),  % make default times
    times = 0:frames-1;
    srate = 1000*(length(times)-1)/(times(length(times))-times(1));
    fprintf('Using sampling rate %g Hz.\n',srate);
elseif length(times) == 3
    mintime = times(1);
    frames = times(2);
    srate = times(3);
    times = mintime:1000/srate:mintime+(frames-1)*1000/srate;
    fprintf('Using sampling rate %g Hz.\n',srate);
else
    % Note: might use default srate read from icadefs here...
    srate = 1000*(length(times)-1)/(times(end)-times(1));
end
if length(times) ~= frames
    fprintf(...
        '\nerpimage(): length(data)(%d) ~= length(sortvar)(%d) * length(times)(%d).\n\n',...
        framestot,              length(sortvar),   length(times));
    return
end

if decfactor == 0,
    decfactor = DEFAULT_DECFACTOR;
elseif decfactor > ntrials
    fprintf('Setting variable decfactor to max %d.\n',ntrials)
    decfactor = ntrials;
end
%
%% %%%%%%%%%%%%%%% Collect optional args %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin > 6
    flagargs = [];
    
    a = 6;
    while a < nargin % for each remaining Arg
        a = a + 1;
        
        Arg = eval(['arg' int2str(a-6)]);
        if Caxflag == YES
            if size(Arg,1) ~= 1 || size(Arg,2) > 2
                help erpimage
                fprintf('\nerpimage(): caxis arg must be a scalar or (1,2) vector.\n');
                return
            end
            if size(Arg,2) == 2
                Caxis = Arg;
            else
                caxfraction = Arg;
            end
            Caxflag = NO;
            
        elseif timestretchflag == YES % Added -JH
            timeStretchMarks = Arg{1};
            timeStretchMarks = round(1+(timeStretchMarks-times(1))*srate/1000); % convert from ms to frames -sm
            [smc smr] = find(diff(timeStretchMarks') < 0);
            if ~isempty(smr)
                fprintf('\nerpimage(): Timewarp event latencies not in ascending order in trial %d.\n',smr)
                return
            end
            
            timeStretchMarks = [ ...
                repmat(1, [size(timeStretchMarks,1), 1]), ...% Epoch begins
                timeStretchMarks, ...
                repmat(length(times), [size(timeStretchMarks,1), 1])]; % Epoch ends
            if length(Arg) < 2 || isempty(Arg{2})
                timeStretchRef = median(timeStretchMarks);
            else
                timeStretchRef = Arg{2};
                timeStretchRef = round(1+(timeStretchRef-times(1))*srate/1000); % convert from ms to frames -sm
                timeStretchRef = [1 timeStretchRef length(times)];  % add epoch beginning, end
            end
            if length(Arg) < 3 || isempty(Arg{3})
                timeStretchColors = {};
            else
                timeStretchColors = Arg{3};
            end
            fprintf('The %d events specified in each trial will be time warped to latencies:',length(timeStretchRef)-2);
            fprintf(' %.0f', times(1)+1000*(timeStretchRef(2:end-1)-1)/srate); % converted from frames to ms -sm
            fprintf(' ms\n');
            timestretchflag = NO;
        elseif Coherflag == YES
            if length(Arg) > 3 || length(Arg) < 1
                help erpimage
                fprintf('\nerpimage(): coher arg must be size <= 3.\n');
                return
            end
            coherfreq = Arg(1);
            if size(Arg,2) == 1
                coherfreq = Arg(1);
            else
                coherfreq = Arg(1:2);
            end;
            if size(Arg,2) == 3
                Cohsigflag = YES;
                alpha  = Arg(3);
                if alpha < 0 || alpha > 0.1
                    fprintf('\nerpimage(): alpha value %g out of bounds.\n',alpha);
                    return
                end
            end
            Coherflag = NO;
            Erpflag = YES;  % plot amp, coher below erp time series
        elseif Topoflag == YES;
            if length(Arg) < 2
                help erpimage
                fprintf('\nerpimage(): topo arg must be a list of length 2 or 3.\n');
                return
            end
            topomap = Arg{1};
            eloc_file = Arg{2};
            if length(Arg) > 2, eloc_info = Arg{3};
            else                eloc_info = [];
            end;
            Topoflag = NO;
        elseif Specflag == YES;
            if length(Arg) ~= 2 
                error('\nerpimage(): ''spec'' flag argument must be a numeric array of length 2.\n');
                return
            end
            lospecHz = Arg(1);
            hispecHz = Arg(2);
            Specflag = NO;
        elseif SpecAxisflag == YES;
            SpecAxis = Arg;
            if ~strcmpi(SpecAxis,'lin') && ~strcmpi(SpecAxis,'log')
              error('\nerpimage(): spectrum axis type must be ''lin'' or ''log''.\n');
              return
            end
            if strcmpi(SpecAxis,'lin')
               SpecAxis = 'linear';   % convert to MATLAB Xscale keyword
            end
            SpecAxisflag = NO;
        elseif Renormflag == YES
            renorm = Arg;
            Renormflag = NO;
        elseif NoShowflag == YES
            NoShow = Arg;
            if strcmpi(NoShow, 'off'), NoShow = 'no'; end;
            NoShowflag = NO;
        elseif Alignflag == YES
            aligntime = Arg;
            Alignflag = NO;
        elseif percentileflag == YES
            percentiles = Arg;
            percentileflag = NO;
        elseif Limitflag == YES
            %  [lotime hitime loerp hierp loamp hiamp locoher hicoher]
            if size(Arg,1) ~= 1 || size(Arg,2) < 2 || size(Arg,2) > 9
                help erpimage
                fprintf('\nerpimage(): limits arg must be a vector sized (1,2<->9).\n');
                return
            end
            if  ~isnan(Arg(1)) & (Arg(2) <= Arg(1))
                help erpimage
                fprintf('\nerpimage(): time limits out of order or out of range.\n');
                return
            end
            if Arg(1) < min(times)
                Arg(1) = min(times);
                fprintf('Adjusting mintime limit to first data value %g\n',min(times));
            end
            if Arg(2) > max(times)
                Arg(2) = max(times);
                fprintf('Adjusting maxtime limit to last data value %g\n',max(times));
            end
            timelimits = Arg(1:2);
            if length(Arg)> 2
                minerp = Arg(3);
            end
            if length(Arg)> 3
                maxerp = Arg(4);
            end
            if ~isnan(maxerp) & maxerp <= minerp
                help erpimage
                fprintf('\nerpimage(): erp limits args out of order.\n');
                return
            end
            if length(Arg)> 4
                minamp = Arg(5);
            end
            if length(Arg)> 5
                maxamp = Arg(6);
            end
            if maxamp <= minamp
                help erpimage
                fprintf('\nerpimage(): amp limits args out of order.\n');
                return
            end
            if length(Arg)> 6
                mincoh = Arg(7);
            end
            if length(Arg)> 7
                maxcoh = Arg(8);
            end
            if maxcoh <= mincoh
                help erpimage
                fprintf('\nerpimage(): coh limits args out of order.\n');
                return
            end
            if length(Arg)>8
                baseamp = Arg(9);    % for 'allamps'
            end
            Limitflag = NO;
            
        elseif Srateflag == YES
            srate = Arg(1);
            Srateflag = NO;
        elseif Cycleflag == YES
            cycles = Arg;
            Cycleflag = NO;
        elseif Auxvarflag == YES;
            if isa(Arg,'cell')==YES && length(Arg)==2
                auxvar = Arg{1};
                auxcolors = Arg{2};
            elseif isa(Arg,'cell')==YES
                fprintf('\nerpimage(): auxvars argument must be a matrix or length-2 cell array.\n');
                return
            else
                auxvar = Arg; % no auxcolors specified
            end
            [xr,xc] = size(auxvar);
            lns = length(sortvar);
            if xr ~= lns && xc ~= lns
                error('\nerpimage(): auxvar columns different from the number of epochs in data');
            elseif xr == lns && xc ~= lns
                auxvar = auxvar';   % exchange rows/cols
            end
            Auxvarflag = NO;
        elseif Vertflag == YES
            verttimes = Arg;
            Vertflag = NO;
        elseif Horzflag == YES
            horzepochs = Arg;
            Horzflag = NO;
        elseif yerplabelflag == YES
            yerplabel = Arg;
            yerplabelflag = NO;
        elseif Signifflag == YES
            signifs = Arg; % [low_amp hi_amp coher]
            if length(signifs) ~= 3
                fprintf('\nerpimage(): signif arg [%g] must have 3 values\n',Arg);
                return
            end
            Signifflag = NO;
        elseif Allcohersflag == YES
            data2=Arg;
            if size(data2) ~= size(data)
                fprintf('\nerpimage(): allcohers data matrix must be the same size as data.\n');
                return
            end
            Allcohersflag = NO;
        elseif Phaseflag == YES
            n = length(Arg);
            if n > 5
                error('\nerpimage(): Too many arguments for keyword ''phasesort''');
            end
            phargs = Arg;
            
            if phargs(3) < 0
                error('\nerpimage(): Invalid negative frequency argument for keyword ''phasesort''');
            end
            if n>=4
                if phargs(4) < 0
                    error('\nerpimage(): Invalid negative argument for keyword ''phasesort''');
                end
            end
            if min(phargs(1)) < times(1) | max(phargs(1)) > times(end)
                error('\nerpimage(): time for phase sorting filter out of bound.');
            end
            if phargs(2) >= 100 | phargs(2) < -100
                error('\nerpimage(): %-argument for keyword ''phasesort'' must be (-100;100)');
            end
            if length(phargs) >= 4 & phargs(3) > phargs(4)
                error('\nerpimage(): Phase sorting frequency range must be increasing.');
            end
            if length(phargs) == 5
                topphase = phargs(5);
            end
            Phaseflag = NO;
        elseif Sortwinflag == YES % 'ampsort' mean amplitude over a time window
            n = length(Arg);
            sortwinarg = Arg;
            if n > 2
                error('\nerpimage(): Too many arguments for keyword ''sortwin''');
            end
            if min(sortwinarg(1)) < times(1) | max(sortwinarg(1)) > times(end)
                error('\nerpimage(): start time for value sorting out of bounds.');
            end
            if n > 1
                if min(sortwinarg(2)) < times(1) | max(sortwinarg(2)) > times(end)
                    error('\nerpimage(): end time for value sorting out of bounds.');
                end
            end
            if n > 1 & sortwinarg(1) > sortwinarg(2)
                error('\nerpimage(): Value sorting time range must be increasing.');
            end
            Sortwinflag = NO;
        elseif Ampflag == YES % 'ampsort',[center_time,prcnt_reject,minfreq,maxfreq]
            n = length(Arg);
            if n > 4
                error('\nerpimage(): Too many arguments for keyword ''ampsort''');
            end
            ampargs = Arg;
            
            % if ampargs(3) < 0
            %    error('\nerpimage(): Invalid negative argument for keyword ''ampsort''');
            % end
            if n>=4
                if ampargs(4) < 0
                    error('\nerpimage(): Invalid negative argument for keyword ''ampsort''');
                end
            end
            
            if ~isinf(ampargs(1))
                if min(ampargs(1)) < times(1) | max(ampargs(1)) > times(end)
                    error('\nerpimage(): time for amplitude sorting filter out of bounds.');
                end
            end
            
            if ampargs(2) >= 100 | ampargs(2) < -100
                error('\nerpimage(): percentile argument for keyword ''ampsort'' must be (-100;100)');
            end
            
            if length(ampargs) == 4 & abs(ampargs(3)) > abs(ampargs(4))
                error('\nerpimage(): Amplitude sorting frequency range must be increasing.');
            end
            Ampflag = NO;
            
        elseif Valflag == YES % sort by potential value in a given window
            % Usage: 'valsort',[mintime,maxtime,direction]
            n = length(Arg);
            if n > 3
                error('\nerpimage(): Too many arguments for keyword ''valsort''');
            end
            valargs = Arg;
            
            if min(valargs(1)) < times(1) | max(valargs(1)) > times(end)
                error('\nerpimage(): start time for value sorting out of bounds.');
            end
            if n > 1
                if min(valargs(2)) < times(1) | max(valargs(2)) > times(end)
                    error('\nerpimage(): end time for value sorting out of bounds.');
                end
            end
            if n > 1 & valargs(1) > valargs(2)
                error('\nerpimage(): Value sorting time range must be increasing.');
            end
            if n==3 & (~isnumeric(valargs(3)) | valargs(3)==0)
                error('\nerpimage(): Value sorting direction must be +1 or -1.');
            end
            Valflag = NO;
        elseif plotmodeflag == YES
            plotmode = Arg; plotmodeflag = NO;
        elseif titleflag == YES
            titl = Arg; titleflag = NO;
        elseif Erpalphaflag == YES
            erpalpha = Arg(1);
            if erpalpha < MIN_ERPALPHA | erpalpha > MAX_ERPALPHA
                fprintf('\nerpimage(): erpalpha value is out of bounds [%g, %g]\n',...
                    MIN_ERPALPHA,MAX_ERPALPHA);
                return
            end
            Erpalphaflag = NO;
            % -----------------------------------------------------------------------
            % -----------------------------------------------------------------------
            % -----------------------------------------------------------------------
        elseif strcmpi(Arg,'avg_type')
            if a < nargin,
                a=a+1;
                Arg = eval(['arg' int2str(a-6)]);
                if strcmpi(Arg, 'Gaussian'), mvavg_type='gaussian';
                elseif strcmpi(Arg, 'Boxcar'), mvavg_type='boxcar';
                else error('\nerpimage(): Invalid value for optional argument ''avg_type''.');
                end;
            else
                error('\nerpimage(): Optional argument ''avg_type'' needs to be assigned a value.');
            end
        elseif strcmp(Arg,'nosort')
            Nosort = YES;
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     Nosort = YES; a = a+1;
                elseif strcmpi(Arg, 'off') Nosort = NO;  a = a+1;
                end;
            end;
        elseif strcmp(Arg,'showwin')
            Showwin = YES;
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     Showwin = YES; a = a+1;
                elseif strcmpi(Arg, 'off') Showwin = NO;  a = a+1;
                end;
            end;
        elseif strcmp(Arg,'noplot') % elseif strcmp(Arg,'NoShow') % by Luca & Ramon
            NoShowVar = YES;
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     NoShowVar = YES; a = a+1;
                elseif strcmpi(Arg, 'off'), NoShowVar = NO;  a = a+1;
                end;
            end;
        elseif strcmpi(Arg,'replace_ties')
            if a < nargin,
                a = a+1;
                temp = eval(['arg' int2str(a-6)]);
                if strcmpi(temp,'on'),
                    replace_ties = YES;
                elseif strcmpi(temp,'off') replace_ties = NO;
                else
                    error('\nerpimage(): Argument ''replace_ties'' needs to be followed by the string ''on'' or ''off''.');
                end
            else
                error('\nerpimage(): Argument ''replace_ties'' needs to be followed by the string ''on'' or ''off''.');
            end
        elseif strcmpi(Arg,'sortvar_limits')
            if a < nargin,
                a = a+1;
                sortvar_limits = eval(['arg' int2str(a-6)]);
                if ischar(sortvar_limits) || length(sortvar_limits)~=2
                    error('\nerpimage(): Argument ''sortvar_limits'' needs to be followed by a two element vector.');
                end
            else
                error('\nerpimage(): Argument ''sortvar_limits'' needs to be followed by a two element vector.');
            end
        elseif strcmpi(Arg,'erp')
            Erpflag = YES;
            erp_ptiles=1;
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     Erpflag = YES; erp_ptiles=1; a = a+1;
                elseif strcmpi(Arg, 'off') Erpflag = NO;  a = a+1;
                elseif strcmpi(Arg,'1') | (Arg==1) Erplag = YES; erp_ptiles=1; a=a+1;
                elseif strcmpi(Arg,'2') | (Arg==2) Erplag = YES; erp_ptiles=2; a=a+1;
                elseif strcmpi(Arg,'3') | (Arg==3) Erplag = YES; erp_ptiles=3; a=a+1;
                elseif strcmpi(Arg,'4') | (Arg==4) Erplag = YES; erp_ptiles=4; a=a+1;
                end;
            end;
        elseif strcmpi(Arg,'rmerp')
            Rmerp = 'yes';
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     Rmerp = 'yes'; a = a+1;
                elseif strcmpi(Arg, 'off') Rmerp = 'no';  a = a+1;
                end;
            end;
        elseif strcmp(Arg,'cbar') | strcmp(Arg,'colorbar')
            Colorbar = YES;
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     Colorbar = YES; a = a+1;
                elseif strcmpi(Arg, 'off') Colorbar = NO;  a = a+1;
                end;
            end;
        elseif (strcmp(Arg,'allamps') | strcmp(Arg,'plotamps'))
            Allampsflag = YES;
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     Allampsflag = YES; a = a+1;
                elseif strcmpi(Arg, 'off') Allampsflag = NO;  a = a+1;
                end;
            end;
        elseif strcmpi(Arg,'erpstd')
            Erpstdflag = YES;
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     Erpstdflag = YES; a = a+1;
                elseif strcmpi(Arg, 'off') Erpstdflag = NO;  a = a+1;
                end;
            end;
        elseif strcmp(Arg,'noxlabel') | strcmp(Arg,'noxlabels') | strcmp(Arg,'nox')
            NoTimeflag = YES;
            if a < nargin,
                Arg = eval(['arg' int2str(a+1-6)]);
                if strcmpi(Arg, 'on'),     NoTimeflag = YES; a = a+1;
                elseif strcmpi(Arg, 'off') NoTimeflag = NO;  a = a+1;
                end;
            end;
        elseif strcmp(Arg,'plotmode')
            plotmodeflag = YES;
        elseif strcmp(Arg,'sortvarpercent')
            percentileflag = YES;
        elseif strcmp(Arg,'renorm')
            Renormflag = YES;
        elseif strcmp(Arg,'NoShow')
            NoShowflag = YES;
        elseif strcmp(Arg,'caxis')
            Caxflag = YES;
        elseif strcmp(Arg,'title')
            titleflag = YES;
        elseif strcmp(Arg,'coher')
            Coherflag = YES;
        elseif strcmp(Arg,'timestretch') | strcmp(Arg,'timewarp') % Added -JH
            timestretchflag = YES;
        elseif strcmp(Arg,'allcohers')
            Allcohersflag = YES;
        elseif strcmp(Arg,'topo') | strcmp(Arg,'topoplot')
            Topoflag = YES;
        elseif strcmp(Arg,'spec') | strcmp(Arg,'spectrum')
            Specflag = YES;
        elseif strcmp(Arg,'Specaxis') || strcmp(Arg,'specaxis') || strcmp(Arg,'SpecAxis')
            SpecAxisflag = YES;
        elseif strcmpi(Arg,'erpalpha')
            Erpalphaflag = YES;
        elseif strcmp(Arg,'align')
            Alignflag = YES;
        elseif strcmp(Arg,'limits')
            Limitflag = YES;
        elseif (strcmp(Arg,'phase') | strcmp(Arg,'phasesort'))
            Phaseflag = YES;
        elseif strcmp(Arg,'ampsort')
            Ampflag = YES;
        elseif strcmp(Arg,'sortwin')
            Sortwinflag = YES;
        elseif strcmp(Arg,'valsort')
            Valflag = YES;
        elseif strcmp(Arg,'auxvar')
            Auxvarflag = YES;
        elseif strcmp(Arg,'cycles')
            Cycleflag = YES;
        elseif strcmpi(Arg,'yerplabel')
            yerplabelflag = YES;
        elseif strcmpi(Arg,'srate')
            Srateflag = YES;
        elseif strcmpi(Arg,'erp_grid')
            erp_grid = YES;
        elseif strcmpi(Arg,'baseline')
            if a < nargin,
                a = a+1;
                baseline = eval(['arg' int2str(a-6)]);
            else
                error('\nerpimage(): Argument ''baseline'' needs to be followed by a two element vector.');
            end
        elseif strcmpi(Arg,'baselinedb')
            if a < nargin,
                a = a+1;
                baselinedb = eval(['arg' int2str(a-6)]);
                if length(baselinedb) > 2, error('''baselinedb'' need to be a 2 argument vector'); end;
            else
                error('\nerpimage(): Argument ''baselinedb'' needs to be followed by a two element vector.');
            end
        elseif strcmpi(Arg,'filt')
            if a < nargin,
                a = a+1;
                flt = eval(['arg' int2str(a-6)]);
            else
                error('\nerpimage(): Argument ''filt'' needs to be followed by a two element vector.');
            end
        elseif strcmpi(Arg,'erp_vltg_ticks')
            if a < nargin,
                a = a+1;
                erp_vltg_ticks=eval(['arg' int2str(a-6)]);
            else
                error('\nerpimage(): Argument ''erp_vltg_ticks'' needs to be followed by a vector.');
            end
        elseif strcmpi(Arg,'img_trialax_label')
            if a < nargin,
                a = a+1;
                img_ylab = eval(['arg' int2str(a-6)]);
            else
                error('\nerpimage(): Argument ''img_trialax_label'' needs to be followed by a string.');
            end;
        elseif strcmpi(Arg,'img_trialax_ticks')
            if a < nargin,
                a = a+1;
                img_ytick_lab = eval(['arg' int2str(a-6)]);
            else
                error('\nerpimage(): Argument ''img_trialax_ticks'' needs to be followed by a vector of values at which tick marks will appear.');
            end;
        elseif strcmpi(Arg,'cbar_title')
            if a < nargin,
                a = a+1;
                cbar_title = eval(['arg' int2str(a-6)]);
            else
                error('\nerpimage(): Argument ''cbar_title'' needs to be followed by a string.');
            end;
        elseif strcmp(Arg,'vert') ||  strcmp(Arg,'verttimes')
            Vertflag = YES;
        elseif strcmp(Arg,'horz') ||  strcmp(Arg,'horiz') || strcmp(Arg,'horizontal')
            Horzflag = YES;
        elseif strcmp(Arg,'signif') || strcmp(Arg,'signifs') || strcmp(Arg,'sig') || strcmp(Arg,'sigs')
            Signifflag = YES;
        else
            help erpimage
            if isstr(Arg)
                fprintf('\nerpimage(): unknown arg %s\n',Arg);
            else
                fprintf('\nerpimage(): unknown arg %d, size(%d,%d)\n',a,size(Arg,1),size(Arg,2));
            end
            return
        end
    end % Arg
end

if exist('img_ylab','var') || exist('img_ytick_lab','var'),
    oops=0;
    if exist('phargs','var'),
        fprintf('********* Warning *********\n');
        fprintf('Options ''img_ylab'' and ''img_ytick_lab'' have no effect when sorting by phase.\n');
        oops=0;
    elseif exist('valargs','var'),
        fprintf('********* Warning *********\n');
        fprintf('Options ''img_ylab'' and ''img_ytick_lab'' have no effect when sorting by EEG voltage.\n');
        oops=0;
    elseif exist('ampargs','var'),
        fprintf('********* Warning *********\n');
        fprintf('Options ''img_ylab'' and ''img_ytick_lab'' have no effect when sorting by frequency amplitude.\n');
        oops=0;
    end
    if oops
        img_ylab=[];
        img_ytick_lab=[];
    end
end

if   Caxflag == YES ...
        |Coherflag == YES ...
        |Alignflag == YES ...
        |Limitflag == YES
    help erpimage
    fprintf('\nerpimage(): missing option arg.\n')
    return
end
if (Allampsflag | exist('data2')) & ( any(isnan(coherfreq)) | ~Cohsigflag )
    fprintf('\nerpimage(): allamps and allcohers flags require coher freq, srate, and cohsig.\n');
    return
end
if Allampsflag & exist('data2')
    fprintf('\nerpimage(): cannot image both allamps and allcohers.\n');
    return
end
if ~exist('srate') | srate <= 0
    fprintf('\nerpimage(): Data srate must be specified and > 0.\n');
    return
end
if ~isempty(auxvar)
    % whos auxvar
    if size(auxvar,1) ~= ntrials & size(auxvar,2) ~= ntrials
        fprintf('\nerpimage(): auxvar size should be (N,ntrials), e.g., (N,%d)\n',...
            ntrials);
        return
    end
    if size(auxvar,1) == ntrials & size(auxvar,2) ~= ntrials  % make (N,frames)
        auxvar = auxvar';
    end
    if size(auxvar,2) ~= ntrials
        fprintf('\nerpimage(): auxvar size should be (N,ntrials), e.g., (N,%d)\n',...
            ntrials);
        return
    end
    if exist('auxcolors')==YES % if specified
        if isa(auxcolors,'cell')==NO % if auxcolors is not a cell array
            fprintf(...
                '\nerpimage(): auxcolors argument to auxvar flag must be a cell array.\n');
            return
        end
    end
elseif exist('timeStretchRef') & ~isempty(timeStretchRef)
    if ~isnan(aligntime)
        fprintf(['\nerpimage(): options "align" and ' ...
            '"timewarp" are not compatiable.\n']);
        return;
    end
    
    if ~isempty(timeStretchColors)
        if length(timeStretchColors) < length(timeStretchRef)
            nColors = length(timeStretchColors);
            for k=nColors+1:length(timeStretchRef)-2
                timeStretchColors = { timeStretchColors{:} timeStretchColors{1+rem(k-1,nColors)}};
            end
        end
        timeStretchColors = {'' timeStretchColors{:} ''};
    else
        timeStretchColors = { 'k--'};
        for k=2:length(timeStretchRef)
            timeStretchColors = { timeStretchColors{:} 'k--'};
        end
    end
    
    
    auxvarInd = 1-strcmp('',timeStretchColors); % indicate which lines to draw
    newauxvars = ((timeStretchRef(find(auxvarInd))-1)/srate+times(1)/1000) * 1000; % convert back to ms
    fprintf('Overwriting vert with auxvar\n');
    verttimes = [newauxvars'];
    verttimesColors = {timeStretchColors{find(auxvarInd)}};
    newauxvars = repmat(newauxvars, [1 ntrials]);
    
    if isempty(auxvar) % Initialize auxvar & auxcolors
        %  auxvar = newauxvars;
        auxcolors = {timeStretchColors{find(auxvarInd)}};
    else % Append auxvar & auxcolors
        if ~exist('auxcolors')
            % Fill with default color (k-- for now)
            auxcolors = {};
            for j=1:size(auxvar,1)
                auxcolors{end+1} = 'k--';
            end
        end
        for j=find(auxvarInd)
            auxcolors{end+1} = timeStretchColors{j};
        end
        auxvar = [auxvar; newauxvars];
    end
end

% No need to turn off ERP anymore; we now time-stretch potentials
% separately (see tsurdata below) so the (warped) ERP is actually meaningful

if exist('phargs')
    if phargs(3) > srate/2
        fprintf(...
            '\nerpimage(): Phase-sorting frequency (%g Hz) must be less than Nyquist rate (%g Hz).',...
            phargs(3),srate/2);
    end
    
    if frames < cycles*srate/phargs(3)
        fprintf('\nerpimage(): phase-sorting freq. (%g) too low: epoch length < %d cycles.\n',...
            phargs(3),cycles);
        return
    end
    if length(phargs)==4 & phargs(4) > srate/2
        phargs(4) = srate/2;
    end
    if length(phargs)==5 & (phargs(5)>180 | phargs(5) < -180)
        fprintf('\nerpimage(): coher topphase (%g) out of range.\n',topphase);
        return
    end
end
if exist('ampargs')
    if abs(ampargs(3)) > srate/2
        fprintf(...
            '\nerpimage(): amplitude-sorting frequency (%g Hz) must be less than Nyquist rate (%g Hz).',...
            abs(ampargs(3)),srate/2);
    end
    
    if frames < cycles*srate/abs(ampargs(3))
        fprintf('\nerpimage(): amplitude-sorting freq. (%g) too low: epoch length < %d cycles.\n',...
            abs(ampargs(3)),cycles);
        return
    end
    if length(ampargs)==4 & abs(ampargs(4)) > srate/2
        ampargs(4) = srate/2;
        fprintf('> Reducing max ''ampsort'' frequency to Nyquist rate (%g Hz)\n',srate/2)
    end
end
if ~any(isnan(coherfreq))
    if coherfreq(1) <= 0 | srate <= 0
        fprintf('\nerpimage(): coher frequency (%g) out of range.\n',coherfreq(1));
        return
    end
    if coherfreq(end) > srate/2 | srate <= 0
        fprintf('\nerpimage(): coher frequency (%g) out of range.\n',coherfreq(end));
        return
    end
    if frames < cycles*srate/coherfreq(1)
        fprintf('\nerpimage(): coher freq. (%g) too low:  epoch length < %d cycles.\n',...
            coherfreq(1),cycles);
        return
    end
end

if isnan(timelimits)
    timelimits = [min(times) max(times)];
end
if ~isstr(aligntime) & ~isnan(aligntime)
    if ~isinf(aligntime) ...
            & (aligntime < timelimits(1) | aligntime > timelimits(2))
        help erpimage
        fprintf('\nerpimage(): requested align time outside of time limits.\n');
        return
    end
end
%

%% %%%%%%%%%%%%%%  Replace nan's with 0s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
nans = find(isnan(data));
if length(nans)
    fprintf('Replaced %d nan in data with 0s.\n');
    data(nans) = 0;
end
%
%% %%%%%%%%%%%% Reshape data to (frames,ntrials) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if size(data,2) ~= ntrials
    if size(data,1)>1
        % fprintf('frames %d, ntrials %d length(data) %d\n',frames,ntrials,length(data));
        data=reshape(data,1,frames*ntrials);
    end
    data=reshape(data,frames,ntrials);
end
% fprintf('Plotting input data as %d epochs of %d frames sampled at %3.1f Hz.\n',...
%     ntrials,frames,srate);
%
%% %%%%%%%%%%%% Reshape data2 to (frames,ntrials) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if exist('data2') == 1
    if size(data2,2) ~= ntrials
        if size(data2,1)>1
            data2=reshape(data2,1,frames*ntrials);
        end
        data2=reshape(data2,frames,ntrials);
    end
end
%
%% %%%%%%%%%%%%% if sortvar=NaN, remove lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -ad
%
if any(isnan(sortvar))
    nanlocs = find(isnan(sortvar));
    fprintf('Removing %d trials with NaN sortvar values.\n', length(nanlocs));
    data(:,nanlocs) = [];
    sortvar(nanlocs) = [];
    if exist('data2') == 1
        data2(:,nanlocs) = [];
    end;
    if ~isempty(auxvar)
        auxvar(:,nanlocs) = [];
    end
    if ~isempty(verttimes)
        if size(verttimes,1) == ntrials
            verttimes(nanlocs,:) = [];
        end;
    end;
    ntrials = size(data,2);
    if ntrials <= 1, close(gcf); error('\nerpimage(): Too few trials'); end;
end;

%% Create moving average window %%
if strcmpi(mvavg_type,'Gaussian'),
    %construct Gaussian window to weight trials
    if avewidth == 0,
        avewidth = DEFAULT_SDEV;
    elseif avewidth < 1,
        help erpimage
        fprintf('\nerpimage(): Variable avewidth cannot be < 1.\n')
        fprintf('\nerpimage(): avewidth needs to be a positive integer.\n')
        return
    end
    wt_wind=exp(-0.5*([-3*avewidth:3*avewidth]/avewidth).^2)';
    wt_wind=wt_wind/sum(wt_wind); %normalize to unit sum
    avewidth=length(wt_wind);
    
    if avewidth > ntrials
        avewidth=floor((ntrials-1)/6);
        if avewidth==0,
            avewidth=DEFAULT_SDEV; %should be a window with one time point (smallest possible)
        end
        wt_wind=exp(-0.5*([-3*avewidth:3*avewidth]/avewidth).^2)';
        wt_wind=wt_wind/sum(wt_wind);
        fprintf('avewidth is too big for this number of trials.\n');
        fprintf('Changing avewidth to maximum possible size: %d\n',avewidth);
        avewidth=length(wt_wind);
    end
else
    %construct rectangular "boxcar" window to equally weight trials within
    %window
    if avewidth == 0,
        avewidth = DEFAULT_AVEWIDTH;
    elseif avewidth < 1
        help erpimage
        fprintf('\nerpimage(): Variable avewidth cannot be < 1.\n')
        return
    elseif avewidth > ntrials
        fprintf('Setting variable avewidth to max %d.\n',ntrials)
        avewidth = ntrials;
    end
    wt_wind=ones(1,avewidth)/avewidth;
end


%% Filter data with Butterworth filter (if requested) %%%%%%%%%%%%
%
if ~isempty(flt)
    %error check
    if length(flt)~=2,
        error('\nerpimage(): ''filt'' parameter argument should be a two element vector.');
    elseif max(flt)>(srate/2),
        error('\nerpimage(): ''filt'' parameters need to be less than or equal to sampling rate/2 (i.e., %f).',srate/2);
    elseif (flt(2)==(srate/2)) && (flt(1)==0),
        error('\nerpimage(): If second element of ''filt'' parameter is srate/2, then the first element must be greater than 0.');
    elseif abs(flt(2))<=abs(flt(1)),
        error('\nerpimage(): Second element of ''filt'' parameters must be greater than first in absolute value.');
    elseif (flt(1)<0) || (flt(2)<0),
        if (flt(1)>=0) || (flt(2)>=0),
            error('\nerpimage(): BOTH parameters of ''filt'' need to be greater than or equal to zero OR need to be negative.');
        end
        if min(flt)<=(-srate/2),
            error('\nerpimage(): ''filt'' parameters need to be greater than sampling rate/2 (i.e., -%f) when creating a stop band.',srate/2);
        end
    end
    
    fprintf('\nFiltering data with 3rd order Butterworth filter: ');
    if (flt(1)==0),
        %lowpass filter the data
        [B A]=butter(3,flt(2)*2/srate,'low');
        fprintf('lowpass at %.0f Hz\n',flt(2));
    elseif (flt(2)==(srate/2)),
        %highpass filter the data
        [B A]=butter(3,flt(1)*2/srate,'high');
        fprintf('highpass at %.0f Hz\n',flt(1));
    elseif (flt(1)<0)
        %bandstop filter the data
        flt=-flt;
        [B A]=butter(3,flt*2/srate,'stop');
        fprintf('stopband from %.0f to %.0f Hz\n',flt(1),flt(2));
    else
        %bandpass filter the data
        [B A]=butter(3,flt*2/srate);
        fprintf('bandpass from %.0f to %.0f Hz\n',flt(1),flt(2));
    end
    s=size(data);
    for trial=1:s(2),
        data(:,trial)=filtfilt(B,A,double(data(:,trial)));
    end
    if isempty(baseline)
        fprintf('Note, you might want to re-baseline the data using the erpimage ''baseline'' option.\n\n');
    end
end

%% Mean Baseline Each Trial (if requested) %%
if ~isempty(baseline),
    %check argument values for errors
    if baseline(2)<baseline(1),
        error('\nerpimage(): First element of ''baseline'' argument needs to be less than or equal to second argument.');
    elseif baseline(2)<times(1),
        error('\nerpimage(): Second element of ''baseline'' argument needs to be greater than or equal to epoch start time %.1f.',times(1));
    elseif baseline(1)>times(end),
        error('\nerpimage(): First element of ''baseline'' argument needs to be less than or equal to epoch end time %.1f.',times(end));
    end
    
    %convert msec into time points
    if baseline(1)<times(1),
        strt_pt=1;
    else
        strt_pt=find_crspnd_pt(baseline(1),times,1:length(times));
        strt_pt=ceil(strt_pt);
    end
    if baseline(end)>times(end),
        end_pt=length(times);
    else
        end_pt=find_crspnd_pt(baseline(2),times,1:length(times));
        end_pt=floor(end_pt);
    end
    fprintf('\nRemoving pre-stimulus mean baseline from %.1f to %.1f msec.\n\n',times(strt_pt),times(end_pt));
    bsln_mn=mean(data(strt_pt:end_pt,:),1);
    data=data-repmat(bsln_mn,length(times),1);
end


%
%% %%%%%%%%%%%%%%%%% Renormalize sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
switch lower(renorm)
    case 'yes',
%         disp('\nerpimage warning: *** sorting variable renormalized ***');
        sortvar = (sortvar-min(sortvar)) / (max(sortvar) - min(sortvar)) * ...
            0.5 * (max(times) - min(times)) + min(times) + 0.4*(max(times) - min(times));
    case 'no',;
    otherwise,
        if ~isempty(renorm)
            locx = findstr('x', lower(renorm));
            if length(locx) ~= 1, error('\nerpimage: unrecognized renormalizing formula'); end;
            eval( [ 'sortvar =' renorm(1:locx-1) 'sortvar' renorm(locx+1:end) ';'] );
        end;
end;
%
%% %%%%%%%%%%%%%%%%% Align data to sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if isstr(aligntime) | ~isnan(aligntime)
    if ~isstr(aligntime) & isinf(aligntime)
        aligntime= median(sortvar);
        %fprintf('Aligning data to median sortvar.\n');
        % Alternative below: trimmed median - ignore top/bottom 5%
        %   ssv = sort(sortvar); % ssv = 'sorted sortvar'
        %   aligntime= median(ssv(ceil(ntrials/20)):floor(19*ntrials/20));
    end
    
    if ~isstr(aligntime)
        %fprintf('Realigned sortvar plotted at %g ms.\n',aligntime);
        aligndata=zeros(frames,ntrials); % begin with matrix of zeros()
        shifts = zeros(1,ntrials);
        for t=1:ntrials, %%%%%%%%% foreach trial %%%%%%%%%
            shft = round((aligntime-sortvar(t))*srate/1000);
            shifts(t) = shft;
            if shft>0, % shift right
                if frames-shft > 0
                    aligndata(shft+1:frames,t)=data(1:frames-shft,t);
                else
                    %fprintf('No aligned data for epoch %d - shift (%d frames) too large.\n',t,shft);
                end
            elseif shft < 0 % shift left
                if frames+shft > 0
                    aligndata(1:frames+shft,t)=data(1-shft:frames,t);
                else
                    %fprintf('No aligned data for epoch %d - shift (%d frames) too large.\n',t,shft);
                end
            else % shft == 0
                aligndata(:,t) = data(:,t);
            end
        end % end trial
        if ~isempty(auxvar)
            auxvar = auxvar+shifts;
        end;
        %fprintf('Shifted epochs by %d to %d frames.\n',min(shifts),max(shifts));
        data = aligndata;                       % now data is aligned to sortvar
    else
        aligntime = str2num(aligntime);
        if isinf(aligntime),  aligntime= median(sortvar); end;
    end;
end

%
%% %%%%%%%%%%%%%%%%%%%% Remove the ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(Rmerp, 'yes')
    data = data - repmat(nan_mean(data')', [1 size(data,2)]);
end;

%
%% %%%%%%%%%%%%% Sort the data trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if exist('phargs') == 1 % if phase-sort the data trials
    if length(phargs) >= 4 & phargs(3) ~= phargs(4) % find max frequency
        % in specified band
        if exist('psd') == 2 % requires Signal Processing Toolbox
            %fprintf('Computing data spectrum using psd().\n');
            [pxx,freqs] = psd(data(:),max(1024, pow2(ceil(log2(frames)))),srate,frames,0);
        else % EEGLAB native work-around
            %fprintf('Computing data spectrum using spec().\n');
            [pxx,freqs] = spec(data(:),max(1024, pow2(ceil(log2(frames)))),srate,frames,0);
        end;
        % gf = gcf; % figure;plot(freqs,pxx); %xx=axis; %axis([phargs(3) phargs(4) xx(3) xx(4)]); %figure(gf);
        pxx = 10*log10(pxx);
        n = find(freqs >= phargs(3) & freqs <= phargs(4));
        if ~length(n)
            freq = (phargs(3)+phargs(4))/2;
        end
        [dummy maxx] = max(pxx(n));
        freq = freqs(n(maxx));
    else
        freq = phargs(3); % else use specified frequency
    end
%     fprintf('Sorting trials on phase at %.2g Hz.\n',freq);
    
    [amps, cohers, cohsig, ampsig, allamps, allphs] = ...
        phasecoher(data,length(times),srate,freq,cycles,0, ...
        [], [], timeStretchRef, timeStretchMarks);
    
    phwin = phargs(1);
    [dummy minx] = min(abs(times-phwin)); % closest time to requested
    winlen = floor(cycles*srate/freq);
    winloc = minx-linspace(floor(winlen/2), floor(-winlen/2), winlen+1);
    tmprange = find(winloc>0 & winloc<=frames);
    winloc = winloc(tmprange); % sorting window times
    winlocs = winloc;
    %
    %%%%%%%%%%%%%%%%%%%% Compute phsamp and phaseangles %%%%%%%%%%%%%%%%%%%%
    %
    % $$$     [phaseangles phsamp] = phasedet(data,frames,srate,winloc,freq);
    phaseangles = allphs(minx,:);
    phsamp = allamps(minx,:);
    % $$$     % hist(phaseangles,50);
    % $$$     % return
    % $$$
    % $$$     %
    % $$$     % Print facts on commandline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % $$$     %
    % $$$     if length(tmprange) ~=  winlen+1
    % $$$       filtersize = cycles * length(tmprange) / (winlen+1);
    % $$$       timecenter = median(winloc)/srate*1000+times(1); % center of window in ms
    % $$$       phaseangles = phaseangles + 2*pi*(timecenter-phargs(1))*freq;
    % $$$       fprintf('Sorting data epochs by phase at frequency %2.1f Hz: \n', freq);
    % $$$       fprintf(...
    % $$$           '    Data time limits reached -> now uses a %1.1f cycles (%1.0f ms) window centered at %1.0f ms\n', ...
    % $$$           filtersize, 1000/freq*filtersize, timecenter);
    % $$$       fprintf(...
    % $$$           '    Filter length is %d; Phase has been linearly interpolated to latency at %1.0f ms.\n', ...
    % $$$           length(winloc), phargs(1));
    % $$$     else
    % $$$       fprintf(...
    % $$$           'Sorting data epochs by phase at %2.1f Hz in a %1.1f-cycle (%1.0f ms) window centered at %1.0f ms.\n',...
    % $$$           freq,cycles,1000/freq*cycles,times(minx));
    % $$$       fprintf('Phase is computed using a wavelet of %d frames.\n',length(winloc));
    % $$$     end;
    %
    % Reject small (or large) phsamp trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    phargs(2) = phargs(2)/100; % convert rejection rate from % to fraction
    amprej = phargs(2);
    [tmp ampsortidx] = sort(phsamp); % sort amplitudes
    if amprej>=0
        ampsortidx = ampsortidx(ceil(amprej*length(ampsortidx))+1:end); % if amprej==0, select all trials
%         %fprintf('Retaining %d epochs (%g percent) with largest power at the analysis frequency,\n',...
%             length(ampsortidx),100*(1-amprej));
    else % amprej < 0
        amprej = 1+amprej; % subtract from end
        ampsortidx = ampsortidx(1:floor(amprej*length(ampsortidx)));
%         %fprintf('Retaining %d epochs (%g percent) with smallest power at the analysis frequency,\n',...
%             length(ampsortidx),amprej*100);
    end
    %
    % Remove low|high-amplitude trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    data = data(:,ampsortidx); % amp-sort the data, removing rejected-amp trials
    phsamp = phsamp(ampsortidx);           % amp-sort the amps
    phaseangles = phaseangles(ampsortidx); % amp-sort the phaseangles
    sortvar = sortvar(ampsortidx);         % amp-sort the trial indices
    ntrials = length(ampsortidx);          % number of trials retained
    if ~isempty(auxvar)
        auxvar = auxvar(:,ampsortidx);
    end
    if ~isempty(timeStretchMarks)
        timeStretchMarks =  timeStretchMarks(:,ampsortidx);
    end
    
    %
    % Sort remaining data by phase angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    phaseangles = -phaseangles;
    topphase = (topphase/360)*2*pi; % convert from degrees to radians
    ip = find(phaseangles>topphase);
    phaseangles(ip) = phaseangles(ip)-2*pi; % rotate so topphase at top of plot
    
    [phaseangles sortidx] = sort(phaseangles); % sort trials on (rotated) phase
    data    =  data(:,sortidx);                % sort data by phase
    phsamp  =  phsamp(sortidx);                % sort amps by phase
    sortvar = sortvar(sortidx);                % sort input sortvar by phase
    if ~isempty(auxvar)
        auxvar = auxvar(:,sortidx);
    end
    if ~isempty(timeStretchMarks)
        timeStretchMarks =  timeStretchMarks(:,sortidx);
    end
    phaseangles = -phaseangles; % Note: phsangles now descend from pi
    % TEST auxvar = 360 + (1000/256)*(256/5)*phaseangles/(2*pi); % plot phase+360 in ms for test
    
    %fprintf('Size of data = [%d,%d]\n',size(data,1),size(data,2));
    sortidx = ampsortidx(sortidx); % return original trial indices in final sorted order
    %
    % %%%%%%%%%%%%%%% Sort data by amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
elseif exist('ampargs') == 1 % if amplitude-sort
    if length(ampargs) == 4 % find max frequency in specified band
        if exist('psd') == 2
            %fprintf('Computing data spectrum using psd().\n');
            [pxx,freqs] = psd(data(:),max(1024, pow2(ceil(log2(frames)))),srate,frames,0);
        else
            %fprintf('Computing data spectrum using spec().\n');
            [pxx,freqs] = spec(data(:),max(1024, pow2(ceil(log2(frames)))),srate,frames,0);
        end;
        pxx = 10*log10(pxx);
        if ampargs(3) == ampargs(4)
            [freq n] = min(abs(freqs - ampargs(3)));
        else
            n = find(freqs >= abs(ampargs(3)) & freqs <= abs(ampargs(4)));
        end;
        if ~length(n)
            freq = mean([abs(ampargs(3)),abs(ampargs(4))]);
        end
        if ampargs(3)>=0
            [dummy maxx] = max(pxx(n));
            freq = freqs(n(maxx)); % use the highest-power frequency
        else
            freq = freqs(n);  % use all frequencies in the specified range
        end
    else
        freq = abs(ampargs(3)); % else use specified frequency
    end
    if length(freq) == 1
        %fprintf('Sorting data epochs by amplitude at frequency %2.1f Hz \n', freq);
    else
        %fprintf('Sorting data epochs by amplitude at %d frequencies (%2.1f Hz to %.1f Hz) \n',...
%             length(freq),freq(1),freq(end));
    end
    SPECWININCR = 10;	% make spectral sorting time windows increment by 10 ms
    if isinf(ampargs(1))
        ampwins = sortwinarg(1):SPECWININCR:sortwinarg(2);
    else
        ampwins = ampargs(1);
    end
    if ~isinf(ampargs(1)) % single time given
        if length(freq) == 1
            %fprintf('   in a %1.1f-cycle (%1.0f ms) time window centered at %1.0f ms.\n',...
%                 cycles,1000/freq(1)*cycles,ampargs(1));
        else
            %fprintf('   in %1.1f-cycle (%1.0f-%1.0f ms) time windows centered at %1.0f ms.\n',...
%                 cycles,1000/freq(1)*cycles,1000/freq(end)*cycles,ampargs(1));
        end
    else % range of times
        [dummy sortwin_st ] = min(abs(times-ampwins(1)));
        [dummy sortwin_end] = min(abs(times-ampwins(end)));
        if length(freq) == 1
            %fprintf('   in %d %1.1f-cycle (%1.0f ms) time windows centered from %1.0f to  %1.0f ms.\n',...
%                 length(ampwins),cycles,1000/freq(1)*cycles,times(sortwin_st),times(sortwin_end));
        else
            %fprintf('   in %d %1.1f-cycle (%1.0f-%1.0f ms) time windows centered from %1.0f to %1.0f ms.\n',...
%                 length(ampwins),cycles,1000/freq(1)*cycles,1000/freq(end)*cycles,times(sortwin_st),times(sortwin_end));
        end
    end
    
    phsamps = 0; %%%%%%%%%%%%%%%%%%%%%%%%%% sort by (mean) amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%
    minxs = [];
    for f = 1:length(freq)  % use one or range of frequencies
        frq = freq(f);
        [amps, cohers, cohsig, ampsig, allamps, allphs] = ...
            phasecoher(data,length(times),srate,frq,cycles,0, ...
            [], [], timeStretchRef, timeStretchMarks);
        
        for ampwin = ampwins
            [dummy minx] = min(abs(times-ampwin)); % find nearest time point to requested
            minxs = [minxs minx];
            winlen = floor(cycles*srate/frq);
            % winloc = minx-[winlen:-1:0]; % ending time version
            winloc = minx-linspace(floor(winlen/2), floor(-winlen/2), winlen+1);
            tmprange = find(winloc>0 & winloc<=frames);
            winloc = winloc(tmprange); % sorting window frames
            if f==1
                winlocs = [winlocs;winloc];  % store tme windows
            end
            % $$$         [phaseangles phsamp] = phasedet(data,frames,srate,winloc,frq);
            % $$$         phsamps = phsamps+phsamps;  % accumulate amplitudes across 'sortwin'
            phaseangles = allphs(minx,:);
            phsamps = phsamps+allamps(minx,:);  % accumulate amplitudes across 'sortwin'
        end
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(tmprange) ~=  winlen+1 % ?????????
        filtersize = cycles * length(tmprange) / (winlen+1);
        timecenter = median(winloc)/srate*1000+times(1); % center of window in ms
        phaseangles = phaseangles + 2*pi*(timecenter-ampargs(1))*freq(end);
        %fprintf(...
%             '    Data time limits reached -> now uses a %1.1f cycles (%1.0f ms) window centered at %1.0f ms\n', ...
%             filtersize, 1000/freq(1)*filtersize, timecenter);
        %fprintf(...
%             '    Wavelet length is %d; Phase has been linearly interpolated to latency et %1.0f ms.\n', ...
%             length(winloc(1,:)), ampargs(1));
    end
    if length(freq) == 1
        %fprintf('Amplitudes are computed using a wavelet of %d frames.\n',length(winloc(1,:)));
    end
    
    %
    % Reject small (or large) phsamp trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    ampargs(2) = ampargs(2)/100; % convert rejection rate from % to fraction
    [tmp n] = sort(phsamps); % sort amplitudes
    if ampargs(2)>=0
        n = n(ceil(ampargs(2)*length(n))+1:end); % if rej 0, select all trials
%         fprintf('Retaining %d epochs (%g percent) with largest power at the analysis frequency,\n',...
%             length(n),100*(1-ampargs(2)));
    else % ampargs(2) < 0
        ampargs(2) = 1+ampargs(2); % subtract from end
        n = n(1:floor(ampargs(2)*length(n)));
%         fprintf(...
%             'Retaining %d epochs (%g percent) with smallest power at the analysis frequency,\n',...
%             length(n),ampargs(2)*100);
    end
    %
    % Remove low|high-amplitude trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    data = data(:,n); % amp-sort the data, removing rejected-amp trials
    phsamps = phsamps(n);           % amp-sort the amps
    phaseangles = phaseangles(n); % amp-sort the phaseangles
    sortvar = sortvar(n);         % amp-sort the trial indices
    ntrials = length(n);          % number of trials retained
    if ~isempty(auxvar)
        auxvar = auxvar(:,n);
    end
    %
    % Sort remaining data by amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    [phsamps sortidx] = sort(phsamps); % sort trials on amplitude
    data    =  data(:,sortidx);                % sort data by amp
    phaseangles  =  phaseangles(sortidx);      % sort angles by amp
    sortvar = sortvar(sortidx);                % sort input sortvar by amp
    if ~isempty(auxvar)
        auxvar = auxvar(:,sortidx);
    end
    %fprintf('Size of data = [%d,%d]\n',size(data,1),size(data,2));
    
    sortidx = n(sortidx); % return original trial indices in final sorted order
    %
    %%%%%%%%%%%%%%%%%%%%%% Don't Sort trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
elseif Nosort == YES
    %fprintf('Not sorting data on input sortvar.\n');
    sortidx = 1:ntrials;
    %
    %%%%%%%%%%%%%%%%%%%%%% Sort trials on (mean) value %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
elseif exist('valargs')
    [sttime stframe] = min(abs(times-valargs(1)));
    sttime = times(stframe);
    if length(valargs)>1
        [endtime endframe] = min(abs(times-valargs(2)));
        endtime = times(endframe);
    else
        endframe = stframe;
        endtime = times(endframe);
    end
    if length(valargs)==1 || sttime == endtime
%         fprintf('Sorting data on value at time %4.0f ms.\n',sttime);
    elseif length(valargs)>1
%         fprintf('Sorting data on mean value between %4.0f and %4.0f ms.\n',...
%             sttime,endtime);
    end
    if endframe>stframe
        sortval = mean(data(stframe:endframe,:));
    else
        sortval = data(stframe,:);
    end
    [sortval,sortidx] = sort(sortval);
    if length(valargs)>2
        if valargs(3) <0
            sortidx = sortidx(end:-1:1); % plot largest values on top
            % if direction < 0
        end
    end
    data = data(:,sortidx);
    sortvar = sortvar(sortidx);                % sort input sortvar by amp
    if ~isempty(auxvar)
        auxvar = auxvar(:,sortidx);
    end
    if ~isempty(phaseangles)
        phaseangles  =  phaseangles(sortidx);      % sort angles by amp
    end
    winloc = [stframe,endframe];
    winlocs = winloc;
    %
    %%%%%%%%%%%%%%%%%%%%%% Sort trials on sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
else
%     fprintf('Sorting data on input sortvar.\n');
    [sortvar,sortidx] = sort(sortvar);
    data = data(:,sortidx);
    if ~isempty(auxvar)
        auxvar = auxvar(:,sortidx);
    end
    uni_svar=unique_bc(sortvar);
    n_ties=0;
    tie_dist=zeros(1,length(uni_svar));
    loop_ct=0;
    for tie_loop=uni_svar,
        ids=find(sortvar==tie_loop);
        n_ids=length(ids);
        if n_ids>1,
            if replace_ties==YES,
                mn=mean(data(:,ids),2);
                data(:,ids)=repmat(mn,1,n_ids);
                if ~isempty(auxvar)
                    mn=mean(auxvar(:,ids),2);
                    auxvar(:,ids) = repmat(mn,1,n_ids);
                end
            end
            n_ties=n_ties+n_ids;
            loop_ct=loop_ct+1;
            tie_dist(loop_ct)=n_ids;
        end
    end
%     fprintf('%.2f%c of the trials (i.e., %d out of %d) have the same sortvar value as at least one other trial.\n', ...
%         100*n_ties/length(sortvar),37,n_ties,length(sortvar));
%     fprintf('Distribution of number ties per unique value of sortvar:\n');
    if exist('prctile')
        try
%             fprintf('Min: %d, 25th ptile: %d, Median: %d, 75th ptile: %d, Max: %d\n',min(tie_dist),round(prctile(tie_dist,25)), ...
%                 round(median(tie_dist)),round(prctile(tie_dist,75)),max(tie_dist));
        catch
        end;
    end;
    if replace_ties==YES,
%         fprintf('Trials with tied sorting values will be replaced by their mean.\n');
    end
    %fprintf('\n');
end

%if max(sortvar)<0
%   fprintf('Changing the sign of sortvar: making it positive.\n');
%   sortvar = -sortvar;
%end
%
%% %%%%%%%%%%%%%%%%% Adjust decfactor if phargs or ampargs %%%%%%%%%%%%%%%%%%%%%
%
if decfactor < 0
    decfactor = -decfactor;
    invdec = 1;
else
    invdec = 0;
end;
if decfactor > sqrt(ntrials) % if large, output this many trials
    n = 1:ntrials;
    if exist('phargs') & length(phargs)>1
        if phargs(2)>0
            n = n(ceil(phargs(2)*ntrials)+1:end); % trials after rejection
        elseif phargs(2)<0
            n = n(1:floor(phargs(2)*length(n)));  % trials after rejection
        end
    elseif exist('ampargs') & length(ampargs)>1
        if ampargs(2)>0
            n = n(ceil(ampargs(2)*ntrials)+1:end); % trials after rejection
        elseif ampargs(2)<0
            n = n(1:floor(ampargs(2)*length(n)));  % trials after rejection
        end
    end
    if invdec
        decfactor = (length(n)-avewidth)/decfactor;
    else
        decfactor = length(n)/decfactor;
    end;
end
if ~isreal(decfactor), decfactor = imag(decfactor); end;

%
%% %%%%%%%%%%%%%%%% Smooth data using moving average %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
urdata = data; % save data to compute amp, coher on unsmoothed data

if ~Allampsflag & ~exist('data2') % if imaging potential,
    if length(timeStretchRef) > 0 & length(timeStretchMarks) > 0
        %
        % Perform time-stretching here -JH %%%%%%%%%%%%%%%%
        %
        for t=1:size(data,2)
            M = timewarp(timeStretchMarks(t,:)', timeStretchRef');
            data(:,t) = M*data(:,t);
        end
    end
    tsurdata = data;
    % Time-stretching ends here %%%%%%%%%%%%%%
    
    if avewidth > 1 || decfactor > 1
        if Nosort == YES
%             fprintf('Smoothing the data using a window width of %g epochs ',avewidth);
        else
%             fprintf('Smoothing the sorted epochs with a %g-epoch moving window.',...
%                 avewidth);
        end
%         fprintf('\n');
%         fprintf('  and a decimation factor of %g\n',decfactor);
        
        if ~exist('phargs') % if not phase-sorted trials
            [data,outtrials] = movav(data,1:ntrials,avewidth,decfactor,[],[],wt_wind);
            % Note: movav() here sorts using square window
            [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor,[],[],wt_wind);
            
        else % if phase-sorted trials, use circular / wrap-around smoothing
            backhalf  = floor(avewidth/2);
            fronthalf = floor((avewidth-1)/2);
            if avewidth > 2
                [data,outtrials] = movav([data(:,[(end-backhalf+1):end]),...
                    data,...
                    data(:,[1:fronthalf])],...
                    [1:(ntrials+backhalf+fronthalf)],avewidth,decfactor,[],[],wt_wind);
                [outsort,outtrials] = movav([sortvar((end-backhalf+1):end),...
                    sortvar,...
                    sortvar(1:fronthalf)],...
                    1:(ntrials+backhalf+fronthalf),avewidth,decfactor,[],[],wt_wind);
                % Shift elements of outtrials so the first element is 1
                outtrials = outtrials - outtrials(1) + 1;
            else % avewidth==2
                [data,outtrials] = movav([data(:,end),data],...
                    [1:(ntrials+1)],avewidth,decfactor,[],[],wt_wind);
                % Note: movav() here sorts using square window
                [outsort,outtrials] = movav([sortvar(end) sortvar],...
                    1:(ntrials+1),avewidth,decfactor,[],[],wt_wind);
                % Shift elements of outtrials so the first element is 1
                outtrials = outtrials - outtrials(1) + 1;
            end
        end
        for index=1:length(percentiles)
            outpercent{index} = compute_percentile( sortvar, percentiles(index), outtrials, avewidth);
        end;
        if ~isempty(auxvar)
            if ~exist('phargs') % if not phase-sorted trials
                [auxvar,tmp] = movav(auxvar,1:ntrials,avewidth,decfactor,[],[],wt_wind);
            else % if phase-sorted trials
                if avewidth>2
                    [auxvar,tmp] = movav([auxvar(:,[(end-backhalf+1):end]),...
                        auxvar,...
                        auxvar(:,[1:fronthalf])],...
                        [1:(ntrials+backhalf+fronthalf)],avewidth,decfactor,[],[],wt_wind);
                    % Shift elements of tmp so the first element is 1
                    tmp = tmp - tmp(1) + 1;
                else % avewidth==2
                    [auxvar,tmp] = movav([auxvar(:,end),auxvar],[1:(ntrials+1)],avewidth,decfactor,[],[],wt_wind);
                    % Shift elements of tmp so the first element is 1
                    tmp = tmp - tmp(1) + 1;
                end
            end
        end
        %  if ~isempty(sortvar_limits),
        %      fprintf('Output data will be %d frames by %d smoothed trials.\n',...
        %          frames,length(outtrials));
        %      fprintf('Outtrials: %3.2f to %4.2f\n',min(outtrials),max(outtrials));
        %  end
    else % don't smooth
        outtrials = 1:ntrials;
        outsort = sortvar;
    end
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%% Find color axis limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty(Caxis)
        mindat = Caxis(1);
        maxdat = Caxis(2);
        %fprintf('Using the specified caxis range of [%g,%g].\n', mindat, maxdat);
    else
        mindat = min(min(data));
        maxdat = max(max(data));
        maxdat =  max(abs([mindat maxdat])); % make symmetrical about 0
        mindat = -maxdat;
        if ~isempty(caxfraction)
            adjmax = (1-caxfraction)/2*(maxdat-mindat);
            mindat = mindat+adjmax;
            maxdat = maxdat-adjmax;
%             fprintf(...
%                 'The caxis range will be %g times the sym. abs. data range -> [%g,%g].\n',...
%                 caxfraction,mindat,maxdat);
        else
%             fprintf(...
%                 'The caxis range will be the sym. abs. data range -> [%g,%g].\n',...
%                 mindat,maxdat);
        end
    end
end % if ~Allampsflag & ~exist('data2')

%
%% %%%%%%%%%%%%%%%%%%%%%%%% Set time limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isnan(timelimits(1))
    timelimits = [min(times) max(times)];
end
% fprintf('Data will be plotted between %g and %g ms.\n',timelimits(1),timelimits(2));

return;

