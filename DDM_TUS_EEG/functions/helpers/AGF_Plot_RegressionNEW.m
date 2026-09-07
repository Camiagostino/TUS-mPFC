%%%ENTER HERE CRITERIA
function [  ] = AGF_Plot_RegressionNEW( infold, outfold, settings )
%%
%requires:  AGF_shadedErrorBar
%           AGF_cmap.mat

%Correct folder name
if ~strcmp(infold(end),'/')
    infold=[infold '/'];
end;
if ~strcmp(outfold(end),'/')
    outfold=[outfold '/'];
end;

%%
for HideSettings=1
    %Set colors for multicolor plots
    colors       = {'g' 'r' 'b' 'c' 'r' 'k'};
    colors       = {'m' 'g' 'b' 'c' 'r' 'k'};

    %Set color for regressor values
    colorsSingle = 'b';
    colorsPval   = 'b';

    %Set defaults
    EEGplot         = 0;                        % Plot EEG timecourse only if it is present in output.
    transp          = 1;                        % Transparency for shades [1/0]
    PlotSteps       = 50;                       % plot every 50 ms
    MinMaplimits    = [];                       % define a minimum for maplimits (e.g. 1)
    PlotElec        = 'best';                   % Electrode to plot EEG bins for regressor
    PlotElec2       = {'Fz' 'FCz' 'Cz'};
    Confidence      = 0;                        % CI for shades, 0 = sem, 1 = SD
    ExtrapolS       = 100;
    ERP_Im_plot     = 0;                        % Does not assume an ERP image to be present (overwritten otherwise)
    zScoreR2        = 0;                        % Do not zscore R2 value before plotting.
    IndiviPlot      = 1;
    AllTime         = 0;   
    UseValues       = 'b';                      % Use 'b' or 't' values from regression output? g = group level t stat  
    printstring     = 'pdf';                    % Default fileformat for individual plots.
    IndivTypes      = {'Topo' 'Regressor' 'Pval' 'EEG'};
    MaskPval        = 0.05;
    Groups          = {};
    GroupNames      = {''};
    AddString       = ''; 
    TFstring        = '';
    Groupstring     = '';
    SPstring        = '';
    ReduceString    = '';
    GroupStats      = 0;
    ScalpPower      = 0;
    SPBands         = []; %default: no bands to average
    SPBands         = {[1 4]; [5 8]; [9 12]; [13 30]};
    SPNames         = {}; %default: no names for bands
    SPNames         = {'delta' 'theta' 'alpha' 'beta'};
    NN              = 0;
    FileName        = '';
    Level           = 1; %Default: First level average plot
    PLotP           = 1; %Default: Plots a small inset minimum p value in TF plots
    BinPlotTitleStr = 'ERP';
    %Read in regression output
    CMap=dir;
    if sum(ismember({CMap.name}, 'AGF_cmap.mat')) == 1
        Fcmap       = load('AGF_cmap.mat');
    else
        Fcmap       = load('functions/AGF_regression/AGF_cmap.mat');
    end;
    Fcmap       = Fcmap.AGF_cmap;
    if isfield(settings, 'FileName') && ~isempty(settings.FileName)
        F(1).name = settings.FileName;
    else
        F           = dir([infold '*.mat']);
    end;
    Ninclu      = size(F,1);
    
    %load first subjects and determine fields of output
    VP              = load([infold F(1).name]);
    FN1             = fieldnames(VP);
    FN2             = fieldnames(VP.(FN1{1}));
    
    %Test if more than one model is saved within the structure
    if size(FN2,1) > 2
        disp('Warning: Found more than one model in input. Models are: ')
        KPMN = 0;
        for PMN = 1 : 2 : size(FN2,1)
            KPMN = KPMN + 1;
            disp(['Model no ' num2str(KPMN) ': ' FN2{PMN}])
        end;
        MulM = input('\n\nWhich model do you want to use: ', 's');
        MulM = str2num(MulM);
        if MulM>1
            MulM = MulM * 2 - 1;
        end;
        MulI = MulM + 1;
    else
        MulM = 1; MulI = 2;
    end;
    
    Info            = VP.(FN1{1}).(FN2{MulI});
    Data            = VP.(FN1{1}).(FN2{MulM});
    dims            = size(Data.t_values);
    RegressionTime  = Info.Return_Timewindow_ms;
    Electrodes      = Info.Output_Labels;
    Elocs           = Info.Output_Electrodes;
    PlotImage       = 1;
    PlotRegress     = 1;
    PlotPval        = 1;
    PlotEEG         = 1;
    ERPImage_values = [];
    R2_values       = [];
    
    A=load(['data/ChanLocs.mat']);
    AF=fieldnames(A);
    chanlocs=A.(AF{1});
    ChanLoc = chanlocs(AGF_structfind(chanlocs,'labels',Electrodes)); %this is the reduced channel location file that includes only electrodes with data
    PlotReg         = [];
    UPlotReg       = 1 : length(Info.RegNames); %Default: plot all regressors
    
    %%%%%%%%SHOULD BE REMOVED; FOR COMPATIBILITY ONLY%%%%%%%%
    try Info.IncluReg
    catch
        Info.IncluReg = 1 : length(Info.RegNames);
    end;
    
    %Get input from settings
    if isfield(settings, 'plotsteps');  PlotSteps       = settings.plotsteps;       end;
    if isfield(settings, 'maplimits');  MinMaplimits    = settings.maplimits;       end;
    if isfield(settings, 'plotelect');  PlotElec        = settings.plotelect;       end;
    if isfield(settings, 'confidence'); Confidence      = settings.confidence;      end;
    if isfield(settings, 'plotelec2');  PlotElec2       = settings.plotelec2;       end;
    if isfield(settings, 'plotImage');  plotImage       = settings.plotImage;       end;
    if isfield(settings, 'PlotRegress');PlotRegress     = settings.PlotRegress;     end;
    if isfield(settings, 'plotPval');   plotPval        = settings.plotPval;        end;
    if isfield(settings, 'plotEEG');    plotEEG         = settings.plotEEG;         end;
    if isfield(settings, 'ExtrapolS');  ExtrapolS       = settings.ExtrapolS;       end;
    if isfield(settings, 'ERP_Im_plot');ERP_Im_plot     = settings.ERP_Im_plot;     end;
    if isfield(settings, 'AllTime');    AllTime         = settings.AllTime;         end;
    if isfield(settings, 'IndiviPlot'); IndiviPlot      = settings.IndiviPlot;      end;
    if isfield(settings, 'IndivTypes'); IndivTypes      = settings.IndivTypes;      end;
    if isfield(settings, 'MaskPval');   MaskPval        = settings.MaskPval;        end;
    if isfield(settings, 'printstring');printstring     = settings.printstring;     end;
    if isfield(settings, 'Groups');     Groups          = settings.Groups;          end;
    if isfield(settings, 'GroupNames'); GroupNames      = settings.GroupNames;      end;
    if isfield(settings, 'UseValues');  UseValues       = settings.UseValues;       end;
    if isfield(settings, 'AddString');  AddString       = settings.AddString;       end;
    if isfield(settings, 'GroupStats'); GroupStats      = settings.GroupStats;      end;
    if isfield(settings, 'ScalpPower'); ScalpPower      = settings.ScalpPower;      end;
    if isfield(settings, 'SPBands');    SPBands         = settings.SPBands;         end;
    if isfield(settings, 'SPNames');    SPNames         = settings.SPNames;         end;
    if isfield(settings, 'PlotReg');    PlotReg         = settings.PlotReg;         end;
    if isfield(settings, 'level');      Level           = settings.level;           end;
    if isfield(settings, 'PLotP');      PLotP           = settings.PLotP;           end;
    
    if ~isempty(Info.IncluReg) %not all regressors are included
        UPlotReg = intersect(Info.IncluReg,UPlotReg);
        ReduceString = [ReduceString 'Data contains ' num2str(UPlotReg) ' regressors.\n'];
    end;
    
    if ~isempty(PlotReg)
        UPlotReg = intersect(PlotReg,UPlotReg);
        ReduceString = [ReduceString num2str(length(UPlotReg)) ' regressors will be plotted.\n'];
    end;
    dims(2) = length(UPlotReg);
    
    if GroupStats && Level == 2
        warning('2nd level plotting assumes already group statistics as input. Setting GroupdStats to 0.'); GroupStats = 0;
    end;
        
    if Ninclu == 1 && Level == 1
        warning('Only one dataset read in but assuming first-level instead of second level output. Is this correct? Otherwise, set settings.level = 2.'); pause(1)
    end;
    %%%%%%%%%%%%
    %Pre-allocate arrays, separately for TF or time domain data
    try TF = Info.TF;
    catch
        try TF = Data.TF;
        catch
            TF = 'Time domain';
        end;
    end;
    if isstruct(TF)
        Ind_Plt_elec2 = AGF_structfind(Info,'Output_Labels',PlotElec2);
        regress_values  = nan(Ninclu, length(Ind_Plt_elec2), dims(2), dims(3), dims(4));
        if isfield(Data,'PseudoR2')
            R2_values = nan(Ninclu, length(Ind_Plt_elec2), dims(2)+1, dims(3), dims(4));
        end;
        TFstring = ['Plotting TF decomposed regression results over ' num2str(length(Ind_Plt_elec2)) ' electrodes.\n'];
        TFstring = [TFstring '      ---> Frequency range from ' num2str(TF.frequencies(1)) ' to ' num2str(TF.frequencies(2)) ' Hz with ' num2str(TF.stepnumber) ' ' TF.space ' steps.\n'];
                
        if ScalpPower
            if isempty(SPBands)
                disp(['To plot scalp power topographies, fequency bands have to be provided in field settings.SPBands']); return;
            end;
            BandNumber = length(SPBands); %number of frequencies to average over
            
            if length(length(SPNames))~=BandNumber
                disp('Length of names for frequencies to plot does not match number of frequency ranges. Generating new names.')
                SPNames = {}; NN = 1;
            end;
            
            for CM = 1 : BandNumber % find closest matches to these frequencies in the output
                [~,CloseInd1]       = min(abs(TF.freq-SPBands{CM}(1)));
                usedfreq{CM}(1)     = TF.freq(CloseInd1);
                [~,CloseInd2]       = min(abs(TF.freq-SPBands{CM}(2)));
                usedfreq{CM}(2)     = TF.freq(CloseInd2);
                CollapseIndex(CM,:) = [CloseInd1 CloseInd2];
                if NN %generate new names
                    SPNames(CM) = {['Band ' num2str(round(100*usedfreq{CM}(1))/100) ' - ' num2str(round(100*usedfreq{CM}(2))/100)]};
                end;
            end;

            SPstring = [SPstring 'Averaging over ' num2str(BandNumber) ' frequency ranges for scalp topography plots.']
        end;
    else 
        regress_values  = nan(Ninclu, dims(1), dims(2), dims(3));
        if isfield(Data,'PseudoR2')
            R2_values = nan(Ninclu, dims(1), dims(2)+1, dims(3));
        end;
        Ind_Plt_elec2 = 1 : length(Elocs);
    end;
    %%%%%%%%%%%%%
    Nsubjects=Ninclu;
    if Level==2
        disp(['Model contains 2nd level data of ' num2str(Info.TotalTrials) ' subjects.'])
        if strcmp(Info.UseValues2ndLevel,'t') || strcmp(Info.UseValues2ndLevel,'b')
            BinPlotTitleStr = ['first level ' Info.UseValues2ndLevel];
        else
            BinPlotTitleStr = ['first level ERP'];
        end;
        disp(['2nd level model was run on ' BinPlotTitleStr ' data.'])
        Nsubjects=Info.TotalTrials;
    else
         %%%%%%%%%%%%
        %Test if groups are present
        if max([Groups{:}])>Ninclu
            disp('Error: Groups contain more values than input folder datasets.');return;
        else
            Groupstring = ['Comparing ' num2str(length(Groups)) ' groups:\n'];
            for GC = 1 : length(Groups)
                if ~isempty(GroupNames)
                    Groupstring = [Groupstring '   ' num2str(GC) ' = ' GroupNames{GC} ' with n = ' num2str(length(Groups{GC})) '\n'];
                else
                    Groupstring = [Groupstring '   ' num2str(GC) ' = Group ' num2str(GC) ' with n = ' num2str(length(Groups{GC})) '\n'];
                    if GC == length(Groups)
                        Groupstring = [Groupstring '\tNote: To specify names for groups, just set value settings.GroupNames.\n'];
                    end;
                end;
            end;
            Groupstring = [Groupstring '\n'];
        end;
    end;
    
    %Adjust shades
    if Confidence>0 && Confidence<1
        CString = [num2str(Confidence*100) '% CI'];
        ConfFactor = norminv(1+(0.5-(1-Confidence/2)),0,1);  % if confidence is set to calculate CI, scale factor for SE calculation
    elseif Confidence==0
        ConfFactor = 1;
        CString = 'SE';
    elseif Confidence==1
        ConfFactor = sqrt(Nsubjects-1);
        CString = 'SD';
    end;

    if ~isfield(settings, 'plottime')
        PlotTime = [Info.Return_Timewindow_ms(1) Info.Return_Timewindow_ms(end)];
    else
        PlotTime = settings.plottime;
    end;
    
    if isfield(Data, 'EEG_per_regressor')
        EEGplot       = 1;
        dims2         = size(VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor);
        EEG_values    = nan(Ninclu, length(UPlotReg), dims2(2), dims2(3), dims2(4));
        EEG_SD_values = EEG_values;
    end;
    
    if Level == 1
        if strcmp(UseValues,'t')
                WhatIsIt = 'average first level t values';
                if GroupStats
                    WhatIsIt = 'group level t-stat of within t';
                end;
        elseif strcmp(UseValues,'b')
            WhatIsIt = 'average first level b values';
            if GroupStats
                WhatIsIt = 'group level t-stat of within b';
            end;
        end;
    else
        WhatIsIt = '2nd level t values';
    end;
    
    %Test if time windows are okay?
    if PlotTime(1)<RegressionTime(1) || PlotTime(end)>RegressionTime(end)
        disp(['Plot time not found in regression output. Stopping. Plottime: ' num2str(PlotTime(1)) ' to ' num2str(PlotTime(end)) ' and limits of regression: ' num2str(RegressionTime(1)) ' to ' num2str(RegressionTime(end)) '.']);return;
    end;
    
    %Display information on prompt about input and model
    fprintf('\n\n************Starting Plot Regression Script...***********\n')
    if Level == 1
        disp(['Found ' num2str(Ninclu) ' files in folder ' infold '.'])
    end;
    fprintf(Groupstring)
    disp(['Regresison model ' strrep(FN2{MulM},'_',' ') ' includes ' num2str(length(Info.Output_Electrodes)) ' electrodes and ' num2str(min(Info.Return_Timewindow_ms)) ' until ' num2str(max(Info.Return_Timewindow_ms)) 'ms.'])
    disp(['Number of regressors is ' num2str(length(Info.RegNames)) ' and time window for plotting ' num2str(min(PlotTime)) ' to ' num2str(max(PlotTime)) ' ms.'])
    fprintf(TFstring)
    fprintf(ReduceString)
    
    
    if isfield(Data, 'ERP_Image') && ERP_Im_plot
        dims3       = size(VP.(FN1{1}).(FN2{MulM}).ERP_Image);
        if ExtrapolS < dims3(4) %more data availible than extrapolation should be --> do not extrapolate!
            ExtrapolS = dims3(4);
        end;
        
        if prod([Ninclu length(UPlotReg) dims3(2) ExtrapolS dims3(3)])>10^9 %VERY large array, use slower cumulative averaging instead
            %Dimensions are:     Regressor | Electrode | Image Y | Image X
            ERPImage_values= nan(length(UPlotReg), dims3(2), ExtrapolS, dims3(3));
            disp('********')
            disp(['Warning: ERP images are very large, using cumulative averaging (can be slow)!'])
            disp('********')
        else
            if Level==2 || (Ninclu == 1 && Level == 1)
                ERPImage_values= nan(length(UPlotReg), dims3(2), ExtrapolS, dims3(3));
            else
                %Dimensions are: VP | Regressor | Electrode | Image Y | Image X
                ERPImage_values= nan(Ninclu, length(UPlotReg), dims3(2), ExtrapolS, dims3(3));
            end;
        end;
        
        disp(['Will extrapolate ERP image(s) to ' num2str(ExtrapolS) ' datapoints.' ])
    end;
end;
clear Data %remove dummy Data from first subject

%%
%Load data
if Level==2 || (Ninclu == 1 && Level == 1)
    p_values        = VP.(FN1{1}).(FN2{MulM}).Robust_p(Ind_Plt_elec2,UPlotReg,:,:);
    R2_values       = VP.(FN1{1}).(FN2{MulM}).PseudoR2(Ind_Plt_elec2,1,:,:);
    regress_values  = nan([1, size(p_values)]);
    Ninclu = 1;
    Groups = [];
    
    if ~isstruct(TF)
        EEG_values = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,Ind_Plt_elec2,:);  %4D --> Regressor | bins | electrode | time
        if isfield(VP.(FN1{1}).(FN2{MulM}),'EEG_SD_per_regressor')
            EEG_SD_values = VP.(FN1{1}).(FN2{MulM}).EEG_SD_per_regressor(UPlotReg,:,Ind_Plt_elec2,:);
        end;
    end;
    if Level == 2
        AddString = [AddString ' for ' Info.ModelName]; 
    end;
end;

for p = 1 : Ninclu
    %Get fitting INFO file
    VPname = F(p).name(1:end-4);
    disp(['Processing ' VPname])
    VP = load([infold F(p).name]);
    if isstruct(TF)
        %TF plot discards electrodes that are not included in the plot (save memory)
        if strcmp(UseValues,'t')
            regress_values(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).t_values(Ind_Plt_elec2,UPlotReg,:,:);
        elseif strcmp(UseValues,'b')
            regress_values(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).b_values(Ind_Plt_elec2,UPlotReg,:,:);
        else
            disp(['Unrecognized value for parameter UseValues: ' UseValues]);return;
        end;
        
        if ~isempty(R2_values) && Level == 1
            R2_values(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).PseudoR2(Ind_Plt_elec2,[1 UPlotReg],:,:);
        end;
        
    else %Load Timedomain Data
        %Save to grandtotal matrix
        if strcmp(UseValues,'t')
            regress_values(p,:,:,:) = VP.(FN1{1}).(FN2{MulM}).t_values(:,UPlotReg,:);
        elseif strcmp(UseValues,'b')
            regress_values(p,:,:,:) = VP.(FN1{1}).(FN2{MulM}).b_values(:,UPlotReg,:);
        else
            disp(['Unrecognized value for parameter UseValues: ' UseValues]);return;
        end;
        
        if ~isempty(R2_values)  && Level == 1 && Ninclu ~= 1
            R2_values(p,:,:,:) = VP.(FN1{1}).(FN2{MulM}).PseudoR2(:,[1 UPlotReg],:);
        end;
        
        if exist('dims3') && ERP_Im_plot
            for Ec = 1 : size(Electrodes,2)
                for Rc = 1 : dims3(1)
                    GridInt = griddedInterpolant([squeeze(VP.(FN1{1}).(FN2{MulM}).ERP_Image(UPlotReg(Rc),Ec,:,:))]');
                    vec     = AGF_exact_veclength(1, size([squeeze(VP.(FN1{1}).(FN2{MulM}).ERP_Image(UPlotReg(Rc),Ec,:,:))]',1), ExtrapolS);
                    [Xq,Yq] = ndgrid(vec,1:dims3(3));
                    if prod(size(ERPImage_values),2) > 10^9 %use cumulative averaging
                        if p==1
                            ERPImage_values(Rc,Ec,:,:) = GridInt(Xq,Yq);
                        else
                            ERPImage_values(Rc,Ec,:,:) = (squeeze(ERPImage_values(UPlotReg(Rc),Ec,:,:)).*(p-1) + GridInt(Xq,Yq)) / p;
                        end;
                    else
                        if Level==2 || (Ninclu == 1 && Level == 1)
                            ERPImage_values(Rc,Ec,:,:) = GridInt(Xq,Yq);
                        else
                            ERPImage_values(p,Rc,Ec,:,:) = GridInt(Xq,Yq);
                         end;
                    end;
                end;
            end;
        end;
        
        if EEGplot
            if Level == 1 && Ninclu ~= 1
                EEG_values(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,:,:,:); %5D --> VP | Regressor | bins | electrode | time
                if isfield(VP.(FN1{1}).(FN2{MulM}),'EEG_SD_per_regressor')
                    EEG_SD_values = VP.(FN1{1}).(FN2{MulM}).EEG_SD_per_regressor(UPlotReg,:,:,:,:);
                end;
                CountCat = 0;
                for c = 1 : dims(2)
                    Ireg = UPlotReg(c);
                    dimsEEG = size([VP.(FN1{1}).(FN2{MulM}).EEG_Values{Ireg}]);
                    if dimsEEG(1)<dimsEEG(2) %parametric regressor
                        CountCat = CountCat +1;
                        for c2 = 1 : dimsEEG(2)
                            CatReg(p,CountCat,c2) = VP.(FN1{1}).(FN2{MulM}).EEG_Values{Ireg}(c2); % cannot easily be pre-allocated as the size depends on the maximum number of splits in EEG data
                        end;
                    end;
                    if p == length(Ninclu) %average values for parametric regressors
                        if dimsEEG(1)<dimsEEG(2)
                            ALL_VP_Regressor_Categories{c} = squeeze(mean(CatReg(:,CountCat,:),1));
                        else
                            ALL_VP_Regressor_Categories{c} = VP.(FN1{1}).(FN2{MulM}).EEG_Values{Ireg};
                        end;
                    end;
                end;
            end;
        end;
    end;
end;
%%
%Average ERP images if not done before...
if prod(size(ERPImage_values),2) < 10^9 && Level == 1 && Ninclu ~= 1
    ERPImage_values = squeeze(mean(ERPImage_values,1));
end;
%Adjust for inverted Regression
if strcmp(Info.RegNames{end},'EEG')
    Info.RegLables{end+1}={'numeric'; 'numeric'};
end;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check if multiple groups are specified? 
GL = length(Groups);
if GL == 0 
    GL = 1;
else
    disp(['Plotting for ' num2str(GL) ' different groups.'])
end;

if ~exist(outfold)
    mkdir(outfold);
end

close all; %close other plots
figure;
for GroupCount = 1 : GL
    %Addgroupname if missing
    Gstring = ''; UseGroupName='';
    if isempty(Groups)
        GroupVP = 1 : Ninclu;
    else
        GroupVP = Groups{GroupCount};
        try 
            if isempty(GroupNames{GroupCount})
                UseGroupName = ['Group' num2str(GroupCount)];
                Gstring = [' in group ' UseGroupName];
            elseif ~strcmp(GroupNames{GroupCount}(1),' ')
                UseGroupName = [' ' GroupNames{GroupCount}];
                Gstring = [' in group ' UseGroupName];
            else
                UseGroupName = GroupNames{GroupCount};
                Gstring = [' in group ' UseGroupName];
            end;
        catch
            UseGroupName = ['Group' num2str(GroupCount)];
            Gstring = '';
        end;
    end;
    
    clear AllPlots
    if length(PlotSteps>1)
        Ptimes=PlotSteps;
    else
        Ptimes=[PlotTime(1) : PlotSteps : PlotTime(2)];
    end;
    DataTime = find(ismember(RegressionTime,Ptimes)); %Create index for regression datapoints
    AnyDataTime = find(ismember(RegressionTime,PlotTime(1))):find(ismember(RegressionTime,PlotTime(2)));
    %Time Frequency Plots
    if isstruct(TF)
        f_xsize= 2000;
        f_ysize= 1400;
        pHeight = 0.105;
        pWidth = 0.05;
        NumberOfYTicks = 4; 
        XTickSteps = 200; %distance between X Tick Labels
        clc
        ShrinkFactor    = 0.8;      %shrinks the plot a little bit
        OffSet          = 0.1;      %offest from lower left corner  
        %determine small string with short model summary
        ti3 = sprintf([strrep(FN2{MulM},'_',' ')]);
        AbsMinP = 1;
        for T = 1 : length(Info.RegNames)
            if ismember(T,UPlotReg) %this regressor will be plotted
                ti3=[ti3 sprintf(strrep(['\nReg: ' Info.RegNames{T} ' - ' Info.RegLables{T,1} '(' num2str(Info.RegValues{T}(1)) ') & ' Info.RegLables{T,2} '(' num2str(Info.RegValues{T}(2)) ')'  ], '_', ' '))];
            else
                ti3=[ti3 sprintf(strrep(['\nReg (not plotted): ' Info.RegNames{T} ' - ' Info.RegLables{T,1} '(' num2str(Info.RegValues{T}(1)) ') & ' Info.RegLables{T,2} '(' num2str(Info.RegValues{T}(2)) ')'  ], '_', ' '))];
            end;
        end;
        
        for R =  1 : length(UPlotReg) + 1
            if ishandle(1)
                clf; %clear figure 
            end;
            hFig = figure(1);
            set(hFig, 'Position', [100 100 f_xsize f_ysize])%
            handaxes2 = axes('Position', [OffSet OffSet 1 1]);
            set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');
            
            text(-pWidth, 0.88,ti3,'HorizontalAlignment','left','VerticalAlignment', 'top', 'FontSize', 11)
            if R <= length(UPlotReg) %plot normal effects
                %Determine maplimits overall
                
                if GroupStats
                    [~,~,~,stats]=ttest(regress_values(GroupVP,:,R,AnyDataTime,:));
                    Limi = [-max(max(max(abs(squeeze(stats.tstat))))) max(max(max(abs(squeeze(stats.tstat)))))];
                    ti3 = [ti3 sprintf(['\n\nT limits in p = ' num2str(1-tcdf(Limi(2), stats.df(1)))])];
                else
                    Limi = [-max(max(max(abs(squeeze(nanmean(regress_values(GroupVP,:,R,AnyDataTime,:))))))) max(max(max(abs(squeeze(nanmean(regress_values(GroupVP,:,R,AnyDataTime,:)))))))];
                end;
                
                if Level == 1
                    ti=['TF Regressor: ' Info.RegNames{UPlotReg(R)} ' - values: ' Info.RegLables{UPlotReg(R),1} '(' num2str(Info.RegValues{UPlotReg(R)}(1)) ') and '...
                        Info.RegLables{UPlotReg(R),2} '(' num2str(Info.RegValues{UPlotReg(R)}(2)) ') ' UseGroupName]; 
                    PSavename = [outfold  Info.RegNames{UPlotReg(R)} UseGroupName AddString];
                else
                    ti=['2nd Level TF on ' FN2{MulM} ' for: ' Info.RegNames{UPlotReg(R)} ' - values: ' Info.RegLables{UPlotReg(R),1} '(' num2str(Info.RegValues{UPlotReg(R)}(1)) ') and '...
                        Info.RegLables{UPlotReg(R),2} '(' num2str(Info.RegValues{UPlotReg(R)}(2)) ') ' UseGroupName]; 
                    PSavename = [outfold  Info.RegNames{UPlotReg(R)} ' on ' FN2{MulM} ' ' UseGroupName AddString];
                    %Overwrite limits
                    Limi = [-max(max(max(abs(regress_values(GroupVP,:,R,AnyDataTime,:))))) max(max(max(abs(regress_values(GroupVP,:,R,AnyDataTime,:)))))];
                end;
                ti = strrep(ti, '_', ' ');
                ti2 = ['Maplimits: ' num2str(Limi(1)) ' to ' num2str(Limi(2)) ' (' WhatIsIt ')'];
                if MaskPval
                    ti2 = [ti2 ' & masked with p < ' num2str(MaskPval)]; 
                end;
                text(0.5-pWidth, 0.9,ti,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)
                text(0.5-pWidth, 0.87,ti2,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 15)
                for E = 1 : length(Ind_Plt_elec2)
                    ax =  [ChanLoc(Ind_Plt_elec2(E)).X]; ax(abs(ax)<0.01)=0;
                    ay =  [ChanLoc(Ind_Plt_elec2(E)).Y]*ShrinkFactor*-1; ay(abs(ay)<0.01)=0;
                    if abs(ax)>5
                        ax=ax/100; ay=ay/100; %az=az/100;
                    end;
                    ax = (ax / 2 + 0.5 - pHeight / 2) * ShrinkFactor+OffSet;
                    ay = (ay / 2 + 0.5 - pWidth  / 2) * ShrinkFactor+OffSet;

                    handaxes2 = axes('Position', [ay ax pHeight pWidth]);
                    
                    %Determine p-value across subjects
                    if Level == 1
                        [~,p,~,stats]=ttest(squeeze(regress_values(GroupVP,E,R,AnyDataTime,TF.stepnumber:-1:1)));
                        p=squeeze(p)';
                    else
                        p=squeeze(p_values(E,R,:,TF.stepnumber:-1:1))';
                    end;

                    if GroupStats %Plot second level t-stat
                        ThePlot = squeeze(stats.tstat)';
                    else %Plot average of first level stats
                        if Level == 1
                            ThePlot = [squeeze(nanmean(regress_values(GroupVP,E,R,AnyDataTime,TF.stepnumber:-1:1)))]';
                        else
                            ThePlot = [squeeze(regress_values(1,E,R,:,TF.stepnumber:-1:1))]';
                        end;
                    end;
                    
                    if MaskPval
                        ThePlot(p>MaskPval)=0;
                    end;
                    imagesc(ThePlot, Limi); 
                    set(handaxes2, 'Box', 'off')
                    handaxes2.XTick = find(ismember(RegressionTime(AnyDataTime),intersect([-1e+4:XTickSteps:1e+4], [PlotTime+1:PlotTime(end)-1])));
                    handaxes2.XTickLabel = intersect([-1e+4:XTickSteps:1e+4], [PlotTime+1:PlotTime(end)-1]);
                    handaxes2.YTick = round(quantile([1:diff(handaxes2.YLim)],4));
                    handaxes2.YTickLabel = fliplr(round(TF.freq(round(quantile([1:diff(handaxes2.YLim)],NumberOfYTicks)))*10)/10);
                    title([PlotElec2{E}]);
                    if PLotP 
                        text(0.01, 0.1,[' min p = ' num2str(min(min(p)))],'HorizontalAlignment','left','VerticalAlignment', 'top', 'FontSize', 9);
                    end;
                    colormap(Fcmap);  
                    pause(0.001)
                end;
                %create dummy axis with colorbar
                handaxes2 = axes('Position', [0.95 0.93 0.01 pWidth ]);
                imagesc([0 0], Limi);
                set(gca,'Visible','off'); %white image without axes
                colorbar; handaxes2.XTickLabel = []; handaxes2.YTickLabel = []; 
                saveas (hFig,PSavename,printstring);
                
                if ScalpPower %Topographies for collapsed scalp powers for this regressor 
                    Spalte = length(DataTime) + 3; 
                    Zeile = BandNumber;
                    clear AllPlots;
                    AllPlots.Handle=[];

                    if ishandle(1)
                        clf; %clear figure 
                    end; 
                    hFig = figure(1); set(gcf,'Visible', 'on'); pc=1; HC=1;
                    set(hFig, 'Position', [100 100 1800 (Zeile)*200])%
                    set(hFig,'PaperPositionMode','Auto')
                    figureSize = get(gcf,'Position');
                    ti = strrep(['Topographies for collapsed band power for reg ' Info.RegNames{UPlotReg(R)} '\n in ' FN2{MulM} Gstring AddString],'_',' ');
                    DatPlot = squeeze(nanmean(regress_values(GroupVP,:,R,DataTime,:),1));
                    UseLim = [-max(abs([min(min(min(DatPlot))) max(max(max(DatPlot)))])) max(abs([min(min(min(DatPlot))) max(max(max(DatPlot)))]))];
                    %Determine topoplot maplimits for this band frequency
                    for BP = 1 : BandNumber
                        hand(HC)=subplot(Zeile,Spalte,[pc pc+2]); hold(gca, 'on');HC=HC+1;pc=pc+2;
                        PSTR = sprintf(['Collapsed band-power for ' SPNames{BP} '\nRanging from ' num2str(round(100*usedfreq{BP}(1))/100) ' to ' num2str(round(100*usedfreq{BP}(2))/100) ' Hz\nPlotlimits = ' num2str(round(100*UseLim(1))/100) ' to ' num2str(round(100*UseLim(2))/100) ]);
                        text(0,0.5,PSTR,'FontSize',12);
                        set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');

                        for TopoCount = 1 : length(DataTime)
                            hand(HC)=subplot(Zeile,Spalte,TopoCount+pc); hold(gca, 'on');
                            topoplot(squeeze(mean(DatPlot(:,TopoCount,CollapseIndex(BP,1):CollapseIndex(BP,2)),3)), ChanLoc, 'maplimits', UseLim);colormap(Fcmap);        
                            title([num2str(Ptimes(TopoCount)) ' ms']);
                            AllPlots(length(AllPlots(1).Handle)+1).Type    = 'TFTopo';
                            AllPlots(length(AllPlots(1).Handle)+1).Title   = ['Band' num2str(BP) ' TF Topo ' Info.RegNames{UPlotReg(R)} ' ' num2str(Ptimes(TopoCount)) ' ms'];
                            AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
                            pause(0.00001)
                        end;
                        pc = pc + TopoCount + 1;
                    end;
                    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                    text(0.5, 0.99,sprintf(ti),'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)
                    saveas (hFig,[outfold 'Power Topo ' Info.RegNames{UPlotReg(R)} UseGroupName AddString],printstring)
                end;
                
            else %Plot R2 if present
                if ~isempty(R2_values)
                    %Determine maplimits overall
                    if Level == 1
                        if zScoreR2
                            Limi = [-max(max(max(abs(squeeze(zscore(mean(R2_values(GroupVP,:,1,AnyDataTime,:)))))))) max(max(max(abs(squeeze(zscore(mean(R2_values(GroupVP,:,1,AnyDataTime,:))))))))];
                            R2Str = 'z-scored R2';
                        else
                            Limi = [0 max(max(max(abs(squeeze(mean(R2_values(GroupVP,:,1,AnyDataTime,:)))))))];
                            R2Str = 'actual R2';
                        end;
                    else
                        if zScoreR2
                            Limi = [-max(max(max(abs(squeeze(zscore(R2_values(:,1,AnyDataTime,:))))))) max(max(max(abs(squeeze(zscore(R2_values(:,1,AnyDataTime,:)))))))];
                            R2Str = 'z-scored R2';
                        else
                            Limi = [0 max(max(max(abs(squeeze(R2_values(:,1,AnyDataTime,:))))))];
                            R2Str = 'actual R2';
                        end;
                    end;
                    if Level == 1
                        ti=['R2 values for the whole model including ' num2str(R-1) ' regressors ' UseGroupName]; 
                        PSavename = [outfold  ' R2 ' UseGroupName AddString];
                    else
                        ti=['2nd level R2 values for the whole model including ' num2str(R-1) ' regressors ' UseGroupName]; 
                        PSavename = [outfold  ' R2 on ' FN2{MulM} UseGroupName AddString];
                    end;
                    ti = strrep(ti, '_', ' ');
                    ti2 = ['Maplimits: ' num2str(Limi(1)) ' to ' num2str(Limi(2)) ' ' R2Str];
                    text(0.5-pWidth, 0.9,ti,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)
                    text(0.5-pWidth, 0.87,ti2,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 15)
                    for E = 1 : length(Ind_Plt_elec2)
                        ax =  [ChanLoc(Ind_Plt_elec2(E)).X]; ax(abs(ax)<0.01)=0;
                        ay =  [ChanLoc(Ind_Plt_elec2(E)).Y]*ShrinkFactor*-1; ay(abs(ay)<0.01)=0;
                        if abs(ax)>5
                            ax=ax/100; ay=ay/100; %az=az/100;
                        end;
                        ax = (ax / 2 + 0.5 - pHeight / 2) * ShrinkFactor+OffSet;
                        ay = (ay / 2 + 0.5 - pWidth  / 2) * ShrinkFactor+OffSet;

                        handaxes2 = axes('Position', [ay ax pHeight pWidth]);
                        set(handaxes2, 'Box', 'off')
                        if Level == 1
                            ThePlot = [squeeze(mean(R2_values(GroupVP,E,1,AnyDataTime,TF.stepnumber:-1:1)))]';
                        else
                            ThePlot = [squeeze(R2_values(E,1,:,TF.stepnumber:-1:1))]';
                        end;
                        
                        if zScoreR2
                            imagesc(zscore(ThePlot), Limi);
                        else
                            imagesc(ThePlot, Limi);
                        end;
                        handaxes2.XTick = find(ismember(RegressionTime(AnyDataTime),intersect([-1e+4:XTickSteps:1e+4], [PlotTime+1:PlotTime(end)-1])));
                        handaxes2.XTickLabel = intersect([-1e+4:XTickSteps:1e+4], [PlotTime+1:PlotTime(end)-1]);
                        handaxes2.YTick = round(quantile([1:diff(handaxes2.YLim)],4));
                        handaxes2.YTickLabel = fliplr(round(TF.freq(round(quantile([1:diff(handaxes2.YLim)],NumberOfYTicks)))*10)/10);
                        title(PlotElec2{E});
                        colormap(Fcmap);  
                        pause(0.001)
                    end;
                     %create dummy axis with colorbar
                    handaxes2 = axes('Position', [0.95 0.93 0.01 pWidth ]);
                    imagesc([0 0], Limi);
                    set(gca,'Visible','off'); %white image without axes
                    colorbar; handaxes2.XTickLabel = []; handaxes2.YTickLabel = []; 
                    saveas (hFig,PSavename,printstring);
                end;
            end;                
        end;

    else %ERP-like Plots
        for HideErpPlot = 1
            %Image has dimensions: 3 (Info) | n plots | ? 3 Regressor | ? 3 Pval | ? 5 EEG     | ? 5 ERPimage
            Spalte = length(DataTime)        +4       + 3*PlotRegress + 3*PlotPval + 5*PlotEEG + 5*PlotImage; 
            AllPlots.Handle=[];
            Zeile = dims(2)+1;
            if Zeile == 2 %slightly larger plot if only one line of regressors (avoids overlap)
                HeightPerRow = 300;
            else
                HeightPerRow = 200;
            end;
            disp(['Generating ' num2str(length(Ptimes)) ' topo plots...'])

            if ishandle(1)
                clf; %clear figure 
            end;
            hFig = figure(1); set(gcf,'Visible', 'on'); 
            set(hFig, 'Position', [100 100 2400 (dims(2)+1)*HeightPerRow])%
            set(hFig,'PaperPositionMode','Auto')
            figureSize = get(gcf,'Position');
            ti = strrep(['Overview at ' PlotElec ' for ' strrep(FN2{MulM},'_',' ') Gstring AddString],'_',' ');
            HC = 1; %Handle Counter to export subplots

            for c = 1 : Zeile
                disp(['Plotting row ' num2str(c) ' / ' num2str(Zeile)])
                %%%%%%%%%%%R2 Plot and General Information%%%%%%%%%%%%%
                if c == 1 %First row is general R2 of full model
                    if isempty(find(isnan(R2_values)==1)) && ~isempty(R2_values)
                        %Plot first field with general information
                        hand(HC)=subplot(Zeile,Spalte,[1 2]); hold(gca, 'on');HC=HC+1;

                        %Make R2 usable for plot
                        if Level == 1 && Ninclu > 1
                            All_Plot_R2 = squeeze(nanmean(R2_values(GroupVP,:,1,DataTime),1));
                        else
                            All_Plot_R2 = squeeze(R2_values(:,1,DataTime));
                        end;
                        if zScoreR2
                            All_Plot_R2 = zscore(All_Plot_R2);
                        else
                            All_Plot_R2 = All_Plot_R2 - min(min(All_Plot_R2)); %Subtract minimum R2 (otherwise plot looks bad)                    
                        end;
                        maxR2 = max(max(All_Plot_R2));
                        minR2 = min(min(All_Plot_R2));

                        %Generate string
                        if Level == 1
                            PSTR=[ 'Overrall R2 \nn subjects:  ' num2str(length(GroupVP)) '\nn Regressors  ' num2str(dims(2)) '\nMin  ' num2str(round(minR2*100)/100) '\nMax  ' num2str(round(maxR2*100)/100) ];
                        else
                            PSTR=[ 'Overrall R2 \nn subjects:  ' num2str(Info.TotalTrials) '\nn Regressors  ' num2str(dims(2)) '\nMin  ' num2str(round(minR2*100)/100) '\nMax  ' num2str(round(maxR2*100)/100) ];
                        end;
                        str = sprintf(PSTR);  
                        text(0.1,0.3,str,'FontSize',12);
                        set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');
                        
                        %Plot R2 topographies
                        
                        for pc = 4 : length(DataTime)+3
                            hand(HC)=subplot(Zeile,Spalte,pc); hold(gca, 'on');
                            set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);
                            topoplot(All_Plot_R2(:,pc-3), ChanLoc, 'maplimits', [-maxR2 maxR2]);colormap(Fcmap); 
                            title([num2str(Ptimes(pc-3)) ' ms']);
                            AllPlots(length(AllPlots(1).Handle)+1).Type    = 'R2';
                            AllPlots(length(AllPlots(1).Handle)+1).Title   = ['R2 ' num2str(Ptimes(pc-3)) ' ms'];
                            AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
                            pause(0.00001);
                        end;
                        %Determine electrode for R2 plot
                        if strcmp(PlotElec, 'best')
                            [x,y]=find(All_Plot_R2 == maxR2);
                            UseElectrode = x;               
                        else
                            UseElectrode = find(strcmpi({ChanLoc.labels},PlotElec));
                        end;

                        %Plot R2 time course
                        hand(HC)=subplot(Zeile,Spalte,[pc+2 pc+6]); hold(gca, 'on');
                        titel=['R2 '];

                        if Ninclu ~= 1
                            This_t = mean(squeeze(R2_values(GroupVP,UseElectrode,1,:)),1);
                            This_SD = ConfFactor*(std(squeeze(R2_values(GroupVP,UseElectrode,1,:)))./sqrt(length(GroupVP)-1));
                            AGF_shadedErrorBar(RegressionTime, This_t, This_SD, colorsSingle, transp);
                        else
                            if length(size(R2_values)) == 4
                                This_t = squeeze(R2_values(1,UseElectrode,1,:));
                            else
                                This_t = squeeze(R2_values(UseElectrode,1,:));
                            end;
                            This_SD = zeros(length(This_t),1);
                            plot(RegressionTime, This_t, colorsSingle);
                        end;

                        if Confidence
                            title([titel ' shade = ' num2str(Confidence*100) '% CI @' ChanLoc(UseElectrode).labels]);
                        else
                            title([titel ' shade = sem']);
                        end;
                        set(gca, 'FontSize', 9); 
                        axe = gca;
                        axe.XLim = [min(RegressionTime) max(RegressionTime)];
                        AllPlots(length(AllPlots(1).Handle)+1).Type    = 'R2';
                        AllPlots(length(AllPlots(1).Handle)+1).Title   = ['R2 overall time course'];
                        AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
                   end;

                %%%%%%%%%%%%plot single regressor effects%%%%%%%%%%%%%
                else 
                    startPOS = (c-1)*Spalte;
                    DatPlot = squeeze(nanmean(regress_values(GroupVP,:,c-1,DataTime),1));
                    if Ninclu ~= 1
                        [~,DatPlotP]=ttest(squeeze(regress_values(GroupVP,:,c-1,DataTime)));
                    else
                        DatPlotP = squeeze(p_values(:,c-1,DataTime));
                    end;

                    UseLim = [-max(abs([min(min(DatPlot)) max(max(DatPlot))])) max(abs([min(min(DatPlot)) max(max(DatPlot))]))];
                    if MinMaplimits
                        if UseLim(1)>-MinMaplimits
                            UseLim(1)=-MinMaplimits;
                        end;
                        if UseLim(2)<MinMaplimits
                            UseLim(2)=MinMaplimits;
                        end;
                    end;
                    UseLim = round(UseLim*100)./100;
                    [BestElectrode,BestTime]=find([abs(DatPlot)]==max(max(abs(DatPlot(:,:))))); % maximum within the plotted time-frames
                    [BestElectrode_PosAllTime,BestTime_PosAllTime]=find([squeeze(nanmean(regress_values(GroupVP,:,c-1,AnyDataTime),1))]==max(max(squeeze(nanmean(regress_values(GroupVP,:,c-1,AnyDataTime),1))))); % maximum across the whole time
                    [BestElectrode_NegAllTime,BestTime_NegAllTime]=find([squeeze(nanmean(regress_values(GroupVP,:,c-1,AnyDataTime),1))]==min(min(squeeze(nanmean(regress_values(GroupVP,:,c-1,AnyDataTime),1))))); % maximum across the whole time


                    SaveBest(c-1)=BestElectrode;
                    hand(HC)=subplot(Zeile,Spalte,[startPOS+1 startPOS+2]); hold(gca, 'on');HC=HC+1;
                    if ~strcmp(Info.RegLables{UPlotReg(c-1),2},'numeric')
                        try PSTR=[ 'Regressor: ' Info.RegNames{UPlotReg(c-1)}...
                            '\n' num2str(Info.RegValues{UPlotReg(c-1)}(1)) ' = ' num2str(Info.RegLables{UPlotReg(c-1),1})...
                            '\n' num2str(Info.RegValues{UPlotReg(c-1)}(end)) ' = ' num2str(Info.RegLables{UPlotReg(c-1),2}) '\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))];
%                         catch PSTR=[ 'Regressor: ' Info.RegNames{UPlotReg(c-1)}...
%                             '\n' num2str(Info.RegValues{UPlotReg(c-1)}(1)) ' = ' num2str(Info.RegLables{UPlotReg(c-1),1})...
%                             '\n' num2str(Info.RegValues{UPlotReg(c-1),2}) ' = ' num2str(Info.RegLables{UPlotReg(c-1),2}) '\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))];
                        end;
                    elseif ~strcmp(Info.RegNames{UPlotReg(c-1)},'EEG')
                        try PSTR=[ 'Regressor: ' Info.RegNames{UPlotReg(c-1)}...
                                '\nminsplit: ' num2str(round(Info.RegValues{UPlotReg(c-1)}(1)*100)/100)...
                                '\nmaxsplit: ' num2str(round(Info.RegValues{UPlotReg(c-1)}(end)*100)/100)...
                                '\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))];
                        catch PSTR=[ 'Regressor: ' Info.RegNames{UPlotReg(c-1)}...
                                '\nminsplit: ' num2str(round(Info.RegValues{UPlotReg(c-1)}(1)*100)/100)...
                                '\nmaxsplit: ' num2str(round(Info.RegValues{UPlotReg(c-1)}(end)*100)/100)...
                                '\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))];
                        end;
                    else
                        PSTR='EEG predictor';
                    end;
                    if AllTime == 1
                        PSTR = [PSTR '\nMax: ' ChanLoc(BestElectrode_PosAllTime).labels ' ' num2str(RegressionTime(AnyDataTime(BestTime_PosAllTime)))...
                            'ms\nMin: ' ChanLoc(BestElectrode_NegAllTime).labels ' ' num2str(RegressionTime(AnyDataTime(BestTime_NegAllTime))) 'ms'];
                    else

                    end;

                    %If set, mask all above threshold p-values in plots in white
                    if MaskPval<1
                        DatPlot(DatPlotP>MaskPval)=0;
                    end;

                    PSTR=strrep(PSTR, '_', ' ');
                    str = sprintf(PSTR);  
                    text(0.1,0.3,str,'FontSize',12);
                    set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');
                    %%%%%%%%%%%%%PLOT TOPOs%%%%%%%%%%%%
                    for pc = 4 : length(DataTime)+3
                        hand(HC)=subplot(Zeile,Spalte,startPOS+pc); hold(gca, 'on');
                        topoplot(DatPlot(:,pc-3), ChanLoc, 'maplimits', UseLim);colormap(Fcmap);        
                        title([num2str(Ptimes(pc-3)) ' ms']);
                        AllPlots(length(AllPlots(1).Handle)+1).Type    = 'Topo';
                        AllPlots(length(AllPlots(1).Handle)+1).Title   = ['Topo ' Info.RegNames{UPlotReg(c-1)} ' ' num2str(Ptimes(pc-3)) ' ms'];
                        AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1; pause(0.00001);
                    end;

                    CurrPosPlot = pc+1;
                    if strcmp(PlotElec, 'best')
                        UseElectrode = BestElectrode;               
                    else
                        UseElectrode = find(strcmpi({ChanLoc.labels},PlotElec));
                    end;

                    if PlotRegress
                        %Plot regressor time-course
                        hand(HC)=subplot(Zeile,Spalte,[startPOS+CurrPosPlot+1 startPOS+CurrPosPlot+2]); hold(gca, 'on');
                        CurrPosPlot = CurrPosPlot + 3;
                        if Ninclu ~= 1
                            This_t = mean(squeeze(regress_values(GroupVP,UseElectrode,c-1,:)),1);
                            This_SD = ConfFactor*(std(squeeze(regress_values(GroupVP,UseElectrode,c-1,:)))./sqrt(length(GroupVP)-1));
                            AGF_shadedErrorBar(RegressionTime, This_t, This_SD, colorsSingle, transp);
                        else
                            This_t = squeeze(regress_values(Ninclu,UseElectrode,c-1,:));
                            This_SD = zeros(length(This_t),1);
                            plot(RegressionTime, This_t, colorsSingle);
                        end;

                        set(gca,'YDir','reverse', 'FontSize', 9); 
                        axe = gca; axe.XLim = [min(RegressionTime) max(RegressionTime)];
                        title(['Regression weight@' ChanLoc(UseElectrode).labels])
                        AllPlots(length(AllPlots(1).Handle)+1).Type    = 'Regressor';
                        AllPlots(length(AllPlots(1).Handle)+1).Title   = ['Regressor ' Info.RegNames{UPlotReg(c-1)} ' time course'];
                        AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
                    end;

                    if PlotPval
                        %Plot p-value time-course
                        hand(HC)=subplot(Zeile,Spalte,[startPOS+CurrPosPlot+1 startPOS+CurrPosPlot+2]); hold(gca, 'on');
                        CurrPosPlot = CurrPosPlot + 3;
                        if Ninclu ~= 1
                            [~,PP]=ttest(squeeze(regress_values(GroupVP,UseElectrode,c-1,:)));
                        else
                            PP = squeeze(p_values(UseElectrode,c-1,:));
                        end;
                        plot(RegressionTime, PP, colorsPval)
                        axe = gca; axe.XLim = [min(RegressionTime) max(RegressionTime)];
                        set(axe,'yscale','log');
                        title(['P-values @' ChanLoc(UseElectrode).labels]);
                        AllPlots(length(AllPlots(1).Handle)+1).Type    = 'Pval';
                        AllPlots(length(AllPlots(1).Handle)+1).Title   = ['Pval ' Info.RegNames{UPlotReg(c-1)} ' time course'];
                        AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
                    end;

                    %Plot EEG over regressor
                    if EEGplot
                        hand(HC)=subplot(Zeile,Spalte,[startPOS+CurrPosPlot+1 startPOS+CurrPosPlot+4]); hold(gca, 'on');
                        CurrPosPlot = CurrPosPlot + 5;
                        titel='';
                        
                        
                        lejenda={};
                        try EP = length(Info.RegNumbers{UPlotReg(c-1)});
                        catch EP = length(VP.VP_MODEL.Model_1_Resp_Alltrials_Overall.EEG_Values{UPlotReg(c-1)});
                        end;
                        size(EEG_values)
                        for BinCount = 1:EP
%                             if ~isnan(EEG_values(GroupVP,c-1,BinCount,UseElectrode,:,:)
                                if Level == 1 && Ninclu ~= 1
                                    if Ninclu~=1
                                        This_EEG(BinCount,:)         = nanmean(squeeze(EEG_values(GroupVP,c-1,BinCount,UseElectrode,:,:)),1);
                                        This_EEG_SD(BinCount,:)      = ConfFactor*(nanstd(squeeze(EEG_values(GroupVP,c-1,BinCount,UseElectrode,:,:)))./sqrt(length(GroupVP)-1));
                                    else
                                        This_EEG(BinCount,:)         = squeeze(EEG_values(GroupVP,c-1,BinCount,UseElectrode,:,:));
                                        This_EEG_SD(BinCount,:)      = (ConfFactor.*squeeze(EEG_SD_values(c-1,BinCount,UseElectrode,:,:)))./sqrt(Info.RegNumbers{UPlotReg(c-1)}(BinCount)-1);
                                    end;
                                else
                                    This_EEG(BinCount,:)         = squeeze(EEG_values(c-1,BinCount,UseElectrode,:));
                                    This_EEG_SD(BinCount,:)      = (ConfFactor.*squeeze(EEG_SD_values(c-1,BinCount,UseElectrode,:)))./sqrt(Info.RegNumbers{UPlotReg(c-1)}(BinCount)-1);
                                end;
                                AGF_shadedErrorBar(RegressionTime, This_EEG(BinCount,:), This_EEG_SD(BinCount,:), colors{BinCount}, transp);
                                if Info.RegValues{UPlotReg(c-1)}(BinCount) > 10
                                    lejenda(BinCount)       = {[' m = ' num2str(round(Info.RegValues{UPlotReg(c-1)}(BinCount)))]};
                                else
                                    lejenda(BinCount)       = {[' m = ' num2str(round(100*Info.RegValues{UPlotReg(c-1)}(BinCount))/100)]};
                                end;
%                                 if size(EEG_values,3)>2
%                                     if ~isnan(EEG_values(1, c-1, 3, 1, 1))
%                                         
%                                     else
%                                         lejenda(1)       = {num2str(Info.RegLables{UPlotReg(c-1),1})}; lejenda(2)       = {num2str(Info.RegLables{UPlotReg(c-1),2})};
%                                     end;
%                                 else
%                                     lejenda(1)       = {num2str(Info.RegLables{UPlotReg(c-1),1})}; lejenda(2)       = {num2str(Info.RegLables{UPlotReg(c-1),2})};
%                                 end;

%                             end;
                        end;
                        legend(lejenda,'Location','southeast');
                        if Confidence
                            title([titel ' shade = ' num2str(Confidence*100) '% CI @' ChanLoc(UseElectrode).labels]);
                        else
                            title([titel ' shade = sem']);
                        end;
                        set(gca,'YDir','reverse', 'FontSize', 9); 
                        axe = gca;
                        axe.XLim = [min(RegressionTime) max(RegressionTime)];
                        if isfield(Info,'EEG_Factor_Level')
                            axe.YLabel.String=[BinPlotTitleStr ' Factor ' Info.EEG_Factor_Level ' '];
                        else
                            axe.YLabel.String=[BinPlotTitleStr ' '];
                        end;
                        AllPlots(length(AllPlots(1).Handle)+1).Type    = 'EEG';
                        AllPlots(length(AllPlots(1).Handle)+1).Title   = ['EEG ' Info.RegNames{UPlotReg(c-1)} ' time course'];
                        AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
                    end;    

                    %Plot ERP image for regressor (new feature!)
                    if ERP_Im_plot && exist('dims3')
                        hand(HC)=subplot(Zeile,Spalte,[startPOS+CurrPosPlot+1 startPOS+CurrPosPlot+4]); hold(gca, 'on');
                        CurrPosPlot = CurrPosPlot + 5;
                        imagesc(RegressionTime,[1:size(ERPImage_values,3)],squeeze(ERPImage_values(UPlotReg(c-1), UseElectrode,:,:)));colormap(Fcmap);
                        set(gca,'YDir','normal', 'ytick',[10 90],'YTickLabel', {'min' 'max', 'xcolor','w','ycolor','w','box','off'})
                        axe = gca; 
                        axe.XLim = [min(RegressionTime) max(RegressionTime)];
                        axe.YLim = [1 100];
                        title('ERP image');
                        AllPlots(length(AllPlots(1).Handle)+1).Type    = 'ERPimage';
                        AllPlots(length(AllPlots(1).Handle)+1).Title   = ['ERPimage ' Info.RegNames{UPlotReg(c-1)}];
                        AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
                    end;
                end;
            end;

            %Add top title
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 1,ti,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)
            saveas (hFig,[outfold ti ],printstring);

            %Save individual images (may be better for publications etc.
            if IndiviPlot
                %Set order for Items to plot
                IPlots = [];
                for c = 1 : length(IndivTypes)
                    IPlots = [IPlots find(strcmpi({AllPlots.Type},IndivTypes{c}))]
                end;

                %Test if 'Singles' subfolder exists, otherwise create
                if ~exist([outfold  '/singles/'])
                    mkdir([outfold  '/singles/']);
                end;

                for c = 1 : length(IPlots)
                    pfig = figure;
                    hax_new = copyobj(hand(AllPlots(1).Handle(IPlots(c))), pfig);
                    set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
                    colormap(Fcmap);
                    saveas (pfig,[outfold  '/singles/' AllPlots(IPlots(c)).Title], printstring);
                    close(pfig);
                end;
            end;
        end;
    end;
end;
return;
%%
%settings
%.SPBands         = {[1 4]; [5 8]; [9 12]; [13 30]};
%.SPNames         = {'delta' 'theta' 'alpha' 'beta'};




































%output very simple statistics for 2nd level
TimePoint = 310
Ele = 'FCz'
DataTime = find(ismember(RegressionTime,TimePoint));
UseElectrode = find(strcmpi({ChanLoc.labels},Ele));
RegNumber = 1

P_Val = squeeze(p_values(UseElectrode,RegNumber,DataTime))
T_Val = squeeze(regress_values(1,UseElectrode,RegNumber,DataTime))


%%
%Simple plot for SVM
close all;
figure
plot(Info.Return_Timewindow_ms,squeeze(mean(regress_values(:,1,:,:))));hold
plot(Info.Return_Timewindow_ms,squeeze(mean(regress_values(:,2,:,:))))
plot(Info.Return_Timewindow_ms,squeeze(mean(regress_values(:,3,:,:))))

