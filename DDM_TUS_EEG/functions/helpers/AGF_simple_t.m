%%%%JUST PERFORM VERY QUICK AND SIMPLE t-Test over open ALLEEG documents
%Enter ALLEEG numbers for sets
n1=1;        %EEG dataset
n2=2;        %EEG dataset
paired=1;    %Paired (1) ttest oder independet (alles andere) sample ttest
tw=350:450;  %Timewindow in milliseconds
chan=22;     %22 = FCz
mt='max'       %min = minimum, max = maximum, alles andere = mean

dw=[tw(1)/(1000/ALLEEG(n1).srate):tw(end)/(1000/ALLEEG(n1).srate)]+abs(ALLEEG(n1).xmin*ALLEEG(n1).srate);
if paired==1
    if mt == 'max'
        [h,p,ci,stats]=ttest(squeeze(max(ALLEEG(n1).data(chan,dw,:))),squeeze(max(ALLEEG(n2).data(chan,dw,:))));
    elseif mt == 'min'
        [h,p,ci,stats]=ttest(squeeze(min(ALLEEG(n1).data(chan,dw,:))),squeeze(min(ALLEEG(n2).data(chan,dw,:))));
    else
        [h,p,ci,stats]=ttest(squeeze(mean(ALLEEG(n1).data(chan,dw,:))),squeeze(mean(ALLEEG(n2).data(chan,dw,:))));
    end;
else
    if mt == 'max'
        [h,p,ci,stats]=ttest2(squeeze(max(ALLEEG(n1).data(chan,dw,:))),squeeze(max(ALLEEG(n2).data(chan,dw,:))));
    elseif mt == 'min'
        [h,p,ci,stats]=ttest2(squeeze(min(ALLEEG(n1).data(chan,dw,:))),squeeze(min(ALLEEG(n2).data(chan,dw,:))));
    else
        [h,p,ci,stats]=ttest2(squeeze(mean(ALLEEG(n1).data(chan,dw,:))),squeeze(mean(ALLEEG(n2).data(chan,dw,:))));
    end;
end;
if mt == 'max'
    mean1=mean(squeeze(max(ALLEEG(n1).data(chan,dw,:))));
    sem1=std(squeeze(max(ALLEEG(n1).data(chan,dw,:))))/sqrt(size(ALLEEG(n1).data,3)-1);
    max2=max(squeeze(max(ALLEEG(n2).data(chan,dw,:))));
    sem2=std(squeeze(max(ALLEEG(n2).data(chan,dw,:))))/sqrt(size(ALLEEG(n2).data,3)-1);
elseif mt == 'min'
    mean1=mean(squeeze(min(ALLEEG(n1).data(chan,dw,:))));
    sem1=std(squeeze(min(ALLEEG(n1).data(chan,dw,:))))/sqrt(size(ALLEEG(n1).data,3)-1);
    mean2=mean(squeeze(min(ALLEEG(n2).data(chan,dw,:))));
    sem2=std(squeeze(min(ALLEEG(n2).data(chan,dw,:))))/sqrt(size(ALLEEG(n2).data,3)-1);
else
    mean1=mean(squeeze(mean(ALLEEG(n1).data(chan,dw,:))));
    sem1=std(squeeze(mean(ALLEEG(n1).data(chan,dw,:))))/sqrt(size(ALLEEG(n1).data,3)-1);
    mean2=mean(squeeze(mean(ALLEEG(n2).data(chan,dw,:))));
    sem2=std(squeeze(mean(ALLEEG(n2).data(chan,dw,:))))/sqrt(size(ALLEEG(n2).data,3)-1);
end;

%output
sprintf(['Results of test:'...
        '\n p:            ' num2str(p)...
        '\n t:            ' num2str(stats.tstat)...
        '\n ci:           ' num2str(ci(1)) ' to ' num2str(ci(2))...
        '\n Name set 1:   ' ALLEEG(n1).setname...
        '\n Name set 2:   ' ALLEEG(n2).setname...
        '\n mean set 1:   ' num2str(mean1) ' and SE ' num2str(sem1) ' over ' num2str(size(ALLEEG(n1).data,3)) ' subjects'...
        '\n mean set 2:   ' num2str(mean2) ' and SE ' num2str(sem2) ' over ' num2str(size(ALLEEG(n1).data,3)) ' subjects'...
        '\n time-window:  ' num2str(tw(1)) ' to ' num2str(tw(end)) 'ms'...
        '\n electrode:    ' ALLEEG(n1).chanlocs(chan).labels])

    
%%
%create quick differencewave of dataset 1 - 2 in ALLEEG(3)

ALLEEG(3).data=ALLEEG(2).data-ALLEEG(1).data;
EEG.data=ALLEEG(2).data-ALLEEG(1).data;
EEG = eeg_checkset( EEG );
eeglab redraw
%%
figure;pop_topoplot(EEG,1, 214,'Differencewave unfavorable - favorable',[1 1] ,0,'electrodes','on');
figure;pop_topoplot(EEG,1, 260,'Differencewave unfavorable - favorable',[1 1] ,0,'electrodes','on');
figure;pop_topoplot(EEG,1, 300,'Differencewave unfavorable - favorable',[1 1] ,0,'electrodes','on');
figure;pop_topoplot(EEG,1, 450,'Differencewave unfavorable - favorable',[1 1] ,0,'electrodes','on');
figure;pop_topoplot(EEG,1, 500,'Differencewave unfavorable - favorable',[1 1] ,0,'electrodes','on');

