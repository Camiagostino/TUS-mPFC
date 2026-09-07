function [] = DDM_plot(D)
global VP
%%
if ~isempty(findall(0,'Type','Figure'))
    h=clf;
else
    h=figure;
end
set(h, 'Position', [1919 209 1231 677])%
set(h,'PaperPositionMode','Auto')

subplot(3,2,1)
histogram(nanrem(reshape(VP.RT(VP.error==1),[],1)),VP.EqBin_Edge_Group_cor); hold
histogram(D.rt(D.error==1),VP.EqBin_Edge_Group_cor);
legend({'VP' 'DDM'})
title('correct trials')

subplot(3,2,2)
histogram(nanrem(reshape(VP.RT(VP.error==2),[],1)),VP.EqBin_Edge_Group_err); hold
histogram(D.rt(D.error==2),VP.EqBin_Edge_Group_err);
legend({'VP' 'DDM'})
title('error trials')

subplot(3,2,3)
histogram(nanrem(reshape(VP.RT(VP.error==1 & VP.congr==2),[],1)),VP.EqBin_Edge_Group_cor); hold
histogram(D.rt(D.error==1 & VP.congr==2),VP.EqBin_Edge_Group_cor);
legend({'VP' 'DDM'})
title('correct incongruent trials')

subplot(3,2,4)
histogram(nanrem(reshape(VP.RT(VP.error==1 & VP.congr==1),[],1)),VP.EqBin_Edge_Group_cor); hold
histogram(D.rt(D.error==1 & VP.congr==1),VP.EqBin_Edge_Group_cor);
legend({'VP' 'DDM'})
title('correct congruent trials')

subplot(3,2,5)
PSTR = ['Err rate DDM: ' num2str(round2(sum(D.error==2)/length(D.error),0.01)) '\n'];
PSTR = [PSTR 'Err rate VP: ' num2str(round2(sum(VP.error==2)/length(VP.error),0.01)) '\n'];
PSTR = [PSTR 'perc err in inc: ' num2str(round2(sum(D.error==2 & VP.congr==2)/sum(D.error==2),0.01)) '\n'];
PSTR = [PSTR 'perc err in inc: ' num2str(round2(sum(VP.error==2 & VP.congr==2)/sum(VP.error==2),0.01)) '\n'];
PSTR = [PSTR 'posterior: ' num2str(round2(D.posterior,0.1)) '\n\n'];

FN=fieldnames(D.params);
for c = 1 : length(FN)
    PSTR = [PSTR FN{c} ': ' num2str(round2(D.params.(FN{c}),0.01)) '\n'];
end

str = sprintf(PSTR);  
text(0,0.4,str,'FontSize',8);
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');

if isfield(D,'RecoveryParams') %this is the model recovery mode, plot true parameters
    PSTR = 'True model parameters:\n';
    for c = 1 : length(D.RecoveryParams)
        PSTR = [PSTR FN{c} ': ' num2str(round2(D.RecoveryParams(c),0.01)) '\n'];
    end
    str = sprintf(PSTR);  
    text(0.3,0.35,str,'FontSize',8);
    set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');
end

subplot(3,2,6)
PSTR = ['Correct:\nmean RT inc DDM: ' num2str(round2(mean(D.rt(VP.congr==2 & D.error==1)),0.1)) '\n'];
PSTR = [PSTR 'mean RT inc VP : ' num2str(round2(nanmean(VP.RT(VP.congr==2 & VP.error==1)),0.1)) '\n'];
PSTR = [PSTR 'mean RT con DDM: ' num2str(round2(mean(D.rt(VP.congr==1 & D.error==1)),0.1)) '\n'];
PSTR = [PSTR 'mean RT con VP : ' num2str(round2(nanmean(VP.RT(VP.congr==1 & VP.error==1)),0.1)) '\n\nError:\n'];

PSTR = [PSTR 'mean RT inc DDM: ' num2str(round2(mean(D.rt(VP.congr==2 & D.error==2)),0.1)) '\n'];
PSTR = [PSTR 'mean RT inc VP : ' num2str(round2(nanmean(VP.RT(VP.congr==2 & VP.error==2)),0.1)) '\n'];
PSTR = [PSTR 'mean RT con DDM: ' num2str(round2(mean(D.rt(VP.congr==1 & D.error==2)),0.1)) '\n'];
PSTR = [PSTR 'mean RT con VP : ' num2str(round2(nanmean(VP.RT(VP.congr==1 & VP.error==2)),0.1)) '\n\n'];

str = sprintf(PSTR);  
text(0,0.4,str,'FontSize',8);
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');

%%
return