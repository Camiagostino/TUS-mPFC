function [O] = DDM_llf(rt,Error)
%returns the likelihood value of DDM distributions and observations.
global VP
%%
persistent ModelN FitMethod VP_n wP
if isempty(ModelN)
    ModelN = nan(length(VP.n_Group_Quantile_cor),2);         %this is the number of observations in the respective bin
    if ~VP.FitGroup
        if VP.fitmode == 1
           FitMethod = {'EqBin_Edge_Group_cor' 'EqBin_Edge_Group_err'};
           VP_n(:,1) = VP.n_EqBin_Edge_Group_cor; VP_n(:,2) = VP.n_EqBin_Edge_Group_err;
        elseif VP.fitmode == 2
           FitMethod = {'Quant_Edge_Group_cor' 'Quant_Edge_Group_err'};
           VP_n(:,1) = VP.n_Group_Quantile_cor; VP_n(:,2) = VP.n_Group_Quantile_err;
        elseif VP.fitmode == 3
            FitMethod = {'Quant_Edge_VP_cor' 'Quant_Edge_VP_err'};
            VP_n(:,1) = VP.n_VP_Quant_Edge_VP_cor; VP_n(:,2) = VP.n_VP_Quant_Edge_VP_err;
        end
    else
        if VP.fitmode == 1
           FitMethod = {'EqBin_Edge_Group_cor' 'EqBin_Edge_Group_err'};
           VP_n(:,1) = VP.n_EqBin_Edge_Group_cor; VP_n(:,2) = VP.n_EqBin_Edge_Group_err;
        elseif VP.fitmode == 2
           FitMethod = {'Quant_Edge_Group_cor' 'Quant_Edge_Group_err'};
           VP_n(:,1) = VP.n_Group_Quantile_cor; VP_n(:,2) = VP.n_Group_Quantile_err;
        else
            error('When fittig to group data, you cannot use individual quantiles. Set VP.fitmode to 1 or 2.'); return
        end
    end
    VP_n(VP_n==0)=VP.CutLog;
    wP = (VP_n./sum(sum(VP_n))).*numel(VP_n);
    if VP.MixturePercent
        %% determine LL for quantiles in the mixture model
        LikperSample = 1/(length(VP.MixtureBounds(1):VP.MixtureBounds(2))*2);
        FramesCor = VP.(FitMethod{1}); FramesCor(FramesCor==inf) = VP.MixtureBounds(2);
        FramesErr = VP.(FitMethod{2}); FramesErr(FramesErr==inf) = VP.MixtureBounds(2);
        for c = length(VP.(FitMethod{1})):-1:2
            SamplesPerQuantile(c-1,1) = FramesCor(c)-FramesCor(c-1);
            SamplesPerQuantile(c-1,2) = FramesErr(c)-FramesErr(c-1);
        end
        TC = sum(VP_n);
        VP.MixtureLL = [SamplesPerQuantile(:,1).*(LikperSample*(VP.MixturePercent/100*(TC(1)/sum(TC)))) SamplesPerQuantile(:,2).*(LikperSample*(VP.MixturePercent/100*(TC(2)/sum(TC))))];
        if TC(2)==0;VP.MixtureLL(:,2) = VP.CutLog; end %no errors, likelihood would be zero
    end
end
%%
 %keyboard
%%
for c = 1 : length(ModelN)
    ModelN(c,1) = numel(rt(Error==1 & rt>=VP.(FitMethod{1})(c) & rt<VP.(FitMethod{1})(c+1)));
    ModelN(c,2) = numel(rt(Error==2 & rt>=VP.(FitMethod{2})(c) & rt<VP.(FitMethod{2})(c+1)));
end
VP.x2df = prod(size(ModelN)-1);
% [ModelN VP_n]
%%
if VP.FitMethod==2 %calculate log likelihood of observed (participant) bin frequency given data from model
    if VP.FitDirection == 2
        P = ModelN ./ repmat(sum(ModelN(:)),size(ModelN,1),2);
    else
        P = VP_n ./ repmat(sum(VP_n(:)),size(VP_n,1),2);
    end
    if VP.MixturePercent
        P(isnan(P))=0;
        P=P.*(1-VP.MixturePercent/100)+VP.MixtureLL;
    else %no mixture model, just cut P at a low threshold
        P(P<VP.CutLog)=VP.CutLog;
        P(isnan(P))=VP.CutLog; %can happen if model / participant has no errors at all
    end
%     if VP.MCMC ~= 0 %for mcmc sampling, we should scale the LL to allow proper comaprison of extremely small LL differences
%         if VP.MCMC > 0 
%             O.chi_2 = O.chi_2 + VP.MCMC;
%         else
%             O.chi_2 = O.chi_2 + VP.MCMC;
%         end
%         if O.chi_2 >= 0 
%             O.chi_2 = 0-VP.CutLog;
%         end
%     end
    %first multiply to number of obs, then take log
%     VP.LLbyCor = sum(log(P(:,1).*numel(P)))*VP.Maximize;
%     VP.LLbyErr = sum(log(P(:,2).*numel(P)))*VP.Maximize;
    VP.LL_quant = sum(log(P(:).*wP(:)*numel(P)))*VP.Maximize;
    VP.LL_tot = sum(log(P(:)).*(VP_n(:)./VP.Inflation))*VP.Maximize;
                
    %take log first, then multiply with observations
%     VP.LLbyCor = sum(log(P(:,1)).*(VP_n(:)./VP.Inflation))*VP.Maximize;
%     VP.LLbyErr = sum(log(P(:,2)).*(VP_n(:)./VP.Inflation))*VP.Maximize;
%     VP.LL_tot = sum(log(P(:)).*(VP_n(:)./VP.Inflation))*VP.Maximize;

%     %perfect fit would have LL:
%     P2 = repmat(0.05,10,2);
%     VP.LLbyCor
%     VP.LLbyErr
%     VP.LL_tot
%     
%     PerLLCor = sum(log(P2(:,1).*VP_n(:,1)))*VP.Maximize
%     PerLLErr = sum(log(P2(:,2).*VP_n(:,2)))*VP.Maximize
%     sum(log(P2(:).*VP_n(:)))*VP.Maximize
%     ert
else
    ModelN = ModelN / sum(ModelN(:)); %get frequencies of model
    VP_n = VP_n / sum(VP_n(:)); %get frequencies of VP
    O.log = sum(sum((ModelN-VP_n).^2./VP_n))*VP.Maximize*-1;
end
O.LL_tot = VP.LL_tot;
return
%%
ModelN = ModelN / sum(ModelN(:)); %get frequencies of model
VP_n = VP_n / sum(VP_n(:)); %get frequencies of VP
O.chi_2 = sum(sum((ModelN-VP_n).^2./VP_n))*VP.Maximize*-1;
% ModelN(ModelN<VP.CutLog)=VP.CutLog;
% O.chi_2 = sum(sum((VP_n-ModelN).^2./ModelN));

if VP.fit_ratio_congruent_error>= 1 %the ratio of (only) congruent errors should also matter to the model fit
    O.chi_2 = O.chi_2 + (((sum(Error==2 & VP.congr==1) - sum(VP.error==2 & VP.congr==1))^2)/sum(VP.error==2 & VP.congr==1))*VP.weight_congruent_errors;
    if VP.fit_ratio_congruent_error== 2 %also add number of incongruent errors to fit
        O.chi_2 = O.chi_2 + ((sum(Error==2 & VP.congr==2) - sum(VP.error==2 & VP.congr==2))^2)/sum(VP.error==2 & VP.congr==2);
    end
end
% O.chi_2 = O.chi_2 / 100;
% 
% if isfield(VP, 'symbmath') && VP.symbmath && VP.useExp == 1 %symbolic math avoids precision problem, but is cosiderably slower
%     O.L_value = vpa(exp(-sym(O.chi_2)))*VP.Maximize;
% elseif VP.useExp == 1
%     O.L_value = exp(-O.chi_2)*VP.Maximize;
% else
%     O.L_value = -O.chi_2*VP.Maximize;
% end
% if L_value==0
%     L_value=-inf
% end
return