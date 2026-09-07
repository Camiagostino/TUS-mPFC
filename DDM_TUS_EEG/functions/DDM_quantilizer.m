function [] = DDM_quantilizer(rt,error)
%sets the required quantile options for a VP
global VP
% keyboard
%%
VP.n_VP_EqBin_Edge_Group_cor = [];
VP.n_VP_EqBin_Edge_Group_err = [];
VP.n_VP_Quant_Edge_Group_cor = [];
VP.n_VP_Quant_Edge_Group_err = [];
VP.n_VP_Quant_Edge_VP_cor    = [];
VP.n_VP_Quant_Edge_VP_err    = [];

%split the global rt data into n quantiles (to calculate X2 statistics)
Quantile_Spacing = [0:VP.quant_space:1];
Quantile_Number = length(Quantile_Spacing)-1;

QuantileSet.space = Quantile_Spacing;
QuantileSet.extremes = [0 inf];

Group_cor_RT = nanrem(reshape(rt(error==1),[],1));
Group_err_RT = nanrem(reshape(rt(error==2),[],1));

VP_cor_RT = nanrem(reshape(VP.RT(VP.error==1),[],1));
VP_err_RT = nanrem(reshape(VP.RT(VP.error==2),[],1));

VP.NoErr = 0;
if isempty(VP_err_RT); VP.NoErr = 1; end

[VP.Quant_Edge_Group_cor,VP.n_Group_Quantile_cor] = agf_quantile(Group_cor_RT, QuantileSet);
[VP.Quant_Edge_Group_err,VP.n_Group_Quantile_err] = agf_quantile(Group_err_RT, QuantileSet);

[~,VP.EqBin_Edge_Group_cor]     = histcounts(Group_cor_RT,Quantile_Number); VP.EqBin_Edge_Group_cor(end)    = inf;

%Individual VP quantiles and N
VP.Quant_Edge_VP_cor = agf_quantile(VP_cor_RT, QuantileSet); 

%All errror trials
if ~VP.NoErr && numel(VP_err_RT)/VP.Inflation > Quantile_Number
    VP.Quant_Edge_VP_err = agf_quantile(VP_err_RT, QuantileSet); 
    %very particular case: sometimes quantiles can not be evenly space, we have to use correct quantiles then
    if sum(hist(VP.Quant_Edge_VP_err(1:end-1),unique(VP.Quant_Edge_VP_err(1:end-1)))>1)>0 %one value is more than once a quantile edge
        VP.Quant_Edge_VP_err        = VP.Quant_Edge_Group_err; VP.UseGroupErrQuant = 1;
    end
    [~,VP.EqBin_Edge_Group_err]     = histcounts(Group_err_RT,Quantile_Number); VP.EqBin_Edge_Group_cor(end)    = inf;
else %VP has no errors, quantiles are the same as group error quantiles but count will be zero
    VP.UseGroupErrQuant = 1;
    VP.Quant_Edge_VP_err        = VP.Quant_Edge_Group_err;
    VP.EqBin_Edge_Group_err     = VP.EqBin_Edge_Group_err;
end


%VP quantiles based on group edges must be calcualted
for c =  1 : Quantile_Number
    VP.n_VP_EqBin_Edge_Group_cor(c)      = numel(VP.RT(VP.error==1 & VP.RT>=VP.EqBin_Edge_Group_cor(c) & VP.RT<VP.EqBin_Edge_Group_cor(c+1)));
    VP.n_VP_EqBin_Edge_Group_err(c)      = numel(VP.RT(VP.error==2 & VP.RT>=VP.EqBin_Edge_Group_err(c) & VP.RT<VP.EqBin_Edge_Group_err(c+1)));
    VP.n_VP_Quant_Edge_Group_cor(c)      = numel(VP.RT(VP.error==1 & VP.RT>=VP.Quant_Edge_Group_cor(c) & VP.RT<VP.Quant_Edge_Group_cor(c+1)));
    VP.n_VP_Quant_Edge_Group_err(c)      = numel(VP.RT(VP.error==2 & VP.RT>=VP.Quant_Edge_Group_err(c) & VP.RT<VP.Quant_Edge_Group_err(c+1)));
    VP.n_VP_Quant_Edge_VP_cor(c)         = numel(VP.RT(VP.error==1 & VP.RT>=VP.Quant_Edge_VP_cor(c) & VP.RT<VP.Quant_Edge_VP_cor(c+1)));
    VP.n_VP_Quant_Edge_VP_err(c)         = numel(VP.RT(VP.error==2 & VP.RT>=VP.Quant_Edge_VP_err(c) & VP.RT<VP.Quant_Edge_VP_err(c+1)));
end



