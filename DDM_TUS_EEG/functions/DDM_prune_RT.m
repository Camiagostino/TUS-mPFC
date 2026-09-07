function [PrunedRT] = DDM_prune_RT(rt, q, hardcut)
%%
%this function prunes RT of inidvidual subjects by discarding (replace with
%nan) the highest and lower q quantile of the data. If rt is a matrix, it
%is assumed that columns are subjects. hardcut can be an additional cutoff
%point that is applied first to the data.
%%
s = size(rt);
if ~isempty(hardcut)
    rep_hard1 = sum(sum(rt<hardcut(1)));
    rep_hard2 = sum(sum(rt>hardcut(2)));
    rt(rt<hardcut(1))=nan;
    rt(rt>hardcut(2))=nan;
end

if s(2)>1 %multiple subjects
    for c = 1 : s(2)
       [qVP]=lequantile(nanrem(rt(:,c)),'qq', [0 q/2 1-q/2 inf], 'eQuant', [0 inf]); 
       rt(rt(:,c)<qVP(2),c)=nan;
       rt(rt(:,c)>qVP(3),c)=nan;
    end
else
    [qVP]=lequantile(nanrem(rt),'qq', [0 q/2 1-q/2 inf], 'eQuant', [0 inf]); 
    rt(rt<qVP(2))=nan;
    rt(rt>qVP(3))=nan;
end
PrunedRT=rt;
display(['n removed by low hardcut: ' num2str(rep_hard1) ' and higher: ' num2str(rep_hard2)])
return;











