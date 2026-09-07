function [e,N]=agf_quantile(data,s)
% keyboard
%%

if ~isfield(s,'space');     s.space = [.25 .5 .75];     end
if ~isfield(s,'extremes');  s.extremes = 0;             end
%%
%get quantile edges
q = quantile(data, s.space);
e = [s.extremes(1) q(2:end)];

if ~isempty(s.extremes) && numel(s.extremes) > 1
    e(end) = s.extremes(2);
end

for j=2:length(q)
    N(j-1) = sum(data<=e(j)) - sum(data<=e(j-1));
end
return;

