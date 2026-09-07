function [I,rA,rB] = AGF_match_vectors(A,B)
%function that returns the closest possible match between two vectors. It discards entries until no significant difference is observed between both vectors or the means are almost equal.
%Settings:
%%
p_val = 0.3; %above 1 = target divergence in percent

%Idea: sort both vectors, always use largest of smaller mean and lowest of larger mean

if mean(A) < mean(B)
    [sA,iA] = sort(A,'descend');
    [sB,iB] = sort(B);
else
    [sA,iA] = sort(A);
    [sB,iB] = sort(B,'descend');
end

for c = min([length(A) length(B)]) : -1 : 1
    [~,p] = ttest(sA(1:c), sB(1:c));
    if p>p_val
        I(:,1) = sort(iA(1:c));
        I(:,2) = sort(iB(1:c));
        rA = A(I(:,1));
        rB = B(I(:,2));
        disp(['Both vectors are not significantly different and reduced to a size of ' num2str(c) ' entries.'])
        disp(['Mean A = ' num2str(mean(rA)) ' and mean B = ' num2str(mean(rB)) ', p for difference = ' num2str(p) '.'])
        break
    elseif c == 2
        disp('Error: both arrays cannot be made as similar as to reach the target criterion...')
        I(:,1) = NaN;
        I(:,2) = NaN;
        rA = NaN;
        rB = NaN;
        break
    end
end
return