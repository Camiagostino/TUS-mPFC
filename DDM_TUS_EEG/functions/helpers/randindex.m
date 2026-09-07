function x = randindex(n,s)
%returns n random indices between 1 and s. 
r = randperm(s);
x = r(1:n);