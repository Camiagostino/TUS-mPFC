function [par,logL,exF]=QMLE(data,distr, varargin)
%QMLE, quantile maximum likelihood estimation, as expleined in 
%QMLE: fast, robust, and efficient estimation of distribution functions based on quantiles.
%           -Brown and Heathcote, 2003
%
%QMLE: estimating Lognormal, Wald, and Weibull RT distributions with a parameter-dependent lower bound.
%           -Heathcote, Brown, Cousineau, 2004
%see also: 
%A comment on Heathcote, Brown, and Mewhort's QMLE method for response time distributions.
%           -Speckman, Rouder, 2004
%
% [par,logL, exF]=QMLE(data, distr)
%    par= return the parameter of the distribution distr estimated using the
%    QUANTILES in data.  DATA already contains QUANTILE. Use the function
%    "quantile" in matlab or "lequantile" provided with this code.
%    Heathcote et al. propose to use as first and last quantile the extreme
%   of the assumed distribution (for example, 0 and +inf for exp. distr). 
%
%   logL= return the minimum found log-likelihood
%   
%   exF= exit flag. 1=everything went fine. 2=fminsearch did not converge
%   
%   [...]=QMLE(data, distr, 'startPoint',[0 10], 'plotF',1,'N',10)
%   StartPoint= specift the initial "good guess" for the search. If provided
%   it as to be consistent with the number of parameter of the assumed distribution.
%   If not  specified, some classic Reaction Time's good guess will be used. 
%   
%   plotF= plot the parameter estimation while it's been calculated
%
%   N=number of elements in each quantile. If it's all the same, do not
%   specify.
%
%   NOTICE that any custom distribution function can be used, as far as you
%   have a function that describe his pdf. Notice that this function has to
%   end with pdf. (SEE EXAMPLE 2)
%   
%
%   **EXAMPLE 1
%
%   a=normrnd(10,1,10000,1);
%   %we will have a different amount of data in the quantile (like in this example, and we will
%   also specify the extremes of the assumed distribution 
%   [q,~,~,N]=lequantile(a, 'qq',[0.0:0.05:0.8 1], 'eQuant',[-inf inf])
%   [par]=QMLE(q,'norm', 'startPoint',[1 1 ],'plotF',1,'N',N);
%
%
%   **EXAMPLE 2
%   
%   a=exprnd(5, 1000,1); 
%   %this time we will use the equally distributed quantile, so we don't
%   have to specify N
%   [q]=lequantile(a, 'qq', [0:0.2:1], 'eQuant',[0 Inf]);
%   [par]=QMLE(q,'exp', 'startPoint',[200],'plotF',1);
%
%   Created by Valerio Biscione, 11 04 2014

if nargin<2
    distr='norm';
end
p=inputParser;
addParamValue(p, 'startPoint', []); 
addParamValue(p, 'plotF', []);
addParamValue(p, 'N', []); 
parse(p, varargin{:});
startPoint=p.Results.startPoint;
plotF=p.Results.plotF; 
N=p.Results.N;

dpdf=[distr 'pdf']; 

if isempty(N)
    N(1:length(data))=1; 
end 

if isempty(startPoint) %some good starting point for reaction times  distributions
    switch distr        
        case 'norm'
            startPoint=[0.003 0.003];
        case 'reciinvg'
            startPoint=[300 500];
        case 'recinorm'
            startPoint=[0.001 0.0003];
        case 'exgauss'
            startPoint=[200 10 100];
        case 'reciexgauss'
            startPoint=[240 34 78];
        case 'invg'
            startPoint=[300 5000];
        case 'logn'
            startPoint=[4.70 0.5];
    end
end 

if plotF==1
    plotF=@optimplotx;
else
    plotF=[];
end 

[par, logL, exF]= fminsearch(@(x) objFun(x, dpdf,data,N) , startPoint,optimset('PlotFcns',plotF));%'TolX',0.00000001));

end


function QMLE=objFun(x,distr,q,N)
str=[];
for i=1:nargin(eval(['@' distr]))-1
    if i>1
        str=[str ','];
    end
    str=[str 'x(' num2str(i) ')'];
    
end

Lvalue=[];
for j=2:length(q)
  Lvalue(j-1)=[N(j)*log(integral(@(t) eval([distr '(t,' str ')']), q(j-1), q(j)))];
end

Lvalue(isnan(Lvalue))=log(typecast(uint64(1),'double'));
%Deal with negative inf (when pdf=0)
Lvalue(isinf(Lvalue) & (Lvalue<Inf))=log(typecast(uint64(1),'double'));


QMLE=-sum(Lvalue);

end

