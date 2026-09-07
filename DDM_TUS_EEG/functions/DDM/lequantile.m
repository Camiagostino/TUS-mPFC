%calculate the linearly interpolated estimated quantile for the data set,
%according to qq. For details, see  Heathcote, Brown and Mewhort, 2002. 

%More or less identical to quantiles in matlab, a part for the extremes.
%In matlab they are the min and the max of the data, in this function they
%are either like matlab (if nothing is specified) or can be specified by
%the user. 

%data is a COLUMN VECTOR. 
%N is the number of data for each quantile. If the quantile are equally
%spaced, N is equal for each one. Otherwise is not. 

%binCtr is self explenatory.
%binHeight is the heigth of the bin of each quantile assumed that the total
%area is 1. If you want to plot it. 

%Valerio Biscione 30-03-14

function [q,binCtr,binHeigth,N]=lequantile(data,varargin)
p=inputParser;
addParamValue(p, 'eQuant', []); %extreme Quantiles q(1) and q(end)
addParamValue(p, 'qq', [0.0:0.05:1]);
parse(p, varargin{:});
eQ=p.Results.eQuant;
qq=p.Results.qq;
%data=normrnd(1,1,1000,1);
%qq=[0 0.3 0.7 1]; 
%data=[2 4 6]; 
%qq=[0 0.3 0.7 1];

data=sort(data); 
n=length(data);
N(1)=NaN; 
for j=2:length(qq)
    N(j)=(qq(j)-(qq(j-1)))*n; 
end 

for j=2:length(qq)-1
    I=@(j) qq(j).*n+1/2;
    mI= @(j) floor(I(j));
    pI=@(j) ceil(I(j));
    
    q(j)=data(mI(j))+(data(pI(j))-data(mI(j)))*(I(j)-mI(j));
    
    
end
%this is equal to the domain of the distribution. You should be really
%carefull here. Change this according to the distribution you are trying to
%fit. 
if isempty(eQ)
  q(1)=min(data); q(length(qq))=max(data); 
else
 q(1)=eQ(1);   q(length(qq))=eQ(2);
end 

%---calculate the binCntrs
x_l=q(1:end-1); width=q(2:end)-q(1:end-1);
binCtr=x_l+width/2;
area=1./(length(q)-1);
binHeigth=area./width; 
end 