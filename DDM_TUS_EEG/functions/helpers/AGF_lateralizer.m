function [outdata] = lateralizer( data, mode, middle)
%function that creates lateralized data by subtracting the timecourse of
%each electrode on the right hemisphere from the left hemisphere (mode=1)
%and vice versa (mode=2). Returns result of this calculation. 
%Middle sets the value for central electrodes, NaN if they sould be left
%normal. 
%
%AGF, 2012

%defaults
flipind=[3,2,1,8,7,6,5,4,17,16,15,14,13,12,11,10,9,26,25,24,23,22,21,20,19,18,33,32,31,30,29,28,27,41,40,39,38,37,36,35,34,50,49,48,47,46,45,44,43,42,55,54,53,52,51,58,57,56,59,60,62,61,63,64];

for c=1:length(flipind)
    if c~=2 && c~=6 && c~=13 && c~=22 && c~=30 && c~=46 && c~=53 && c~=57 %Not a central electrode
        if mode==1
            if mod(c,2) %left electrode
                outdata(c,:,:)=data(c,:,:)-data(flipind(c),:,:);
            else %right electrode
                outdata(c,:,:)=data(flipind(c),:,:)-data(c,:,:);
            end;
        else
            if mod(c,2) %left electrode
                outdata(c,:,:)=data(flipind(c),:,:)-data(c,:,:);
            else %right electrode
                outdata(c,:,:)=data(c,:,:)-data(flipind(c),:,:);
            end;
        end;
    else
        if ~isnan(middle)
            outdata(c,:,:)=0;
        else
            outdata(c,:,:)=data(c,:,:);
        end;
    end;
end;