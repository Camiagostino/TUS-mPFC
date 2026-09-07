function [Sortstr] = get_sorthandle
%%
%just returns a function to sort the input to fDDM to the DDM parameters.
global VP
FixVsFree = ismember(1:9, VP.Params);
Sortstr='@(x) deal(';k=1;
for c = 1 : 9
    if c == 1 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) '), '];k=k+1;
    elseif c == 1 && ~FixVsFree(c)
        Sortstr = [Sortstr 'DriftRate, '];
    elseif c == 2 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) ')^2, '];k=k+1;
    elseif c == 2 && ~FixVsFree(c)
        Sortstr = [Sortstr 'DRvariance, '];
    elseif c == 3 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) '), '];k=k+1;
    elseif c == 3 && ~FixVsFree(c)
        Sortstr = [Sortstr 'Boundary, '];
    elseif c == 4 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) '), '];k=k+1;
    elseif c == 4 && ~FixVsFree(c)
        Sortstr = [Sortstr 'StartVar, '];
    elseif c == 5 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) ')*1000, '];k=k+1;
    elseif c == 5 && ~FixVsFree(c)
        Sortstr = [Sortstr 'NondecisionTime, '];
    elseif c == 6 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) ')*100, '];k=k+1;
    elseif c == 6 && ~FixVsFree(c)
        Sortstr = [Sortstr 'st, '];
    elseif c == 7 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) '), '];k=k+1;
    elseif c == 7 && ~FixVsFree(c)
        Sortstr = [Sortstr 'k, '];
    elseif c == 8 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) '), '];k=k+1;
    elseif c == 8 && ~FixVsFree(c)
        Sortstr = [Sortstr 'FlankerValue, '];
    elseif c == 9 && FixVsFree(c)
        Sortstr = [Sortstr 'x(' num2str(k) ')^2'];k=k+1;
    elseif c == 9 && ~FixVsFree(c)
        Sortstr = [Sortstr 'sf '];
    end
end
Sortstr = [Sortstr ')'];
return


