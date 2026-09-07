function X = AGF_easy_CSD(input, lookup)      
        
%for i = 1:input.nbchan
%    E{i} = input.chanlocs(i).labels;
%end
%E = E';
%E(59:62)=[];
%M = ExtractMontage('/Volumes/AGF work/SMAC/1auswertung/0myfunctions/data/10-5-System_Mastoids_EGI129.csd',E);

if lookup==0
    [G,H] = GetGH(M);
else
    load( '/Volumes/AGF work/SMAC/1auswertung/0myfunctions/data/CSD_MATRIX_G.mat');
    load( '/Volumes/AGF work/SMAC/1auswertung/0myfunctions/data/CSD_MATRIX_H.mat');
end

%D = input.data([1:58 63:64],:,:);
D = input.data(:,:,:);
X = CSD(D, G, H);