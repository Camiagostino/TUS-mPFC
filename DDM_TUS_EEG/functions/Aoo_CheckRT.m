

datapath = 'data/TUS_EEG/';

F = dir ([datapath '*.mat']);

for i = 1:length(F)

    % load('data/All Behavior4.mat')
    load([datapath F(i).name]);

    figure(1);clf;
    histogram([logg.RT])
    title(['min RT: ' num2str(min([logg.RT])) '; max RT: ' num2str(max([logg.RT]))])
    pause

end