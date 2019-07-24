%% Prepare data
clear; clc; load('data/Evolutionry155.mat');


%% add moving overlap 5 clumns

% Parameter Initialization
% movingWindow = 5;

for ii = 1: size(Evolutionry155,2)
    
    snapshots_xord = Evolutionry155(ii).snapshots_xord;
    N_snapshot = Evolutionry155(ii).N_snapshot;
    
    % moving window size
    movingWindow = size(snapshots_xord(1).WW,1);
    
    % concatenate all sanpshots together
    snapshots_xord_mo5_t = [];
    
    for i = 1:N_snapshot
        
        snapshots_xord_mo5_t = cat(1,snapshots_xord_mo5_t,snapshots_xord(i).WW);
    end
    
    % reshape overlapping matrix
    snapshots_mo5_size = size(snapshots_xord_mo5_t);
    times = snapshots_mo5_size(1) - movingWindow +1;
    
    % overlapping snapshots
    % snapshots_xord_mo5_t1.t = 1:times;
    
    % video sequence
    snapshots_xord_mo5_t1 = [];
    for i = 1: times
        
        snapshots_xord_mo5_t1(i).t = i;
        
        snapshots_xord_mo5_t1(i).WW = snapshots_xord_mo5_t(i:i+movingWindow-1,:);
    end
    
    Evolutionry155(ii).snapshots_xord_mo5 = snapshots_xord_mo5_t1;
    
    Evolutionry155(ii).snapshots_mo5 = times;
end

%% save data
save ('data/Evolutionry155.mat')

%%
% load('data/Evolutionry155.mat');

