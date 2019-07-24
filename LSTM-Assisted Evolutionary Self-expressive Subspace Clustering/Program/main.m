%% Prepare data
clear; clc; load('data/Evolutionry155.mat');

%% Processing sequences

N_snapshot_max = 0;

for i = 1:length(Evolutionry155)
    
    N_snapshot = Evolutionry155(i).N_snapshot;
    
    if i == 1
        
        N_snapshot_max = N_snapshot;
    end
    N_snapshot_max = max(N_snapshot_max,N_snapshot);
end

% parameter initiali_errzation

% parameter initiali_errzation
[WW1, WW2, C1] = deal(cell(length(Evolutionry155),1)); % Initial predictor and target Sequences

[C2, A, spectralC] = deal(cell(length(Evolutionry155),N_snapshot_max)); % initial clustering results

results = struct([]); % initial final results

num_algs = 1;

Tot = length(Evolutionry155);

for i = 1:Tot
    
    %% data initiali_errzation
    fprintf('Sequences: %i out of %i\n',i,length(Evolutionry155));
    % extract out the ith video sequence
    
    % s = Evolutionry155(i).s;
    N_snapshot = Evolutionry155(i).N_snapshot;
    ngroups = Evolutionry155(i).N_motion;
    % F = Evolutionry155(i);
    N = Evolutionry155(i).N;
    
    % xord = Evolutionry155(i).xord;
    % x1ord = Evolutionry155(i).x1ord;
    % yord = Evolutionry155(i).yord;
    
    snapshots_xord = Evolutionry155(i).snapshots_xord;
    % snapshots_x1ord = Evolutionry155(i).snapshots_x1ord;
    % snapshots_yord = Evolutionry155(i).snapshots_yord;
    
    % Reshape the ith sequence
    
    % temporary parameter initiali_errzation
    [WWt1,WWt2] = deal([]);
    
    errorss = zeros(num_algs,N_snapshot);
%     timess = zeros(num_algs,N_snapshot);
    %% data precessing
    nKeypoints = 0;
    for ii = 1:N_snapshot
        
        WW = snapshots_xord(ii).WW;
        kappa = 2e-7;
          
        
        % Dimension deduction
        [U,S,V] = svd(WW',0);
        
                
        % column normali_errzation
        WW = cnormalize(U(:,1:4*ngroups)');
        
        
        % predictors
        nKeypoints = size(WW,2); % # of key points selected
                
        WWt1(:,ii) = WW(:);
             
        
        % Target
        WWt =WW'*WW;
        WWt2(:,ii) = WWt(:);
    end
    
    WW1{i} = WWt1;
    WW2{i} = WWt2;
     
    
    %% Define LSTM network architecture
    featureDimension = size(WWt1,1);
    
    numHiddenUnits = ceil(nKeypoints/5); % tunable hyperparam
    
    numResponses = nKeypoints^2 - nKeypoints;
    
    paddingSize = nKeypoints;
    
    lambda = 0.1; % tunable hyperparam
    
    
    layers = [ ...
        sequenceInputLayer(featureDimension)
        lstmLayer(numHiddenUnits,'OutputMode','sequence')
        fullyConnectedLayer(numResponses)
        myPaddingLayer(paddingSize)
        myRegressionLayer('Evolving', lambda)];
    
    
    maxEpochs = 50;  %60
    %     miniBatchSize = 20;
    
    
    options = trainingOptions('adam', ...
        'MaxEpochs',maxEpochs, ...
        'InitialLearnRate',0.001, ...
        'GradientThreshold',1, ...
        'Shuffle','never', ...
        'Verbose',0); %'Plots','training-progress',...
    %% Network training
    ts=cputime;
    
    net = trainNetwork(WWt1,WWt2,layers,options);
    
    C1{i} = double(predict(net,WWt1));
    
    %% sepctral clustering
    for iii = 1:N_snapshot
        
        % reshape Cl
        C2{i,iii} = reshape(C1{i}(:,iii),nKeypoints,nKeypoints);
        
        % create affinity matrix
        A{i,iii} = abs(C2{i,iii})+abs(C2{i,iii}');
        
        % Spectral Clustering
        [~,~,~,~,spectralC{i,iii},~]=spectralcluster(A{i,iii},ngroups,ngroups);
        
        % error calculation
        errorss(1,iii) = missclass(spectralC{i,iii},N,ngroups)/sum(N)*100;
    end

    results(i).time = (cputime-ts);
    results(i).error = errorss;
end

%% Save and average

Avg_err = 0;
Avg_tim = 0;
RES = 0;

for i = 1:2 %Tot
    
    ali_err = results(i).error;
    ali_tim = results(i).time;
    
    sz = size(ali_err,2);
    
    Avg_err = Avg_err + sum(ali_err,2)/sz;
    Avg_tim = Avg_tim + ali_tim/sz;
    
    ali_err = results(i).error(:,2:end);
    sz = size(ali_err,2);
    RES = RES + sum(ali_err,2)/sz;
end

Avg_err = Avg_err/Tot
RES = RES/Tot
Avg_tim = Avg_tim/Tot

clearvars -except results
save('Allresults.mat');
