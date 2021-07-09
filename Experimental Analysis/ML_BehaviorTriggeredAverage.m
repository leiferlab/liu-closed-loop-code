function [BTA, behaviorCounts, BTA_RMSD, BTA_stats] = ML_BehaviorTriggeredAverage(Behaviors, Stim, bootstrap)
    %finds the behavior triggered average, optionally determine the
    %significance of the BTA by shuffling transitions randomly and finding
    %the BTA; not used in paper
    fps = 14;
    BTA_seconds_before_and_after = 10;
    number_of_random_shuffles = 10;
    trainfrac = .8;  % fraction of data to use for training

    number_of_behaviors = size(Behaviors{1},1);
    seconds_before = BTA_seconds_before_and_after;
    seconds_after = BTA_seconds_before_and_after;
    BTA_length = (fps*seconds_before)+(fps*seconds_after)+1;
    
    if nargin < 3
        %default for bootstrap is true
        bootstrap = true;
    end
    
    %getting the design matrix
    StimRows = cell(size(Stim)); 
    for track_index = 1:length(Stim)
        StimRows{track_index} = makeStimRows(transpose(Stim{track_index}),BTA_length);
    end
    X = vertcat(StimRows{:}); %build design matrix
    nThi = size(X,1);
    clear StimRows
    
    %get triggers accounting for edges
    no_edge_Behaviors = circshift_triggers(Behaviors, BTA_seconds_before_and_after, false, true);
    Y = horzcat(no_edge_Behaviors{:})';
    clear no_edge_Behaviors

    % Divide data into "training" and "test" sets for cross-validation
    ntrain = ceil(nThi*trainfrac);  % number of training samples
    ntest = nThi-ntrain; % number of test samples
    randsequence = randperm(nThi);
    iitest = randsequence(1:ntest); % time indices for test
    iitrain = randsequence(ntest+1:end);   % time indices for training
    Xtrain = X(iitrain,:); % training stimulus
    Xtest = X(iitest,:); % test stimulus
    Ytrain = Y(iitrain,:);
    Ytest = Y(iitest,:);
    clear X randsequence iitest iitrain randsequence  %get rid of original design matrix to save memory

    [BTA, behaviorCounts, BTA_RMSD] = ML_BehaviorTriggeredAverage_whitened_ridge(Xtrain,Xtest,Ytrain,Ytest);
    
    if bootstrap
        BTA_norm = sqrt(sum(BTA.^2, 2));

        %perform bootstrapping
        shuffle_norms = zeros(number_of_behaviors,number_of_random_shuffles);
        for shuffle_index = 1:number_of_random_shuffles
            % do the BTA many many times
            
            %getting the design matrix
            StimRows = cell(size(Stim)); 
            for track_index = 1:length(Stim)
                StimRows{track_index} = makeStimRows(transpose(Stim{track_index}),BTA_length);
            end
            X = vertcat(StimRows{:}); %build design matrix
            clear StimRows

            % Divide data into "training" and "test" sets for cross-validation
            ntrain = ceil(nThi*trainfrac);  % number of training samples
            ntest = nThi-ntrain; % number of test samples
            randsequence = randperm(nThi);
            iitest = randsequence(1:ntest); % time indices for test
            iitrain = randsequence(ntest+1:end);   % time indices for training
            Xtrain = X(iitrain,:); % training stimulus
            Xtest = X(iitest,:); % test stimulus
            Ytrain = Y(iitrain,:);
            Ytest = Y(iitest,:);
            clear X randsequence iitest iitrain randsequence  %get rid of original design matrix to save memory

            [shuffle_BTA, ~, ~] = ML_BehaviorTriggeredAverage_whitened_ridge(Xtrain,Xtest,Ytrain,Ytest);
            
            %get the L2 norm of the shuffled_BTAs
            shuffle_norm = sqrt(sum(shuffle_BTA.^2, 2));
            shuffle_norms(:,shuffle_index) = shuffle_norm;

            disp(num2str(shuffle_index));
        end

        allnorms = [BTA_norm, shuffle_norms];
        BTA_percentile = zeros(number_of_behaviors,1);
        for behavior_index = 1:number_of_behaviors
            norm_ranks = tiedrank(allnorms(behavior_index,:)) / (number_of_random_shuffles+1);
            BTA_percentile(behavior_index) = norm_ranks(1);
        end
        BTA_stats.BTA_norm = BTA_norm;
        BTA_stats.shuffle_norms = shuffle_norms;
        BTA_stats.BTA_percentile = BTA_percentile;
    end
    
    BTA_stats.mean_subtracted = true;  
end
