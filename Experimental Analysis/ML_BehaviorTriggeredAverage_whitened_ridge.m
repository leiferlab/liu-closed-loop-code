function [BTA, behaviorCounts, BTA_RMSD] = ML_BehaviorTriggeredAverage_whitened_ridge(Xtrain,Xtest,Ytrain,Ytest)
%This function performs the core BTA caculation using maximum likelihood
%with whitening and ridge regression; not used in paper

%nThi = size(X,1);
ntfilt = size(Xtrain,2);
number_of_behaviors = size(Ytrain,2);

% Set up grid of lambda values (ridge parameters)
lamvals = 2.^(15:46); % it's common to use a log-spaced set of values
nlam = length(lamvals);

BTA = zeros(number_of_behaviors,ntfilt);
BTA_RMSD = zeros(number_of_behaviors,ntfilt);
behaviorCounts = sum(Ytrain,1);

for behavior_index = 1:number_of_behaviors
    % Allocate space for train and test errors
    msetrain = zeros(nlam,1);  % training error
    msetest = zeros(nlam,1);   % test error
    w_ridge = zeros(ntfilt,nlam); % filters for each lambda

    spstrain = Ytrain(:,behavior_index);
    spstest =  Ytest(:,behavior_index);

    % Precompute some quantities (X'X and X'*y) for training and test data
    XXtr = Xtrain'*Xtrain; %can be replaced by precomputed stimulus covariance
    
    XYtr = Xtrain'*spstrain;  % spike-triggered average, training data
    Imat = eye(ntfilt); % identity matrix of size of filter
    
%     figure; hold on;
    for jj = 1:nlam
        % Compute ridge regression estimate
        w = (XXtr+lamvals(jj)*Imat) \ XYtr; 

        % Compute MSE
        msetrain(jj) = mean((spstrain-Xtrain*w).^2); % training error
        msetest(jj) = mean((spstest-Xtest*w).^2); % test error

        % store the filter
        w_ridge(:,jj) = w;

%         % plot it
%         plot(w);
%         title(['ridge estimate: lambda = ', num2str(lamvals(jj))]);
%         xlabel('time before spike (s)'); drawnow;
    end
%     hold off;

    %the filter is the the one with the min MSE for the testing set
    [min_msetest, min_idx] = min(msetest);
    w_ridge_min = w_ridge(:,min_idx);
%     plot(w_ridge_min)
    BTA(behavior_index,:) = w_ridge_min';
    %the error bar on the BTA is the sqrt of the mean squared error
    BTA_RMSD(behavior_index,:) = sqrt(min_msetest);
end

