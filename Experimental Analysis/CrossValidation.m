% function [K_LNP, K_Shuffle] = CrossValidation()
%CrossValidation compares the predicted reversal rate with a null model;
%not used in paper
% %   Detailed explanation goes here
    fps = 14;
    number_of_trials = 100;
    dt = 1/(fps*60);
    training_ratio = .8;
    load('reference_embedding.mat');
    number_of_behaviors = max(L(:))-1;
    relevant_track_fields = {'BehavioralTransition','Frames','LEDPower','LEDVoltage2Power'};    
    %folders = getfoldersGUI();
    
    %load the tracks
    %[allTracks, folder_indecies, track_indecies] = loadtracks(folders,relevant_track_fields);
    
    % Get binary array of when behavior transitions are
    allTracks = get_behavior_triggers(allTracks,false);

    LNPScore = zeros(number_of_behaviors, number_of_trials);
    ML_LNPScore = zeros(number_of_behaviors, number_of_trials);
    ShuffleScore = zeros(number_of_behaviors, number_of_trials);

    for trial_index = 1:number_of_trials
        %take half the tracks and use them to fit the model while use the other
        %half to cross validate
        PermutatedTracks = randperm(length(allTracks));
        fit_indecies = PermutatedTracks(1:floor(length(allTracks)*training_ratio));
        validation_indecies = PermutatedTracks(floor(length(allTracks)*training_ratio)+1:end);

        fitTracks = allTracks(fit_indecies);
        validationTracks = allTracks(validation_indecies);
        fit_folder_indecies = folder_indecies(fit_indecies);

        %fit the LNP
        [LNPStats, meanLEDPower, stdLEDPower] = FitLNP(fitTracks, fit_folder_indecies, folders, true);
        validationTracks(1).PredictedRate = []; %preallocate memory
        %calculate the predicted rate for each validation track
        for validation_track_index = 1:length(validationTracks)
			validationTracks(validation_track_index).PredictedRate = zeros(number_of_behaviors, length(validationTracks(validation_track_index).LEDPower));
        	for behavior_index = 1:number_of_behaviors
            	validationTracks(validation_track_index).PredictedRate(behavior_index,:) = PredictLNP(validationTracks(validation_track_index).LEDPower, LNPStats(behavior_index).linear_kernel, LNPStats(behavior_index).non_linearity_fit);
        	end
        end
        Behaviors = [validationTracks.Behaviors];
        PredictedRate = [validationTracks.PredictedRate];
        for behavior_index = 1:number_of_behaviors
            LNPScore(behavior_index,trial_index) = compare_predicted_and_actual_rates(PredictedRate(behavior_index,:),Behaviors(behavior_index,:));
        end
        
        %fit the LNP with ridge regression and ML
        [ML_LNPStats, meanLEDPower, stdLEDPower] = ML_FitLNP(fitTracks, fit_folder_indecies, folders, true);
        %calculate the predicted rate for each validation track
        for validation_track_index = 1:length(validationTracks)
			validationTracks(validation_track_index).PredictedRate = zeros(number_of_behaviors, length(validationTracks(validation_track_index).LEDPower));
        	for behavior_index = 1:number_of_behaviors
            	validationTracks(validation_track_index).PredictedRate(behavior_index,:) = PredictLNP(validationTracks(validation_track_index).LEDPower, ML_LNPStats(behavior_index).linear_kernel, ML_LNPStats(behavior_index).non_linearity_fit);
        	end
        end
        Behaviors = [validationTracks.Behaviors];
        PredictedRate = [validationTracks.PredictedRate];
        for behavior_index = 1:number_of_behaviors
            ML_LNPScore(behavior_index,trial_index) = compare_predicted_and_actual_rates(PredictedRate(behavior_index,:),Behaviors(behavior_index,:));
        end
        
        
        %shuffle and refit
        fitTracks = get_behavior_triggers(fitTracks,true);

        %fit the shuffled LNP
        [Shuffled_LNPStats, meanLEDPower, stdLEDPower] = FitLNP(fitTracks, fit_folder_indecies, folders, true);
        %calculate the predicted rate for each validation track
        for validation_track_index = 1:length(validationTracks)
			validationTracks(validation_track_index).PredictedRate = zeros(number_of_behaviors, length(validationTracks(validation_track_index).LEDPower));
        	for behavior_index = 1:number_of_behaviors
            	validationTracks(validation_track_index).PredictedRate(behavior_index,:) = PredictLNP(validationTracks(validation_track_index).LEDPower, Shuffled_LNPStats(behavior_index).linear_kernel, Shuffled_LNPStats(behavior_index).non_linearity_fit);
        	end
        end
        Behaviors = [validationTracks.Behaviors];
        PredictedRate = [validationTracks.PredictedRate];
        for behavior_index = 1:number_of_behaviors
            ShuffleScore(behavior_index,trial_index) = compare_predicted_and_actual_rates(PredictedRate(behavior_index,:),Behaviors(behavior_index,:));
        end
    end
    

    
%     figure
%     hist(LNPScore)
    
%     parfor trial_index = 1:trial_number
%         %take half the tracks and use them to fit the model while use the other
%         %half to cross validate
%         PermutatedTracks = randperm(length(allTracks));
%         fitTracks = allTracks(PermutatedTracks(1:floor(length(allTracks)/2)));
%         validationTracks = allTracks(PermutatedTracks(floor(length(allTracks)/2)+1:end));
% 
%         %fit the LNP
%         [linear_kernel, non_linearity_fit, ~, ~, ~, ~, ~, ~] = FitLNP(fitTracks);
% 
%         exp_fit_a = non_linearity_fit.a;
%         exp_fit_b = non_linearity_fit.b;
% 
%         %calculate the predicted rate for each validation track
%         for validation_track_index = 1:length(validationTracks)
%             validationTracks(validation_track_index).PredictedRate = PredictLNP(validationTracks(validation_track_index).LEDVoltages, linear_kernel, exp_fit_a, exp_fit_b);
%         end
%         Behaviors = double(~[validationTracks.Behaviors]);
%         %Behaviors(Behaviors == 0) = -1;
%         PredictedRate = [validationTracks.PredictedRate];
%         rdt = PredictedRate.*dt;
%         PredictedProbability = rdt.*exp(rdt);
%         ShuffledBehaviors = Behaviors(randperm(length(Behaviors)));
%         
%         %K_LNP(trial_index) = dot(Behaviors, PredictedRate)/norm(Behaviors)/norm(PredictedRate);
%         ShuffleScore(trial_index) = dot(ShuffledBehaviors, PredictedProbability)/norm(ShuffledBehaviors)/norm(PredictedProbability);
%     end
%     figure
%     hist(ShuffleScore)


    p = CompareTwoHistograms(LNPScore(1,:), ML_LNPScore(1,:), 'LNP Score', 'MLLNP Score')
    % end