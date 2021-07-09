function success = calculate_behaviors(folder_name)
% given tracks that has embeddings, classify every time point into a
% behavior
    addpath(genpath(pwd))
    %set up parameters
    parameters = load_parameters(folder_name);
    if parameters.TrackOnly
        success = true;
        return
    end
    load('reference_embedding.mat')
    number_of_behaviors = max(L(:)-1);
    num_velocity_behaviors = numel(velocity_based_behavior_names);

    relevant_track_fields = {'Embeddings','Velocity'};

    %% Load tracks
    Tracks = load_single_folder(folder_name, relevant_track_fields);
    if isempty(Tracks)
        error('Empty Tracks');
    end


    %get the stereotyped behaviors
    Tracks = find_stereotyped_behaviors(Tracks, L, xx, velocity_based_behavior_edges, parameters);

    % Get binary array of when certain behaviors start
    Tracks(1).Behaviors = [];

    for track_index = 1:length(Tracks)
        % for behavioral mapping
        triggers = false(number_of_behaviors, size(Tracks(track_index).Embeddings,1)); %a binary array of when behaviors occur
        for behavior_index = 1:number_of_behaviors
            transition_indecies = Tracks(track_index).BehavioralTransition(:,1) == behavior_index;
            %transition into of
            transition_start_frames = Tracks(track_index).BehavioralTransition(transition_indecies,2);
            triggers(behavior_index,transition_start_frames) = true;
%                 %transition out of
%                 transition_end_frames = Tracks(track_index).BehavioralTransition(transition_indecies,3);
%                 triggers(behavior_index,transition_end_frames) = true;
        end
        Tracks(track_index).Behaviors = triggers(:,1:size(Tracks(track_index).Embeddings,1));
        % for velocity based classification
        
        triggers = false(num_velocity_behaviors, numel(Tracks(track_index).Velocity)); %a binary array of when behaviors occur
        transition_indecies = diff(Tracks(track_index).VelocityBehavior) ~= 0;
        transition_indecies = find([false, transition_indecies]);
        
        for transition_index = transition_indecies
            %transition into of
            triggers(Tracks(track_index).VelocityBehavior(transition_index),transition_index) = true;
        end
        Tracks(track_index).VelocityBehaviorTriggers = triggers(:,1:size(Tracks(track_index).Embeddings,1));

    end
    
    %save
    savetracks(Tracks, folder_name);
    success = true;    
 end