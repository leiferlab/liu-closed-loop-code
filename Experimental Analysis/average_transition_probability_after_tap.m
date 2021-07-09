%function [tap_transition_rate_of_interest,control_transition_rate_of_interest,tap_transition_rate_of_interest_std,control_transition_rate_of_interest_std,h,p,tap_transition_total_count,control_transition_total_count,tap_observation_total_count,control_observation_total_count] = average_transition_probability_after_tap(folders_platetap, behavior_from, behavior_to)
% this function looks at the conditional probability after a platetap and compares
% it to the control of the time point in between platetaps. If the
% behavior_from is 0, it is ignored
    load('reference_embedding.mat')
    %load tracks
    relevant_track_fields = {'BehavioralTransition','Frames'};

    %load stimuli.txt from the first experiment
    normalized_stimuli = 1; %delta function
    time_window_before = 0;
%     time_window_after = 14; %transition rate average for 1 seconds after tap
    time_window_after = 28; %transition rate average for 2 seconds after tap

    number_of_behaviors = max(L(:)-1);
    all_edge_pairs = get_edge_pairs(number_of_behaviors);
    number_of_edges = length(all_edge_pairs);

    tap_observation_total_count = 0;
    control_observation_total_count = 0;  
    tap_transition_counts = zeros(number_of_behaviors,number_of_behaviors);
    control_transition_counts = zeros(number_of_behaviors,number_of_behaviors);
    rows_per_page = 9;
    fps = 14;
    folders_platetap = getfoldersGUI();
    
    %% behavioral rate compare
    for folder_index = 1:length(folders_platetap)
        %for each experiment, search for the occurance of each stimulus after
        %normalizing to 1
        LEDVoltages = load([folders_platetap{folder_index}, filesep, 'LEDVoltages.txt']);
        % LEDVoltages = LEDVoltages(randperm(length(LEDVoltages))); %optional, randomly permuate the taps
        %LEDVoltages(LEDVoltages>0) = 1; %optional, make the stimulus on/off binary

        %find when each stimuli is played back by convolving the time
        %reversed stimulus (cross-correlation)
        xcorr_ledvoltages_stimulus = padded_conv(LEDVoltages, normalized_stimuli);
        peak_thresh = 0.99.*max(xcorr_ledvoltages_stimulus); %the peak threshold is 99% of the max (because edge effects)
        [~, tap_frames] = findpeaks(xcorr_ledvoltages_stimulus, 'MinPeakHeight', peak_thresh,'MinPeakDistance',14);
        
        %generate a series of control taps
        control_frame_shift = round((tap_frames(2)-tap_frames(1))/2); %the control taps are exactly in between taps
        control_LEDVoltages = circshift(LEDVoltages,[0,control_frame_shift]);
        xcorr_ledvoltages_stimulus = padded_conv(control_LEDVoltages, normalized_stimuli);
        [~, control_frames] = findpeaks(xcorr_ledvoltages_stimulus, 'MinPeakHeight', peak_thresh,'MinPeakDistance',14);
        
        %load the tracks for this folder
        [current_tracks, ~, ~] = loadtracks(folders_platetap(folder_index),relevant_track_fields);
        %generate the Behavior matricies
        current_tracks = get_directional_behavior_triggers(current_tracks);
        current_tracks = BehavioralTransitionToBehavioralAnnotation(current_tracks);

        
        %get the transitions probability for tap condition
        for critical_frame_index = 1:length(tap_frames)
            %for every time a stimulus is delivered, look at a certain range of
            %frames only if the track fits certain criteria
            current_critical_frame = tap_frames(critical_frame_index);
            if current_critical_frame + time_window_after <= length(LEDVoltages) && current_critical_frame - time_window_before >= 1
                %get tracks that last through the entire duration of the window
                tracks_within_critical_window = FilterTracksByTime(current_tracks,current_critical_frame - time_window_before, current_critical_frame + time_window_after, true);
                if ~isempty(tracks_within_critical_window)
                    tracks_on_current_critical_frame = FilterTracksByTime(tracks_within_critical_window,current_critical_frame, current_critical_frame);
                    tap_observation_total_count = tap_observation_total_count + length(tracks_on_current_critical_frame); % keep track of how many observations we take
                    
                    for tracks_within_critical_window_index = 1:length(tracks_within_critical_window)
                        %loop through all the selected tracks, find
                        %behavior_from and behavior_to, and add to each of
                        %their counts
                        behavioral_annotation = tracks_within_critical_window(tracks_within_critical_window_index).BehavioralAnnotation;
                        behavior_from = behavioral_annotation(1);
                        behavioral_annotation(behavioral_annotation == behavior_from) = 0;
                        behavior_to = behavioral_annotation(find(behavioral_annotation,1));
                        if ~isempty(behavior_from) &&  ~isempty(behavior_to) && (behavior_from > 0) && (behavior_to > 0) && (behavior_from <= number_of_behaviors) && (behavior_to <= number_of_behaviors)
                            tap_transition_counts(behavior_from, behavior_to) = tap_transition_counts(behavior_from, behavior_to) + 1;
                        elseif ~isempty(behavior_from) && (behavior_from > 0)
                            %no transition out of behavior_from obseved
                            tap_transition_counts(behavior_from, behavior_from) = tap_transition_counts(behavior_from, behavior_from) + 1;
                        end
                    end                    
                end
            end
        end
        
        %get the transitions rates for control condition
        for critical_frame_index = 1:length(control_frames)
            %for every time a stimulus is delivered, look at a certain range of
            %frames only if the track fits certain criteria
            current_critical_frame = control_frames(critical_frame_index);
            if current_critical_frame + time_window_after <= length(LEDVoltages) && current_critical_frame - time_window_before >= 1
                %get tracks that last through the entire duration of the window
                tracks_within_critical_window = FilterTracksByTime(current_tracks,current_critical_frame - time_window_before, current_critical_frame + time_window_after, true);
                if ~isempty(tracks_within_critical_window)
                    tracks_on_current_critical_frame = FilterTracksByTime(tracks_within_critical_window,current_critical_frame, current_critical_frame);
                    control_observation_total_count = control_observation_total_count + length(tracks_on_current_critical_frame); % keep track of how many observations we take
                    
                    for tracks_within_critical_window_index = 1:length(tracks_within_critical_window)
                        %loop through all the selected tracks, find
                        %behavior_from and behavior_to, and add to each of
                        %their counts
                        behavioral_annotation = tracks_within_critical_window(tracks_within_critical_window_index).BehavioralAnnotation;
                        behavior_from = behavioral_annotation(1);
                        behavioral_annotation(behavioral_annotation == behavior_from) = 0;
                        behavior_to = behavioral_annotation(find(behavioral_annotation,1));
                        if ~isempty(behavior_from) &&  ~isempty(behavior_to) && (behavior_from > 0) && (behavior_to > 0)
                            control_transition_counts(behavior_from, behavior_to) = control_transition_counts(behavior_from, behavior_to) + 1;
                        elseif ~isempty(behavior_from) && (behavior_from > 0)
                            %no transition out of behavior_from obseved
                            control_transition_counts(behavior_from, behavior_from) = control_transition_counts(behavior_from, behavior_from) + 1;
                        end
                    end                    
                end
            end
        end        
    end

%calculate condidtional prob and significance testing
tap_difference_significant = false(number_of_behaviors, number_of_behaviors);
tap_pvalue = eye(number_of_behaviors);
tap_transition_prob = tap_transition_counts;
control_transition_prob = control_transition_counts;
tap_transition_prob_std = tap_transition_counts;
control_transition_prob_std = control_transition_counts;

for behavior_from = 1:number_of_behaviors
    tap_transition_prob(behavior_from,:) = tap_transition_counts(behavior_from,:) ./ sum(tap_transition_counts(behavior_from,:));
    tap_transition_prob_std(behavior_from,:) = sqrt(tap_transition_counts(behavior_from,:)) ./ sum(tap_transition_counts(behavior_from,:));
    control_transition_prob(behavior_from,:) = control_transition_counts(behavior_from,:) ./ sum(control_transition_counts(behavior_from,:));
    control_transition_prob_std(behavior_from,:) = sqrt(control_transition_counts(behavior_from,:)) ./ sum(control_transition_counts(behavior_from,:));
%     tap_transition_prob(behavior_from,:) = tap_transition_counts(behavior_from,:) ./ sum(tap_transition_counts(behavior_from,:)) .* fps .* 60;
%     tap_transition_prob_std(behavior_from,:) = sqrt(tap_transition_counts(behavior_from,:)) ./ sum(tap_transition_counts(behavior_from,:)) .* fps .* 60;
%     control_transition_prob(behavior_from,:) = control_transition_counts(behavior_from,:) ./ sum(control_transition_counts(behavior_from,:)) .* fps .* 60;
%     control_transition_prob_std(behavior_from,:) = sqrt(control_transition_counts(behavior_from,:)) ./ sum(control_transition_counts(behavior_from,:)) .* fps .* 60;
    for behavior_to = 1:number_of_behaviors
        tap_pvalue(behavior_from,behavior_to) = testPoissonSignificance(tap_transition_counts(behavior_from,behavior_to),control_transition_counts(behavior_from,behavior_to),sum(tap_transition_counts(behavior_from,:)),sum(control_transition_counts(behavior_from,:)),0,2);
        tap_difference_significant(behavior_from,behavior_to) = tap_pvalue(behavior_from,behavior_to) < 0.05./ double(number_of_behaviors.*number_of_behaviors);
    end
end

%plot it
figure
for behavior_from = 1:number_of_behaviors
    for behavior_to = 1:number_of_behaviors
        if control_transition_counts(behavior_from,behavior_to) == 0 && tap_transition_counts(behavior_from,behavior_to) == 0
        else
            scrollsubplot(rows_per_page,double(number_of_behaviors),double((behavior_from-1)*number_of_behaviors+behavior_to))
            barwitherr([control_transition_prob_std(behavior_from, behavior_to); tap_transition_prob_std(behavior_from, behavior_to)],[control_transition_prob(behavior_from, behavior_to); tap_transition_prob(behavior_from, behavior_to)],'FaceColor',behavior_colors(behavior_to,:))
            ax = gca;
            if tap_difference_significant(behavior_from,behavior_to)
                sigstar({[1,2]},0.05);
            end

            title(['n=', num2str(control_transition_counts(behavior_from, behavior_to)),', ',num2str(tap_transition_counts(behavior_from, behavior_to))],'Color', 'k', 'FontWeight', 'normal', 'Fontsize', 14)
            box('off')
            set(gca,'XTick','')
            set(gca,'fontsize',14)
            axis([0 3 0 1]);
            set(gca,'XTick','')
            set(gca,'YTick','')
        end
    end
end