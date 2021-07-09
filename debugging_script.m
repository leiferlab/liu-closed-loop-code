%% use behavioral annotations to determine decision making +
% plot response probability heatmap given head and tail intensities 
% (this is only compatible without delayed stimulus)

n_tracks = zeros(1,n_sti); %number of tracks in each sti intensities
transition_count = 0;

%record all transitions in a 3D matrix indexed from_behavior, to_behavior, stimulus_index
behavioral_transition_matrix_for_all_stimuli = zeros(numel(my_behavior_names),numel(my_behavior_names),n_sti);
behavioral_transition_stats = [];

for stimulus_index = 1:n_sti
    % loop through the variou stimulation conditions
    n_tracks(stimulus_index) = numel(video_possible_frames{stimulus_index});
    if n_tracks <= 5
        %skip analysis for conditions with few tracks
       continue 
    end

    %calculate total stimulus duration, which encompasses the delayed
    %second stimulus. The delay is not applicapble in msot cases
    if saved_head_or_tail_first(stimulus_index)
        %head is first
        total_stimulus_duration_for_stim = saved_stimulus_delays(stimulus_index) + saved_tail_stimulus_durations(stimulus_index);
        head_or_tail_string = 'HeadStimulusFirst';
    else
        %tail is first
        total_stimulus_duration_for_stim = saved_stimulus_delays(stimulus_index) + saved_head_stimulus_durations(stimulus_index);
        head_or_tail_string = 'TailStimulusFirst';
    end

    %determine the behavior before stimulus and then subsequent transitions
    starting_count_for_stimulus = transition_count; 
    for stimulation_event_index = 1:numel(video_possible_frames{stimulus_index})
        %loop through all the animal-stimulation events for the stimulus condition
        stimulation_track_index = video_possible_tracks{stimulus_index}(stimulation_event_index);
        stimulation_frame_index = video_possible_frames{stimulus_index}(stimulation_event_index);
        behavioral_annotation_for_track = allTracks(stimulation_track_index).VelocityBehavior;
        %behavior before is the behavior a smoothing window before the start of the stimulus
        starting_behavior_annotations_for_stim = behavioral_annotation_for_track(stimulation_frame_index - (parameters.TrackingSmoothingWindow*parameters.SampleRate));
        next_behavior_annotations_for_stim = starting_behavior_annotations_for_stim;
        next_behavior_onset_frame = 0;
        next_behavior_duration = 0;

        for frame_index = stimulation_frame_index+1:stimulation_frame_index+total_stimulus_duration_for_stim
            %scan forwards until the behavior changes up until the end of
            %the stimulus, record the duration of the next behavior
            behavioral_annotation_for_frame = behavioral_annotation_for_track(frame_index);

            if behavioral_annotation_for_frame < 1 || behavioral_annotation_for_frame > numel(my_behavior_names) || behavioral_annotation_for_frame == starting_behavior_annotations_for_stim
                continue
            else
                %change detected
                change_detected = false;
                %keep track of the most extreme change instead of any
                %change, since the velocity behaviors are continuous
                if abs(starting_behavior_annotations_for_stim-behavioral_annotation_for_frame) > abs(starting_behavior_annotations_for_stim-next_behavior_annotations_for_stim)
                    next_behavior_annotations_for_stim = behavioral_annotation_for_frame;
                    change_detected = true;
                end

                if change_detected
                    %find the duration of this behavior, count the transition if it is longer than a threshold (0.5s)
                    temp_next_behavior_duration = find(behavioral_annotation_for_frame(frame_index:end) ~= next_behavior_annotations_for_stim, 1, 'first');
                    if isempty(temp_next_behavior_duration)
                        temp_next_behavior_duration = numel(allTracks(stimulation_track_index).Frames) - frame_index;
                    end
                    if temp_next_behavior_duration > parameters.StereotypedBehaviorMinTime * parameters.SampleRate
                        next_behavior_onset_frame = frame_index;
                        next_behavior_duration = temp_next_behavior_duration;
                    end
                end
            end
        end
        
        if starting_behavior_annotations_for_stim > 0 && next_behavior_annotations_for_stim > 0 && next_behavior_onset_frame > 0
            % save the transition stats only if conditions are satisfied
            transition_count = transition_count+1;
            behavioral_transition_stats(transition_count).worm_index = stimulation_track_index;
            behavioral_transition_stats(transition_count).stimulus_index = stimulus_index;
            behavioral_transition_stats(transition_count).starting_behavior = starting_behavior_annotations_for_stim;
            behavioral_transition_stats(transition_count).next_behavior = next_behavior_annotations_for_stim;
            behavioral_transition_stats(transition_count).time_elapsed_since_stim_onset = (next_behavior_onset_frame - stimulation_frame_index) / parameters.SampleRate; % latency between stim and next behavior
            behavioral_transition_stats(transition_count).next_behavior_duration = next_behavior_duration;
        end
    end
    
    %create the behavioral transition matrix
    behavioral_transition_matrix = zeros(numel(my_behavior_names),numel(my_behavior_names));
    behavioral_transition_next_behavior_durations = cell(numel(my_behavior_names),numel(my_behavior_names));
    for transition_index = starting_count_for_stimulus+1:transition_count
        behavioral_transition_next_behavior_durations{behavioral_transition_stats(transition_index).starting_behavior, behavioral_transition_stats(transition_index).next_behavior} = ...
        [behavioral_transition_next_behavior_durations{behavioral_transition_stats(transition_index).starting_behavior, behavioral_transition_stats(transition_index).next_behavior}, 
        behavioral_transition_stats(transition_index).next_behavior_duration];
        
        behavioral_transition_matrix(behavioral_transition_stats(transition_index).starting_behavior,behavioral_transition_stats(transition_index).next_behavior) = ...
            behavioral_transition_matrix(behavioral_transition_stats(transition_index).starting_behavior,behavioral_transition_stats(transition_index).next_behavior) + 1;
    end
    behavioral_transition_matrix_for_all_stimuli(:,:,stimulus_index) = behavioral_transition_matrix;
end

%loop through the behaviors and plot the probability of behavioral
%response with respect to each stimulus intensity, i.e. heatmap
behavioral_response_for_all_stimuli = sum(behavioral_transition_matrix_for_all_stimuli,1); %display only the result probability
for behavioral_index = 1:numel(my_behavior_names)
    response_heatmap_for_behavior = zeros(numel(parameters.RailsIntensities),numel(parameters.RailsIntensities));
    for stimulus_index = 1:n_sti
        head_intensity_index = find(parameters.RailsIntensities == saved_head_stimulus_intensities(stimulus_index));
        tail_intensity_index = find(parameters.RailsIntensities == saved_tail_stimulus_intensities(stimulus_index));            
        response_heatmap_for_behavior(head_intensity_index,tail_intensity_index) = behavioral_response_for_all_stimuli(1,behavioral_index,stimulus_index) ./ squeeze(sum(behavioral_response_for_all_stimuli(:,:,stimulus_index)));
    end

    figure
    plot(1:numel(parameters.RailsIntensities))
    hold on
    imagesc(response_heatmap_for_behavior)
    for behavior_from_intensity_index = 1:numel(parameters.RailsIntensities)
        for behavior_to_intensity_index = 1:numel(parameters.RailsIntensities)
            text(behavior_from_intensity_index,behavior_to_intensity_index, ...
                num2str(round(response_heatmap_for_behavior(behavior_to_intensity_index, behavior_from_intensity_index), 2)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
        end
    end
    xlabel('Tail Stimulus Intensity (uW/mm2)')
    ylabel('Head Stimulus Intensity (uW/mm2)')
    set(gca, 'XTick', 1:numel(parameters.RailsIntensities))
    set(gca, 'XTickLabel', abs(round(parameters.RailsIntensities)))
    set(gca, 'YTick', 1:numel(parameters.RailsIntensities))
    set(gca, 'YTickLabel', abs(round(parameters.RailsIntensities)))
    axis tight
    caxis([0 max(response_heatmap_for_behavior(:))])
%         caxis([0 0.8])
    colorbar
    colormap(othercolor('OrRd9'))
    chart_title = ['Probability of ', my_behavior_names{behavioral_index}, ' response'];

    title(chart_title)
end

%create a an array of pi charts
x = zeros(n_sti, 1);
y = zeros(n_sti, 1);
s = ones(n_sti,1);
graph_data = zeros(n_sti,numel(my_behavior_names));
for stimulus_index = 1:n_sti
    head_intensity_index = find(parameters.RailsIntensities == saved_head_stimulus_intensities(stimulus_index));
    tail_intensity_index = find(parameters.RailsIntensities == saved_tail_stimulus_intensities(stimulus_index));            
    x(stimulus_index) = tail_intensity_index; %saved_tail_stimulus_intensities(stimulus_index) + axis_offset;
    y(stimulus_index) = head_intensity_index; %saved_head_stimulus_intensities(stimulus_index) + axis_offset;
    graph_data(stimulus_index,:) = squeeze(behavioral_response_for_all_stimuli(1,:,stimulus_index));
end
graph_labels = [];
graph_legend = my_behavior_names;
xlab = 'Tail Stimulus Intensity (uW/mm2)';

