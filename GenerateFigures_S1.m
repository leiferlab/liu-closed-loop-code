%%% Please use this code to generate the following figures presented in the
%%% publication Liu et al.: Fig. 2, Fig. 3, supplementary fig. S1,
%%% supplementary fig. S4, supplementary fig. S6, information related to
%%% anterior vs posterior stimulation in supplementary table S1,
%%% Supplementary Video S2, Supplementary Video S3, Supplementary Video S4, 
%%% and Supplementary Video S5.

clear
clc
close all

cd('/projects/LEIFER/Mochi/APIData')

%% ask for folders and initialize
%%% addpath(genpath(pwd))
load('reference_embedding.mat')
relevant_track_fields = {'Centerlines','Path','Frames','AlignedStimulus', 'Velocity','Length','VelocityBehavior'};
%%%select folders
folders = getfoldersGUI();
parameters = load_parameters(folders{1}); % load parameters for the first folder

%% optional videoplotting
plot_video = false;
N_rows = 4;
N_columns = 5;

%% constants
max_traces = N_rows*N_columns;
min_avg_worm_length = 0.7 * parameters.CameraPixeltommConversion; %in mm (Please note: set this to zero while calculating parameter presented in supp. table S1)
use_max_velocity = true; %use max velocity for 2d histograms, otherwise mean velocity is used
worm_region_count = 2;
worm_region_head = 1; 
worm_region_tail = worm_region_count;

num_velocity_behaviors = numel(velocity_based_behavior_names);

BOOTSTRAP_N = 1000; %number of bootstraps for estimating CI

normalized_stimuli = 1; %delta function
time_window_before = 8*parameters.SampleRate;
time_window_after = 8*parameters.SampleRate;
total_window_frames = time_window_before+time_window_after+1;
fps = parameters.SampleRate;
stim_similarity_thresh = 1.2;  
stim_durations = parameters.RailsDurations;
stim_durations = stim_durations.*parameters.SampleRate; %put durations in frames

%%%velocity analysis
velocity_time_window_before = time_window_before;
velocity_time_window_after = time_window_after;
n_bins = 40;

edges = linspace(-0.2,0.2,n_bins);
display_field_name = 'Velocity';
top_percentile_cutoff = 95;
if use_max_velocity
    boxcar_window = ones(1,time_window_before) ./ (time_window_before);
else
	boxcar_window = ones(1,time_window_before+time_window_after+1) ./ (time_window_before+time_window_after+1);
end

number_of_behaviors = num_velocity_behaviors;
saved_stimulus_count = 0;
saved_head_stimulus_intensities = [];
saved_head_stimulus_durations = [];
saved_tail_stimulus_intensities = [];
saved_tail_stimulus_durations = [];
saved_stimulus_delays = [];
saved_head_or_tail_first = [];

%%%for every combination of possible stimulation category, make an entry
for head_intensity_index = 1:numel(parameters.RailsIntensities)
   for head_duration_index = 1:numel(parameters.RailsDurations)
       for tail_intensity_index = 1:numel(parameters.RailsIntensities)
           for tail_duration_index = 1:numel(parameters.RailsDurations)
               for delay_index = 1:numel(parameters.RailsDelays)
                   for head_or_tail_first = [true false]
                       if parameters.RailsDelays(delay_index) == 0 && ~head_or_tail_first
                           %%% head always goes first if there is no delay
                           continue
                       end
                       saved_stimulus_count = saved_stimulus_count + 1;
                       saved_head_stimulus_intensities(saved_stimulus_count) = parameters.RailsIntensities(head_intensity_index);
                       saved_head_stimulus_durations(saved_stimulus_count) = parameters.RailsDurations(head_duration_index) * parameters.SampleRate;
                       saved_tail_stimulus_intensities(saved_stimulus_count) = parameters.RailsIntensities(tail_intensity_index);
                       saved_tail_stimulus_durations(saved_stimulus_count) = parameters.RailsDurations(tail_duration_index) * parameters.SampleRate;
                       saved_stimulus_delays(saved_stimulus_count) = parameters.RailsDelays(delay_index) * parameters.SampleRate;
                       saved_head_or_tail_first(saved_stimulus_count) = head_or_tail_first;

                       all_behavior_transitions_for_frame{saved_stimulus_count} = cell(1,total_window_frames);
                       all_behavior_annotations_for_frame{saved_stimulus_count} = cell(1,total_window_frames);
                       velocities{saved_stimulus_count} = [];
                       video_possible_tracks{saved_stimulus_count} = [];
                       video_possible_frames{saved_stimulus_count} = [];
                   end
               end
           end
       end
   end
end
video_duration = total_window_frames; % total duration in seconds


%% allow user to select the folder to save as
if plot_video
    pathname = uigetdir('', 'Select Video Output Folder')
    if isequal(pathname,0)
        %cancel
       return
    end
end

%% loop through every stimulus event
allTracks = [];
stim_peaks_folder=[];  %%% to count the number of all stimulus events

for folder_index = 1:length(folders)
    sprintf('Analyzing folder %d out of %d folders', folder_index, length(folders))
    %%%load the tracks for this folder
    [current_tracks, ~, ~] = loadtracks(folders{folder_index},relevant_track_fields);
    
    %%%delete worms that are not long enough, also add folder info
    track_indecies_to_delete = [];
    for track_index = 1:length(current_tracks)
        if mean(current_tracks(track_index).Length) < min_avg_worm_length
            track_indecies_to_delete = [track_indecies_to_delete, track_index];
        else
            current_tracks(track_index).folder_index = folder_index;
            current_tracks(track_index).within_folder_track_index = track_index;
        end
    end
    current_tracks(track_indecies_to_delete) = [];
    if isempty(current_tracks)
        continue
    end
    
    current_param = load_parameters(folders{folder_index});
    for track_index = 1:length(current_tracks)
        %%%process the stimulus going from n centerline points to a single
        %%%point by taking the median
        nd_sampled_stimulus = centerline_sampled_stimulus_to_nd_sampled_stimulus(current_tracks(track_index).AlignedStimulus, worm_region_count)';
        head_stimulus = nd_sampled_stimulus(worm_region_head,:);
        tail_stimulus = nd_sampled_stimulus(worm_region_tail,:);
        
        single_dimension_stimulus = abs(head_stimulus) + abs(tail_stimulus);
        single_dimension_stimulus = double((single_dimension_stimulus - median(single_dimension_stimulus)) > 0);
        [~, stim_peaks, widths, ~] = findpeaks(single_dimension_stimulus, 'MinPeakDistance',round(min(current_param.InterTriggerInterval*current_param.SampleRate/stim_similarity_thresh,numel(current_tracks(track_index).Frames)/2-1)));
        stim_peaks = [stim_peaks, numel(single_dimension_stimulus)];
        stim_peaks_folder=[stim_peaks_folder stim_peaks];
        
        for peak_index = 1:numel(stim_peaks)-1
            if stim_peaks(peak_index)-time_window_before+1 > 0 && stim_peaks(peak_index)+time_window_after <= numel(current_tracks(track_index).Frames)
                %%% make sure the track is entirely in the window
                
                %%%scan the head and tail stimuli from the element before the peak for 
                [head_stim_frame, head_stim_intensity, head_stim_stim_duration] = read_next_rails_stim(head_stimulus, stim_peaks(peak_index)-1, stim_peaks(peak_index+1));
                [tail_stim_frame, tail_stim_intensity, tail_stim_stim_duration] = read_next_rails_stim(tail_stimulus, stim_peaks(peak_index)-1, stim_peaks(peak_index+1));
                
                %%%find the delay between head and tail stimuli
                head_or_tail_first = tail_stim_frame >= head_stim_frame;
                stimulus_delay = abs(head_stim_frame - tail_stim_frame);

                %%% match this instance with the stimulus possibilities we expect
                matching_head_intensities = or(and(saved_head_stimulus_intensities <= head_stim_intensity*stim_similarity_thresh, saved_head_stimulus_intensities >= head_stim_intensity/stim_similarity_thresh), ...
                    and(-saved_head_stimulus_intensities <= head_stim_intensity*-stim_similarity_thresh, -saved_head_stimulus_intensities >= head_stim_intensity/-stim_similarity_thresh));
                matching_head_durations = and(head_stim_stim_duration >= saved_head_stimulus_durations./stim_similarity_thresh, head_stim_stim_duration <= saved_head_stimulus_durations.*stim_similarity_thresh);
                matching_tail_intensities = or(and(saved_tail_stimulus_intensities <= tail_stim_intensity*stim_similarity_thresh, saved_tail_stimulus_intensities >= tail_stim_intensity/stim_similarity_thresh), ...
                    and(-saved_tail_stimulus_intensities <= tail_stim_intensity*-stim_similarity_thresh, -saved_tail_stimulus_intensities >= tail_stim_intensity/-stim_similarity_thresh));
                matching_tail_durations = and(tail_stim_stim_duration >= saved_tail_stimulus_durations./stim_similarity_thresh, tail_stim_stim_duration <= saved_tail_stimulus_durations.*stim_similarity_thresh);

                matching_stimulus_delays = and(stimulus_delay >= saved_stimulus_delays./stim_similarity_thresh, stimulus_delay <= saved_stimulus_delays.*stim_similarity_thresh);
                matching_head_or_tail_first = head_or_tail_first == saved_head_or_tail_first;

                if saved_stimulus_count > 0
                    current_stim_index = find(matching_head_intensities & matching_head_durations & ...
                        matching_tail_intensities & matching_tail_durations & ...
                        matching_stimulus_delays & matching_head_or_tail_first,1);
                else
                    current_stim_index = [];
                end
                if isempty(current_stim_index)
                    %%%no entry, skip it
                    continue
                end
                
                %%%keep track of velocities
                velocities{current_stim_index} = [velocities{current_stim_index}; current_tracks(track_index).Velocity(stim_peaks(peak_index)-velocity_time_window_before+1:stim_peaks(peak_index)+velocity_time_window_after)'];
                
                for frame_shift = -time_window_before:time_window_after
                    current_frame = stim_peaks(peak_index) + frame_shift;
                    if current_frame <= length(current_tracks(track_index).Frames) && current_frame >= 1
                        %%%make sure the current frame is in range
                        %%%cut up tracks to each frame
                        all_behavior_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1} = [all_behavior_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1}, current_tracks(track_index).VelocityBehavior(current_frame)];
                    end
                end
                %%%add it to the possible video plotting with conditions
                video_possible_frames{current_stim_index} = [video_possible_frames{current_stim_index}, stim_peaks(peak_index)];
                video_possible_tracks{current_stim_index} = [video_possible_tracks{current_stim_index}, (length(allTracks) + track_index)];
            end
        end
    end
    allTracks = [allTracks, current_tracks];
end

n_sti=saved_stimulus_count;
behavior_counts_for_frame = zeros(number_of_behaviors,n_sti,total_window_frames);
behavior_ratios_for_frame = zeros(number_of_behaviors,n_sti,total_window_frames);
bootstrapped_behavioral_ratios_for_frame = zeros(number_of_behaviors,n_sti,total_window_frames, BOOTSTRAP_N);

for stimulus_index = 1:n_sti
    %%% plot the transition rates centered on stim delivery
    total_counts_for_frame = zeros(1,total_window_frames);
    for frame_index = 1:total_window_frames
        for behavior_index = 1:number_of_behaviors
            behavior_counts_for_frame(behavior_index,stimulus_index,frame_index) = sum(all_behavior_annotations_for_frame{stimulus_index}{frame_index}==behavior_index);
        end
        %%%get ratio
        behavior_ratios_for_frame(:,stimulus_index,frame_index) = behavior_counts_for_frame(:,stimulus_index,frame_index)./sum(behavior_counts_for_frame(:,stimulus_index,frame_index)); 
        
        %%%get bootstrapped ratios
        for bootstrap_index = 1:BOOTSTRAP_N
            bootstrapped_behavioral_observations = datasample(all_behavior_annotations_for_frame{stimulus_index}{frame_index},numel(all_behavior_annotations_for_frame{stimulus_index}{frame_index}));
            for behavior_index = 1:number_of_behaviors
                bootstrapped_behavioral_ratios_for_frame(behavior_index,stimulus_index,frame_index,bootstrap_index) = sum(bootstrapped_behavioral_observations==behavior_index) / numel(bootstrapped_behavioral_observations);
            end
        end
    end
end
%% 2 plot the behavioral ratios as a function of time
n_tracks=zeros(1,n_sti); %number of tracks in each sti intensities
my_colors = velocity_based_behavior_colors;
my_behavior_names = velocity_based_behavior_names;
for stimulus_index = 1:n_sti
    n_tracks(stimulus_index) = round(mean(arrayfun(@(x) size(x{1},2), [all_behavior_transitions_for_frame{1,stimulus_index}])));
    if n_tracks <= 5
       continue 
    end

    figure
    hold on;grid on
    head_stim_y_location = 1.1;
    tail_stim_y_location = 1.0;
    %%%head stimulus color
    if saved_head_stimulus_intensities(stimulus_index) < 0
        %%%blue
        head_stimulus_color = 1 - (abs(saved_head_stimulus_intensities(stimulus_index)) / parameters.avgPowerBlue);
        head_stimulus_color = [head_stimulus_color head_stimulus_color 1];
    else
        %%%red
        head_stimulus_color = 1 - (saved_head_stimulus_intensities(stimulus_index) / parameters.avgPowerRed);
        head_stimulus_color = [1 head_stimulus_color head_stimulus_color];
    end
    %%%tail stimulus color
    if saved_tail_stimulus_intensities(stimulus_index) < 0
        %%%blue
        tail_stimulus_color = 1 - (abs(saved_tail_stimulus_intensities(stimulus_index)) / parameters.avgPowerBlue);
        tail_stimulus_color = [tail_stimulus_color tail_stimulus_color 1];
    else
        %%%red
        tail_stimulus_color = 1 - (saved_tail_stimulus_intensities(stimulus_index) / parameters.avgPowerRed);
        tail_stimulus_color = [1 tail_stimulus_color tail_stimulus_color];
    end
    text(0, head_stim_y_location, 'Head ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    text(0, tail_stim_y_location, 'Tail ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    if saved_head_or_tail_first(stimulus_index)
        %%%head is first
        rectangle('Position',[0 head_stim_y_location saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',head_stimulus_color)
        rectangle('Position',[saved_stimulus_delays(stimulus_index)/parameters.SampleRate tail_stim_y_location saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',tail_stimulus_color)
        head_or_tail_string = 'Head Stimulus First';
    else
        %%%tail is first
        rectangle('Position',[0 tail_stim_y_location saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',tail_stimulus_color)
        rectangle('Position',[saved_stimulus_delays(stimulus_index)/parameters.SampleRate head_stim_y_location saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',head_stimulus_color)
        head_or_tail_string = 'Tail Stimulus First';
    end
    
    for behavior_index = 1:number_of_behaviors
        plot(-time_window_before/fps:1/fps:time_window_after/fps, squeeze(behavior_ratios_for_frame(behavior_index,stimulus_index,:)), '-', 'color', my_colors(behavior_index,:),'Linewidth', 3,'DisplayName',my_behavior_names{behavior_index});
    end
    hold off
    xlabel('Time (s)') % x-axis label
    ylabel('Behavioral Ratio') % y-axis label
% % %     title({head_or_tail_string, ...
% % %          ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          ['Head Stimulus Duration = ', num2str(saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          [' Tail Stimulus Duration = ', num2str(saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Delay Duration = ', num2str(saved_stimulus_delays(stimulus_index)/parameters.SampleRate), 's'] ...
% % %          ['(n = ', num2str(n_tracks(stimulus_index)), ' events)']});

    title({head_or_tail_string, ...
         ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
         ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2']});
    ax = gca;
    ax.FontSize = 10;
    axis([-time_window_before/fps time_window_after/fps 0  head_stim_y_location*2-tail_stim_y_location])

end

%% Plot the before/after bar charts
my_colors = velocity_based_behavior_colors;
my_behavior_names = velocity_based_behavior_names;
for stimulus_index = 1:n_sti
    %%% plot the change in fraction of animals before or after the stimulus
    if saved_head_or_tail_first(stimulus_index)
        %%%head is first
        total_stim_duration = max(saved_head_stimulus_durations(stimulus_index), saved_tail_stimulus_durations(stimulus_index)+saved_stimulus_delays(stimulus_index));
        head_or_tail_string = 'Head Stimulus First';
    else
        %%%tail is first
        total_stim_duration = max(saved_head_stimulus_durations(stimulus_index)+saved_stimulus_delays(stimulus_index), saved_tail_stimulus_durations(stimulus_index));
        head_or_tail_string = 'Tail Stimulus First';
    end
    
% % % %     %%% if plotting only the four plots shown in fig. 2m-p
% % % %     if stimulus_index~=[1 5 21 25]
% % % %         continue
% % % %     end
         
    figure('position', [0 0 400 400])
    hold on

    ax = gca;
    ax.FontSize = 10;
    
    %%%keep track the bars
    p = zeros(1, number_of_behaviors); % the p-values
    bar_vals = zeros(number_of_behaviors,2);
    error_vals = zeros(number_of_behaviors*2,2);
    for behavior_index = 1:number_of_behaviors
        %%%before
        bar_vals(behavior_index,1) = behavior_ratios_for_frame(behavior_index,stimulus_index,velocity_time_window_before-(2*parameters.SampleRate));  
        %%%after
        bar_vals(behavior_index,2) = behavior_ratios_for_frame(behavior_index,stimulus_index,velocity_time_window_before+total_stim_duration);
        
        %%%construct 95% confidence interval for the before and after from
        %%%bootstrapped values
        error_vals_index = (behavior_index-1)*2+1;
        error_vals(error_vals_index,1) = abs(bar_vals(behavior_index,1) - prctile(squeeze(bootstrapped_behavioral_ratios_for_frame(behavior_index,stimulus_index,velocity_time_window_before-(2*parameters.SampleRate),:)),0.025));
        error_vals(error_vals_index,2) = abs(bar_vals(behavior_index,2) - prctile(squeeze(bootstrapped_behavioral_ratios_for_frame(behavior_index,stimulus_index,velocity_time_window_before+total_stim_duration,:)),0.025));
        error_vals(error_vals_index+1,1) = abs(prctile(squeeze(bootstrapped_behavioral_ratios_for_frame(behavior_index,stimulus_index,velocity_time_window_before-(2*parameters.SampleRate),:)), 0.975) - bar_vals(behavior_index,1));
        error_vals(error_vals_index+1,2) = abs(prctile(squeeze(bootstrapped_behavioral_ratios_for_frame(behavior_index,stimulus_index,velocity_time_window_before+total_stim_duration,:)), 0.975) - bar_vals(behavior_index,2));
        
        %%%determine if the difference is significant using wilcoxon ranksum
        binary_observations_for_behavior_before = double(all_behavior_annotations_for_frame{stimulus_index}{velocity_time_window_before-(2*parameters.SampleRate)}==behavior_index);
        binary_observations_for_behavior_after = double(all_behavior_annotations_for_frame{stimulus_index}{velocity_time_window_before+total_stim_duration}==behavior_index);
        p(behavior_index) = ranksum(binary_observations_for_behavior_before,binary_observations_for_behavior_after);
    end
    
    b = bar(bar_vals, 'group','linewidth', 2);
    %%% Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(bar_vals);
    %%% Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    %%% Set the position of each error bar in the centre of the main bar
    %%% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    all_x = zeros(2,ngroups);
    for bar_i = 1:2
        %%%before and after
        x = (1:ngroups) - groupwidth/2 + (2*bar_i-1) * groupwidth / (2*nbars);
        all_x(bar_i,:) = x;
        %%%make the bar graph appropriate colors
        if bar_i == 1
            %%%before
            b(bar_i).FaceColor = [1 1 1];
        else
            %%%after
            b(bar_i).FaceColor = [0 0 0];
        end
        
        for behavior_index = 1:number_of_behaviors
            hErr = errorbar(x(behavior_index), bar_vals(behavior_index,bar_i), error_vals((bar_i-1)*2+behavior_index,1), error_vals((bar_i-1)*2+behavior_index,2),...
            'color',my_colors(behavior_index,:), 'Linewidth', 2);
            
            %%%resize the error bars
            myerrorbar = hErr.Bar;                           % hidden property/handle
            drawnow                                 % populate b's properties
            vd = myerrorbar.VertexData;
            N = 1;                           % number of error bars
            newLength = 0.05;
            leftInds = N*2+1:2:N*6;
            rightInds = N*2+2:2:N*6;
            vd(1,leftInds,1) = [x(behavior_index)-newLength, x(behavior_index)-newLength];
            vd(1,rightInds,1) = [x(behavior_index)+newLength, x(behavior_index)+newLength];
            myerrorbar.VertexData = vd;
        end
    end
    
    siggroups = cell(1,number_of_behaviors);
    for behavior_index = 1:number_of_behaviors
       %%%plot star for significance
       siggroups{behavior_index} = squeeze(all_x(:,behavior_index));
    end
    p(p>0.05) = nan;
    sigstar(siggroups, p) %plot the significance
   
    axis([0 number_of_behaviors+1 0 0.8])
    set(gca, 'XTick', 1:number_of_behaviors)
    set(gca, 'XTickLabel', my_behavior_names);
    ylabel('Fraction of Animals') % y-axis label
% % %     title({head_or_tail_string, ...
% % %          ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          ['Head Stimulus Duration = ', num2str(saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          [' Tail Stimulus Duration = ', num2str(saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Delay Duration = ', num2str(saved_stimulus_delays(stimulus_index)/parameters.SampleRate), 's'] ...
% % %          ['(n = ', num2str(n_tracks(stimulus_index)), ' events)']});

    title({['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
         ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2']});
    hold off 
end

%% velocity histogram vs time
time_axis = -velocity_time_window_before+1:velocity_time_window_after;
time_axis = time_axis/parameters.SampleRate;
velocity_density = zeros(numel(edges)-1, numel(time_axis));
for stimulus_index = 1:n_sti
    if isempty(velocities{stimulus_index})
        continue
    end
    for time_index = 1:length(time_axis)
        velocity_density(:,time_index) = histcounts(velocities{stimulus_index}(:,time_index), edges);
    end
    mean_velocities = mean(velocities{stimulus_index},1);

    figure
    hold on
    imagesc(time_axis, edges, velocity_density)
    set(gca, 'YDir', 'normal')
    plot(time_axis,mean_velocities, '-k', 'linewidth', 3)
    head_stim_y_location = edges(end) + (abs(edges(end)-edges(1))/10);
    tail_stim_y_location = edges(end);
    %%%head stimulus color
    if saved_head_stimulus_intensities(stimulus_index) < 0
        %%%blue
        head_stimulus_color = 1 - (abs(saved_head_stimulus_intensities(stimulus_index)) / parameters.avgPowerBlue);
        head_stimulus_color = [head_stimulus_color head_stimulus_color 1];
    else
        %%%red
        head_stimulus_color = 1 - (saved_head_stimulus_intensities(stimulus_index) / parameters.avgPowerRed);
        head_stimulus_color = [1 head_stimulus_color head_stimulus_color];
    end
    %%%tail stimulus color
    if saved_tail_stimulus_intensities(stimulus_index) < 0
        %%%blue
        tail_stimulus_color = 1 - (abs(saved_tail_stimulus_intensities(stimulus_index)) / parameters.avgPowerBlue);
        tail_stimulus_color = [tail_stimulus_color tail_stimulus_color 1];
% % %         tail_stimulus_color = [1 1 1];
    else
        %%%red
        tail_stimulus_color = 1 - (saved_tail_stimulus_intensities(stimulus_index) / parameters.avgPowerRed);
        tail_stimulus_color = [1 tail_stimulus_color tail_stimulus_color];
    end
    text(0, head_stim_y_location, 'Head ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    text(0, tail_stim_y_location, 'Tail ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    if saved_head_or_tail_first(stimulus_index)
        %head is first
        rectangle('Position',[0 head_stim_y_location saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',head_stimulus_color)
        rectangle('Position',[saved_stimulus_delays(stimulus_index)/parameters.SampleRate tail_stim_y_location saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',tail_stimulus_color)
        head_or_tail_string = 'Head Stimulus First';
    else
        %tail is first
        rectangle('Position',[0 tail_stim_y_location saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',tail_stimulus_color)
        rectangle('Position',[saved_stimulus_delays(stimulus_index)/parameters.SampleRate head_stim_y_location saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',head_stimulus_color)
        head_or_tail_string = 'Tail Stimulus First';
    end
    axis([-time_window_before/fps time_window_after/fps edges(1) head_stim_y_location*2-tail_stim_y_location])
    xlabel('Time (s)')
    ylabel('Velocity (mm/s)')
% % %     title({head_or_tail_string, ...
% % %          ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          ['Head Stimulus Duration = ', num2str(saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          [' Tail Stimulus Duration = ', num2str(saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Delay Duration = ', num2str(saved_stimulus_delays(stimulus_index)/parameters.SampleRate), 's'] ...
% % %          ['(n = ', num2str(n_tracks(stimulus_index)), ' events)']});

    title({head_or_tail_string, ...
         ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
         ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2']});
end

%% make sorted velocity traces and plot the movies
my_colors = velocity_based_behavior_colors;
my_behavior_names = velocity_based_behavior_names;
for stimulus_index = 1:n_sti
    
% % %     %%% uncomment the following to just plot the figures shown in paper
% % %     if stimulus_index~=[1 5 21 25]
% % %         continue
% % %     end

    if isempty(velocities{stimulus_index})
        continue
    end

    current_max_traces = min(N_rows*N_columns, size(velocities{stimulus_index},1));
    random_order = randperm(size(velocities{stimulus_index},1));
    velocities_to_display = velocities{stimulus_index}(random_order(1:current_max_traces),:);
    
    selected_index = 4; %indicates which individual worm trace to plot
    %%%change the selected index for special cases of 0 and tail stimulation
    if round(saved_head_stimulus_intensities(stimulus_index)) == 0 && round(saved_tail_stimulus_intensities(stimulus_index)) == 0
        %%% zero head and zeo tail stimulus
        selected_index = floor(current_max_traces/2);
    elseif round(saved_head_stimulus_intensities(stimulus_index)) == 0 && round(saved_tail_stimulus_intensities(stimulus_index)) > 0
        %%% tail only
        selected_index = current_max_traces-3;
    end
    
    %%%find total stimulus duration
    if saved_head_or_tail_first(stimulus_index)
        %%%head goes first
        total_stimulus_duration = max(saved_head_stimulus_durations(stimulus_index), saved_tail_stimulus_durations(stimulus_index) + saved_stimulus_delays(stimulus_index));
    else
        %%%tail goes first
        total_stimulus_duration = max(saved_tail_stimulus_durations(stimulus_index), saved_head_stimulus_durations(stimulus_index) + saved_stimulus_delays(stimulus_index));
    end
    
    mean_velocity_before_stimulus = mean(velocities_to_display(:,1:velocity_time_window_before),2);
    mean_velocity_during_stimulus = mean(velocities_to_display(:,velocity_time_window_before+1:total_stimulus_duration+velocity_time_window_before+1),2);
   
    %%%sort by the velocity during stimuation
    [~,sorted_index] = sort(mean_velocity_during_stimulus);

    velocities_to_display = velocities_to_display(sorted_index,:);
    time_axis = -velocity_time_window_before:velocity_time_window_after;
    time_axis = time_axis/parameters.SampleRate;
    figure
    ax = gca;
    ax.Clipping = 'off';
    hold on
    h = imagesc(flipud(velocities_to_display));
    set(h, 'XData', time_axis);
    colormap redblue
    caxis([edges(1) edges(end)])
    
    head_stim_y_location = current_max_traces + 2;
    tail_stim_y_location = current_max_traces + 1;    %head stimulus color
    if saved_head_stimulus_intensities(stimulus_index) < 0
        %%%blue
        head_stimulus_color = 1 - (abs(saved_head_stimulus_intensities(stimulus_index)) / parameters.avgPowerBlue);
        head_stimulus_color = [head_stimulus_color head_stimulus_color 1];
    else
        %%%red
        head_stimulus_color = 1 - (saved_head_stimulus_intensities(stimulus_index) / parameters.avgPowerRed);
        head_stimulus_color = [1 head_stimulus_color head_stimulus_color];
    end
    %%%tail stimulus color
    if saved_tail_stimulus_intensities(stimulus_index) < 0
        %%%blue
        tail_stimulus_color = 1 - (abs(saved_tail_stimulus_intensities(stimulus_index)) / parameters.avgPowerBlue);
        tail_stimulus_color = [tail_stimulus_color tail_stimulus_color 1];
% % %         tail_stimulus_color = [1 1 1];
    else
        %%%red
        tail_stimulus_color = 1 - (saved_tail_stimulus_intensities(stimulus_index) / parameters.avgPowerRed);
        tail_stimulus_color = [1 tail_stimulus_color tail_stimulus_color];
    end
    text(0, head_stim_y_location, 'Head ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    text(0, tail_stim_y_location, 'Tail ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    if saved_head_or_tail_first(stimulus_index)
        %%%head is first
        rectangle('Position',[0 head_stim_y_location saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',head_stimulus_color)
        rectangle('Position',[saved_stimulus_delays(stimulus_index)/parameters.SampleRate tail_stim_y_location saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',tail_stimulus_color)
        head_or_tail_string = 'HeadStimulusFirst';
    else
        %%%tail is first
        rectangle('Position',[0 tail_stim_y_location saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',tail_stimulus_color)
        rectangle('Position',[saved_stimulus_delays(stimulus_index)/parameters.SampleRate head_stim_y_location saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',head_stimulus_color)
        head_or_tail_string = 'TailStimulusFirst';
    end
    %%%plot an arrow next to the selected trace
    anarrow = annotation('arrow');
    anarrow.Parent = gca;
    anarrow.Position = [-velocity_time_window_before/parameters.SampleRate-0.5, current_max_traces-selected_index+1, 0.5, 0];
    
    axis tight
    set(gca, 'YTick', fliplr(current_max_traces:-5:1))
    set(gca, 'YTickLabel', fliplr(1:5:current_max_traces))

    xlabel('Time (s)')
    ylabel('Velocity Traces (mm/s)')
% % %     title({head_or_tail_string, ...
% % %          ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          ['Head Stimulus Duration = ', num2str(saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          [' Tail Stimulus Duration = ', num2str(saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Delay Duration = ', num2str(saved_stimulus_delays(stimulus_index)/parameters.SampleRate), 's'] ...
% % %          ['(n = ', num2str(n_tracks(stimulus_index)), ' events)']});    

    title({head_or_tail_string, ...
         ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
         ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2']});  

     if plot_video
        possible_frames = video_possible_frames{stimulus_index}(random_order(sorted_index));
        possible_tracks = video_possible_tracks{stimulus_index}(random_order(sorted_index));
        saveFileName = fullfile(pathname,[head_or_tail_string, '_HeadStimulusIntensity_', num2str(round(saved_head_stimulus_intensities(stimulus_index))), ...
         '_HeadStimulusDuration_', num2str(saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate), 's', ...
         '_TailStimulusIntensity_', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), ...
         '_TailStimulusDuration_', num2str(saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate), 's', ...
         '_DelayDuration_', num2str(saved_stimulus_delays(stimulus_index)/parameters.SampleRate), 's', ...
         '_n_', num2str(n_tracks(stimulus_index))]);

        make_tiled_movies_given_instances(allTracks, folders, saveFileName, possible_tracks, possible_frames, N_rows, N_columns, video_duration, display_field_name, selected_index)
     end
    
    %%% generate single animal velocity trace
    figure
    ax = gca;
    ax.Clipping = 'off';
    time_axis = -velocity_time_window_before+1:velocity_time_window_after;
    time_axis = time_axis/parameters.SampleRate;
    z = zeros(size(time_axis));
    velocity_timeseires = velocities_to_display(selected_index,:);
    %%% Plot the line with width 8 so we can see the colors well.
    surface([time_axis;time_axis], [velocity_timeseires;velocity_timeseires], [z;z], [velocity_timeseires;velocity_timeseires],...
        'FaceColor', 'no',...
        'EdgeColor', 'interp',...
        'LineWidth', 8);
    colormap redblue
    caxis([edges(1) edges(end)])
    hold on
    plot(time_axis, velocity_timeseires, 'k', 'LineWidth', 2)
    grid on;
    head_stim_y_location = edges(end) + (abs(edges(end)-edges(1))/10);
    tail_stim_y_location = edges(end);
    %%%head stimulus color
    if saved_head_stimulus_intensities(stimulus_index) < 0
        %%%blue
        head_stimulus_color = 1 - (abs(saved_head_stimulus_intensities(stimulus_index)) / parameters.avgPowerBlue);
        head_stimulus_color = [head_stimulus_color head_stimulus_color 1];
    else
        %%%red
        head_stimulus_color = 1 - (saved_head_stimulus_intensities(stimulus_index) / parameters.avgPowerRed);
        head_stimulus_color = [1 head_stimulus_color head_stimulus_color];
    end
    %%%tail stimulus color
    if saved_tail_stimulus_intensities(stimulus_index) < 0
        %%%blue
        tail_stimulus_color = 1 - (abs(saved_tail_stimulus_intensities(stimulus_index)) / parameters.avgPowerBlue);
        tail_stimulus_color = [tail_stimulus_color tail_stimulus_color 1];
%         tail_stimulus_color = [1 1 1];
    else
        %%%red
        tail_stimulus_color = 1 - (saved_tail_stimulus_intensities(stimulus_index) / parameters.avgPowerRed);
        tail_stimulus_color = [1 tail_stimulus_color tail_stimulus_color];
    end
    text(0, head_stim_y_location, 'Head ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    text(0, tail_stim_y_location, 'Tail ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    if saved_head_or_tail_first(stimulus_index)
        %%%head is first
        rectangle('Position',[0 head_stim_y_location saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',head_stimulus_color)
        rectangle('Position',[saved_stimulus_delays(stimulus_index)/parameters.SampleRate tail_stim_y_location saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',tail_stimulus_color)
        head_or_tail_string = 'Head Stimulus First';
    else
        %%%tail is first
        rectangle('Position',[0 tail_stim_y_location saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',tail_stimulus_color)
        rectangle('Position',[saved_stimulus_delays(stimulus_index)/parameters.SampleRate head_stim_y_location saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate abs(head_stim_y_location-tail_stim_y_location)],'FaceColor',head_stimulus_color)
        head_or_tail_string = 'Tail Stimulus First';
    end
    axis([-time_window_before/fps time_window_after/fps edges(1)-0.1 head_stim_y_location*2-tail_stim_y_location])
    xlabel('Time (s)')
    ylabel('Velocity (mm/s)')
% % %     set(gca, 'XTick', -time_window_before/fps:5:time_window_after/fps)
% % %     title({head_or_tail_string, ...
% % %          ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          ['Head Stimulus Duration = ', num2str(saved_head_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
% % %          [' Tail Stimulus Duration = ', num2str(saved_tail_stimulus_durations(stimulus_index)/parameters.SampleRate), 's'], ...
% % %          ['Delay Duration = ', num2str(saved_stimulus_delays(stimulus_index)/parameters.SampleRate), 's'] ...
% % %          ['Selected Index = ', num2str(selected_index)]});
    title({head_or_tail_string, ...
         ['Head Stimulus Intensity = ', num2str(round(saved_head_stimulus_intensities(stimulus_index))), 'uW/mm2'], ...
         ['Tail Stimulus Intensity = ', num2str(round(saved_tail_stimulus_intensities(stimulus_index))), 'uW/mm2']});
     set(gca, 'color', [0.5 0.5 0.5])
     hold off;
end

%% use behavioral annotations to determine decision making +
% % % plot response probability heatmap given head and tail intensities 
% % % (this is only compatible without delayed stimulus)

n_tracks = zeros(1,n_sti); %number of tracks in each sti intensities
transition_count = 0;

%%% record all transitions in a 3D matrix indexed from_behavior, to_behavior, stimulus_index
behavioral_transition_matrix_for_all_stimuli = zeros(numel(my_behavior_names),numel(my_behavior_names),n_sti);
behavioral_transition_stats = [];

for stimulus_index = 1:n_sti
    %%% loop through the variou stimulation conditions
    n_tracks(stimulus_index) = numel(video_possible_frames{stimulus_index});
    if n_tracks <= 5
        %%%skip analysis for conditions with few tracks
       continue 
    end

    %%%calculate total stimulus duration, which encompasses the delayed
    %%%second stimulus. The delay is not applicapble in msot cases
    if saved_head_or_tail_first(stimulus_index)
        %%%head is first
        total_stimulus_duration_for_stim = saved_stimulus_delays(stimulus_index) + saved_tail_stimulus_durations(stimulus_index);
        head_or_tail_string = 'HeadStimulusFirst';
    else
        %%%tail is first
        total_stimulus_duration_for_stim = saved_stimulus_delays(stimulus_index) + saved_head_stimulus_durations(stimulus_index);
        head_or_tail_string = 'TailStimulusFirst';
    end

    %%%determine the behavior before stimulus and then subsequent transitions
    starting_count_for_stimulus = transition_count; 
    for stimulation_event_index = 1:numel(video_possible_frames{stimulus_index})
        %%%loop through all the animal-stimulation events for the stimulus condition
        stimulation_track_index = video_possible_tracks{stimulus_index}(stimulation_event_index);
        stimulation_frame_index = video_possible_frames{stimulus_index}(stimulation_event_index);
        behavioral_annotation_for_track = allTracks(stimulation_track_index).VelocityBehavior;
        %%%behavior before is the behavior a smoothing window before the start of the stimulus
        
        %%% code when I just want to know the starting and ending behavior state
        starting_behavior_annotations_for_stim = behavioral_annotation_for_track(stimulation_frame_index - (2*parameters.TrackingSmoothingWindow*parameters.SampleRate)); %%% we are looking at beh. state 2 sec before stim onset
        next_behavior_annotations_for_stim = behavioral_annotation_for_track(stimulation_frame_index + (parameters.TrackingSmoothingWindow*parameters.SampleRate));

        if starting_behavior_annotations_for_stim > 0 && next_behavior_annotations_for_stim > 0
            %%% save the transition stats only if conditions are satisfied
            transition_count = transition_count+1;
            behavioral_transition_stats(transition_count).worm_index = stimulation_track_index;
            behavioral_transition_stats(transition_count).stimulus_index = stimulus_index;
            behavioral_transition_stats(transition_count).starting_behavior = starting_behavior_annotations_for_stim;
            behavioral_transition_stats(transition_count).next_behavior = next_behavior_annotations_for_stim;
        end
    end
    
    %%%create the behavioral transition matrix
    behavioral_transition_matrix = zeros(numel(my_behavior_names),numel(my_behavior_names));
    for transition_index = starting_count_for_stimulus+1:transition_count
        behavioral_transition_matrix(behavioral_transition_stats(transition_index).starting_behavior,behavioral_transition_stats(transition_index).next_behavior) = ...
            behavioral_transition_matrix(behavioral_transition_stats(transition_index).starting_behavior,behavioral_transition_stats(transition_index).next_behavior) + 1;
    end
    behavioral_transition_matrix_for_all_stimuli(:,:,stimulus_index) = behavioral_transition_matrix;
end

% % %loop through the behaviors and plot the probability of behavioral
% % %response with respect to each stimulus intensity, i.e. heatmap
behavioral_response_for_all_stimuli = sum(behavioral_transition_matrix_for_all_stimuli,1); %display only the result probability
ax31=figure('Renderer', 'painters', 'Position', [440 290 798 600]);
for behavioral_index = 1:numel(my_behavior_names)
    response_heatmap_for_behavior = zeros(numel(parameters.RailsIntensities),numel(parameters.RailsIntensities));
    for stimulus_index = 1:n_sti
        head_intensity_index = find(parameters.RailsIntensities == saved_head_stimulus_intensities(stimulus_index));
        tail_intensity_index = find(parameters.RailsIntensities == saved_tail_stimulus_intensities(stimulus_index));            
        response_heatmap_for_behavior(head_intensity_index,tail_intensity_index) = behavioral_response_for_all_stimuli(1,behavioral_index,stimulus_index) ./ squeeze(sum(behavioral_response_for_all_stimuli(:,:,stimulus_index)));
    end

% %     figure
    subplot(2,2,behavioral_index)
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
% %     xlabel('Tail Stimulus Intensity (uW/mm2)')
% %     ylabel('Head Stimulus Intensity (uW/mm2)')
    set(gca, 'XTick', 1:numel(parameters.RailsIntensities))
    set(gca, 'XTickLabel', abs(round(parameters.RailsIntensities)))
    set(gca, 'YTick', 1:numel(parameters.RailsIntensities))
    set(gca, 'YTickLabel', abs(round(parameters.RailsIntensities)))
    axis tight
    caxis([0 max(response_heatmap_for_behavior(:))])
%         caxis([0 0.8])
    colorbar
    colormap(othercolor('OrRd9'))
% %     caxis([0 0.9])
% %     chart_title = ['Probability of ', my_behavior_names{behavioral_index}, ' response'];
    chart_title = [my_behavior_names{behavioral_index}];
    title(chart_title)
    hold off;
    
    %% Code for the fit of probabilities
    z=response_heatmap_for_behavior; %Probability
    x=repmat([0, 20, 40, 60, 80],[ 5,1]);
    y=x';
% %     disp('plane fits for behavior index:')
    behavioral_index;
    [sf, gof] = fit([reshape(x.',1,[])', reshape(y.',1,[])'],reshape(z.',1,[])','poly11');
    z_intercept(behavioral_index) = sf.p00;
    xcoeff(behavioral_index) = sf.p10;
    ycoeff(behavioral_index) = sf.p01;
    ci(:,:,behavioral_index) = confint(sf); % three columns: for z_intercept, xcoeff, ycoeff
    
% %     figure
% %     plot(sf)
% %     hold on;
% %     surf(x,y,z)
% %     title(chart_title)
end

%%% Plot plane fit coefficients
B=1; %backward
P=2;
F=3;
S=4; %sprint
VECTORS = [xcoeff(B),  xcoeff(S), xcoeff(F), xcoeff(P);
          ycoeff(B), ycoeff(S), ycoeff(F), ycoeff(P)];
figure;
hold on;

% % xlim(1.2*[-max(max(abs(VECTORS))), max(max(abs(VECTORS)))])
% % ylim(1.2*[-max(max(abs(VECTORS))), max(max(abs(VECTORS)))])
xlim([-0.012 0.012])
ylim([-0.012 0.012])

yneg= VECTORS(2,:) - [ci(1,3,B), ci(1,3,S), ci(1,3,F), ci(1,3,P)]; 
ypos = VECTORS(2,:) - [ci(2,3,B), ci(2,3,S), ci(2,3,F), ci(2,3,P)];
xneg = VECTORS(1,:) - [ci(1,2,B), ci(1,2,S), ci(1,2,F), ci(1,2,P)];
xpos = VECTORS(1,:) - [ci(2,2,B), ci(2,2,S), ci(2,2,F), ci(2,2,P)];
errorbar(VECTORS(1,:),VECTORS(2,:),yneg,ypos,xneg,xpos,'ok', 'LineWidth',3,'HandleVisibility','off')

hold on;
plot([0,  xcoeff(B)], [0, ycoeff(B)], 'color', my_colors(1,:),'LineWidth', 3)
plot([0,  xcoeff(P)], [0, ycoeff(P)], 'color', my_colors(2,:),'LineWidth', 3)
plot([0,  xcoeff(F)], [0, ycoeff(F)], 'color', my_colors(3,:),'LineWidth', 3)
plot([0,  xcoeff(S)], [0, ycoeff(S)], 'color', my_colors(4,:),'LineWidth', 3)

axis square
grid on
xticks([-.01, 0, .01])
yticks([-.01, 0, .01])
set(gca,'fontsize', 16);
legend({'Reversals','Pause','Forward','Sprint'})
xlabel('\partialP/\partialI Tail Stim (mm^2/uW)')
ylabel('\partialP/\partialI Head Stim (mm^2/uW)')
% % title(['Initial behavior state: ', my_behavior_names{initial_behavior_of_interest}])
%%%

%% create a an array of pi charts
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
ylab = 'Head Stimulus Intensity (uW/mm2)';
lab = 1;
bubblepie(x,y,s,graph_data,graph_labels,graph_legend,xlab,ylab,lab)
colormap(my_colors)
ax = gca;
for intensity_index = 1:numel(parameters.RailsIntensities)
    ax_labels{intensity_index} = num2str(round(parameters.RailsIntensities(intensity_index)));
end
set(gca,'xtick',1:numel(parameters.RailsIntensities),'XTickLabel',ax_labels)
set(gca,'ytick',1:numel(parameters.RailsIntensities),'YTickLabel',ax_labels)

%%%%% This portion of the code is used to generate figure S1 and figure S2.
%%%%% For thsi portion of the code use the complete dataset (select "All tags" at the very first step)
%% get closed loop lags
lags = [];
experiment_drifts = [];
for folder_index = 1:length(folders)
    folder_name = folders{folder_index};
    loaded_variable = load([folder_name, filesep, 'timestamps.mat']);
    processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;
    camera_frames = 1:numel(processed_decoded_camera_frames);
    lags = [lags, camera_frames - double(processed_decoded_camera_frames)];
    experiment_drifts = [experiment_drifts, loaded_variable.experimentmaxdrift];
end
bin_edges = 0:10;
lag_prob = histcounts(lags,bin_edges) / numel(lags);

%%
ax1=figure('Renderer', 'painters', 'Position', [440 290 694 520]);
plot(bin_edges(1:end-1), lag_prob , '-o','LineWidth',2,'MarkerFaceColor',[0 0.45,0.74],'MarkerSize',10,'LineWidth',2)
for label_index = 1:numel(bin_edges)-2
    text(bin_edges(1+label_index),lag_prob(1+label_index), [num2str(round(lag_prob(1+label_index)*100,3)), '%'],'VerticalAlignment','bottom','HorizontalAlignment','left')
end
% % % set(gca,'YScale','log')
xlabel('Round-trip Latency (frames)')
ylabel('Probability')

ax1=gca;
set(ax1, 'XTick', bin_edges(1:end-1))

ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% % % set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
XOppTickLabels = round(bin_edges(1:end-1) * (1000/parameters.SampleRate));
% % % Set the x-tick and y-tick  labels for the second axes
set(ax2, 'XTickLabel', XOppTickLabels);
set(ax2, 'YTickLabel', []);
% % % ax2.FontSize = 14;
% % % ax1.FontSize=14;
% % % 
% % % xlim([0 9])
xlabel(ax2,'Round-trip Latency (ms)')
%%
%%%% plot the drifts
figure3=figure;
axes3 = axes('Parent',figure3);
bin_edges = 0:10:200;
experiment_drift_prob = histcounts(experiment_drifts*1000,bin_edges);
bar(bin_edges(1:end-1), experiment_drift_prob , 'hist')
yticks([0:25:50])
set(axes3,'YTickLabel',...
    {'0','25','50'});
% % % for label_index = 1:numel(bin_edges)-2
% % %     text(bin_edges(1+label_index),experiment_drift_prob(1+label_index), [num2str(round(experiment_drift_prob(1+label_index)*100,3)), '%'],'VerticalAlignment','top','HorizontalAlignment','left')
% % % end
xlim([0 150])
ylim([0 50])
xlabel('Spatial Drift (um)')
ylabel('Experiment Plate Count')
ax = gca;
ax.FontSize = 14;
box off;

%% tracking frame rates
MAX_LAG = 15;
track_count_vs_tracking_lag = zeros(255,MAX_LAG); %max 255 tracks, 100 frames lag


for folder_index = 1:length(folders)
    folder_name = folders{folder_index};
    loaded_variable = struct2cell(load([folder_name, filesep, 'labview_tracks.mat']));
    Tracks = [loaded_variable{:}];
    if isempty(Tracks)
        return
    end
    TrackFramesDiff = [];
    TrackFrames = [];
    [~,sort_idx] = sort([Tracks.WormIndex]);
    Tracks = Tracks(sort_idx);
    tracks_to_throw_out = []; %throw out tracks with less than 3 datapoints
    for track_index = 1:length(Tracks)
        if length(Tracks(track_index).Frames) < 3
            tracks_to_throw_out = [tracks_to_throw_out, track_index];
        end
    end
    Tracks(tracks_to_throw_out) = [];
    for track_index = 1:numel(Tracks)
        TrackFramesDiff = [TrackFramesDiff, diff(Tracks(track_index).Frames)'];
        TrackFrames = [TrackFrames, (Tracks(track_index).Frames(2:end))'];
    end
    TrackFramesDiff = double(TrackFramesDiff);
    TrackFramesDiff(TrackFramesDiff<0) = 1;
    
    %%% plot # of worms vs tracking frame rate
    [sorted_all_track_frames,sorted_indecies] = sort(TrackFrames);
    sorted_all_track_frames_diff = TrackFramesDiff(sorted_indecies);
    
    current_frame = sorted_all_track_frames(1);
    starting_all_frame_index = 1;
    for all_frame_index = 1:numel(sorted_all_track_frames)
        if sorted_all_track_frames(all_frame_index) > current_frame
            track_count_vs_tracking_lag(all_frame_index-starting_all_frame_index, min(MAX_LAG,sorted_all_track_frames_diff(starting_all_frame_index))) = ...
                track_count_vs_tracking_lag(all_frame_index-starting_all_frame_index, min(MAX_LAG,sorted_all_track_frames_diff(starting_all_frame_index))) + 1;
            current_frame = sorted_all_track_frames(all_frame_index);
            starting_all_frame_index = all_frame_index;
        else
            continue
        end
    end
end

bin_edges = 0:MAX_LAG;
lag_counts = sum(track_count_vs_tracking_lag,1);
lag_prob = lag_counts / sum(lag_counts);
%%
figure('Renderer', 'painters', 'Position', [440 290 694 520])
plot(bin_edges(2:end), lag_prob , '-o','MarkerFaceColor',[0 0.45,0.74],'MarkerSize',10,'LineWidth',2)
% % % label the first label_n datapoints
% % % label_n = 4;
label_n = numel(lag_prob);
for label_index = 1:label_n
    text(bin_edges(label_index),lag_prob(label_index), [num2str(round(lag_prob(label_index)*100,2)), '%'],'VerticalAlignment','top','HorizontalAlignment','left')
end
set(gca, 'XTick', bin_edges(1:end-1))
set(gca,'YScale','log')
xlabel('Tracking Update Time (frames)')
ylabel('Log Probability')

ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% % % set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
XOppTickLabels = round(parameters.SampleRate./bin_edges(1:end-1),2);
% % % Set the x-tick and y-tick  labels for the second axes
set(ax2, 'XTickLabel', XOppTickLabels);
set(ax2, 'YTickLabel', []);
xlabel(ax2,'Tracking Frame Rate (Hz)')

%% plot the probability of worm count at any given time
worm_counts = sum(track_count_vs_tracking_lag,2);
worm_counts_prob = worm_counts ./ sum(worm_counts);

worm_counts_enumerated = [];
for count_index = 1:numel(worm_counts)
    worm_counts_enumerated = [worm_counts_enumerated, ones(1, worm_counts(count_index)).* count_index];
end
mean_worm_counts = mean(worm_counts_enumerated);
std_worm_counts = std(worm_counts_enumerated);

%%
figure
bar(worm_counts_prob , 'hist')
axis([0 80 0 0.1])
xlabel('Tracked Worm in Frame')
ylabel('Probability')
title(['mean : ', num2str(round(mean_worm_counts)), '    ', 'std: ', num2str(round(std_worm_counts))])
% % % ylim([0 0.1])

%% prep for track duration and worm lengths distribution
relevant_track_fields = {'Path','Frames','Length','BehavioralTransition'};        

all_track_durations = [];
all_track_lengths = [];
for folder_index = 1:numel(folders)
    folder_name = folders{folder_index};
    Tracks = load_single_folder(folder_name, relevant_track_fields);
    parameters = load_parameters(folder_name);
    track_durations = zeros(1,numel(Tracks));
    
    track_lengths = zeros(1,numel(Tracks));
    for track_index = 1:numel(Tracks)
        track_durations(track_index) = numel(Tracks(track_index).Frames);
        track_lengths(track_index) = mean(Tracks(track_index).Length);
    end

    all_track_durations = [all_track_durations, track_durations./parameters.SampleRate];
    all_track_lengths = [all_track_lengths, track_lengths./parameters.CameraPixeltommConversion*1000];
end

%% plot track durations
figure5=figure;
axes5 = axes('Parent',figure5);
hold(axes5,'on');
[track_duration_prob, bin_edges,~] = histcounts(all_track_durations, 100,'Normalization','probability');
bin_centers = (bin_edges(2:end)+bin_edges(1:end-1))./2;
bar(bin_centers, track_duration_prob, 'hist')
xticks([0:600:1800])
set(axes5,'XTickLabel',...
    {'0','10','20','30'});
ylim([0 0.12])
xlim([0 1800])
yticks([0:0.06:0.12])
xlabel('Track Durations (minutes)')
ylabel('Probability')
ax = gca;
ax.FontSize = 14;
box off;
%% plot worm lengths
figure4=figure;
axes4 = axes('Parent',figure4);
hold(axes4,'on');
[track_lengths_prob, bin_edges,~] = histcounts(all_track_lengths, 'Normalization','probability');
bin_centers = (bin_edges(2:end)+bin_edges(1:end-1))./2;
% % % bin_centers=bin_centers(1:2:end);
% % % track_lengths_prob=track_lengths_prob(1:2:end);
hold on
b=bar(bin_centers, track_lengths_prob, 'hist');
% % % b.CData(1,50:end)=[0 1 0];
h=line([700 700], [0 0.1]);
h.Color= [1 0 0];
h.LineWidth= 2;
axis([400 1200 0 0.1]) 
yticks([0:0.05:0.1])
% % % set(gca, 'XTick', bin_edges)
xlabel('Worm Lengths (um)')
ylabel('Probability')
ax = gca;
ax.FontSize = 14;
box off;

disp('Analysis finished!')
%% print the total duration of single animal recordings in hrs and stimulus events
sprintf(['total number of plates= ', num2str(numel(folders))])
sprintf(['total stim events = ', num2str(numel(stim_peaks_folder)), ' events'])
sprintf(['cumulative recording duration = ', num2str(numel(vertcat(allTracks.Frames))/parameters.SampleRate/3600), ' animal hours'])
sprintf(['animals per frame: mean = ', num2str(round(mean_worm_counts)), ' std = ', num2str(round(std_worm_counts))])
sprintf(['total stim events on valid worms= ', num2str(numel(behavioral_transition_stats)), ' events'])