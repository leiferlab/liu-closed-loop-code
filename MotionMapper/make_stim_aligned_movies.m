%This script allows the user to pick a watershed region in a t-SNE map and
%finds when the animal stays in that region and plots their behaivor, not used in paper


%% STEP 1: establish plotting constants 
N_rows = 10;
N_columns = 11;
N = N_rows*N_columns;
fps = 14;
frames_before = 5*fps-1;
frames_after = 10*fps;
normalized_stimuli = 1; %delta function
search_window = 2*fps;

duration = frames_before+frames_after+1;
parameters = load_parameters();
load('reference_embedding.mat')
relevant_track_fields = {'Frames','Centerlines','BehavioralTransition'};%'Direction','Velocity';


folders = getfoldersGUI();
% allow user to select the folder to save as
%pathname = uigetdir('', 'Select Save Folder')
if isequal(pathname,0)
    %cancel
   return
end
%% STEP 2: load tracks
[allTracks, folder_indecies, track_indecies] = loadtracks(folders,relevant_track_fields );
allTracks = BehavioralTransitionToBehavioralAnnotation(allTracks);

for track_index = 1:length(allTracks)
    %we need to remember exactly when the track started
    allTracks(track_index).FirstFrame = allTracks(track_index).Frames(1);
end

LEDVoltages = load([folders{1}, filesep, 'LEDVoltages.txt']);

%find when each stimuli is played back by convolving the time
%reversed stimulus (cross-correlation)
xcorr_ledvoltages_stimulus = padded_conv(LEDVoltages, normalized_stimuli);
peak_thresh = 0.99.*max(xcorr_ledvoltages_stimulus); %the peak threshold is 99% of the max (because edge effects)
[~, critical_frames] = findpeaks(xcorr_ledvoltages_stimulus, 'MinPeakHeight', peak_thresh,'MinPeakDistance',14);

%% STEP 3: filter tracks
saveFileName = fullfile(pathname,'sample.mp4');
time_aligned_tracks = [];
time_aligned_folder_indecies = [];
time_aligned_track_indecies = [];

for critical_frame_index = 1:length(critical_frames)
    %for every time a stimulus is delivered, look at a certain range of
    %frames
    critical_window_start = critical_frames(critical_frame_index) - frames_before;
    critical_window_end = critical_frames(critical_frame_index) + frames_after;
    if critical_window_end <= length(LEDVoltages) && critical_window_start >= 1 
        %make sure the current window is in range
        [tracks_in_critical_window, track_indecies_preserved] = FilterTracksByTime(allTracks, critical_window_start, critical_window_end, true);
        time_aligned_tracks = [time_aligned_tracks, tracks_in_critical_window];
        time_aligned_folder_indecies = [time_aligned_folder_indecies, folder_indecies(track_indecies_preserved)];
        time_aligned_track_indecies = [time_aligned_track_indecies, track_indecies(track_indecies_preserved)];
    end
end

%% STEP 5: get N tracks randomly
perm = randperm(length(time_aligned_tracks));
selected_indecies = perm(1:N);
selected_tracks = time_aligned_tracks(selected_indecies);
selected_folder_indecies = time_aligned_folder_indecies(selected_indecies);
selected_track_indecies = time_aligned_track_indecies(selected_indecies);

%% STEP 6: sort them by their resulting behaviors
% first by transition into (if at all)
% then by when 

sorted_indecies = [];
sorting_struct = cell(1,length(behavior_names)+1);

%loop through the tracks, looking for the next transitional behavior
for track_index = 1:N
    %starting from tap, scan forward.
    %cases:
    %take the next behavior if it is within search_window (current behavior does not count unless we don't change for entire 2s duration),
    %if still transitioning by the end of search_window, take the next behavior (if it exists, if doesn't exist, save it for last)
    behavior_annotation = selected_tracks(track_index).BehavioralAnnotation(frames_before:end);
    initial_behavior = behavior_annotation(1);
    if initial_behavior == 0
        transition_happened = true;
    else
        transition_happened = false;
    end
    for frame_index = 2:length(behavior_annotation)+1
        if frame_index == length(behavior_annotation)+1
            %we run out of frames, and the animal is still "transitioning"
            sorting_struct{length(behavior_names)+1} = [sorting_struct{length(behavior_names)+1}; track_index, frame_index];
        elseif transition_happened && behavior_annotation(frame_index) > 0
            %next behavior found, in the window
            sorting_struct{behavior_annotation(frame_index)} = [sorting_struct{behavior_annotation(frame_index)}; track_index, frame_index];
            break
        elseif frame_index >= search_window && ~transition_happened
            % the worm remained in this state the entire duration
            sorting_struct{initial_behavior} = [sorting_struct{initial_behavior}; track_index, frame_index];
            break
        elseif behavior_annotation(frame_index) ~= initial_behavior
            %transition happened
            transition_happened = true;
        end
    end
end

%loop through the sorted struct, sort by time
for sorting_struct_index = 1:length(behavior_names)+1
    if ~isempty(sorting_struct{sorting_struct_index})
        [~,I] = sort(sorting_struct{sorting_struct_index}(:,2));
        current_sorted_indecies = sorting_struct{sorting_struct_index}(:,1);
        current_sorted_indecies = current_sorted_indecies(I);
        sorted_indecies = [sorted_indecies, current_sorted_indecies'];
    end
end

selected_tracks = selected_tracks(sorted_indecies);
selected_folder_indecies = selected_folder_indecies(sorted_indecies);
selected_track_indecies = selected_track_indecies(sorted_indecies);

%% STEP 7: plot the behaviors

%load the worm images
clear required_worm_images
required_worm_images(N).worm_images = [];
for worm_images_index = 1:N
    image_file = fullfile([folders{selected_folder_indecies(worm_images_index)},filesep,'individual_worm_imgs',filesep,'worm_', num2str(selected_track_indecies(worm_images_index)), '.mat']);
    required_worm_images(worm_images_index) = load(image_file);
end


behavior_figure = figure('Position', [0, 0, N_columns*100, N_rows*100]);
outputVideo = VideoWriter(saveFileName,'MPEG-4');
outputVideo.FrameRate = fps;
open(outputVideo)
time_series = -frames_before:frames_after;

for relative_frame_index = 1:frames_before+frames_after+1
    for subplot_index = 1:N
        worm_frame_index = selected_tracks(subplot_index).Frames(relative_frame_index) - selected_tracks(subplot_index).FirstFrame + 1;
        if worm_frame_index < 1 || worm_frame_index > size(required_worm_images(subplot_index).worm_images,3)
            %the video does not exist, skip
            continue
        else
            behavior_annotation = selected_tracks(subplot_index).BehavioralAnnotation(relative_frame_index);
            if behavior_annotation > 0
                current_color = behavior_colors(behavior_annotation,:);
            else
                current_color = [];
            end
                
            subplot_tight(N_rows,N_columns,subplot_index,0);
            plot_worm_frame(required_worm_images(subplot_index).worm_images(:,:,worm_frame_index), squeeze(selected_tracks(subplot_index).Centerlines(:,:,relative_frame_index)), ...
            current_color, [], [], [], [], [], 0);
            if subplot_index == N-floor(N_columns/2)
                time_text = [datestr(abs(time_series(relative_frame_index))/24/3600/fps,'SS.FFF'), ' s'];
                if time_series(relative_frame_index) < 0
                    time_text = ['-', time_text];
                else
                    time_text = [' ', time_text];
                end
                text(size(required_worm_images(subplot_index).worm_images,1)/2,size(required_worm_images(subplot_index).worm_images,2)*.7,time_text,'color','red','fontsize',10,'VerticalAlignment','top','HorizontalAlignment','center')
            end
        end
        
        
    end


    %pause
    writeVideo(outputVideo, getframe(gcf));
    clf
end
close(outputVideo)
close(behavior_figure)
