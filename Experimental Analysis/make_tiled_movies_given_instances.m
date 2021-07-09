function [] = make_tiled_movies_given_instances(Tracks, folders, saveFileName, possible_tracks, possible_frames, N_rows, N_columns, video_duration,display_field_name, highlight_index)
    parameters = load_parameters(folders{1});
    fps = parameters.SampleRate;
    load('reference_embedding.mat')

    %% STEP 1: establish plotting constants 
    if nargin < 4
        N_rows = 4;
        N_columns = 4;
        frames_before = 3*fps-1;
        frames_after = 3*fps;        
        display_field_name = 'Velocity';
        highlight_index = 0;
    else
        frames_before = round(video_duration/2)-1;
        frames_after = round(video_duration/2);        
    end
    N = N_rows*N_columns;


%%
    if length(possible_tracks) < N
        display(['Output video is truncated. There are only ' num2str(length(possible_frames)) ' observed behaviors, fewer than the required ' num2str(N)])
%         return
    end


    %% STEP 5: get N points that are from different tracks and fits the criteria
    selected_tracks = [];
    selected_frames = [];
    current_index = 1;
    while length(selected_tracks) < N
        current_track_number = possible_tracks(current_index);
        current_track_length = numel(Tracks(current_track_number).Frames);
        current_frame_number = possible_frames(current_index);
        if current_frame_number - frames_before < 1 || current_frame_number + frames_after > current_track_length
            %this point will be cut out at some point, throw it out
        else
            data_point_accepted = false;
            if ismember(current_track_number, selected_tracks)
                previously_found_indecies = find(selected_tracks==current_track_number);
                previously_found_frames = selected_frames(previously_found_indecies);
                covered_frames = current_frame_number-frames_before:current_frame_number+frames_after;
                if ~sum(ismember(previously_found_frames,covered_frames))               
                    %this track has been represented, but this behavior is out
                    %of range
                    data_point_accepted = true;
                end
            else
                %this track has not been represented
                data_point_accepted = true;
            end
            if data_point_accepted
                selected_tracks = [selected_tracks, current_track_number];
                selected_frames = [selected_frames, current_frame_number];  
            end
        end
        current_index = current_index + 1;
    end
    N = current_index-1;

    %% STEP 7: plot the behaviors
    %load the relevant worm images and timestamps
    clear required_worm_images
    required_worm_images(N).worm_images = [];
    required_worm_images(N).possible_decoded_camera_frames = [];
    required_worm_images(N).processed_decoded_camera_frames = [];
    required_worm_images(N).projector_image_files = [];
    required_worm_images(N).image_size = [];
    for sample_index = 1:N
        track_index = selected_tracks(sample_index);
        folder_index = Tracks(track_index).folder_index;
        image_file = fullfile([folders{folder_index},filesep,'individual_worm_imgs',filesep,'worm_', num2str(Tracks(track_index).within_folder_track_index), '.mat']);
        loaded_variable = load(image_file);
        required_worm_images(sample_index).worm_images = loaded_variable.worm_images;
        required_worm_images(sample_index).image_size = [size(loaded_variable.worm_images,1), size(loaded_variable.worm_images,2)];
        
        %determine if this folder has been processed before
        previously_processed_index = [];
        for index = 1:sample_index-1
            if Tracks(selected_tracks(index)).folder_index == folder_index
                previously_processed_index = index;
                break
            end
        end
        
        if ~isempty(previously_processed_index)
            %we have information on this folder before, copy it
            required_worm_images(sample_index).possible_decoded_camera_frames = required_worm_images(previously_processed_index).possible_decoded_camera_frames;
            required_worm_images(sample_index).processed_decoded_camera_frames = required_worm_images(previously_processed_index).processed_decoded_camera_frames;
            required_worm_images(sample_index).projector_image_files = required_worm_images(previously_processed_index).projector_image_files;
        else
            %get stimulus timing
            projector_image_directory = [folders{folder_index}, filesep, 'ConvertedProjectorFrames', filesep];
            projector_image_files = dir([projector_image_directory, '*.png']); %get all the image files
            for image_file_index = 1:length(projector_image_files)
                file_name = strsplit(projector_image_files(image_file_index).name,'.');
                file_name = file_name{1};
                decoded_camera_frame = strsplit(file_name,'_');
                decoded_camera_frame = str2double(decoded_camera_frame{2});
                projector_image_files(image_file_index).CameraFrame = decoded_camera_frame;
            end
            required_worm_images(sample_index).possible_decoded_camera_frames = [projector_image_files.CameraFrame];

            loaded_variable = load([folders{folder_index}, filesep, 'timestamps.mat']);
            required_worm_images(sample_index).processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;
            required_worm_images(sample_index).projector_image_files = projector_image_files;
        end
    end
    
    behavior_figure = figure('Position', [0, 0, N_columns*150, N_rows*150]);
    outputVideo = VideoWriter(saveFileName,'MPEG-4');
    outputVideo.FrameRate = fps;
    open(outputVideo)

    for relative_frame_index = -frames_before:frames_after
        for subplot_index = 1:N
            in_track_index = selected_frames(subplot_index) + relative_frame_index;
            track_index = selected_tracks(subplot_index);
            if in_track_index < 1 || in_track_index > size(required_worm_images(subplot_index).worm_images,3)
                %the video does not exist, skip
                continue
            else
                subplot_tight(N_rows,N_columns,subplot_index,0);
                
                %get a stimulus image for this frame, crop it if necessary
                current_stimulus_image_decoded_camera_frame = required_worm_images(subplot_index).processed_decoded_camera_frames(Tracks(track_index).Frames(in_track_index));
                current_projector_image_directory = [folders{Tracks(track_index).folder_index}, filesep, 'ConvertedProjectorFrames', filesep];
                curProjImage = imread([current_projector_image_directory, required_worm_images(subplot_index).projector_image_files(find(required_worm_images(subplot_index).possible_decoded_camera_frames == current_stimulus_image_decoded_camera_frame,1,'first')).name]);
                centroid_x = double(round(Tracks(track_index).Path(in_track_index,1)));
                centroid_y = double(round(Tracks(track_index).Path(in_track_index,2)));
                image_top_left_corner_x = centroid_x-required_worm_images(subplot_index).image_size(1)/2;
                image_top_left_corner_y = centroid_y-required_worm_images(subplot_index).image_size(2)/2;
                image_bottom_right_corner_x = image_top_left_corner_x+required_worm_images(subplot_index).image_size(1);
                image_bottom_right_corner_y = image_top_left_corner_y+required_worm_images(subplot_index).image_size(2);
                cropped_image = imcrop(curProjImage, [image_top_left_corner_x, image_top_left_corner_y, (required_worm_images(subplot_index).image_size-1)]);
                %pad the image if necessary
                if image_top_left_corner_x < 1 || image_top_left_corner_y < 1
                    %pad the front
                    cropped_image = padarray(cropped_image, [max(1-image_top_left_corner_y,0), max(1-image_top_left_corner_x,0), 0], 0, 'pre');
                end
                if image_bottom_right_corner_x > size(curProjImage,2) || image_bottom_right_corner_y > size(curProjImage,1)
                    %pad the end
                    cropped_image = padarray(cropped_image, [max(image_bottom_right_corner_y-size(curProjImage,1)-1,0), max(image_bottom_right_corner_x-size(curProjImage,2)-1,0), 0], 0, 'post');
                end
                
                I = required_worm_images(subplot_index).worm_images(:,:,in_track_index);
                current_behavior = Tracks(track_index).BehavioralAnnotation(in_track_index);
                if current_behavior < 1 || current_behavior > length(behavior_colors)
                    current_behavior_color = [1 1 1];
                else
                    current_behavior_color = behavior_colors(Tracks(track_index).BehavioralAnnotation(in_track_index),:);
                end
                %convert the stimulus in uW/mm^2 to a color
                current_stimulus_along_centerline = Tracks(track_index).AlignedStimulus(in_track_index,:);
                for centerline_point_index = 1:length(current_stimulus_along_centerline) %convert power intensity back to pixel intensity
                    if current_stimulus_along_centerline(centerline_point_index) < 0
                        % blue
                        current_stimulus_along_centerline(centerline_point_index) = max(current_stimulus_along_centerline(centerline_point_index) / parameters.avgPowerBlue * 255, -255);
                    else
                        current_stimulus_along_centerline(centerline_point_index) = min(current_stimulus_along_centerline(centerline_point_index) / parameters.avgPowerRed * 255, 255);
                    end
                end
%                 plot_worm_frame(I, squeeze(Tracks(track_index).Centerlines(:,:,in_track_index)), current_behavior_color, ...
%                     [], current_stimulus_along_centerline, ...
%                     cropped_image, num2str(Tracks(track_index).(display_field_name)(in_track_index)));
                
                plot_worm_frame(I, squeeze(Tracks(track_index).Centerlines(:,:,in_track_index)), [0 1 0], ...
                    [], current_stimulus_along_centerline, ...
                    cropped_image, num2str(subplot_index));
                
                if ismember(subplot_index,highlight_index)
                    hold on
                    rectangle('position', [2,2, (size(I)-2)], 'edgecolor', [1 1 0], 'linewidth', 2);
                    hold off
                end
                
                ga = axes('Position',[0,0,1,1],'Xlim',[0,400],'Ylim',[0,400],'tag','ga');
                % set print margins
                topm = 400; botm = 0;
                rgtm = 400; lftm = 0;
                ctrm = (rgtm-lftm)/2;

                time_text = datestr(abs(relative_frame_index)/24/3600/fps,'SS.FFF');
                if relative_frame_index < 0
                    time_text = ['-', time_text];
                end
                text(ctrm,botm+35,time_text,'color','red','fontsize',20,'VerticalAlignment','top','HorizontalAlignment','center')

                % make sure the plot is visible
                set(ga,'vis','off');

            end
        end

%        pause
        writeVideo(outputVideo, getframe(behavior_figure));
        clf
    end
    close(outputVideo)
    close(behavior_figure)
end