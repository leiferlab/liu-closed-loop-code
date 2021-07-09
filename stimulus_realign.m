function success = stimulus_realign(folder_name)
% tracks and saves individual worms for all the images in a directory
    addpath(genpath(pwd))
    parameters = load_parameters(folder_name); %load experiment parameters
    relevant_fields = {'Path', 'Frames', 'Centerlines'};
    
    %% STEP 1: initialize %%
    ImageSize = 2*round(parameters.TrackingWormBoxHalfSize * parameters.CameraPixeltommConversion);
    image_size = [ImageSize, ImageSize];
    
    %% STEP 2: Find stimulus images and time stamps %%
    projector_image_directory = [folder_name, filesep, 'ConvertedProjectorFrames', filesep];
    projector_image_files = dir([projector_image_directory, '*.png']); %get all the image files
    for image_file_index = 1:length(projector_image_files)
        file_name = strsplit(projector_image_files(image_file_index).name,'.');
        file_name = file_name{1};
        decoded_camera_frame = strsplit(file_name,'_');
        decoded_camera_frame = str2double(decoded_camera_frame{2});
        projector_image_files(image_file_index).CameraFrame = decoded_camera_frame;
    end
    possible_decoded_camera_frames = [projector_image_files.CameraFrame];
    
    Tracks = load_single_folder(folder_name, relevant_fields);
    if isempty(Tracks)
        success = false;
        return
    end
    %get timestamps for alignment
    loaded_variable = load([folder_name, filesep, 'timestamps.mat']);
    processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;
    
    %% STEP 3: loop through camera frames, align the stimulus to the centerline%%
    frame_count = length(processed_decoded_camera_frames);
    %get where each track begins and ends in terms of frames and put them
    %in a sparse binary matrix
    tracks_start_in_frame = logical(sparse(length(Tracks), frame_count));
    tracks_end_in_frame = logical(sparse(length(Tracks), frame_count));
    for track_index = 1:length(Tracks)
        tracks_start_in_frame(track_index, Tracks(track_index).Frames(1)) = true;
        tracks_end_in_frame(track_index, Tracks(track_index).Frames(end)) = true;
    end

    current_track_indecies = [];
    current_stimulus_image_decoded_camera_frame = -1;
    curProjImage = [];
    
    Tracks(1).AlignedStimulus = [];
    centerline_number_of_points = size(Tracks(1).Centerlines, 1);
    for frame_index = 1:frame_count
        tracks_that_start_in_this_frame = find(tracks_start_in_frame(:,frame_index));
        if ~isempty(tracks_that_start_in_this_frame)
            %%%there are tracks that start in this frame%%%
            previous_length = numel(current_track_indecies);
            for new_track_index = 1:length(tracks_that_start_in_this_frame)
                current_track_indecies(previous_length+new_track_index) = tracks_that_start_in_this_frame(new_track_index);
                if parameters.MultiColorStimulus
                    Tracks(tracks_that_start_in_this_frame(new_track_index)).AlignedStimulus = zeros(numel(Tracks(tracks_that_start_in_this_frame(new_track_index)).Frames), centerline_number_of_points, 3);
                else
                    Tracks(tracks_that_start_in_this_frame(new_track_index)).AlignedStimulus = zeros(numel(Tracks(tracks_that_start_in_this_frame(new_track_index)).Frames), centerline_number_of_points);
                end
            end
        end

        %%%image processing%%%
        if ~isempty(current_track_indecies) && processed_decoded_camera_frames(frame_index) > current_stimulus_image_decoded_camera_frame
            current_stimulus_image_decoded_camera_frame = processed_decoded_camera_frames(frame_index);
            curProjImage = imread([projector_image_directory, projector_image_files(find(possible_decoded_camera_frames == current_stimulus_image_decoded_camera_frame,1,'first')).name]);
        end
            
        for current_track_indecies_index = 1:length(current_track_indecies)
            %for each track in this frame, get the image
            track_index = current_track_indecies(current_track_indecies_index);
            in_track_index = frame_index - Tracks(track_index).Frames(1) + 1;
            centroid_x = double(round(Tracks(track_index).Path(in_track_index,1)));
            centroid_y = double(round(Tracks(track_index).Path(in_track_index,2)));
            image_top_left_corner_x = centroid_x-image_size(1)/2;
            image_top_left_corner_y = centroid_y-image_size(2)/2;
            image_bottom_right_corner_x = image_top_left_corner_x+image_size(1);
            image_bottom_right_corner_y = image_top_left_corner_y+image_size(2);

            cropped_image = imcrop(curProjImage, [image_top_left_corner_x, image_top_left_corner_y, (image_size-1)]);
            %pad the image if necessary
            if image_top_left_corner_x < 1 || image_top_left_corner_y < 1
                %pad the front
                cropped_image = padarray(cropped_image, [max(1-image_top_left_corner_y,0), max(1-image_top_left_corner_x,0), 0], 0, 'pre');
            end
            if image_bottom_right_corner_x > size(curProjImage,2) || image_bottom_right_corner_y > size(curProjImage,1)
                %pad the end
                cropped_image = padarray(cropped_image, [max(image_bottom_right_corner_y-size(curProjImage,1)-1,0), max(image_bottom_right_corner_x-size(curProjImage,2)-1,0), 0], 0, 'post');
            end
            
            %calculate stim for every centerline point
            red_cropped_image = squeeze(cropped_image(:,:,1));
            green_cropped_image = squeeze(cropped_image(:,:,2));
            blue_cropped_image = squeeze(cropped_image(:,:,3));
            current_centerline = round(squeeze(Tracks(track_index).Centerlines(:,:,in_track_index)));
            current_centerline(current_centerline < 1) = 1;
            current_centerline(current_centerline > ImageSize) = ImageSize;
            for centerline_point_index = 1:centerline_number_of_points
                red_signal = double(red_cropped_image(current_centerline(centerline_point_index,1),current_centerline(centerline_point_index,2)));
                green_signal = double(green_cropped_image(current_centerline(centerline_point_index,1),current_centerline(centerline_point_index,2)));
                blue_signal = double(blue_cropped_image(current_centerline(centerline_point_index,1),current_centerline(centerline_point_index,2)));
                %power conversion
                red_signal = red_signal * parameters.avgPowerRed / 255; 
                green_signal = red_signal * parameters.avgPowerGreen / 255;
                blue_signal = blue_signal * parameters.avgPowerBlue / 255;
                if parameters.MultiColorStimulus
                    Tracks(track_index).AlignedStimulus(in_track_index, centerline_point_index,:) = [red_signal, green_signal, blue_signal];
                else
                    if red_signal > 0 && blue_signal > 0
                        %conflicting signal, make it nan
                        Tracks(track_index).AlignedStimulus(in_track_index, centerline_point_index) = nan;
                    else
                        Tracks(track_index).AlignedStimulus(in_track_index, centerline_point_index) = red_signal - blue_signal;
                    end
                end
            end
        end

        tracks_that_end_in_this_frame = find(tracks_end_in_frame(:,frame_index));
        if ~isempty(tracks_that_end_in_this_frame)
            %%%there are tracks that end in this frame%%%
            image_stack_indecies = [];
            for ending_track_index = 1:length(tracks_that_end_in_this_frame)
                track_index = tracks_that_end_in_this_frame(ending_track_index);
                current_track_indecies_index = find(current_track_indecies == track_index);
                image_stack_indecies = [image_stack_indecies, current_track_indecies_index];
            end
            current_track_indecies(image_stack_indecies) = []; %stop tracking these tracks
        end
    end
    
%% STEP 4: Save the tracks %%
    savetracks(Tracks,folder_name);
    
%% STEP FINAL: return 
    success = true;
end
