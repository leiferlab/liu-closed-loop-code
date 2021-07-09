function success = save_individual_images(folder_name)
% tracks and saves individual worms for all the images in a directory
    addpath(genpath(pwd))
    parameters = load_parameters(folder_name); %load experiment parameters
    
    %% STEP 1: initialize %%
    number_of_images_for_median_projection = parameters.PostAnalysisNumberofMedianFilterImages;
    ImageSize = 2*round(parameters.TrackingWormBoxHalfSize * parameters.CameraPixeltommConversion);
    image_size = [ImageSize, ImageSize];
    MinTrackingWormAreaPixels = parameters.CameraPixeltommConversion*parameters.CameraPixeltommConversion*parameters.TrackingMinWormArea;

    %% STEP 2: Load images and labview tracks from the directory %%
    camera_image_directory = [folder_name, filesep, 'DecodedCameraFrames', filesep];
    image_files = dir([camera_image_directory, '*.png']); %get all the image files
    Tracks = load_labview_tracks(folder_name, parameters);
    if isempty(Tracks)
        error('Empty Tracks');
    end
    
    %% STEP 3: Get the median projection for use as background%%
    medianProj = imread([camera_image_directory, image_files(1).name]);
    medianProjCount = min(number_of_images_for_median_projection, length(image_files) - 1); 
    medianProj = zeros(size(medianProj,1), size(medianProj,2), medianProjCount);
    for frame_index = 1:medianProjCount
        curImage = imread([camera_image_directory, image_files(floor((length(image_files)-1)*frame_index/medianProjCount)).name]);
        medianProj(:,:,frame_index) = curImage;
    end
    medianProj = median(medianProj, 3);
    medianProj = uint8(medianProj);
    
    
    %% STEP 4: Post-Track Filtering to get rid of invalid tracks, recaculate some parameters, and collect statistics %%
    DeleteTracks = false(1,length(Tracks));
    first_frames = zeros(1,length(Tracks));
    last_frames = zeros(1,length(Tracks));
    min_pixel_displacement = parameters.MinDisplacement * parameters.CameraPixeltommConversion;
    
    if ~isempty(Tracks)
        Tracks(1).DeletionReason = [];
    end
    for track_index = 1:length(Tracks)
        first_frames(track_index) = Tracks(track_index).Frames(1);
        last_frames(track_index) = Tracks(track_index).Frames(end);
    end
    for track_index = 1:length(Tracks)
        if length(Tracks(track_index).Frames) < parameters.MinTrackLength*parameters.SampleRate
            %get rid of tracks that are too short
            Tracks(track_index).DeletionReason = [Tracks(track_index).DeletionReason, 'Short Track length '];
            DeleteTracks(track_index) = true;
        elseif mean(Tracks(track_index).Size) < parameters.MinAverageWormArea
            %get rid of worms that are too small
            Tracks(track_index).DeletionReason = [Tracks(track_index).DeletionReason, 'Small Average Worm Size '];
            DeleteTracks(track_index) = true;
        else
            %find the maximum displacement from the first time point.
            %correct for dirts that don't move
            position_relative_to_start = transpose(Tracks(track_index).Path - repmat(Tracks(track_index).Path(1,:),size(Tracks(track_index).Path,1),1));
            euclideian_distances_relative_to_start = sqrt(sum(position_relative_to_start.^2,1)); %# The two-norm of each column
            if max(euclideian_distances_relative_to_start) < min_pixel_displacement
                Tracks(track_index).DeletionReason = [Tracks(track_index).DeletionReason, 'Short Displacement '];
                DeleteTracks(track_index) = true;
            end
        end
    end
    
    %process and save deleted tracks
    if parameters.TrackingDebugMode
        deleted_tracks = Tracks(DeleteTracks);
        if ~isempty(deleted_tracks)
            fields_to_keep = {'WormIndex','Path','Size','Frames','DeletionReason'};
            track_field_names = fieldnames(deleted_tracks);
            fields_for_removal = setdiff(track_field_names,fields_to_keep);
            deleted_tracks = rmfield(deleted_tracks,fields_for_removal);
            saveFileName = [folder_name, filesep, 'tracking_deleted_tracks.mat'];
            save(saveFileName, 'deleted_tracks', '-v7.3');
        end
    end
    
    Tracks(DeleteTracks) = [];
    Tracks = rmfield(Tracks,'DeletionReason');        

    %recalculate speed and direction
    for track_index = 1:length(Tracks)
        % Smooth track data by rectangular sliding window of size WinSize;
        Tracks(track_index).SmoothX = RecSlidingWindow(Tracks(track_index).Path(:,1)', parameters.TrackingSmoothingWindow*parameters.SampleRate);
        Tracks(track_index).SmoothY = RecSlidingWindow(Tracks(track_index).Path(:,2)', parameters.TrackingSmoothingWindow*parameters.SampleRate);
        % Calculate Direction & Speed
        Xdif = CalcDif(Tracks(track_index).SmoothX, parameters.TrackingSmoothingWindow*parameters.SampleRate) * parameters.SampleRate;
        Ydif = -CalcDif(Tracks(track_index).SmoothY, parameters.TrackingSmoothingWindow*parameters.SampleRate) * parameters.SampleRate;    % Negative sign allows "correct" direction
        % cacluation (i.e. 0 = Up/North)
        Ydif(Ydif == 0) = eps;     % Avoid division by zero in direction calculation
        Tracks(track_index).Direction = atan(Xdif./Ydif) * 360/(2*pi);	    % In degrees, 0 = Up ("North")
        NegYdifIndexes = find(Ydif < 0);
        Index1 = find(Tracks(track_index).Direction(NegYdifIndexes) <= 0);
        Index2 = find(Tracks(track_index).Direction(NegYdifIndexes) > 0);
        Tracks(track_index).Direction(NegYdifIndexes(Index1)) = Tracks(track_index).Direction(NegYdifIndexes(Index1)) + 180;
        Tracks(track_index).Direction(NegYdifIndexes(Index2)) = Tracks(track_index).Direction(NegYdifIndexes(Index2)) - 180;
        Tracks(track_index).Speed = sqrt(Xdif.^2 + Ydif.^2) / parameters.CameraPixeltommConversion;		% In mm/sec
    end    
    
%% STEP 6: save each worms images, recalculate select stats%%
    savePath = [folder_name, filesep, 'individual_worm_imgs', filesep];
    if ~exist(savePath, 'dir')
        mkdir(savePath)
    end
    
    delete_extra_individual_worm_images(folder_name, 0); %delete previous .mat files
    
    frame_count = length(image_files)-1;
    %get where each track begins and ends in terms of frames and put them
    %in a sparse binary matrix
    tracks_start_in_frame = logical(sparse(length(Tracks), frame_count));
    tracks_end_in_frame = logical(sparse(length(Tracks), frame_count));
    for track_index = 1:length(Tracks)
        Tracks(track_index).Eccentricity = zeros(numel(Tracks(track_index).Frames),1); %preallocate for Eccentricity
        tracks_start_in_frame(track_index, Tracks(track_index).Frames(1)) = true;
        tracks_end_in_frame(track_index, Tracks(track_index).Frames(end)) = true;
    end

    current_image_stacks = [];
    starting_binary_threshold = parameters.TrackingBinaryThresholdLevel / 255;

    for frame_index = 1:frame_count
        tracks_that_start_in_this_frame = find(tracks_start_in_frame(:,frame_index));
        if ~isempty(tracks_that_start_in_this_frame)
            %%%there are tracks that start in this frame%%%
            previous_length = length(current_image_stacks);
            current_image_stacks(previous_length+length(tracks_that_start_in_this_frame)).TrackIndex = []; %preallocate memory
            for new_track_index = 1:length(tracks_that_start_in_this_frame)
                track_index = tracks_that_start_in_this_frame(new_track_index);
                current_image_stacks(previous_length+new_track_index).TrackIndex = track_index;
                current_image_stacks(previous_length+new_track_index).Images = zeros([image_size, length(Tracks(track_index).Frames)], 'uint8');
            end
        end

        %%%image processing%%%
        curImage = imread([camera_image_directory, image_files(frame_index).name]);
        subtractedImage = curImage - medianProj; %subtract median projection  - imageBackground
        % Convert frame to a binary image 
            
        for image_stack_index = 1:length(current_image_stacks)
            %for each track in this frame, get the cropped image
            track_index = current_image_stacks(image_stack_index).TrackIndex;
            in_track_index = frame_index - Tracks(track_index).Frames(1) + 1;
            centroid_x = double(round(Tracks(track_index).Path(in_track_index,1)));
            centroid_y = double(round(Tracks(track_index).Path(in_track_index,2)));
            image_top_left_corner_x = centroid_x-image_size(1)/2;
            image_top_left_corner_y = centroid_y-image_size(2)/2;
            image_bottom_right_corner_x = image_top_left_corner_x+image_size(1);
            image_bottom_right_corner_y = image_top_left_corner_y+image_size(2);

            %lets binary threshold here. make sure that our blob is at
            %least the min size requirement. lower thresh if we need it
            cropped_image = imcrop(subtractedImage, [image_top_left_corner_x, image_top_left_corner_y, (image_size-1)]);
            current_thresholded_area = 0;
            current_threshold = starting_binary_threshold;
            while current_thresholded_area < MinTrackingWormAreaPixels && current_threshold > 0
                cropped_mask_image = im2bw(cropped_image, current_threshold);  % For tracking bright objects on a dark background
                [cropped_mask_image_L,numComponents] = bwlabel(cropped_mask_image);
                % take the region closest to the worm image center
                if numComponents > 0
                    STATS = regionprops(cropped_mask_image_L, {'Area','Eccentricity','Centroid'});
                    [~,closest_index] = pdist2(vertcat(STATS.Centroid),image_size/2,'euclidean','Smallest',1);
                    current_thresholded_area = STATS(closest_index).Area;
                else
                    current_thresholded_area = 0;
                end
                current_threshold = current_threshold - (1/255); %lower threshold to keep going
            end
            %replace the area, eccentricity calculations for this track
            Tracks(track_index).Size(in_track_index) = STATS(closest_index).Area;
            Tracks(track_index).Eccentricity(in_track_index) = STATS(closest_index).Eccentricity;
            single_worm = cropped_mask_image_L == closest_index; %get an binary mask of largest blob in the cropped image
            single_worm = bwmorph(single_worm, 'fill');
            worm_frame = imcrop(subtractedImage, [image_top_left_corner_x, image_top_left_corner_y, (image_size-1)]);
            worm_frame(~single_worm) = 0; %mask

            %pad the image if necessary
            if image_top_left_corner_x < 1 || image_top_left_corner_y < 1
                %pad the front
                worm_frame = padarray(worm_frame, [max(1-image_top_left_corner_y,0), max(1-image_top_left_corner_x,0)], 0, 'pre');
            end
            if image_bottom_right_corner_x > size(subtractedImage,2) || image_bottom_right_corner_y > size(subtractedImage,1)
                %pad the end
                worm_frame = padarray(worm_frame, [max(image_bottom_right_corner_y-size(subtractedImage,1)-1,0), max(image_bottom_right_corner_x-size(subtractedImage,2)-1,0)], 0, 'post');
            end

            current_image_stacks(image_stack_index).Images(:,:,in_track_index) = worm_frame;
            
        end

        tracks_that_end_in_this_frame = find(tracks_end_in_frame(:,frame_index));
        if ~isempty(tracks_that_end_in_this_frame)
            %%%there are tracks that end in this frame, do the computation%%%
            image_stack_indecies = [];
            for ending_track_index = 1:length(tracks_that_end_in_this_frame)
                track_index = tracks_that_end_in_this_frame(ending_track_index);
                image_stack_index = find([current_image_stacks.TrackIndex] == track_index);
                image_stack_indecies = [image_stack_indecies, image_stack_index];
                worm_images = current_image_stacks(image_stack_index).Images;
                save([folder_name, filesep, 'individual_worm_imgs', filesep, 'worm_', num2str(track_index), '.mat'], 'worm_images', '-v7.3');
            end
            current_image_stacks(image_stack_indecies) = []; %clear the memory of these images
        end
    end
    
%% STEP 5: Save the tracks %%
    savetracks(Tracks,folder_name);
    
    
%% STEP FINAL: return 
    success = true;
end
