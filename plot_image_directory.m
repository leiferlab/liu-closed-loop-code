function success = plot_image_directory(folder_name)
% plots the individual worm videos and all the track videos
    addpath(genpath(pwd))
    %% STEP 1: initialize %%
    parameters = load_parameters(folder_name); %load experiment parameters
    load('reference_embedding.mat')
    number_of_images_for_median_projection = parameters.PostAnalysisNumberofMedianFilterImages;
    inset_magification = 10;
    if parameters.IndividualVideoPlottingFrameRate > 0
        %plot individual worms
        relevant_track_fields = {'Centerlines','UncertainTips','Eccentricity','BehavioralTransition'...
            'Direction','Speed','TotalScore','AlignedStimulus','Path','Frames','Size'};
    else
        %plot whole field only
        relevant_track_fields = {'Centerlines','Path','Frames','Size','BehavioralTransition'};        
    end
    %% Load tracks
    Tracks = load_single_folder(folder_name, relevant_track_fields);
    if isempty(Tracks)
        error('Empty Tracks');
    end
    Tracks = BehavioralTransitionToBehavioralAnnotation(Tracks);
    
    %% STEP 2: plot individual worms
    if parameters.IndividualVideoPlottingFrameRate > 0
        individual_figure = figure;
        individual_worm_videos(Tracks, folder_name, parameters, behavior_colors);
        close(individual_figure)
    end
    
    %% STEP 3: Load images and other properties from the directory %%
    % check if preferences indicate not to plot
    if parameters.PlottingFrameRate <= 0
        success = true;
        return
    end
    
    % Get all the image files
    camera_image_directory = [folder_name, filesep, 'DecodedCameraFrames', filesep];
    image_files = dir([camera_image_directory, '*.png']); %get all the image files
    
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
    
    loaded_variable = load([folder_name, filesep, 'timestamps.mat']);
    processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;
   
    
    % Load deleted tracks if we are debugging mode
    all_deleted_tracks = [];
    if parameters.TrackingDebugMode
        deleted_track_file_name = [folder_name, filesep, 'tracking_deleted_tracks.mat'];
        if exist(deleted_track_file_name, 'file') == 2
            load(deleted_track_file_name);
        else
            deleted_tracks = [];
        end
        all_deleted_tracks = deleted_tracks;
        deleted_track_file_name = [folder_name, filesep, 'centerline_deleted_tracks.mat'];
        if exist(deleted_track_file_name, 'file') == 2
            load(deleted_track_file_name);
            for track_index = 1:length(deleted_tracks)
                deleted_tracks(track_index).DeletionReason = 'Low Centerline Score';
            end
        else
            deleted_tracks = [];
        end

        all_deleted_tracks = [all_deleted_tracks, deleted_tracks];
        clear deleted_tracks
    end
    
    %% STEP 4: Get the median z projection %%
    medianProj = imread([camera_image_directory, image_files(1).name]);
    medianProjCount = min(number_of_images_for_median_projection, length(image_files) - 1); 
    medianProj = zeros(size(medianProj,1), size(medianProj,2), medianProjCount);
    for frame_index = 1:medianProjCount
        curImage = imread([camera_image_directory, image_files(floor((length(image_files)-1)*frame_index/medianProjCount)).name]);
        medianProj(:,:,frame_index) = curImage;
    end
    medianProj = median(medianProj, 3);
    medianProj = uint8(medianProj);
    
    %% STEP 5: plot all the tracks on top of stimulus and raw images
    % Setup figure for plotting tracker results
    % -----------------------------------------
    WTFigH = findobj('Tag', 'WTFIG');
    if isempty(WTFigH)
        WTFigH = figure('Name', 'Tracking Results', ...
            'NumberTitle', 'off', ...
            'Tag', 'WTFIG','units','normalized','outerposition',[0 0 2 2]);
    else
        figure(WTFigH);
    end

    frames_per_plot_time = round(parameters.SampleRate/parameters.PlottingFrameRate);
    
    %save subtracted avi
%    outputVideo = VideoWriter(fullfile([folder_name, filesep, 'processed']),'MPEG-4'); not supported in linux
    video_file_name = fullfile([folder_name, filesep, 'processed', '.avi']);
    outputVideo = VideoWriter(video_file_name,'Uncompressed AVI');
%     outputVideo.Quality = 100;
    outputVideo.FrameRate = parameters.PlottingFrameRate;
    open(outputVideo)
    warning('off','images:initSize:adjustingMag');
    current_track = 0;
    for frame_index = 1:frames_per_plot_time:length(image_files) - 1
        % Get Camera Frame and Projector Frames
        curImage = imread([camera_image_directory, image_files(frame_index).name]) * parameters.CameraImagePlottingMultiplier;
        current_stimulus_image_decoded_camera_frame = processed_decoded_camera_frames(frame_index);
        curProjImage = imread([projector_image_directory, projector_image_files(find(possible_decoded_camera_frames == current_stimulus_image_decoded_camera_frame,1,'first')).name]);
        nolag_projector_file_index = find(possible_decoded_camera_frames >= frame_index,1,'first');
        if ~isempty(nolag_projector_file_index)
            nolagProjImage = imread([projector_image_directory, projector_image_files(nolag_projector_file_index).name]);
        else
            nolagProjImage = imread([projector_image_directory, projector_image_files(1).name]);
        end
        
        if ~parameters.PlottingRawImage
            curImage = curImage - uint8(medianProj); %subtract median projection  - imageBackground
        end
        
        %combine the stimulus image and the current image 
        curImage = cat(3, curImage, curImage, curImage); % convert greyscale to RGB
        curImage = curImage + curProjImage;
        active_tracks = PlotFrame(WTFigH, curImage, Tracks, frame_index, all_deleted_tracks, behavior_colors);
        if ~isempty(active_tracks)
            %draw inset video
            figure(WTFigH);
            axis manual;
            hold on;
            % get the smallest active track
            if current_track == 0 || frame_index > max(Tracks(current_track).Frames)
                current_track = min(active_tracks);
                loaded_file = load([folder_name, filesep, 'individual_worm_imgs', filesep, 'worm_', num2str(current_track), '.mat']);
                worm_images = loaded_file.worm_images;
                image_size = [size(worm_images,1), size(worm_images,2)];
            end
            in_track_index = frame_index-Tracks(current_track).Frames(1)+1;

            
            %get a cropped stimulus image for this worm
            centroid_x = double(round(Tracks(current_track).Path(in_track_index,1)));
            centroid_y = double(round(Tracks(current_track).Path(in_track_index,2)));
            image_top_left_corner_x = centroid_x-image_size(1)/2;
            image_top_left_corner_y = centroid_y-image_size(2)/2;
            image_bottom_right_corner_x = image_top_left_corner_x+image_size(1);
            image_bottom_right_corner_y = image_top_left_corner_y+image_size(2);
            cropped_curr_image = imcrop(curProjImage, [image_top_left_corner_x, image_top_left_corner_y, (image_size-1)]);
            cropped_nolag_image = imcrop(nolagProjImage, [image_top_left_corner_x, image_top_left_corner_y, (image_size-1)]);
            %pad the image if necessary
            if image_top_left_corner_x < 1 || image_top_left_corner_y < 1
                %pad the front
                cropped_curr_image = padarray(cropped_curr_image, [max(1-image_top_left_corner_y,0), max(1-image_top_left_corner_x,0), 0], 0, 'pre');
                cropped_nolag_image = padarray(cropped_nolag_image, [max(1-image_top_left_corner_y,0), max(1-image_top_left_corner_x,0), 0], 0, 'pre');
            end
            if image_bottom_right_corner_x > size(curProjImage,2) || image_bottom_right_corner_y > size(curProjImage,1)
                %pad the end
                cropped_curr_image = padarray(cropped_curr_image, [max(image_bottom_right_corner_y-size(curProjImage,1)-1,0), max(image_bottom_right_corner_x-size(curProjImage,2)-1,0), 0], 0, 'post');
                cropped_nolag_image = padarray(cropped_nolag_image, [max(image_bottom_right_corner_y-size(curProjImage,1)-1,0), max(image_bottom_right_corner_x-size(curProjImage,2)-1,0), 0], 0, 'post');
            end
            %treat the nolag image
            cropped_nolag_image = uint8(bwmorph(im2bw(rgb2gray(cropped_nolag_image),0.05),'remove') .* 255);
            cropped_nolag_image = imresize(cropped_nolag_image, size(cropped_nolag_image).*inset_magification);
            
            I = squeeze(worm_images(:,:,in_track_index));
            I_resize = imadjust(imresize(I, size(I).*inset_magification));
            cropped_curr_image_resize = imresize(cropped_curr_image, [size(cropped_curr_image,1)*inset_magification, size(cropped_curr_image,2)*inset_magification]);
            
            imshow(cat(3, I_resize, I_resize, I_resize) + cat(3,cropped_nolag_image,cropped_nolag_image,cropped_nolag_image) + cropped_curr_image_resize);
            plot(Tracks(current_track).Centerlines(1,2,in_track_index).*inset_magification, Tracks(current_track).Centerlines(1,1,in_track_index).*inset_magification, '.', 'Color', [0 1 0], 'markersize',50)
            plot(Tracks(current_track).Centerlines(:,2,in_track_index).*inset_magification, Tracks(current_track).Centerlines(:,1,in_track_index).*inset_magification, '-', 'Color', [0 1 0], 'LineWidth',4)

            rectangle('Position',[0,0,size(I_resize)],'EdgeColor', 'g', 'LineWidth',5,'LineStyle','-')
            rectangle('Position',[Tracks(current_track).Path(in_track_index,1)-(size(I,1)/2),Tracks(current_track).Path(in_track_index,2)-(size(I,2)/2),size(I)],'EdgeColor', 'g', 'LineWidth',2,'LineStyle','-')
 
            %add scale bars
            % 1 mm on the inset
            half_distance = 1*parameters.CameraPixeltommConversion*inset_magification/2;
            line([size(I_resize,2)/2-half_distance, size(I_resize,2)/2+half_distance], ...
                0.8*[size(I_resize,1), size(I_resize,1)], 'Color','yellow','LineWidth',3)
            text(size(I_resize,2)/2, 0.8*size(I_resize,1), '1 mm', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment','top','fontsize',20,'color','yellow')
            hold off
        end
        
        %add scale bar to main
        hold on
        % 1 cm on the inset
        half_distance = 10*parameters.CameraPixeltommConversion/2;
        line([size(curImage,2)/2-half_distance, size(curImage,2)/2+half_distance], ...
            0.95*[size(curImage,1), size(curImage,1)], 'Color','yellow','LineWidth',3)
        text(size(curImage,2)/2, 0.95*size(curImage,1), '1 cm', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment','top','fontsize',20,'color','yellow')
        
        hold off

        FigureName = ['Tracking Results for Frame ', num2str(frame_index)];
        set(WTFigH, 'Name', FigureName);
        writeVideo(outputVideo, getframe(WTFigH));
        
    end
    close(outputVideo) 
    close(WTFigH)
    video_transcode(video_file_name);
%     success = true;
    success = false;
end
