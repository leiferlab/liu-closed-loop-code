function [] = individual_worm_videos(Tracks, folder_name, parameters, behavior_colors)
    plotting_fps = min(parameters.IndividualVideoPlottingFrameRate,parameters.SampleRate);
%     ImageSize = 2*round(parameters.TrackingWormBoxHalfSize * parameters.CameraPixeltommConversion);
%     image_size = [ImageSize, ImageSize];

    %get stimulus
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


    if ~exist([folder_name, filesep, 'individual_worm_videos'], 'dir')
        mkdir([folder_name, filesep, 'individual_worm_videos'])
    end
    % Plots a single worm over time along with its centerline and its stimulus
    frames_per_plot_time = round(parameters.SampleRate/plotting_fps);
    for track_index = 1:length(Tracks)
        loaded_file = load([folder_name, filesep, 'individual_worm_imgs', filesep, 'worm_', num2str(track_index), '.mat']);
        worm_images = loaded_file.worm_images;
        image_size = [size(worm_images,1), size(worm_images,2)];
        video_file_name = fullfile([folder_name, filesep, 'individual_worm_videos', filesep, 'worm_', num2str(track_index), '.avi']);
        outputVideo = VideoWriter(video_file_name,'Uncompressed AVI');
        outputVideo.FrameRate = plotting_fps;
        open(outputVideo)

        for in_track_index = 1:frames_per_plot_time:size(worm_images,3)
            %get a stimulus image for this frame
            current_stimulus_image_decoded_camera_frame = processed_decoded_camera_frames(Tracks(track_index).Frames(in_track_index));
            curProjImage = imread([projector_image_directory, projector_image_files(find(possible_decoded_camera_frames == current_stimulus_image_decoded_camera_frame,1,'first')).name]);
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
            I = squeeze(worm_images(:,:,in_track_index));
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

            plot_worm_frame(I, squeeze(Tracks(track_index).Centerlines(:,:,in_track_index)), current_behavior_color, ...
                Tracks(track_index).UncertainTips(in_track_index), current_stimulus_along_centerline, ...
                cropped_image, []);
            
%             IWFig = findobj('Tag', ['IWFig', num2str(plotting_index)]);
%             writeVideo(outputVideo, getframe(IWFig));
             writeVideo(outputVideo, getframe(gcf));
        end
        close(outputVideo) 
        video_transcode(video_file_name)
    end
end