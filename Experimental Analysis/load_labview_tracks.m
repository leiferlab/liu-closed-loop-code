function Tracks = load_labview_tracks(folder_name, parameters)
    if exist([folder_name, filesep, 'labview_tracks.mat'], 'file') == 2 && exist([folder_name, filesep, 'timestamps.mat'], 'file') == 2
        loaded_variable = struct2cell(load([folder_name, filesep, 'labview_tracks.mat']));
        Tracks = [loaded_variable{:}];
        if isempty(Tracks)
            return
        end
        [~,sort_idx] = sort([Tracks.WormIndex]);
        Tracks = Tracks(sort_idx);
        tracks_to_throw_out = []; %throw out tracks with less than 3 datapoints
        for track_index = 1:length(Tracks)
            if length(Tracks(track_index).Frames) < 3
                tracks_to_throw_out = [tracks_to_throw_out, track_index];
            end
        end
        Tracks(tracks_to_throw_out) = [];
        %load timestamps
        loaded_variable = load([folder_name, filesep, 'timestamps.mat']);
        % process decoded camera frames to iron out kinks. This should also
        % give us an idea of how many frames we are dropping while delivering stimulus
        timestamps = loaded_variable.timestamps;
        if any(diff(timestamps) < 0)
            error('Timestamps not monotoncially increasing. The experiment likely terminated before completetion');
        end
        %get a list of avaible stimulus frames
        projector_image_directory = [folder_name, filesep, 'ConvertedProjectorFrames', filesep];
        image_files = dir([projector_image_directory, '*.png']); %get all the image files
        for image_file_index = 1:length(image_files)
            file_name = strsplit(image_files(image_file_index).name,'.');
            file_name = file_name{1};
            decoded_camera_frame = strsplit(file_name,'_');
            decoded_camera_frame = str2double(decoded_camera_frame{2});
            image_files(image_file_index).CameraFrame = decoded_camera_frame;
        end
        possible_decoded_camera_frames = sort([image_files.CameraFrame]);
        
        [processed_decoded_camera_frames, invalid_lags_count, dropped_frame_fraction, unmatched_stimulus_frame_count] = process_decoded_camera_frames(loaded_variable.decodedcameraframes, possible_decoded_camera_frames, parameters);
        save([folder_name, filesep, 'timestamps.mat'], 'processed_decoded_camera_frames', '-append');
        invalid_lags_count
        dropped_frame_fraction
        unmatched_stimulus_frame_count
        %change into matlab tracks by interpolating between dropped track
        %frames
%         Tracks = rmfield(Tracks,'Stimulus'); %get rid of the stimulus field, we will recalculate anyway
        %get rid of tracks that are less than 3 frames long
        indecies_to_remove = [];
        for track_index = 1:length(Tracks)
           if length(Tracks(track_index).Frames) < 3
               indecies_to_remove = [indecies_to_remove, track_index];
           end
        end
        Tracks(indecies_to_remove) = [];
        Tracks(1).Time = [];
        for track_index = 1:length(Tracks)
            Tracks(track_index) = labview_track_to_matlab_track(Tracks(track_index), timestamps, processed_decoded_camera_frames);
        end
    else
        Tracks = [];
    end
end