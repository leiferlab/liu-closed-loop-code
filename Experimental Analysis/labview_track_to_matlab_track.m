function Track = labview_track_to_matlab_track(Track, timestamps, decoded_camera_frames)
    % converts a labview track to a matlab track by interpolating between
    % dropped tracking frames. Ensure the first dimension is matching
    % frames
    
    Track.WormIndex = Track.WormIndex + 1; %change to matlab indexing
    
    %correct for dropped frames during tracking by interpolating
    labview_tracking_frames = double(Track.Frames); 
    if length(labview_tracking_frames) < 3
        Track = [];
        return
    end
    %remove the weird thing where the last few camera buffer index is 0
    anomaly_index = find(diff(labview_tracking_frames) < 0, 1); %locate where trames start to decrease
    if ~isempty(anomaly_index)
        labview_tracking_frames(anomaly_index:end) = (labview_tracking_frames(anomaly_index):labview_tracking_frames(anomaly_index)+numel(labview_tracking_frames(anomaly_index:end-1)))';
    end
    labview_tracking_frames = labview_tracking_frames + 1; %change to matlab indexing
    labview_tracking_frame_count = length(labview_tracking_frames);
    last_track_frame = min(labview_tracking_frames(end),length(timestamps)-1); %ensure last track frame camera frame is within the experiment
    new_track_frames = labview_tracking_frames(1):last_track_frame;
    
%     %tracking stutters bug in labview..
%     non_duplicate_tracking_frames = diff(labview_tracking_frames)~=0;
%     labview_tracking_frames = labview_tracking_frames(non_duplicate_tracking_frames);

    if labview_tracking_frame_count >= 3
        field_names = fieldnames(Track);
        for field_name_index = 1:length(field_names)
            field_name = field_names{field_name_index};
            if strcmp(field_name, 'Frames')
                Track.Frames = new_track_frames'; % keep Frames int32
            else
                current_field_value = single(Track.(field_name));
                current_field_dimension_sizes = size(current_field_value);
                interpolation_dimension = find(current_field_dimension_sizes == labview_tracking_frame_count, 1);
                if ~isempty(interpolation_dimension)
                    current_field_value = shiftdim(current_field_value,interpolation_dimension-1); %make the time dimension first
                    new_field_dimension_sizes = size(current_field_value);
                    new_field_dimension_sizes(1) = length(new_track_frames);
                    current_field_value = reshape(current_field_value, labview_tracking_frame_count, []);
                    %first dimension to be interpolated
                    % x is track_frames
                    % v is current_field_values
                    % xq is evenly spaced frames from track_frames(1) to
                    % labview_tracking_frames(end) or last camera frame, whichever is smaller
                    % iterate over 2nd dim
                    new_current_field_value = zeros(length(new_track_frames),size(current_field_value,2),'single');
                    for idx = 1:size(current_field_value,2)
                        %remove tracking stutters and interpolate
                        new_current_field_value(:,idx) = interp1(double(labview_tracking_frames),current_field_value(:,idx),double(new_track_frames));
                    end
                    % reshape back to original shape (correcting for last frame)
                    new_current_field_value = reshape(new_current_field_value, new_field_dimension_sizes);
                    Track.(field_name) = new_current_field_value;
                end
            end
        end        
    end
    
    %pad the stimulus for when the track is immature, correct for stimulation dropped frames
    if ~isempty(Track.Stimulus)
        Track.Stimulus = Track.Stimulus';
        %pad the immature period with 0 stimulus
        padding_before = zeros(Track.MatureCameraBuffer - int32(new_track_frames(1)) - 1, size(Track.Stimulus,2),'single');
        Track.Stimulus = vertcat(padding_before,Track.Stimulus);
        padding_after = zeros(length(new_track_frames)-size(Track.Stimulus,1), size(Track.Stimulus,2),'single');
        Track.Stimulus = vertcat(Track.Stimulus,padding_after);
        Track.Stimulus = Track.Stimulus(1:length(new_track_frames),:);
        %drop frame correction starts from mature camera buffer
        % process first, filter out weirdness
        frame_indecies_for_stim = decoded_camera_frames(new_track_frames) - int32(new_track_frames(1)) + 1;  %the frame indecies for stimulus, instead of camera
        frame_indecies_for_stim(frame_indecies_for_stim <= 0) = 1;
        frame_indecies_for_stim(frame_indecies_for_stim > size(Track.Stimulus,1)) = size(Track.Stimulus,1);
        Track.Stimulus = Track.Stimulus(frame_indecies_for_stim,:); % correct the stimulus by shifting
    end

    %generate the Time field
    Track.Time = timestamps(Track.Frames);
end