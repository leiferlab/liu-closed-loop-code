function [decoded_camera_frames, invalid_lags_count, dropped_frame_fraction, unmatched_stimulus_frame_count] = process_decoded_camera_frames(decoded_camera_frames, possible_decoded_camera_frames, parameters)
% process the decoded camera frames to ensure no weirdness.
% enforce rules:
% 1. the decodeded camera frame cannot be further ahead in time than its own index (i.e. collected camera frame)
% 2. the decoded camera frame cannot lag x seconds behind than its own index
% 3. the decoded camera frame must be monotonic
% resolution by replacing with NaNs and then finally interpolating through NaNs
% if beginning and end are NaNs, replace with lag of 1 frame / frame from the closest known good frame
    decoded_camera_frames = double(decoded_camera_frames)';
    max_lag_in_frames = parameters.MaxStimulusLag*parameters.SampleRate;
    camera_frames = 0:numel(decoded_camera_frames)-1;
    lag = camera_frames - decoded_camera_frames;
    invalid_lags = or(lag > max_lag_in_frames, lag <= 0);
    %make sure that decoded_camera_frames is monotonically increasing
    current_highest_number = 0;
    for index = 1:length(decoded_camera_frames)
        if ~invalid_lags(index)
            % this frame is a valid lag
            if decoded_camera_frames(index) < current_highest_number
                invalid_lags(index) = true;
                % the decoded camera frame is decreasing, it is not valid
            else
                % the decoded camera frame is increasing, update it 
                current_highest_number = decoded_camera_frames(index);
            end
        end
    end
    invalid_lags_count = sum(invalid_lags);
    decoded_camera_frames(invalid_lags) = nan;
    if isnan(decoded_camera_frames(1))
        %make lag of 1 if first frame is missing
        decoded_camera_frames(1) = 0;
    end
    if isnan(decoded_camera_frames(end))
        %find the lag at the last non-nan element, and use that lag as the last frame lag if it is missing
        reversed_decoded_camera_frames = fliplr(decoded_camera_frames);
        reversed_lag = fliplr(lag);
        decoded_camera_frames(end) = numel(decoded_camera_frames) - reversed_lag(find(~isnan(reversed_decoded_camera_frames),1));
    end   
    
    %interpolate through nans
    nanx = isnan(decoded_camera_frames);
    t    = 1:numel(decoded_camera_frames);
    decoded_camera_frames(nanx) = round(interp1(t(~nanx), decoded_camera_frames(~nanx), t(nanx)));
    
    %correct any stimulus frames that do not exist
    matched_frames = ismember(decoded_camera_frames,possible_decoded_camera_frames);
    unmatched_stimulus_frame_indecies = find(matched_frames == 0);
    unmatched_stimulus_frame_count = numel(unmatched_stimulus_frame_indecies);
    for unmatched_stimulus_frame_indecies_index = 1:unmatched_stimulus_frame_count
        %match the unmatched frame to the nearest matched frame
        current_frame_index = unmatched_stimulus_frame_indecies(unmatched_stimulus_frame_indecies_index);
        closest_backwards = find(matched_frames(1:current_frame_index-1),1,'last');
        closest_forwards = find(matched_frames(current_frame_index+1:end),1,'first') + current_frame_index;
        if isempty(closest_backwards)
            replacement_index = closest_forwards;
        elseif isempty(closest_forwards)
            replacement_index = closest_backwards;
        elseif current_frame_index - closest_backwards > closest_forwards - current_frame_index
            replacement_index = closest_forwards;
        else
            replacement_index = closest_backwards;
        end
        decoded_camera_frames(current_frame_index) = decoded_camera_frames(replacement_index);
    end
    
    
    %calculate dropped frames
    dropped_frames_counting = diff(decoded_camera_frames) - 1;
    dropped_frames_counting(dropped_frames_counting < 0) = 0;
    dropped_frames_count = sum(dropped_frames_counting);
    dropped_frame_fraction = dropped_frames_count / numel(decoded_camera_frames);
    
    %return decoded_camera_frames
    decoded_camera_frames = int32(decoded_camera_frames);
end