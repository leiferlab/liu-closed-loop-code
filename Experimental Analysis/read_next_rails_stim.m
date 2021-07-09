function [first_stim_frame, stim_value, stim_duration] = read_next_rails_stim(stim_time_series, start_frame, end_frame)
%given a rails timeseries and a starting point, extract the beginning and
%duration of the next pulse
    end_frame = min(end_frame, numel(stim_time_series));
    [~, stim_peaks, widths, ~] = findpeaks(abs(stim_time_series(start_frame:end_frame)));
    if ~isempty(stim_peaks)
        first_stim_frame = stim_peaks(1) + start_frame - 1;
        stim_value = stim_time_series(first_stim_frame);
        stim_duration = widths(1);
    else
        first_stim_frame = start_frame + 1;
        stim_value = stim_time_series(first_stim_frame);
        stim_duration = 0;        
    end
end

