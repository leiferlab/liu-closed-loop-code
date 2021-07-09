function subsample_coord = SubSampleTrack(coord, time, fps, interpolate)
% not used in paper

    subsample_coord = zeros(floor(max(time)*fps),2);
    time_index = 2;
    for i = 1:length(subsample_coord)
        curr_time = i / fps;
        while time(time_index) < curr_time
            time_index = time_index + 1;
        end
        if time(time_index) == curr_time
            %exact time
            subsample_coord(i,:) = coord(time_index,:);
        else
            if interpolate
                %interpolate
                subsample_coord(i,:) = coord(time_index-1,:) + ((coord(time_index,:)-coord(time_index-1,:)) ./ (time(time_index)-time(time_index-1))*(curr_time-time(time_index-1)));
            else
                %round
                subsample_coord(i,:) = round(coord(time_index-1,:) + ((coord(time_index,:)-coord(time_index-1,:)) ./ (time(time_index)-time(time_index-1))*(curr_time-time(time_index-1))),0);
            end
        end
    end
end