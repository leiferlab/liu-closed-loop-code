function Tracks = FilterTracksBySize(Tracks, min_size, max_size)

    selected_indecies = [];
    for track_index = 1:length(Tracks)
       worm_size = mean(Tracks(track_index).Size); 
       if worm_size > min_size && worm_size < max_size
          selected_indecies = [selected_indecies, track_index]; 
       end
    end

    Tracks = Tracks(selected_indecies);

end