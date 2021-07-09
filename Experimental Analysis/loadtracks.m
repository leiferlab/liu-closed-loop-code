function [ allTracks, folder_indecies, track_indecies ] = loadtracks(folders, field_names)
%get the tracks given folders
    if nargin < 2
        field_names = {};
    end
    
    allTracks = struct([]);
    folder_indecies = [];
    track_indecies = [];

    if iscell(folders)
        for folder_index = 1:length(folders)
            folder_name = folders{folder_index};
            Tracks = load_single_folder(folder_name, field_names);
            folder_indecies = [folder_indecies, repmat(folder_index,1,length(Tracks))];
            track_indecies = [track_indecies, 1:length(Tracks)];
            
            % add info about the track's folder and within folder index
            Tracks(1).folder_index = [];
            Tracks(1).within_folder_track_index = [];
            for track_index = 1:length(Tracks)
                Tracks(track_index).folder_index = folder_index;
                Tracks(track_index).within_folder_track_index = track_index;
            end
            
            allTracks = [allTracks, Tracks];
        end
    else
        allTracks = load_single_folder(folders, field_names);
        folder_indecies = ones(1,length(allTracks));
        track_indecies = 1:length(allTracks);
    end
    
    
end

