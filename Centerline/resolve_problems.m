function [newTracks] = resolve_problems(folder_name, parameters)
% displays a GUI to help user decide how to resolve centerline problems, not used in paper

    min_track_length = parameters.MinTrackLength;
%     plotting_index = 1;
    
    Tracks = load_single_folder(folder_name);
    if isempty(Tracks)
        return
    end
    modifications_index = 1;
    Modifications = [];
    for track_index = 1:length(Tracks)
        Track = Tracks(track_index);
        potential_problems = Track.PotentialProblems;
        if sum(potential_problems) > 0
            loaded_file = load([folder_name, filesep, 'individual_worm_imgs', filesep, 'worm_', num2str(track_index), '.mat']);
            worm_images = loaded_file.worm_images;

            frames_to_show = conv(single(potential_problems), ones(1, 28), 'same'); %show for 2 sec around the problem
            frames_to_show = frames_to_show > 0;
            worm_frame_start_index = 0;
            while worm_frame_start_index <= size(worm_images, 3)
                [worm_frame_start_index, worm_frame_end_index] = find_next_section(frames_to_show, worm_frame_start_index, 'f');
                if isempty(worm_frame_start_index)
                    break
                else
                    %call the gui for resolution
                    h = resolve_problems_gui(worm_images, Track, worm_frame_start_index, worm_frame_end_index, track_index);
                    movegui(h, 'center');
                    uiwait(h);
                    action = h.UserData{7};
                    current_frame = h.UserData{6};
                    close(h);
                    switch action
                        %depending on what the user selected, 
                        case 1
                            %no action: repress the error
                            Track.PotentialProblems(worm_frame_start_index:worm_frame_end_index) = -1;
                        case 2
                            %flip head/tail before
                            Track.Centerlines(:,:,1:current_frame) = flip(Track.Centerlines(:,:,1:current_frame),1);
                        case 3
                            %flip head/tail after
                            Track.Centerlines(:,:,current_frame:end) = flip(Track.Centerlines(:,:,current_frame:end),1);
                        case 4
                            %delete track
                            Modifications(modifications_index).TrackIndex = track_index;
                            Modifications(modifications_index).Action = 1; %delete
                            modifications_index = modifications_index + 1;
                            break
                        case 5
                            %delete section and split track
                            Modifications(modifications_index).TrackIndex = track_index;
                            Modifications(modifications_index).Action = 2;
                            Modifications(modifications_index).StartFrame = Track.Frames(worm_frame_start_index);
                            Modifications(modifications_index).EndFrame = Track.Frames(worm_frame_end_index);
                            modifications_index = modifications_index + 1;
                    end
                end
                worm_frame_start_index = worm_frame_end_index + 1;
            end
        end
    end
    
    if isempty(Modifications)
        newTracks = Tracks;
        return
    end
    
    newTracks = [];
    modification_track_indecies = [Modifications.TrackIndex];
    current_track_index = 1;
    current_end_index = length(Tracks);
    for track_index = 1:length(Tracks)
        Track = Tracks(track_index);
        %get if there are modifications
        if ismember(track_index, modification_track_indecies)
            %there are modifications, get all the modifications
            modifications_to_this_track = Modifications(modification_track_indecies == track_index);
            %see if there is a single delete command
            modifications_to_this_track_actions = [modifications_to_this_track.Action];
            if ismember(1, modifications_to_this_track_actions)
                %there is a delete command, delete the track
                if current_track_index == current_end_index
                    delete([folder_name, filesep, 'individual_worm_imgs', filesep, 'worm_', num2str(current_track_index), '.mat'])
                else
                    rename_individual_worm_images(folder_name,current_track_index+1,current_end_index,-1);
                end
                current_end_index = current_end_index - 1;
            else
                %there are one or many split commands
                split_tracks = [];
                split_modifications = modifications_to_this_track(modifications_to_this_track_actions == 2);
                
                loaded_file = load([folder_name, filesep, 'individual_worm_imgs', filesep, 'worm_', num2str(current_track_index), '.mat']);
                Track.WormImages = loaded_file.worm_images;
                
                for split_modification_index = 1:length(split_modifications)
                    current_split_modification_begin = split_modifications(split_modification_index).StartFrame;
                    current_split_modification_end = split_modifications(split_modification_index).EndFrame;
                    new_subtrack = FilterTracksByTime(Track, Track.Frames(1), current_split_modification_begin);
                    if length(new_subtrack.Frames) >= min_track_length
                        split_tracks = [split_tracks, new_subtrack];
                    end
                    Track = FilterTracksByTime(Track, current_split_modification_end, Track.Frames(end));
                end
                
                if length(Track.Frames) >= min_track_length
                    split_tracks = [split_tracks, Track];
                end                
                
                %shift up by the number of new tracks added
                rename_individual_worm_images(folder_name, current_track_index+1, current_end_index, length(split_tracks)-1);
                current_end_index = current_end_index + length(split_tracks) - 1;
                
                %split the worm frames exactly as before
                for saveindex = 1:length(split_tracks)
                    worm_images = split_tracks(saveindex).WormImages;
                    save([folder_name, filesep, 'individual_worm_imgs', filesep, 'worm_', num2str(current_track_index+saveindex-1), '.mat'], 'worm_images');
                end
                
                split_tracks = rmfield(split_tracks, 'WormImages');
                
                newTracks = [newTracks, split_tracks];
                current_track_index = length(newTracks) + 1;
            end
        else
            %no modifications, add the old track
            newTracks = [newTracks, Track];
            current_track_index = length(newTracks) + 1;
        end
    end
    delete_extra_individual_worm_images(folder_name, length(newTracks))
    
end