function active_track_indecies = PlotFrame(FigH, Frame, Tracks, frame_index, all_deleted_tracks, behavior_colors)

figure(FigH)
clf;
%imshow(Frame,'InitialMagnification',300);
imshow(Frame,'InitialMagnification',40, 'Border','tight'); %this determines the figure size that the mpeg encodes, which has to be a particular multiple
hold on;
active_track_indecies = [];

%plot after analysis
if ~isempty(Tracks)
    track_indecies_in_frame = find(vertcat(Tracks.Frames) == frame_index);
    frameSum = 0;
    currentActiveTrack = 1; %keeps the index of the track_indecies_in_frame
%     myColors = winter(length(track_indecies_in_frame));
    for track_index = 1:length(Tracks)
        if currentActiveTrack > length(track_indecies_in_frame)
            %all active tracks found
            break;
        end
        if track_indecies_in_frame(currentActiveTrack) - frameSum <= length(Tracks(track_index).Frames) 
            %active track found
            in_track_index = track_indecies_in_frame(currentActiveTrack) - frameSum;
            current_behavior = Tracks(track_index).BehavioralAnnotation(in_track_index);
%             if current_behavior < 1 || current_behavior > length(behavior_colors)
                current_behavior_color = [1 1 0];
%             else
%                 current_behavior_color = behavior_colors(Tracks(track_index).BehavioralAnnotation(in_track_index),:);
%             end
            plot(Tracks(track_index).Path(1:in_track_index,1), Tracks(track_index).Path(1:in_track_index,2), 'Color', current_behavior_color);
            plot(Tracks(track_index).Path(in_track_index,1), Tracks(track_index).Path(in_track_index,2),'x' , 'Color', current_behavior_color);
            text(double(Tracks(track_index).Path(in_track_index,1))+10, double(Tracks(track_index).Path(in_track_index,2))+10, num2str(track_index), 'Color', current_behavior_color);
            currentActiveTrack = currentActiveTrack + 1;
            active_track_indecies = [active_track_indecies, track_index];
        end
        frameSum = frameSum + length(Tracks(track_index).Frames);
    end
end

%deleted tracks plotting
if ~isempty(all_deleted_tracks)
    track_indecies_in_frame = find(vertcat(all_deleted_tracks.Frames) == frame_index);
    frameSum = 0;
    currentActiveTrack = 1; %keeps the index of the track_indecies_in_frame
    myColors = autumn(length(track_indecies_in_frame));
    for track_index = 1:length(all_deleted_tracks)
        if currentActiveTrack > length(track_indecies_in_frame)
            %all active tracks found
            break;
        end
        if track_indecies_in_frame(currentActiveTrack) - frameSum <= length(all_deleted_tracks(track_index).Frames) 
            %active track found
            in_track_index = track_indecies_in_frame(currentActiveTrack) - frameSum;
            if size(all_deleted_tracks(track_index).Path,2) <= 1
                %tracks that are only 1 frame
                plot(all_deleted_tracks(track_index).Path(1), all_deleted_tracks(track_index).Path(2),'o' , 'Color', myColors(currentActiveTrack,:));
                if ~isempty(all_deleted_tracks(track_index).DeletionReason)
                    text(double(all_deleted_tracks(track_index).Path(1))+10, double(all_deleted_tracks(track_index).Path(2))+10, all_deleted_tracks(track_index).DeletionReason, 'Color', myColors(currentActiveTrack,:));
                end
            else
                plot(all_deleted_tracks(track_index).Path(1:in_track_index,1), all_deleted_tracks(track_index).Path(1:in_track_index,2), 'Color', myColors(currentActiveTrack,:));
                plot(all_deleted_tracks(track_index).Path(in_track_index,1), all_deleted_tracks(track_index).Path(in_track_index,2),'o' , 'Color', myColors(currentActiveTrack,:));
                if ~isempty(all_deleted_tracks(track_index).DeletionReason)
                    text(double(all_deleted_tracks(track_index).Path(in_track_index,1))+10, double(all_deleted_tracks(track_index).Path(in_track_index,2))+10, all_deleted_tracks(track_index).DeletionReason, 'Color', myColors(currentActiveTrack,:));
                end
            end

            currentActiveTrack = currentActiveTrack + 1;
            active_track_indecies = [active_track_indecies, track_index];
        end
        frameSum = frameSum + length(all_deleted_tracks(track_index).Frames);
    end
end

axis tight
hold off;    % So not to see movie replay