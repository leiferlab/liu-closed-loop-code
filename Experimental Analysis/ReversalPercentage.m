function [Percentage, reversal_counts, frame_count]= ReversalPercentage(folders, Tracks, frames)
%not used in paper
    fps = 14;
    allTracks = [];
    relevant_fields = {'Frames','Pirouettes'};

    if nargin < 1 %no folders are given, ask user to select
        folders = {};
        while true
            folder_name = uigetdir
            if folder_name == 0
                break
            else
                folders{length(folders)+1} = folder_name;
            end
        end
    end
    
    if nargin < 3 %Folders are given but no Tracks are given
        for folder_index = 1:length(folders)
            folder_name = folders{folder_index};
            allTracks = loadtracks({folder_name},relevant_fields);
            try
                load('parameters.txt')
                frames = parameters(length(parameters));
            catch
                parameters = readtable('parameters.txt', 'Delimiter', '\t');
                frames = parameters{1,{'FrameCount'}};
            end
        end
    else
        allTracks = Tracks;
    end
    

    %divide bins by bin size
    reversal_counts = zeros(1, ceil(frames));
    frame_count = zeros(1, ceil(frames));
    tracksCentered = [];
    pirouetteCount = 0;

    for track = 1:length(allTracks)
        pirouettes = allTracks(track).Pirouettes;
        frames = allTracks(track).Frames;
        for pirouette_index = 1:size(pirouettes,1)
            pirouetteStart = pirouettes(pirouette_index,1);
            pirouetteEnd = pirouettes(pirouette_index,2);
            
            if pirouetteStart <= length(frames)
                %the cut tracks still holds information beyond the filtered
                %time points, so check for it
                if pirouetteEnd < length(frames)
                    reversal_counts(frames(pirouetteStart):frames(pirouetteEnd)) = reversal_counts(frames(pirouetteStart):frames(pirouetteEnd)) + 1;
                else
                    reversal_counts(frames(pirouetteStart):length(frames)) = reversal_counts(frames(pirouetteStart):length(frames)) + 1;
                end
            end
        end
        for frame_index = 1:length(frames)
            frame_count(frames(frame_index)) = frame_count(frames(frame_index)) + 1;
        end
    end

    %remove all reversals that are counted outside of the frames
    for frame_count_index = 1:length(frame_count)
        if frame_count(frame_count_index) == 0
            reversal_counts(frame_count_index) = 0;
        end
    end
    
    Percentage = reversal_counts./frame_count;
    Percentage(isnan(Percentage)) = 0; %replace nan with 0

    if nargin < 1
        figure
        plot(Percentage, 'bo-')
        %legend(num2str(tracksByVoltage(voltage_index).voltage));
        xlabel(['minutes (', num2str(sum(reversal_counts)), ' reversals analyzed) average reversal rate = ', num2str(sum(reversal_counts)/sum(frame_count)* fps * 60)]) % x-axis label
        ylabel('reversals per worm per min') % y-axis label
        axis([1 frames/bin_size 0 3])
    end
end