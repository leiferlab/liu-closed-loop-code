function Experiments = TimeSegmentCompare()
%TimeSegmentCompare compares the LNP fits from different parts of
%experiments, not used in paper

%   Detailed explanation goes here
    fps = 14;
    EarlyStartFrame = 1;
    EarlyEndFrame = 10*60*fps;
    LateStartFrame = 20*60*fps+1;
    LateEndFrame = 30*60*fps;
    frames = 0;
    
    folders = {};
    allTracks = [];
    
    [filename,pathname] = uigetfile('*.mat','Select Experiment Group');
    
    if isequal(filename,0) || isequal(pathname,0)
        %cancel
       return
    else
        openFileName = fullfile(pathname,filename);
        if exist(openFileName, 'file')
          % File exists.  Load the folders
          load(openFileName)
          folders = {Experiments(1:end-1).Folder};
        end
    end
    
    Experiments = [];
    
    %load the tracks
    for folder_index = 1:length(folders)
        folder_name = folders{folder_index};
        cd(folder_name) %open the directory of image sequence
        load('tracks.mat')

        if length(allTracks) == 0
            allTracks = Tracks;
        else
            allTracks = [allTracks, Tracks];
        end
        try
            load('parameters.txt')
            frames = max(parameters(length(parameters)),frames);
        catch
            parameters = readtable('parameters.txt', 'Delimiter', '\t');
            frames = max(parameters{1,{'FrameCount'}},frames);
        end
    end
    
    % Get binary array of when certain behaviors start and the predicted
    % behavioral rate
    for track_index = 1:length(allTracks)
        pirouettes = allTracks(track_index).Pirouettes;
        behaviors = zeros(1, length(allTracks(track_index).LEDVoltages)); %a binary array of when behaviors occur
        for pirouette_index = 1:size(pirouettes,1)
            pirouetteStart = pirouettes(pirouette_index,1);
            behaviors(pirouetteStart) = 1;
        end
        allTracks(track_index).Behaviors = logical(behaviors);
    end
    
    EarlyTracks = FilterTracksByTime(allTracks, EarlyStartFrame, EarlyEndFrame);
    LateTracks = FilterTracksByTime(allTracks, LateStartFrame, LateEndFrame);
    
    %fit the LNP
    [Experiments(1).linear_kernel, Experiments(1).non_linearity_fit, Experiments(1).BTA, Experiments(1).meanLEDVoltage, Experiments(1).pirouetteCount, Experiments(1).bin_edges, Experiments(1).filtered_signal_histogram, Experiments(1).filtered_signal_given_reversal_histogram] = FitLNP(EarlyTracks);
    [Experiments(2).linear_kernel, Experiments(2).non_linearity_fit, Experiments(2).BTA, Experiments(2).meanLEDVoltage, Experiments(2).pirouetteCount, Experiments(2).bin_edges, Experiments(2).filtered_signal_histogram, Experiments(2).filtered_signal_given_reversal_histogram] = FitLNP(LateTracks);
    [Experiments(3).linear_kernel, Experiments(3).non_linearity_fit, Experiments(3).BTA, Experiments(3).meanLEDVoltage, Experiments(3).pirouetteCount, Experiments(3).bin_edges, Experiments(3).filtered_signal_histogram, Experiments(3).filtered_signal_given_reversal_histogram] = FitLNP(allTracks);
    
    %get speed
    [Experiments(1).Speed, Experiments(1).speed_sum, Experiments(1).frame_count] = SpeedHistogram([],fps*60,EarlyTracks,frames);
    [Experiments(2).Speed, Experiments(2).speed_sum, Experiments(2).frame_count] = SpeedHistogram([],fps*60,LateTracks,frames);
    [Experiments(3).Speed, Experiments(3).speed_sum, Experiments(3).frame_count] = SpeedHistogram([],fps*60,allTracks,frames);

    %get reversal rate
    [Experiments(1).ReversalRate, Experiments(1).ReversalCounts, Experiments(1).FrameCounts] = ReversalRate([],fps*60,EarlyTracks,frames);
    [Experiments(2).ReversalRate, Experiments(2).ReversalCounts, Experiments(2).FrameCounts] = ReversalRate([],fps*60,LateTracks,frames);
    [Experiments(3).ReversalRate, Experiments(3).ReversalCounts, Experiments(3).FrameCounts] = ReversalRate([],fps*60,allTracks,frames);
    
%     %save LEDVoltages
%     Experiments(folder_index).LEDVoltages = LEDVoltages;
%     %find the filtered signal (not used in LNP fitting)
%     filtered_signal = conv(Experiments(folder_index).LEDVoltages, Experiments(folder_index).linear_kernel);
%     Experiments(folder_index).FilteredSignal = filtered_signal(1:length(Experiments(folder_index).LEDVoltages)); %cut off the tail

    
    PlotExperimentGroup(Experiments);
end