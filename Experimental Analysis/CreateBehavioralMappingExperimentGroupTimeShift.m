% analyzes a group of experiments and saves the properties
% they will be saved inside the first folder; not used in paper
 %function LNPStats = CreateBehavioralMappingExperimentGroup()
    %clear all;
    %set up parameters
    recalculateSpectra = true;

    parameters.numProcessors = 15;
    parameters.numProjections = 19;
    parameters.pcaModes = 5;
    parameters.samplingFreq = 14;
    parameters.minF = 0.3;
    parameters.maxF = parameters.SampleRate ./ 2; %nyquist frequency
    parameters.trainingSetSize = 55000;
    parameters.subsamplingIterations = 10;
    parameters = setRunParameters(parameters);
    
    max_timeshift = round((parameters.omega0 + sqrt(2+parameters.omega0^2))./(4*pi.*parameters.minF)*parameters.samplingFreq);

    
    Prefs = load_excel_prefs;
    load('reference_embedding.mat')
    number_of_behaviors = max(L(:));
  
    [filename,pathname] = uiputfile('*.mat','Save Experiment Group As');
    
    if isequal(filename,0) || isequal(pathname,0)
        %cancel
       return
    else
        saveFileName = fullfile(pathname,filename);
        if exist(saveFileName, 'file')
          % File exists.  Load the folders
          previously_loaded_experiments = load(saveFileName);
          folders = previously_loaded_experiments.folders;
        else
          % File does not exist. Ask for experiment folders
            folders = getfolders();
        end
    end
    
    allTracks = [];
    for folder_index = 1:length(folders)
        %single experiment
        folder_name = folders{folder_index};
        cd(folder_name) %open the directory of image sequence
        load('tracks.mat')
        load('LEDVoltages.txt')
        
        try
            experiment_parameters = load('parameters.txt');
            frames = experiment_parameters(length(experiment_parameters));
        catch
            experiment_parameters = readtable('parameters.txt', 'Delimiter', '\t');
            frames = experiment_parameters{1,{'FrameCount'}};
        end
        
        if recalculateSpectra || ~isfield(Tracks, 'BehavioralTransition')
            %get the spectra
            [Spectra, ~, ~, ~] = generate_spectra(Tracks, parameters, Prefs);
            
            %cut the spectra and Tracks
            for track_index = 1:length(Tracks)
                Spectra{track_index} = Spectra{track_index}(1:end-max_timeshift,:);
                Tracks(track_index) = CutTrackByLocalFrame(Tracks(track_index), 1, size(Spectra{track_index},1));
            end            
            
            data = vertcat(Spectra{:});
            [embeddingValues,~] = findEmbeddings(data,trainingSetData,trainingEmbedding,parameters);
            clear data
     
            % cut the embeddings
            Tracks(1).Embeddings = []; %preallocate memory
            start_index = 1;
            for track_index = 1:length(Spectra)
                end_index = start_index + size(Spectra{track_index},1) - 1;
                Tracks(track_index).Embeddings = embeddingValues(start_index:end_index, :);
                start_index = end_index + 1;
            end

            %get the stereotyped behaviors
            Tracks = find_stereotyped_behaviors(Tracks, L, xx);
        end

        % Get binary array of when certain behaviors start
        Tracks(1).Behaviors = [];
        for track_index = 1:length(Tracks)
            triggers = false(number_of_behaviors, length(Tracks(track_index).LEDVoltages)); %a binary array of when behaviors occur
            for behavior_index = 1:number_of_behaviors
                transition_indecies = Tracks(track_index).BehavioralTransition(:,1) == behavior_index;
                %transition into of
                transition_start_frames = Tracks(track_index).BehavioralTransition(transition_indecies,2);
                triggers(behavior_index,transition_start_frames) = true;
%                 %transition out of
%                 transition_end_frames = Tracks(track_index).BehavioralTransition(transition_indecies,3);
%                 triggers(behavior_index,transition_end_frames) = true;
            end
            Tracks(track_index).Behaviors = triggers;
        end

        if isempty(allTracks)
            allTracks = Tracks;
        else
            allTracks = [allTracks, Tracks];
        end
    end
    
    
    
    %the very last entry in Experiments is the average of all experiments
    %fit the LNP
    [LNPStats, meanLEDPower, stdLEDPower] = FitLNP(allTracks);

    PlotBehavioralMappingExperimentGroup(LNPStats, meanLEDPower, stdLEDPower, L, density, xx);
    save(saveFileName, 'folders', 'LNPStats', 'L', 'density', 'xx');

%  end