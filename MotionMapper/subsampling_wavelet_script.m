% not used in paper

% Projections = {Tracks.ProjectedEigenValues};
% Projections = cellfun(@transpose, Projections, 'UniformOutput', false);

%% get the experiment folders
folders = [];
while true
    if isempty(folders)
        start_path = '';
    else
        start_path = fileparts(fullfile(folders{length(folders)}, '..', 'tracks.mat')); %display the parent folder
    end
    folder_name = uigetdir(start_path, 'Select Experiment Folder')
    if folder_name == 0
        break
    else
        folders{length(folders)+1} = folder_name;
    end
end

%% load the tracks
allTracks = struct([]);
folder_indecies = [];
track_indecies = [];

for folder_index = 1:length(folders)
    curDir = folders{folder_index};
    if exist([curDir, '\tracks.mat'], 'file') == 2
        load([curDir, '\tracks.mat'])
        allTracks = [allTracks, Tracks];
        folder_indecies = [folder_indecies, repmat(folder_index,1,length(Tracks))];
        track_indecies = [track_indecies, 1:length(Tracks)];
    end
end

Projections = {allTracks.ProjectedEigenValues};
clear('Tracks');

%% set up parameters
parameters.numProcessors = 7;
parameters.numProjections = 19;
parameters.pcaModes = 5;
parameters.samplingFreq = 14;
parameters.minF = 0.3;
parameters.maxF = 7;
parameters.trainingSetSize = 25000;
parameters = setRunParameters(parameters);

L = length(Projections);


%% generate spectra
poolobj = gcp('nocreate'); 
if isempty(poolobj)
    parpool(7)
end
Spectra = cell(size(Projections));
for track_index = 1:L
    [Spectra{track_index},f] = findWavelets(Projections{track_index}',parameters.pcaModes,parameters);  
end
poolobj = gcp('nocreate'); 
delete(poolobj);
f = fliplr(f);

% plot_data = flipud(Spectra{2}');
% pcaSpectra = flipud(mat2cell(plot_data, repmat(parameters.numPeriods, 1, parameters.pcaModes)));
% %pcaSpectra{5} = pcaSpectra{2} - pcaSpectra{3};
% figure
% for i = 1:length(pcaSpectra)
%     subplot(length(pcaSpectra), 1, i)
%     imagesc(pcaSpectra{i});
%     ax = gca;
%     ax.YTick = 1:5:parameters.numPeriods;
%     ax.YTickLabel = num2cell(round(f(mod(1:length(f),5) == 1), 1));
%     ylabel({['PCA Mode ', num2str(i)], 'Frequency (Hz)'});
%     
%     ax.XTickLabel = round(ax.XTick/parameters.samplingFreq, 1);
%     
%     if i == length(pcaSpectra)
%         xlabel('Time (s)');
%     end
% end
%% Use subsampled t-SNE to find the training set
data = vertcat(Spectra{:});
fprintf(1, 'Finding Training Set\n');
[trainingSetData, trainingSetAmps, trainingKey] = ...
    runEmbeddingSubSampling(data, parameters);

%% Embed training set
% amps = sum(data,2);
% data(:) = bsxfun(@rdivide,data,amps);
% 
% skipLength = round(length(data(:,1))/parameters.trainingSetSize);
% 
% trainingSetData = data(skipLength:skipLength:end,:);
% trainingAmps = amps(skipLength:skipLength:end);
% trainingKey = 1:size(data,1);
% trainingKey = trainingKey(skipLength:skipLength:end);
% parameters.signalLabels = log10(trainingAmps);

fprintf(1,'Finding t-SNE Embedding for Training Set\n');
[trainingEmbedding,betas,P,errors] = run_tSne(trainingSetData,parameters);


%% Find All Embeddings

fprintf(1,'Finding t-SNE Embedding for all Data\n');
% embeddingValues = cell(L,1);
% i=1;

[embeddingValues,~] = findEmbeddings(data,trainingSetData,trainingEmbedding,parameters); %[embeddingValues{i},~]


%% cut the embeddings and generate trainingTracks and trainingFrames
Embeddings = cell(size(Spectra));
trainingTracks = zeros(length(trainingKey),1);
trainingFrames = zeros(length(trainingKey),1);
start_index = 1;
training_start_index = 1;
for track_index = 1:length(Spectra)
    end_index = start_index + size(Spectra{track_index},1) - 1;
    Embeddings{track_index} = embeddingValues(start_index:end_index, :);
    dataIndecies = trainingKey(trainingKey >= start_index & trainingKey <= end_index);
    training_end_index = training_start_index + length(dataIndecies) - 1;
    trainingTracks(training_start_index:training_end_index) = track_index;
    trainingFrames(training_start_index:training_end_index) = dataIndecies - start_index + 1;
    
    training_start_index = training_end_index + 1;
    start_index = end_index + 1;
end
clear trainingKey

%% Make density plots
maxVal = max(max(abs(combineCells(Embeddings))));
maxVal = round(maxVal * 1.1);

sigma = maxVal / 40;
numPoints = 501;
rangeVals = [-maxVal maxVal];

[xx,density] = findPointDensity(combineCells(Embeddings),sigma,numPoints,rangeVals);
maxDensity = max(density(:));

% densities = zeros(numPoints,numPoints,L);
% for i=1:L
%     [~,densities(:,:,i)] = findPointDensity(Embeddings{i},sigma,numPoints,rangeVals);
% end


for track_index = 1:length(allTracks);
    plot_embedding = Embeddings{track_index};
    image_file = fullfile([folders{folder_indecies(track_index)}, '\individual_worm_imgs\worm_', num2str(track_indecies(track_index)), '.mat']);
    save_file = fullfile([folders{folder_indecies(track_index)}, '\individual_worm_imgs\behaviormap_', num2str(track_indecies(track_index))]);
    load(image_file);

    behavior_figure = figure('Position', [100, 100, 500, 250]);
    outputVideo = VideoWriter(save_file,'MPEG-4');
    outputVideo.FrameRate = 14;
    open(outputVideo)
    for worm_frame_index = 1:size(plot_embedding, 1)
        subplot_tight(1,2,2,0);

        plot_worm_frame(worm_images(:,:,worm_frame_index), squeeze(allTracks(track_index).Centerlines(:,:,worm_frame_index)), ...
        allTracks(track_index).UncertainTips(worm_frame_index), ...
        allTracks(track_index).Eccentricity(worm_frame_index), allTracks(track_index).Direction(worm_frame_index), ...
        allTracks(track_index).Speed(worm_frame_index),  allTracks(track_index).TotalScore(worm_frame_index), 0);

        freezeColors

        subplot_tight(1,2,1,0);
        hold on
        imagesc(xx,xx,density)
        axis equal tight off xy
        caxis([0 maxDensity * .8])
        colormap(jet)
        plot(plot_embedding(worm_frame_index,1), plot_embedding(worm_frame_index,2), '*k')
        hold off

        writeVideo(outputVideo, getframe(gcf));
        clf
        %colorbar
    end
    close(outputVideo)
    close(behavior_figure)
end

% figure
% 
% N = ceil(sqrt(L));
% M = ceil(L/N);
% maxDensity = max(densities(:));
% for i=1:L
%     subplot(M,N,i)
%     imagesc(xx,xx,densities(:,:,i))
%     axis equal tight off xy
%     caxis([0 maxDensity * .8])
%     colormap(jet)
%     title(['Data Set #' num2str(i)],'fontsize',12,'fontweight','bold');
% end
% 
