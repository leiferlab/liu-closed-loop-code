%% STEP 1: set up parameters
subsampling = true;
relevant_track_fields = {'Amps','Spectra'};

%% STEP 2: get the experiment folders
folders = getfoldersGUI();

%% STEP 2: Load the analysis preferences %%
parameters = load_parameters();
edgeEffectTime = round(sqrt(1/parameters.minF)*parameters.SampleRate);

%% STEP 3: load the tracks into memory
[allTracks, folder_indecies, track_indecies] = loadtracks(folders, relevant_track_fields);

%% STEP 4: convert tracks to cell format
L = length(allTracks);
Spectra = cell(1,L);
amps = cell(1,L);
SpectraFrames = cell(1,L);
SpectraTracks = cell(1,L);
for track_index = 1:length(allTracks)
    Spectra{track_index} = allTracks(track_index).Spectra;
    amps{track_index} = allTracks(track_index).Amps;
    SpectraFrames{track_index} = 1:size(Spectra{track_index},1);
    SpectraTracks{track_index} = repmat(track_index,1,size(Spectra{track_index},1));
end

clear allTracks


%% STEP 5: Get a set of "training spectra" without edge effects
TrainingSpectra = cell(1,L);
TrainingSpectraFrames = cell(1,L);
TrainingSpectraTracks = cell(1,L);
TrainingAmps = cell(1,L);

for track_index = 1:L
    TrainingSpectra{track_index} = Spectra{track_index}(edgeEffectTime:end-edgeEffectTime,:);
    TrainingSpectraFrames{track_index} = SpectraFrames{track_index}(edgeEffectTime:end-edgeEffectTime);
    TrainingSpectraTracks{track_index} = SpectraTracks{track_index}(edgeEffectTime:end-edgeEffectTime);  
    TrainingAmps{track_index} = amps{track_index}(edgeEffectTime:end-edgeEffectTime); 
    Spectra{track_index} = []; %optional clearing of memory
end

clear Spectra amps

%% STEP 6A: initialize training input
training_input_data = vertcat(TrainingSpectra{:}); %these timpoints will be randomly sampled from
clear TrainingSpectra
training_input_frames = [TrainingSpectraFrames{:}];
clear TrainingSpectraFrames
training_input_tracks = [TrainingSpectraTracks{:}];
clear TrainingSpectraTracks
training_amps = vertcat(TrainingAmps{:}); 
clear TrainingAmps

%% STEP 6B Option 1: Find training set by sampling uniformly
if ~subsampling
    skipLength = round(length(training_input_data(:,1))/parameters.trainingSetSize);
    trainingSetData = training_input_data(skipLength:skipLength:end,:);
    trainingSetAmps = training_amps(skipLength:skipLength:end);
    trainingSetFrames = training_input_frames(skipLength:skipLength:end);
    trainingSetTracks = training_input_tracks(skipLength:skipLength:end);
    clear training_input_data
else
%% STEP 6B Option 2: Find training set by embedding several iterations

    %save workspace
    save('temp.mat','-regexp','^(?!(training_input_data|training_amps|training_input_frames|training_input_tracks)$).','-v7.3')
    clearvars('-except', 'subsampling', 'training_input_data', 'training_amps', 'training_input_frames', 'training_input_tracks', 'parameters')

    %constants
    N = parameters.trainingSetSize;
    iterations = parameters.subsamplingIterations;
    numPerDataSet = round(N/iterations);

    if iterations*N > size(training_input_data, 1)
        error('too many t-SNE iterations, not enough data')
    end

    trainingSetData = zeros(numPerDataSet*iterations,size(training_input_data,2));
    trainingSetAmps = zeros(numPerDataSet*iterations,1);
    trainingSetTracks = zeros(numPerDataSet*iterations,1);
    trainingSetFrames = zeros(numPerDataSet*iterations,1);

    for iteration_index = 1:iterations
        fprintf(1,['Finding training set contributions from data set #' ...
        num2str(iteration_index) '\n']);

        currentIdx = (1:numPerDataSet) + (iteration_index-1)*numPerDataSet;

        %randomly sample without replacement from our data
        [iteration_data, sampled_indecies] = datasample(training_input_data,N,1,'replace',false);
        iteration_amps = training_amps(sampled_indecies,:);

        %save data
        save('training_input_data.mat', 'training_input_data', 'training_amps', '-v7.3')
        clear training_input_data training_amps

        %run t-SNE embedding
        [iterationEmbedding,~,~,~] = run_tSne(iteration_data,parameters);

        %re-load data
        load('training_input_data.mat');

        %find the templates
        [trainingSetData(currentIdx,:),trainingSetAmps(currentIdx),selectedIndecies] = ...
        findTemplatesFromData(iteration_data,iterationEmbedding,iteration_amps,...
                            numPerDataSet,parameters,sampled_indecies);
        trainingSetFrames(currentIdx) = training_input_frames(selectedIndecies);
        trainingSetTracks(currentIdx) = training_input_tracks(selectedIndecies);

        %delete the points because we are sampling without replacement
        training_input_data(sampled_indecies, :) = [];
        training_amps(sampled_indecies, :) = [];
        training_input_frames(sampled_indecies) = [];
        training_input_tracks(sampled_indecies) = [];
    end
end

%clean memory
%delete('training_input_data.mat')
clear iteration_data iteration_amps sampled_indecies iterationEmbedding amps
clear training_input_data training_input_frames training_input_tracks training_input_amps

%% STEP 7: Embed the training set 
%clear memory
clear TrainingSpectra TrainingSpectraFrames TrainingSpectraTracks

parameters.signalLabels = log10(trainingSetAmps);

fprintf(1,'Finding t-SNE Embedding for Training Set\n');
[trainingEmbedding,betas,P,errors] = run_tSne(trainingSetData,parameters);

if subsampling
    %reload
    load('temp.mat')
    delete('temp.mat')
end

%% STEP 8: Find All Embeddings
fprintf(1,'Finding t-SNE Embedding for all Data\n');
% embeddingValues = cell(L,1);
% i=1;

[allTracks, folder_indecies, track_indecies] = loadtracks(folders, relevant_track_fields);

L = length(allTracks);
Spectra = cell(1,L);
for track_index = 1:length(allTracks)
    Spectra{track_index} = allTracks(track_index).Spectra;
end

data = vertcat(Spectra{:});
clear Spectra allTracks
% phi_dt = data(:,end); %get phase velocity
% % phi_dt = phi_dt - min(phi_dt) + eps; % make all values non-zero positive
% % phi_dt = phi_dt ./ max(phi_dt); %normalize to 1
% phi_dt = phi_dt ./ parameters.pcaModes; % weigh the phase velocity as a PCA mode (1/5)
% 
% % normalize the phase velocity
% data = data(:,1:end-1);
% amps = sum(data,2);
% data(:) = bsxfun(@rdivide,data,amps);
% data = [data, phi_dt];
% 
% amps = sum(data,2);
% data(:) = bsxfun(@rdivide,data,amps);

[embeddingValues,outputStatistics] = findEmbeddings(data,trainingSetData,trainingEmbedding,parameters); %[embeddingValues{i},~]
clear data

[allTracks, folder_indecies, track_indecies] = loadtracks(folders, relevant_track_fields);

%% STEP 9: cut the embeddings
Embeddings = cell(1,length(allTracks));
% trainingTracks = zeros(length(trainingKey),1);
% trainingFrames = zeros(length(trainingKey),1);
start_index = 1;
training_start_index = 1;
for track_index = 1:length(allTracks)
    end_index = start_index + size(allTracks(track_index).Spectra,1) - 1;
    Embeddings{track_index} = embeddingValues(start_index:end_index, :);
    start_index = end_index + 1;
end
clear trainingKey

%% STEP 10: Find watershed regions and make density plot
maxVal = max(max(abs(embeddingValues)));
maxVal = round(maxVal * 1.1);

% sigma = maxVal / 40; %change smoothing factor if necessary
sigma = 4.3; %median(D(:,kdNeighbors+1)); %change smoothing factor if necessary
numPoints = 501;
rangeVals = [-maxVal maxVal];

[xx,density] = findPointDensity(embeddingValues,sigma,501,[-maxVal maxVal]);
maxDensity = max(density(:));
density(density < 10e-6) = 5;
L = watershed(-density,8);
[ii,jj] = find(L==0);

L(L==1) = max(L(:))+1;
L = L - 1;

watershed_centroids = regionprops(L, 'centroid');
watershed_centroids = vertcat(watershed_centroids.Centroid);
watershed_centroids = round(watershed_centroids);


density(density == 5) = 0;
%modify jet map
my_colormap = othercolor('OrRd9');
my_colormap(1,:) = [1 1 1];

figure
hold on
imagesc(xx,xx,density)
caxis([0 maxDensity * .6])
colormap(my_colormap)
plot(xx(jj),xx(ii),'k.')
% for region_index = 1:size(watershed_centroids,1)-1
%     text(xx(watershed_centroids(region_index,1)), ...
%         xx(watershed_centroids(region_index,2)), ...
%         num2str(region_index), 'color', 'k', ...
%         'fontsize', 12, 'horizontalalignment', 'center', ...
%         'verticalalignment', 'middle');
% end
axis equal tight xy
hold off
colorbar
xlimits=round(get(gca,'xlim'));
set(gca,'xtick',xlimits);
ylimits=round(get(gca,'ylim'));
set(gca,'ytick',ylimits);



%% STEP 11: Make density plots for different regions

% sigma = maxVal / 40;
% numPoints = 501;
% rangeVals = [-maxVal maxVal];
% 
% [xx,density] = findPointDensity(embeddingValues,sigma,numPoints,rangeVals);
maxDensity = max(density(:));

%reload the tracks
[allTracks, folder_indecies, track_indecies] = loadtracks(folders);

my_colormap = othercolor('OrRd9');
my_colormap(1,:) = [1 1 1];

for track_index = 1:length(allTracks);
    plot_embedding = allTracks(track_index).Embeddings;
    image_file = fullfile([folders{folder_indecies(track_index)},filesep,'individual_worm_imgs',filesep,'worm_',num2str(track_indecies(track_index)),'.mat']);
    save_file = fullfile([folders{folder_indecies(track_index)},filesep,'individual_worm_imgs',filesep,'behaviormap_', num2str(track_indecies(track_index))]);
    load(image_file);

    behavior_figure = figure('Position', [500, 500, 500, 250]);
    outputVideo = VideoWriter(save_file,'MPEG-4');
    outputVideo.FrameRate = 14;
    open(outputVideo)
    for worm_frame_index = 1:size(plot_embedding, 1)
        subplot_tight(1,2,2,0);

        plot_worm_frame(worm_images(:,:,worm_frame_index), squeeze(allTracks(track_index).Centerlines(:,:,worm_frame_index)), ...
        [], ...
        [], [], ...
        [],  [], 0);

        freezeColors

        subplot_tight(1,2,1,0);
        hold on
        imagesc(xx,xx,density)
        axis equal tight off xy
%         caxis([0 maxDensity * .8])
%         colormap(jet)
        caxis([0 maxDensity * .6])
        colormap(my_colormap)
        plot(xx(jj),xx(ii),'k.') %watershed borders
        plot(plot_embedding(worm_frame_index,1), plot_embedding(worm_frame_index,2), 'ob', 'MarkerSize', 15, 'LineWidth', 3)
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
