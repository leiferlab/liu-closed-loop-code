%%% This matlab code will help the user to re-generate figure S2 as well as
%%% the numbers presented in supplementary table S1. In this table, you can 
%%% get the numbers presented in columns a) Cumulative Recording Length (animal-hours)
%%% and b) Animals per frame, for AML67 and AML470 open loop and closed-loop experiments

clear
clc
close all

%%% users can manually select folders or choose specific folders from below

%%% manual selection of folders
% % select folders 
% % folders = getfoldersGUI();

%%% to plot supplementary figure S2 uncomment the following line of code
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/reviewer_comments/tracking_frame_rate/folder_name_list_for_fig_S2/all_data/folders_list_all_data.mat')

%%% to find the animals per frame and cumulative recording length for AML67 open loop uncomment the following line of code
% load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/reviewer_comments/tracking_frame_rate/folder_name_list_for_fig_S2/AML67_open_loop/folder_list_AML67_open_loop.mat')

%%% to find the animals per frame and cumulative recording length for AML67 closed-loop uncomment the following line of code
% load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/reviewer_comments/tracking_frame_rate/folder_name_list_for_fig_S2/AML67_closed_loop/folder_list_AML67_closed_loop.mat')

%%% to find the animals per frame and cumulative recording length for AML470 open loop uncomment the following line of code
% load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/reviewer_comments/tracking_frame_rate/folder_name_list_for_fig_S2/AML470_open_loop/folder_list_AML470_open_loop.mat')

%%% to find the animals per frame and cumulative recording length for AML470 closed-loop uncomment the following line of code
% load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/reviewer_comments/tracking_frame_rate/folder_name_list_for_fig_S2/AML470_closed_loop/folder_list_AML470_closed_loop.mat')

%%% loading parameters
parameters = load_parameters(folders{1}); 
load('reference_embedding.mat')
relevant_track_fields = {'BehavioralTransition','Frames','AlignedStimulus','Velocity','EllipseRatio','Length'};

%%%%% section 1: To calculate the cumulative recording length
allTracks = [];
duration_each_folder=[];
for folder_index = 1:length(folders)
    
    sprintf('Analyzing folder %d out of %d folders', folder_index, length(folders))
    %%% load the tracks for this folder
    [current_tracks, ~, ~] = loadtracks(folders{folder_index},relevant_track_fields);
    
    allTracks = [allTracks, current_tracks];
    dummy_duration_each_folder=numel(vertcat(current_tracks.Frames))/30/60/60;
    duration_each_folder=[duration_each_folder; dummy_duration_each_folder];
end

%%% print the total duration of single animal recordings in hrs and stimulus events
sprintf(['total duration = ', num2str(numel(vertcat(allTracks.Frames))/parameters.SampleRate/3600), ' hours'])

%% Section 2: To get closed loop lags
lags = [];
experiment_drifts = [];
for folder_index = 1:length(folders)
    sprintf('Analyzing folder %d out of %d folders', folder_index, length(folders))
    folder_name = folders{folder_index};
    loaded_variable = load([folder_name, filesep, 'timestamps.mat']);
    processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;
    camera_frames = 1:numel(processed_decoded_camera_frames);
    lags = [lags, camera_frames - double(processed_decoded_camera_frames)];
    experiment_drifts = [experiment_drifts, loaded_variable.experimentmaxdrift];
end
bin_edges = 0:10;
lag_prob = histcounts(lags,bin_edges) / numel(lags);

%% plot figure S2a
ax1=figure('Renderer', 'painters', 'Position', [440 290 694 520]);
plot(bin_edges(1:end-1), lag_prob , '-o','LineWidth',4,'MarkerFaceColor',[0 0.45,0.74],'MarkerSize',15)
for label_index = 1:numel(bin_edges)-2
    text(bin_edges(1+label_index),lag_prob(1+label_index), [num2str(round(lag_prob(1+label_index)*100,3)), '%'],'VerticalAlignment','bottom','HorizontalAlignment','left')
end
% set(gca,'YScale','log')
xlabel('Round-trip Latency (frames)')
ylabel('Probability')

ax1=gca;
set(ax1, 'XTick', bin_edges(1:end-1))

ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
XOppTickLabels = round(bin_edges(1:end-1) * (1000/parameters.SampleRate));
% Set the x-tick and y-tick  labels for the second axes
set(ax2, 'XTickLabel', XOppTickLabels);
set(ax2, 'YTickLabel', []);
% ax2.FontSize = 14;
% ax1.FontSize=14;
% 
% xlim([0 9])
xlabel(ax2,'Round-trip Latency (ms)')
%%
%% plot the drifts (figure S2c)
figure3=figure;
axes3 = axes('Parent',figure3);
bin_edges = 0:10:200;
experiment_drift_prob = histcounts(experiment_drifts*1000,bin_edges);
bar(bin_edges(1:end-1), experiment_drift_prob , 'hist')
yticks([0:20:40])
set(axes3,'YTickLabel',...
    {'0','20','40'});
% for label_index = 1:numel(bin_edges)-2
%     text(bin_edges(1+label_index),experiment_drift_prob(1+label_index), [num2str(round(experiment_drift_prob(1+label_index)*100,3)), '%'],'VerticalAlignment','top','HorizontalAlignment','left')
% end
axis([0 200 0 40])
xlabel('Spatial Drift (um)')
ylabel('Experiment Plate Count')
ax = gca;
ax.FontSize = 14;
box off;

%% Section 3: To get tracking frame rates
MAX_LAG = 15;
track_count_vs_tracking_lag = zeros(255,MAX_LAG); %max 255 tracks, 100 frames lag

for folder_index = 1:length(folders)
    
    sprintf('Analyzing folder %d out of %d folders', folder_index, length(folders))
    
    folder_name = folders{folder_index};
    loaded_variable = struct2cell(load([folder_name, filesep, 'labview_tracks.mat']));
    Tracks = [loaded_variable{:}];
    if isempty(Tracks)
        return
    end
    TrackFramesDiff = [];
    TrackFrames = [];
    [~,sort_idx] = sort([Tracks.WormIndex]);
    Tracks = Tracks(sort_idx);
    tracks_to_throw_out = []; %throw out tracks with less than 3 datapoints
    for track_index = 1:length(Tracks)
        if length(Tracks(track_index).Frames) < 3
            tracks_to_throw_out = [tracks_to_throw_out, track_index];
        end
    end
    Tracks(tracks_to_throw_out) = [];
    for track_index = 1:numel(Tracks)
        TrackFramesDiff = [TrackFramesDiff, diff(Tracks(track_index).Frames)'];
        TrackFrames = [TrackFrames, (Tracks(track_index).Frames(2:end))'];
    end
    TrackFramesDiff = double(TrackFramesDiff);
    TrackFramesDiff(TrackFramesDiff<0) = 1;
    
    % plot # of worms vs tracking frame rate
    [sorted_all_track_frames,sorted_indecies] = sort(TrackFrames);
    sorted_all_track_frames_diff = TrackFramesDiff(sorted_indecies);
    
    current_frame = sorted_all_track_frames(1);
    starting_all_frame_index = 1;
    for all_frame_index = 1:numel(sorted_all_track_frames)
        if sorted_all_track_frames(all_frame_index) > current_frame
            track_count_vs_tracking_lag(all_frame_index-starting_all_frame_index, min(MAX_LAG,sorted_all_track_frames_diff(starting_all_frame_index))) = ...
                track_count_vs_tracking_lag(all_frame_index-starting_all_frame_index, min(MAX_LAG,sorted_all_track_frames_diff(starting_all_frame_index))) + 1;
            current_frame = sorted_all_track_frames(all_frame_index);
            starting_all_frame_index = all_frame_index;
        else
            continue
        end
    end
end

bin_edges = 0:MAX_LAG;
lag_counts = sum(track_count_vs_tracking_lag,1);
lag_prob = lag_counts / sum(lag_counts);

%% plot figure S2b
figure('Renderer', 'painters', 'Position', [440 290 694 520])
plot(bin_edges(2:end), lag_prob , '-o','Color',[0 0.45,0.74],'MarkerFaceColor',[0 0.45,0.74],'MarkerSize',15,'LineWidth',4)
% label the first label_n datapoints
% label_n = 4;
label_n = numel(lag_prob);
for label_index = 1:label_n
    text(bin_edges(label_index),lag_prob(label_index), [num2str(round(lag_prob(label_index)*100,2)), '%'],'VerticalAlignment','top','HorizontalAlignment','left')
end
set(gca, 'XTick', bin_edges(1:end-1))
set(gca,'YScale','log')
xlabel('Tracking Update Time (frames)')
ylabel('Log Probability')

ax1=gca;
ax1.FontSize = 16;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
XOppTickLabels = round(parameters.SampleRate./bin_edges(1:end-1),2);
% Set the x-tick and y-tick  labels for the second axes
set(ax2, 'XTickLabel', XOppTickLabels);
set(ax2, 'YTickLabel', []);
xlabel(ax2,'Tracking Frame Rate (Hz)')

%% Section 4: plot the probability of worm count at any given time
worm_counts = sum(track_count_vs_tracking_lag,2);
worm_counts_prob = worm_counts ./ sum(worm_counts);

worm_counts_enumerated = [];
for count_index = 1:numel(worm_counts)
    worm_counts_enumerated = [worm_counts_enumerated, ones(1, worm_counts(count_index)).* count_index];
end

mean_worm_counts = mean(worm_counts_enumerated);
std_worm_counts = std(worm_counts_enumerated);

sprintf(['Animals per frame: Mean = ', num2str(mean_worm_counts), ' stdev = ', num2str(std_worm_counts)])

%% plot figure S2f
figure
bar(worm_counts_prob , 'hist','EdgeColor',[1 0 0])
axis([0 100 0 0.1])
xlabel('Tracked Worms in Frame')
ylabel('Probability')
title(['mean : ', num2str(round(mean_worm_counts)), '    ', 'std: ', num2str(round(std_worm_counts))])
ylim([0 0.04])
% xlim([0 100])
yticks([0:0.02:0.04])
ax = gca;
ax.FontSize = 14;
box off;

%% Section 5: Prep for track duration and worm lengths distribution
relevant_track_fields = {'Path','Frames','Length','BehavioralTransition'};        

all_track_durations = [];
all_track_lengths = [];
for folder_index = 1:numel(folders)
    sprintf('Analyzing folder %d out of %d folders', folder_index, length(folders))
    folder_name = folders{folder_index};
    Tracks = load_single_folder(folder_name, relevant_track_fields);
    parameters = load_parameters(folder_name);
    track_durations = zeros(1,numel(Tracks));
    
    track_lengths = zeros(1,numel(Tracks));
    for track_index = 1:numel(Tracks)
        track_durations(track_index) = numel(Tracks(track_index).Frames);
        track_lengths(track_index) = mean(Tracks(track_index).Length);
    end

    all_track_durations = [all_track_durations, track_durations./parameters.SampleRate];
    all_track_lengths = [all_track_lengths, track_lengths./parameters.CameraPixeltommConversion*1000];
end

%% plot track durations (plot figure S2e)
figure5=figure;
axes5 = axes('Parent',figure5);
hold(axes5,'on');
[track_duration_prob, bin_edges,~] = histcounts(all_track_durations, 100,'Normalization','probability');
bin_centers = (bin_edges(2:end)+bin_edges(1:end-1))./2;
bar(bin_centers, track_duration_prob, 'hist')
xticks([0:600:1800])
set(axes5,'XTickLabel',...
    {'0','10','20','30'});
ylim([0 0.18])
xlim([0 1800])
yticks([0:0.06:0.18])
xlabel('Track Durations (minutes)')
ylabel('Probability')
ax = gca;
ax.FontSize = 14;
box off;

%% plot worm lengths (plot figure S2d)
figure4=figure;
axes4 = axes('Parent',figure4);
hold(axes4,'on');
[track_lengths_prob, bin_edges,~] = histcounts(all_track_lengths, 'Normalization','probability');
bin_centers = (bin_edges(2:end)+bin_edges(1:end-1))./2;
% bin_centers=bin_centers(1:2:end);
% track_lengths_prob=track_lengths_prob(1:2:end);
hold on
b=bar(bin_centers, track_lengths_prob, 'hist');
% b.CData(1,50:end)=[0 1 0];
h=line([550 550], [0 0.05]);
h.Color= [1 0 0];
h.LineWidth= 2;
axis([400 1200 0 0.05]) 
yticks([0:0.025:0.05])
% set(gca, 'XTick', bin_edges)
xlabel('Worm Lengths (um)')
ylabel('Probability')
ax = gca;
ax.FontSize = 14;
box off;

disp('Analysis finished')