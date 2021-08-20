%%%% This code is used in the analysis for generating data for table 1 and
%%%% figure 4 for AML67 and AML470. The number of instances in which reversals
%%%% are generated for closed loop and open loop assays are saved in th evariables
%%%% count_number_of_reversal_folder_turns and
%%%% count_number_of_reversal_folder_rails respectively.

close all
clear
clc

main_folder=('/projects/LEIFER/Sandeep/APIData/20200902_RunRailsTriggeredByTurning_Sandeep_AML67_10ulRet');
cd(main_folder)

%%%%%%% user input parameters %%%%%%
test_stimulus_duration=5; %%% if you want to test a specific stim duration
stim_threshold_min=30; %%%% min stim power 
stim_threshold_max=50;   %%%% max stim power

%%%%%%% determine stim color %%%%%%%%
if contains(main_folder, 'blue')
    stim_color=[0 0 1];
elseif contains(main_folder,'purple')
    stim_color=[1 0 1];
elseif contains(main_folder,'red')
    stim_color=[1 0 0]; 
else
    stim_color=[0 0 1];
end

%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%
behavior_to_study=2;  %%%  2= turn, 3=reverse 
two_stim_delivered=11; %%% when both blue and red stim are delivered simultaneously
opto_stim_color=1; %%% red=1, green=2, blue=3

max_allowed_stim_duration_diff=10; %%% this parameter is used for sorting stim of diff widths   
minimum_length_of_worm=0.55;  %%%% 
ellipse_ratio_threshold=3.6;
negative_vel_threshold=-0.1;
frames_in_reversal_threshold=1; %%%% this is the time for which worms must reverse
test_low_ellipse_ratio_worm_threshold=3.9; %%% if the mean ellipse ratio of certian frames below this threshold then that trial is ignored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if contains(main_folder, 'Full')
    testtype='rails';
    analysis_rails=1;
    analysis_stim_while_turning=11;
elseif contains(main_folder,'Turn')
    testtype='turns';
    analysis_rails=11;
    analysis_stim_while_turning=1;
elseif contains(main_folder,'Reversing')
    testtype='turns';
    analysis_rails=11;
    analysis_stim_while_turning=1;
else
    analysis_rails=1;
    analysis_stim_while_turning=11;  
end

%%%% This if loop will pick the correct behavior to select whether it is reversing or turning.
if contains(main_folder, 'Reversing')
    correct_behavior_threshold=negative_vel_threshold;
    analysis_stim_while_turning=1;
    behavior_to_study=3;
elseif contains(main_folder,'Turn')
    correct_behavior_threshold=ellipse_ratio_threshold;
    analysis_stim_while_turning=1;
    behavior_to_study=2;
end

if contains(main_folder,'TwoColors')
    two_stim_delivered=1;
end

%%%% ask for folders
load('reference_embedding.mat')
relevant_track_fields = {'BehavioralTransition','Frames','AlignedStimulus','Velocity','EllipseRatio','Length'};

%%%% select folders
folders = getfoldersGUI();
parameters = load_parameters(folders{1}); % load parameters for the first folder

%%%%% constants
normalized_stimuli = 1; %delta function
fps = parameters.SampleRate;
stim_similarity_thresh = 1.1;  % combine data in case the intensity is within some range (to take care of box to box variation)
number_of_behaviors = max(L(:)-1);
rails_durations=parameters.RailsDurations;
rails_durations=sort(rails_durations);
rails_durations=unique(rails_durations);  %%% in case, in th elabview program same time was entered twice

index_durations=find(rails_durations==test_stimulus_duration);
for rails_dur=index_durations   %%% to debug or analyze particular stim duration

rails_dur_itr=rails_durations(1,rails_dur);
current_rails_dur=rails_dur_itr*parameters.SampleRate;

vel_xlim_turns_1=501; %%% defining the window in which I will look for negative vel for stim on turns
vel_xlim_turns_2=vel_xlim_turns_1+current_rails_dur; %%% defining the window in which I will look for negative vel for stim on turns
vel_xlim_rails_1=501; %%% defining the window in which I will look for negative vel for stim on rails
vel_xlim_rails_2=vel_xlim_rails_1+current_rails_dur; %%% defining the window in which I will look for negative vel for stim on rails

time_window_before =10*parameters.SampleRate;
time_window_after = rails_dur_itr*parameters.SampleRate;
total_window_frames = time_window_before+time_window_after+1;

stimulus_intensities = [];

all_behavior_transitions_for_frame = {};
all_behavior_annotations_for_frame = {};

%%%%% loop through each folder
velocity_at_stim_while_turning_folder=[];
ellipse_ratio_at_stim_while_turning_folder=[];
count_number_of_reversal_folder_turns=[];
count_number_of_reversal_folder_rails=[];
beh_state_on_stim_folder_new_protocol=[];  %%%% to find the turns on stim onset
ellipse_ratio_during_turn_state_folder=[];
just_the_ellipse_ratio_worm_data_folder=[];
all_ellipse_ratio_folder=[];
cumulative_recording_duration_folder=0; %% cumulative reconding length
trial_count_folder=0;       %% to count the number of worm tracks in an assay
worms_tracked_simultaneously_folders=[];   %% to count the number of worms tracked simultaneously
error_code_folder=[]; %% to determine the criteria due to which a worm is not considered for analysis

for folder_index = [1:length(folders)]

    %%%% let us load the parameter file for each folder to determine the worm length threshold
    parameters_length = load_parameters(folders{folder_index});
    length_of_worm_threshold=minimum_length_of_worm*parameters_length.CameraPixeltommConversion;
    
    sprintf('Analyzing data %d out of %d datasets for %d sec duration', folder_index, length(folders),current_rails_dur/parameters.SampleRate)
    %%%%load the tracks for this folder
    [current_tracks, ~, ~] = loadtracks(folders{folder_index},relevant_track_fields);
     
    %%%%% generate the Behavior matricies
    current_tracks = BehavioralTransitionToBehavioralAnnotation(current_tracks);
    current_tracks = get_behavior_triggers(current_tracks);
    
    if size(current_tracks,2)==0  %%%% skip the iteration if it is an empty worm track
        continue
    end
        
    all_behavior_data{:,folder_index} = horzcat(current_tracks.BehavioralAnnotation);
    
    current_param = load_parameters(folders{folder_index});
    
    stim_while_turning_peaks_array=[];
    fullrails_peaks_array=[];
    
    velocity_at_stim_while_turning_trial=[];
    ellipse_ratio_at_stim_while_turning_trial=[];
    count_number_of_reversal_trial_turns=[];
    
    velocity_at_stim_while_turning_trial_rails=[];
    ellipse_ratio_at_stim_on_rails_trial=[];
    count_number_of_reversal_trial_rails=[];
    just_the_ellipse_ratio_worm_data_trial=[];

    beh_state_on_stim_trial_new_protocol=[];  %%%% to find the turns on stim onset
    ellipse_ratio_during_turn_state_trial=[];
    all_ellipse_ratio_trial=[];
    cumulative_recording_duration_trial=0;
    trial_count_trial=0;
    frame_details_trial=[];
    error_code_trial=[];
    
    for track_index = 1:length(current_tracks)
        
        %%%% calculating cumulative recording duration and num of trials for each trial
        cumulative_recording_duration_trial=cumulative_recording_duration_trial+size(current_tracks(track_index).Frames,1);
        trial_count_trial=size(current_tracks,2);
        
        %%%% cancatenating Frame informaion to find out the number of worms tracked simultaneously at a time.
        
        frame_details_trial=[frame_details_trial
            current_tracks(track_index).Frames];
                
        %%%%% if no stimulus was delivered then ignore the iteration
        if nnz(current_tracks(track_index).AlignedStimulus(:,10))==0
            continue
        end
        
        if two_stim_delivered==1
           current_tracks(track_index).AlignedStimulus(:,:,1)=current_tracks(track_index).AlignedStimulus(:,:,opto_stim_color); 
        end
            
        mid_cline_index=size(current_tracks(track_index).AlignedStimulus,2)/2;
        
        %%%% if the stim is blue then multiply stim by -1closestIndex
        if min(min(current_tracks(track_index).AlignedStimulus))<0
            current_tracks(track_index).AlignedStimulus(:,:)=-1*current_tracks(track_index).AlignedStimulus(:,:); 
        end
                
        %%%%%%%%%%%%%%%%% to ignore smaller length worms %%%%%%%%%%%%%%%%%%
        length_of_current_worm=mean(current_tracks(track_index).Length);
        length_of_worm_matrix{folder_index,track_index}=length_of_current_worm;
        
        if length_of_current_worm<length_of_worm_threshold   %%% ignoring the lower 10 percentile worm length
            dummy_error_code(:,1)=folder_index;
            dummy_error_code(:,2)=track_index;
            dummy_error_code(:,3)=1;
            dummy_error_code(:,4)=1;
            error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code for short length =1
            continue
        end

        %%%%%%%%%%%%%%% determining missing behavior indices %%%%%%%%%%%%%%
        
        original_velocity_for_correction=current_tracks(track_index).Velocity;
        smooth_velocity_for_correction=smooth(original_velocity_for_correction,10);   %%%% smoothing the velocity data
        current_tracks(track_index).Velocity=smooth_velocity_for_correction;  %%% redefinig the velocity info with smootherd data
        smooth_velocity_for_correction=smooth_velocity_for_correction';
        original_velocity_matrix_for_correction{track_index,:}=smooth_velocity_for_correction;
        
        original_ellipse_ratio_for_correction=current_tracks(track_index).EllipseRatio;
        smooth_ellipse_ratio_for_correction=smooth(original_ellipse_ratio_for_correction,10)'; %% smotting ellipse ratio
        current_tracks(track_index).EllipseRatio=smooth_ellipse_ratio_for_correction'; %% redefining the ellipse ratio with smoothed data
        smooth_ellipse_ratio_for_correction_matrix{track_index,:}=smooth_ellipse_ratio_for_correction;
        
        %%%%%% Finding the ellipse ratio when the worm is in the turn state
        index_during_turn_state=find(current_tracks(track_index).BehavioralAnnotation==2);
        ellipse_ratio_during_turn_state=current_tracks(track_index).EllipseRatio([index_during_turn_state],1);
        ellipse_ratio_during_turn_state_trial=[ellipse_ratio_during_turn_state_trial
            ellipse_ratio_during_turn_state];
        
        %%%%% cancatenating all the ellipse ratio of the worm
        all_ellipse_ratio_track=current_tracks(track_index).EllipseRatio;
        all_ellipse_ratio_trial=[all_ellipse_ratio_trial 
            all_ellipse_ratio_track];

        %%%% if pipeline does'nt assign behavior to one of the three states then we will take velocity into consideration to assign forward or reverse
        %%%% assigning beh as turn when the ellipse ratio was below the ellipse ration threshold
        for xyz=1:size(current_tracks(track_index).EllipseRatio,1)
            if smooth_ellipse_ratio_for_correction(1,xyz)<ellipse_ratio_threshold
                current_tracks(track_index).BehavioralAnnotation(1,xyz)=2;
            end
        end
        
        original_data_pre_correction=current_tracks(track_index).BehavioralAnnotation;
        original_data_matrix_pre_correction{track_index,:}=original_data_pre_correction;
        
        data_unassigned_pre_correction=find(original_data_pre_correction==0); %%%% when beh data is unassigned
        
        corrected_data_for_missing_behaviors=original_data_pre_correction;
        for xx=[data_unassigned_pre_correction]
            if smooth_velocity_for_correction(1,xx)<=-0.06
                corrected_data_for_missing_behaviors(1,xx)=3;   %%% assigning reversal
            else
                corrected_data_for_missing_behaviors(1,xx)=1;   %%% assigning forward
            end
        end
        
        current_tracks(track_index).BehavioralAnnotation=corrected_data_for_missing_behaviors;  %%%% we are reassigning behavior index for missing frames
        corrected_data_matrix_for_missing_behaviors{track_index,:}=corrected_data_for_missing_behaviors;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        original_stim=current_tracks(track_index).AlignedStimulus(:,mid_cline_index);
        [raw_stim_peaks_amplitude, raw_stim_peaks_locs,raw_stim_widths,~] = findpeaks(current_tracks(track_index).AlignedStimulus(:,mid_cline_index), 'MinPeakDistance',current_param.InterTriggerInterval);
        
        %%%%BEGIN ANDYS REWRITE
        
        %%%start with the original stimulus time series and binarize the
        %%%stimulus into on and off depending on whether it exceeds a stimulu
        %%%sthreshold; then take the derivatives and find instances above zero 
        %%%and add 1 to account for teh derivative time shift effect
        stim_onsets = find(diff([original_stim>stim_threshold_min & original_stim<stim_threshold_max]) >0)+1;
        
        %%%Now we will apply increasingly selective criteria to only end up 
        %%%with stim_onset events that occur exactly during the onset of a turn
        
        %%%require that thes occur during times that mochi's pipeline classifies as turn
        data_turn=find(current_tracks(track_index).BehavioralAnnotation==behavior_to_study); %% when worm is turning
        stims_delivered_during_mochis_turns = intersect(stim_onsets, data_turn);
        
        %%%get the binary timeseries for when ellipse ratio dips below thresohold then take diff        
        er_neg_crossing_timeseries = [0; diff((current_tracks(track_index).EllipseRatio<ellipse_ratio_threshold))>0];
        
        %%%%then find negative values (which are the negative crossing events)
        er_turn_onset_indices = find(er_neg_crossing_timeseries>0);
        
        %%%Then get a list of all indices within +-10 frames of the crossing
        %%%event
        er_turn_onset_zones_indices = find(movmax((er_neg_crossing_timeseries>0),round(parameters.SampleRate/1.5))==1);
        
        %%%then filter valid stims further and impose that they fall in this
        %%%temporal window tied to the ellipsoid ratio crossing event
        valid_stim = intersect(stims_delivered_during_mochis_turns,er_turn_onset_zones_indices);
        valid_stim_matrix{track_index,:}=valid_stim;
        
        final_valid_stim=find_stim_on_turn_onsets(stims_delivered_during_mochis_turns,er_turn_onset_indices);
        final_valid_stim_matrix{track_index,:}=final_valid_stim;
    
        for mn=1:size(stim_onsets,1)
            
            beh_state_on_stim_track_new_protocol(:,1)=folder_index;
            beh_state_on_stim_track_new_protocol(:,2)=track_index;
            beh_state_on_stim_track_new_protocol(:,3)=stim_onsets(mn);
            
            if ismember(stim_onsets(mn),final_valid_stim)
                beh_state_on_stim_track_new_protocol(:,4)=1;  %% these are the valid stims on turn onset
            else
                beh_state_on_stim_track_new_protocol(:,4)=0;  %% these are the in-valid stims. They were not delivered on turn onset
            end
            
            beh_state_on_stim_trial_new_protocol=[beh_state_on_stim_trial_new_protocol
                beh_state_on_stim_track_new_protocol];    
            
        end

        %%% END ANDYS REWRITE
        
        raw_stim_peaks_amplitude=raw_stim_peaks_amplitude';
        raw_stim_peaks_locs=raw_stim_peaks_locs';
        raw_stim_widths=raw_stim_widths';
        
        %%%% determining peaks of the same widths                      
        raw_stim_peaks=[];
        raw_stim_peaks_height=[];
        raw_stim_duration=[];
        for i=1:size(raw_stim_widths,2)
            if (abs(current_rails_dur - raw_stim_widths(1,i)))<max_allowed_stim_duration_diff
            raw_stim_peaks=[raw_stim_peaks raw_stim_peaks_locs(1,i)];  %%% these are all the peaks with same width
            raw_stim_peaks_height=[raw_stim_peaks_height raw_stim_peaks_amplitude(1,i)];
            raw_stim_duration=[raw_stim_duration raw_stim_widths(1,i)];
            end
        end
        
        %%%%% consider stim of particular duration for this iteration
        
        new_stim_data=zeros(size(original_stim));
        for i=1:size(raw_stim_peaks,2)
            new_stim_data([raw_stim_peaks(1,i):(raw_stim_peaks(1,i)+raw_stim_duration(1,i))],1)=raw_stim_peaks_height(1,i);
        end
        
        new_stim_data(new_stim_data<stim_threshold_min)=0;
        new_stim_data(new_stim_data>stim_threshold_max)=0;
        
        if analysis_stim_while_turning==1  
        %%%% determining when the stim in on/off turn
        data_stim=find(new_stim_data>stim_threshold_min & new_stim_data<stim_threshold_max); %%%% when worm is stimulated

        raw_stim_final=zeros(1,size(current_tracks(track_index).AlignedStimulus,1));
        
        raw_data_stim_peaks=[];
        raw_stim_final(data_stim)=max(current_tracks(track_index).AlignedStimulus(:,mid_cline_index));
        
        [~, raw_data_stim_peaks] = findpeaks(raw_stim_final, 'MinPeakDistance',current_param.InterTriggerInterval); 
        
        all_data_stim_peaks_cell{track_index,1}=raw_data_stim_peaks;  %% we will capture all the raw stim peaks
        num_all_data_stim_peaks_array(track_index,1)=size(raw_data_stim_peaks,2);

        stim_on_turn=intersect(data_stim,data_turn); %%% correct stim based on turn 

        stim_on_turn_final=zeros(1,size(current_tracks(track_index).AlignedStimulus,1));
        stim_on_turn_final(stim_on_turn)=max(current_tracks(track_index).AlignedStimulus(:,mid_cline_index)); %%% this is the correct stimulation array
        [~, stim_on_turn_peaks] = findpeaks(stim_on_turn_final, 'MinPeakHeight',0.1, 'MinPeakDistance',current_param.InterTriggerInterval); 
        
        [minValue, closestIndex] = min(abs(stim_on_turn_peaks - raw_data_stim_peaks.'));
        stim_while_turning_peaks_dummy=unique(raw_data_stim_peaks(closestIndex));   %%%% unique is making sure that peaks do not get counted multiple times
        
        if behavior_to_study==2 || behavior_to_study==3
        stim_while_turning_peaks_old=[];
        
        for yz=1:length(stim_while_turning_peaks_dummy)
            if current_tracks(track_index).BehavioralAnnotation(1,stim_while_turning_peaks_dummy(yz))==behavior_to_study
                stim_while_turning_peaks_old=[stim_while_turning_peaks_old stim_while_turning_peaks_dummy(yz)];
            end
        end
        end
        
        stim_while_not_turning_peaks=setdiff(raw_data_stim_peaks,stim_while_turning_peaks_old);
        
        original_stim=current_tracks(track_index).AlignedStimulus(:,mid_cline_index);
        stim_while_turning_final=zeros(size(original_stim));
        
        stim_while_turning_peaks=intersect(final_valid_stim,raw_stim_peaks); %%% this array has the correctly defined stim associated turn peaks and peaks of same width 
 
        for i=1:length(stim_while_turning_peaks)
            
            sprintf('current_track_index: %d and stim_at_peaks: %d',track_index,stim_while_turning_peaks(i))
            
            %%% to avoid error due to 'findpeaks' function when the array begins with peak
            if any(stim_while_turning_peaks(i)<10)
                disp('Error 10: Stim on turn peaks within 10 frames')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=[];
                dummy_error_code(:,4)=9; %% error code is 9 when worm track begins with stim
                error_code_trial=[error_code_trial
                    dummy_error_code];  
                continue
            end
        
            if stim_while_turning_peaks(i)<=500
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=2;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% if the worm is not tracked for 17 sec pre-stimulus onset error code =2
                continue
            end
            
            if size(current_tracks(track_index).AlignedStimulus(:,10),1)<=(stim_while_turning_peaks(i)+current_rails_dur)
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=3;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code=3 if the worm is not tracked for the entirety stimulus duration
                continue
            end
                
            dummy_array_with_velocity_info_turns=[];
            dummy_array_with_ellipse_ratio_info_turns=[];
            dummy_array_with_ellipse_ratio_info_turns=[];
            locs_er=[];
            time_spend_in_turns=[];
            
            index_turn_peak=find(stim_while_turning_peaks(i)==raw_stim_peaks_locs);
            stim_while_turning_final([raw_stim_peaks_locs(1,index_turn_peak):(raw_stim_peaks_locs(1,index_turn_peak)+raw_stim_widths(1,index_turn_peak))],1)=raw_stim_peaks_amplitude(1,index_turn_peak);
            
            stim_while_turning_peaks_array(track_index,i)=stim_while_turning_peaks(i);
            
            dummy_array_with_velocity_info_turns=current_tracks(track_index).Velocity((stim_while_turning_peaks(i)-500):(stim_while_turning_peaks(i)+current_rails_dur),1);
            
            velocity_at_stim_while_turning_trial=[velocity_at_stim_while_turning_trial
                dummy_array_with_velocity_info_turns];
            
            if max(dummy_array_with_velocity_info_turns(1:400,:))<=0.01 && min(dummy_array_with_velocity_info_turns(1:400,:))>=-0.01  %%%% ignore trials when worms hardly moved
                disp('Worm is not moving: Trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=4;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code for not moving worm =4
                continue
            end
            
            dummy_array_with_ellipse_ratio_info_turns=current_tracks(track_index).EllipseRatio((stim_while_turning_peaks(i)-500):(stim_while_turning_peaks(i)+current_rails_dur),1);
       
            just_the_ellipse_ratio_worm_data_trial=[just_the_ellipse_ratio_worm_data_trial
                dummy_array_with_ellipse_ratio_info_turns];
            
            %%%%%%%%%% to get rid of trials when two worms collide
            if any(dummy_array_with_ellipse_ratio_info_turns>8)
                disp('worms collided during trial')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=5;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when worm collide =5
                continue
            end
            
            %%%%% to get rid of trials with very low ellipse ratio in turns analysis 
            randIdcs_turns = randi([1 400],1,50);
            test_low_ellipse_ratio_worm = dummy_array_with_ellipse_ratio_info_turns(randIdcs_turns,:);

            if mean(test_low_ellipse_ratio_worm)<test_low_ellipse_ratio_worm_threshold
                disp('Worm with low mean ellipse ratio')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=6;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when worm's mean ellipse ratio is low =6
                continue
            end
        
            ellipse_ratio_at_stim_while_turning_trial=[ellipse_ratio_at_stim_while_turning_trial
                dummy_array_with_ellipse_ratio_info_turns];
            
            %%% make sure that stim is delivered during turns and not reversals 
            if any(dummy_array_with_velocity_info_turns(495:500,:)<(negative_vel_threshold/2))
                disp('Stim while reversal: trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=7;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when stim delivered to reversing worms =7
                continue
            end
                                            
            %%%%%% detecting reversals in the stimulus window 
            indices_below_set_vel_threshold_turns = find(dummy_array_with_velocity_info_turns(vel_xlim_turns_1:vel_xlim_turns_2,:) < negative_vel_threshold); 

            if size(indices_below_set_vel_threshold_turns,1)>=frames_in_reversal_threshold 
                dummy_array_with_number_of_reversal_turns(:,1)=folder_index;
                dummy_array_with_number_of_reversal_turns(:,2)=track_index;
                dummy_array_with_number_of_reversal_turns(:,3)=stim_while_turning_peaks(i);
                dummy_array_with_number_of_reversal_turns(:,4)=1;
            else
                dummy_array_with_number_of_reversal_turns(:,1)=folder_index;
                dummy_array_with_number_of_reversal_turns(:,2)=track_index;
                dummy_array_with_number_of_reversal_turns(:,3)=stim_while_turning_peaks(i);
                dummy_array_with_number_of_reversal_turns(:,4)=0;
            end

            count_number_of_reversal_trial_turns=[count_number_of_reversal_trial_turns
                dummy_array_with_number_of_reversal_turns];  %%% cancatenating the info if a worm reversed or not
                               
        end
        
        if analysis_stim_while_turning==1
            single_dimension_stimulus=stim_while_turning_final';
        end
        
        end
        
        %%%%%%%%%%%%%%%%%%% analysis for open loop %%%%%%%%%%%%%%%%%%%%%%%%
        
        if analysis_rails==1
            single_dimension_stimulus=new_stim_data';
        end
        
        xcorr_ledvoltages_stimulus = abs(padded_conv(single_dimension_stimulus, normalized_stimuli));
        [~, all_stim_peaks] = findpeaks(xcorr_ledvoltages_stimulus, 'MinPeakDistance',current_param.InterTriggerInterval);

        %%%%%%%%%%%% counting number of reversals for rails %%%%%%%%%%%%%%%
        if analysis_rails==1
            stim_peaks=setdiff(all_stim_peaks,final_valid_stim); %%% get rid of the stims which occured on the onset of turns
            
            for ab=1:length(stim_peaks)
            
            sprintf('current_track_index: %d and stim_at_peaks: %d',track_index,stim_peaks(ab))
                
            if stim_peaks(ab)<=500
                disp('Stim in the first 500 frames: trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=2;   %% worm not tracked for 17 seconds before stimulus onset
                error_code_trial=[error_code_trial 
                dummy_error_code];   
                continue
            end
            
            if size(current_tracks(track_index).AlignedStimulus(:,10),1)<=(stim_peaks(ab)+current_rails_dur)
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=3; %% worm not tracked for the complete stimulus duration
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end
               
            dummy_array_with_velocity_info_rails=[];
            
            stim_on_rails_peaks_array(track_index,ab)=stim_peaks(ab);
            
            dummy_array_with_velocity_info_rails=current_tracks(track_index).Velocity((stim_peaks(ab)-500):(stim_peaks(ab)+current_rails_dur),1);
            
            if max(dummy_array_with_velocity_info_rails(1:400,:))<=0.01 && min(dummy_array_with_velocity_info_rails(1:400,:))>=-0.01  %%%% ignore trials when worms hardly moved
                disp('Worm is not moving: Trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=4;  %% worm is not moving (stationary most of the time)
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end
            
            velocity_at_stim_while_turning_trial_rails=[velocity_at_stim_while_turning_trial_rails
                dummy_array_with_velocity_info_rails];
            
            dummy_array_with_ellipse_ratio_info_rails=current_tracks(track_index).EllipseRatio((stim_peaks(ab)-500):(stim_peaks(ab)+current_rails_dur),1);
            
            %%%%%%%%%% to get rid of trials when two worms collide
            if any(dummy_array_with_ellipse_ratio_info_rails>8)
                disp('worms collided during trial')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=5;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when worm collide =5
                continue
            end
            
            %%%%% to get rid of trials with very low ellipse ratio in turns analysis 
            randIdcs_rails = randi([1 400],1,50);
            test_low_ellipse_ratio_worm_rails = dummy_array_with_ellipse_ratio_info_rails(randIdcs_rails,:);
                     
            if mean(test_low_ellipse_ratio_worm_rails)<test_low_ellipse_ratio_worm_threshold
                disp('Worm with low mean ellipse ratio')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=6;  %% low mean ellipse ratio of worm
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end

            ellipse_ratio_at_stim_on_rails_trial=[ellipse_ratio_at_stim_on_rails_trial
                dummy_array_with_ellipse_ratio_info_rails];
            
            %%% make sure that stim is delivereed during forward motion and not reversals
            if any(dummy_array_with_velocity_info_rails(495:500,:)<negative_vel_threshold/2)
                disp('Stim while reversal: trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=7;  %% stim delivered during reversal
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end
            
            %%% make sure that stim is delivered during forward motion and not turns
            if any(dummy_array_with_ellipse_ratio_info_rails(490:500,:)<=ellipse_ratio_threshold)
                disp('Stim on turn: trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=8; %% stim delivered during turn
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end
            
            indices_below_set_vel_threshold_rails = find(dummy_array_with_velocity_info_rails(vel_xlim_rails_1:vel_xlim_rails_2,:) < negative_vel_threshold); 

            %%% when no reversal is detected
            if isempty(indices_below_set_vel_threshold_rails)

                disp('No reversal detected')
                dummy_array_with_number_of_reversal_rails(:,1)=folder_index;
                dummy_array_with_number_of_reversal_rails(:,2)=track_index;
                dummy_array_with_number_of_reversal_rails(:,3)=stim_peaks(ab);
                dummy_array_with_number_of_reversal_rails(:,4)=0;
                count_number_of_reversal_trial_rails=[count_number_of_reversal_trial_rails
                dummy_array_with_number_of_reversal_rails];
                continue
            end

            %%%% if reversal is detected

            if size(indices_below_set_vel_threshold_rails,1)>=frames_in_reversal_threshold 
                dummy_array_with_number_of_reversal_rails(:,1)=folder_index;
                dummy_array_with_number_of_reversal_rails(:,2)=track_index;
                dummy_array_with_number_of_reversal_rails(:,3)=stim_peaks(ab);
                dummy_array_with_number_of_reversal_rails(:,4)=1;
            else
                dummy_array_with_number_of_reversal_rails(:,1)=folder_index;
                dummy_array_with_number_of_reversal_rails(:,2)=track_index;
                dummy_array_with_number_of_reversal_rails(:,3)=stim_peaks(ab);
                dummy_array_with_number_of_reversal_rails(:,4)=0;
            end

            count_number_of_reversal_trial_rails=[count_number_of_reversal_trial_rails
                dummy_array_with_number_of_reversal_rails];
            
            end
        end

    end  
    
    %%%% counting the number of occurence of unique elements in frame_details_trial
    %%%% This will give me the number of worms tracked at each given frame
    [worms_tracked_simultaneously_trials,GroupR_frames_trials] = groupcounts(frame_details_trial);
    
    worms_tracked_simultaneously_folders=[worms_tracked_simultaneously_folders
        worms_tracked_simultaneously_trials];
    
    %%%% cancatenating the data from each assay (or trial) and generate for all the folders
    just_the_ellipse_ratio_worm_data_folder=[just_the_ellipse_ratio_worm_data_folder
        just_the_ellipse_ratio_worm_data_trial];
        
    velocity_at_stim_while_turning_folder=[velocity_at_stim_while_turning_folder
        velocity_at_stim_while_turning_trial];
    
    ellipse_ratio_at_stim_while_turning_folder=[ellipse_ratio_at_stim_while_turning_folder
    ellipse_ratio_at_stim_while_turning_trial];
    
    count_number_of_reversal_folder_turns=[count_number_of_reversal_folder_turns
        count_number_of_reversal_trial_turns];
    
    count_number_of_reversal_folder_rails=[count_number_of_reversal_folder_rails
    count_number_of_reversal_trial_rails];

    beh_state_on_stim_folder_new_protocol=[beh_state_on_stim_folder_new_protocol
        beh_state_on_stim_trial_new_protocol];
    
    ellipse_ratio_during_turn_state_folder=[ellipse_ratio_during_turn_state_folder
        ellipse_ratio_during_turn_state_trial];
    
    all_ellipse_ratio_folder=[all_ellipse_ratio_folder 
        all_ellipse_ratio_trial];
    
    error_code_folder=[error_code_folder
                error_code_trial];  
    
    cumulative_recording_duration_folder=cumulative_recording_duration_folder+cumulative_recording_duration_trial;
    trial_count_folder=trial_count_folder+trial_count_trial;

end

end

disp('Calculation finished')
