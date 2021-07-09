% This script performs behavioral analysis without a cluster

% analysis options
decode_camera_video = 1;
saving_individual_image = 1;
finding_centerline = 1;
resolving_problems = 1;
calculate_behavior = 1;
stimulus_realignment = 1;
plotting = 1;
parameters = load_parameters(); %load default parameters


%% STEP 1: Get folders to analyze
[folders, folder_count] = getfoldersGUI();

%% STEP 2: Decode Camera Video
if decode_camera_video
    for folder_index = 1:folder_count
        folder_name = folders{folder_index};
        if ~exist([folder_name, filesep, 'DecodedCameraFrames'], 'dir')
            mkdir([folder_name, filesep, 'DecodedCameraFrames'])
        end
        system([pwd, filesep, 'LabviewVIs', filesep, 'ffmpeg -y -i ', folder_name, filesep, 'CameraFrames.mkv -compression_algo raw -pix_fmt gray -start_number 0 ', folder_name, filesep, 'DecodedCameraFrames\Frame_%06d.png']);
    end
end

%% STEP 3: Track and save the individual worm images %%
if saving_individual_image
    'Saving individual worm images...'
    for folder_index = 1:folder_count
        folder_name = folders{folder_index}
        save_individual_images(folder_name);
    end
end

%% STEP 4: Find centerlines %%
if finding_centerline
    'Getting Centerlines...'
    for folder_index = 1:folder_count
        folder_name = folders{folder_index}
        find_centerlines(folder_name);
    end 
end

%% STEP 6: Resolve centerline and tracking problems
if resolving_problems
    'Resolve Issues'
    for folder_index = 1:folder_count
        folder_name = folders{folder_index}
        auto_resolve_problems(folder_name);
    end 
end

%% STEP 7: do behavioral mapping
if calculate_behavior
   'Getting Behaviors'
    for folder_index = 1:folder_count
        folder_name = folders{folder_index}
        calculate_spectra(folder_name);
        calculate_embeddings(folder_name);
        calculate_behaviors(folder_name);
    end
end

%% STEP 8: realign stimulus
if stimulus_realignment
   'Realigning Stimulus'
    for folder_index = 1:folder_count
        folder_name = folders{folder_index}
        stimulus_realign(folder_name);
    end
end

%% STEP 9: Plot
if plotting
    'Plotting...'
    for folder_index = 1:folder_count
        folder_name = folders{folder_index}
        plot_image_directory(folder_name);
    end 
end
