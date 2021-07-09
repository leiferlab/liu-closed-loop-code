function [TestLNPStats] = ValidateLNP(Tracks, folder_indecies, folders, LNPStats)
%Predicts the behavioral rate based on LNP model parameters, not used in paper

%   Detailed explanation goes here
    numbins = 10;

    all_LEDVoltages = cell(1,length(folders));
    meanLEDVoltages = cell(1,length(folders));
    for folder_index = 1:length(folders)
        folder_name = folders{folder_index};
        % Load Voltages
        fid = fopen([folder_name, filesep, 'LEDVoltages.txt']);
        all_LEDVoltages{folder_index} = transpose(cell2mat(textscan(fid,'%f','HeaderLines',0,'Delimiter','\t'))); % Read data skipping header
        meanLEDVoltages{folder_index} = mean(all_LEDVoltages{folder_index});
        fclose(fid);
    end
    
    number_of_behaviors = size(Tracks(1).Behaviors,1);

    allLEDPower = [Tracks.LEDPower];
%     allLEDPower = [Tracks.LEDVoltages];
%     meanLEDPower = mean(allLEDPower);
%     stdLEDPower = std(allLEDPower);
% 
%     allLEDVoltages = [Tracks.LEDVoltages];
%     meanLEDVoltages = mean(allLEDVoltages);

    all_behaviors = horzcat(Tracks.Behaviors);
    linear_kernel = vertcat(LNPStats.linear_kernel);
    
    TestLNPStats(number_of_behaviors).linear_kernel = [];
    TestLNPStats(number_of_behaviors).trigger_count = [];
    
    TestLNPStats(number_of_behaviors).non_linearity_fit = [];
    TestLNPStats(number_of_behaviors).bin_edges = [];
    TestLNPStats(number_of_behaviors).filtered_signal_histogram = [];
    TestLNPStats(number_of_behaviors).filtered_signal_given_reversal_histogram = [];
    TestLNPStats(number_of_behaviors).non_linearity_fit = [];
    TestLNPStats(number_of_behaviors).non_linearity = [];
    TestLNPStats(number_of_behaviors).bin_centers = [];
    TestLNPStats(number_of_behaviors).errors = [];
    
    for behavior_index = 1:number_of_behaviors
        if isempty(nonzeros(linear_kernel(behavior_index,:)))
            %special case: flat kernel
            bin_centers = 0:numbins-1;
            non_linearity = zeros(1,numbins);
            TestLNPStats(behavior_index).bin_edges = bin_centers;
            TestLNPStats(behavior_index).filtered_signal_histogram = [];
            TestLNPStats(behavior_index).filtered_signal_given_reversal_histogram = [];
            TestLNPStats(behavior_index).non_linearity_fit = fit(bin_centers',non_linearity','exp1');
            TestLNPStats(behavior_index).non_linearity = non_linearity;
            TestLNPStats(behavior_index).bin_centers = bin_centers;
            TestLNPStats(behavior_index).errors = zeros(1,numbins);
        else
            %calculate the filtered LEDVoltages for all experiments
            all_filtered_LEDVoltages = cell(1,length(folders));
            for folder_index = 1:length(folders)
                %convolve the linear kernels with the input signal of LED voltages
                all_filtered_LEDVoltages{folder_index} = padded_conv(all_LEDVoltages{folder_index}-meanLEDVoltages{folder_index}, linear_kernel(behavior_index,:));
%                 figure
%                 hold on
%                 plot(all_LEDVoltages{folder_index}-meanLEDVoltages{folder_index})
%                 plot(all_filtered_LEDVoltages{folder_index})
            end

            %get all the filtered signals concatenated together
            all_filtered_signal = zeros(1, length(allLEDPower));
            current_frame_index = 1;
            for track_index = 1:length(Tracks)
                current_LEDVoltages2Power = Tracks(track_index).LEDVoltage2Power;
                filtered_signal = current_LEDVoltages2Power .* all_filtered_LEDVoltages{folder_indecies(track_index)}(Tracks(track_index).Frames);
                all_filtered_signal(current_frame_index:current_frame_index+length(Tracks(track_index).Frames)-1) = filtered_signal;
                current_frame_index = current_frame_index+length(Tracks(track_index).Frames);
%                 plot(filtered_signal)
            end

%             plot(all_filtered_signal)
            
            %make histogram of filtered signal
%             current_bin_edges = LNPStats(behavior_index).bin_edges; %use the previously fitted bin edges
            current_bin_edges = linspace(min(all_filtered_signal), max(all_filtered_signal), numbins+1);
            current_bin_edges(end) = current_bin_edges(end) + 1;
            TestLNPStats(behavior_index).bin_edges = current_bin_edges;
            [current_filtered_signal_histogram, ~] = histc(all_filtered_signal, current_bin_edges);
            current_filtered_signal_histogram = current_filtered_signal_histogram(1:end-1);
            TestLNPStats(behavior_index).filtered_signal_histogram = current_filtered_signal_histogram;
            
            %get histogram of filtered_signal given a reversal
            current_filtered_signal_given_behavior = all_filtered_signal(all_behaviors(behavior_index,:));
            current_filtered_signal_given_behavior_histogram = histc(current_filtered_signal_given_behavior, current_bin_edges);
            current_filtered_signal_given_behavior_histogram = current_filtered_signal_given_behavior_histogram(1:end-1);
            TestLNPStats(behavior_index).filtered_signal_given_reversal_histogram = current_filtered_signal_given_behavior_histogram;
        %     figure
        %     bar(bin_edges(1:end-1), filtered_signal_given_reversal_histogram');
        %     set(gca,'XTick',round(bin_edges*100)/100)
            [TestLNPStats(behavior_index).non_linearity_fit, TestLNPStats(behavior_index).non_linearity, ...
                TestLNPStats(behavior_index).bin_centers, TestLNPStats(behavior_index).errors] = ...
                fit_nonlinearity(current_filtered_signal_given_behavior_histogram, current_filtered_signal_histogram, current_bin_edges);
        end
        TestLNPStats(behavior_index).linear_kernel = LNPStats(behavior_index).linear_kernel;
        TestLNPStats(behavior_index).trigger_count = sum(all_behaviors(behavior_index,:));
    end
end

