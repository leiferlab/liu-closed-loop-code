% not used in paper

%folder_name = uigetdir
cd(folder_name) %open the directory of image sequence
allFiles = dir(); %get all the tif files
fps = 14;
pirouetteCount = 0;
tracksByVoltage = struct('maxVoltage', {}, 'LEDVoltages', {}, 'reversalCounts', {}, 'frameCounts', {});

for file_index = 1:length(allFiles)
    if allFiles(file_index).isdir && ~strcmp(allFiles(file_index).name, '.') && ~strcmp(allFiles(file_index).name, '..')
        folder = strcat(folder_name, filesep, allFiles(file_index).name)
        cd(folder)
        load('LEDVoltages.txt')
        maxVoltage = max(LEDVoltages);
        if isempty(find([tracksByVoltage.maxVoltage] == maxVoltage))
            %no entry with this max voltage before
            trackByVoltageIndex = length(tracksByVoltage) + 1;
            tracksByVoltage(trackByVoltageIndex).maxVoltage = maxVoltage;
            tracksByVoltage(trackByVoltageIndex).LEDVoltages = LEDVoltages;
            tracksByVoltage(trackByVoltageIndex).reversalCounts = zeros(1, length(LEDVoltages));
            tracksByVoltage(trackByVoltageIndex).frameCounts = zeros(1, length(LEDVoltages));            
        else
            trackByVoltageIndex = find([tracksByVoltage.maxVoltage] == maxVoltage);
        end
        [~, reversalCounts, frameCounts] = ReversalRate({folder}, 1);
        tracksByVoltage(trackByVoltageIndex).reversalCounts = tracksByVoltage(trackByVoltageIndex).reversalCounts + reversalCounts;
        tracksByVoltage(trackByVoltageIndex).frameCounts = tracksByVoltage(trackByVoltageIndex).frameCounts + frameCounts;
    end
end 

tracksByVoltage = nestedSortStruct(tracksByVoltage, 'maxVoltage'); %sort it
pirouetteCount = sum([tracksByVoltage.reversalCounts]);

mymarkers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
mycolors = jet(length(tracksByVoltage));
%plot stimuli
figure
scrollsubplot(3,1,1)
hold on;
for voltage_index = 1:length(tracksByVoltage)
%     plot(1/fps:1/fps:length(tracksByVoltage(voltage_index).LEDVoltages)/fps, tracksByVoltage(voltage_index).LEDVoltages, 'color', mycolors(voltage_index,:), 'marker', mymarkers{mod(voltage_index,numel(mymarkers))+1}, 'DisplayName', strcat(num2str(tracksByVoltage(voltage_index).maxVoltage*0.013756+0.000378), ' mW/mm^2'));
    plot(1/fps:1/fps:length(tracksByVoltage(voltage_index).LEDVoltages)/fps, tracksByVoltage(voltage_index).LEDVoltages, 'color', mycolors(voltage_index,:), 'marker', mymarkers{mod(voltage_index,numel(mymarkers))+1}, 'DisplayName', strcat(num2str(tracksByVoltage(voltage_index).maxVoltage), ' V'));
    %legend(num2str(tracksByVoltage(voltage_index).voltage));
end
xlabel('time (s)') % x-axis label
ylabel('voltage (V)') % y-axis label
legend('show');
hold off;

%plot cumulative reversals
scrollsubplot(3,1,2)
hold on;
for voltage_index = 1:length(tracksByVoltage)
    plot(1/fps:1/fps:length(tracksByVoltage(voltage_index).reversalCounts)/fps, cumsum(tracksByVoltage(voltage_index).reversalCounts), 'color', mycolors(voltage_index,:), 'marker', mymarkers{mod(voltage_index,numel(mymarkers))+1}, 'DisplayName', strcat(num2str(tracksByVoltage(voltage_index).maxVoltage), ' V'));
    %legend(num2str(tracksByVoltage(voltage_index).voltage));
end
xlabel(strcat('time in seconds (', num2str(pirouetteCount), ' reversals analyzed)')) % x-axis label
ylabel('cumulative reversal count') % y-axis label
legend('show');
hold off;

%plot reversal probability
scrollsubplot(3,1,3)
hold on;
for voltage_index = 1:length(tracksByVoltage)
    reversal_probability = tracksByVoltage(voltage_index).reversalCounts./tracksByVoltage(voltage_index).frameCounts;
    reversal_probability = smoothts(reversal_probability, 'g', 10*14, 50);
    plot(1/fps:1/fps:length(tracksByVoltage(voltage_index).reversalCounts)/fps, reversal_probability, 'color', mycolors(voltage_index,:), 'marker', mymarkers{mod(voltage_index,numel(mymarkers))+1}, 'DisplayName', strcat(num2str(tracksByVoltage(voltage_index).maxVoltage), ' V'));
    %legend(num2str(tracksByVoltage(voltage_index).voltage));
end
xlabel(strcat('time in seconds (', num2str(pirouetteCount), ' reversals analyzed)')) % x-axis label
ylabel('reversal probability') % y-axis label
legend('show');
hold off;