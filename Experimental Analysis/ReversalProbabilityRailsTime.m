% not used in paper

folder_name = uigetdir
cd(folder_name) %open the directory of image sequence
allFiles = dir(); %get all the tif files
fps = 14;
allTracks = struct('maxVoltage', {}, 'LEDVoltages', {}, 'reversalCounts', {}, 'frameCounts', {}, 'order', {});
order = 1;

for file_index = 1:length(allFiles)
    if allFiles(file_index).isdir && ~strcmp(allFiles(file_index).name, '.') && ~strcmp(allFiles(file_index).name, '..')
        folder = strcat(folder_name, filesep, allFiles(file_index).name);
        cd(folder)
        load('LEDVoltages.txt')
        maxVoltage = max(LEDVoltages);
        allTracks(file_index).maxVoltage = maxVoltage;
        allTracks(file_index).LEDVoltages = LEDVoltages;
        allTracks(file_index).order = order;
        order = order + 1;
        [~, reversalCounts, frameCounts] = ReversalRate({folder}, 1);
        allTracks(file_index).reversalCounts = reversalCounts;
        allTracks(file_index).frameCounts = frameCounts;
    end
end 

%allTracks = nestedSortStruct(allTracks, 'maxVoltage'); %sort it
pirouetteCount = sum([allTracks.reversalCounts]);

figure(1);
hold on;
grid on;
orderInd = cell(10, 1);
counter = 1;
myColors = jet(10);

for track = 1:length(allFiles)
    if (allTracks(track).maxVoltage == 5)
        reversalRate = allTracks(track).reversalCounts ./ allTracks(track).frameCounts;
        plot(1/fps:1/fps:length(allTracks(track).reversalCounts)/fps, cumsum(reversalRate), 'color', myColors(counter, :), 'LineWidth', 2);
        orderInd{counter} = num2str(allTracks(track).order);
        counter = counter + 1;
    end
    if (counter > 10)
        break
    end
end
legend(orderInd);

figure(2);
hold on;
grid on;
orderInd = cell(10, 1);
counter = 1;
myColors = jet(10);

for track2 = track:length(allFiles)
    if (allTracks(track2).maxVoltage == 5)
        reversalRate = allTracks(track2).reversalCounts ./ allTracks(track2).frameCounts;
        plot(1/fps:1/fps:length(allTracks(track2).reversalCounts)/fps, cumsum(reversalRate), 'color', myColors(counter, :), 'LineWidth', 2);
        orderInd{counter} = num2str(allTracks(track2).order);
        counter = counter + 1;
    end
    if (counter > 10)
        break
    end
end
legend(orderInd);