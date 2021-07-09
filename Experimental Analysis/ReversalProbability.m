% not used in paper
tracksByVoltage = struct('voltage', {0}, 'reversalCount', {0}, 'totalCount', {0});

for track = 1:length(allTracks)
    pirouettes = allTracks(track).Pirouettes;
    LEDVoltages = allTracks(track).LEDVoltages;
    maxVoltage = max(LEDVoltages);
    if length(find([tracksByVoltage.voltage] == maxVoltage)) == 0
        %no entry with this max voltage before
        trackByVoltageIndex = length(tracksByVoltage) + 1;
        tracksByVoltage(trackByVoltageIndex).voltage = maxVoltage;
        tracksByVoltage(trackByVoltageIndex).reversalCount = 0;
        tracksByVoltage(trackByVoltageIndex).totalCount = 0;
    else
        trackByVoltageIndex = find([tracksByVoltage.voltage] == maxVoltage);
    end
    tracksByVoltage(trackByVoltageIndex).totalCount = tracksByVoltage(trackByVoltageIndex).totalCount + 1;
    lightPulse = find(LEDVoltages == maxVoltage, 1);
    for pirouette_index = 1:size(pirouettes,1)
        pirouetteStart = pirouettes(pirouette_index,1);
        %get where the voltage transition is, if there is any:
        if lightPulse - pirouetteStart < 2 && lightPulse - pirouetteStart > -10
            %the worm has reversed because of light
            tracksByVoltage(trackByVoltageIndex).reversalCount = tracksByVoltage(trackByVoltageIndex).reversalCount + 1;
        else
            %the worm has reversed spontanously
            
        end
    end
end

tracksByVoltage = nestedSortStruct(tracksByVoltage, 'voltage'); %sort it

tracksByVoltage = transpose(squeeze(cell2mat(struct2cell(tracksByVoltage))));
plot(tracksByVoltage(2:end,1), tracksByVoltage(2:end,2)/tracksByVoltage(2:end,3))
%legend(num2str(tracksByVoltage(voltage_index).voltage));
