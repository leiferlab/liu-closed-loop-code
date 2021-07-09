function parameters = load_parameters(folder_name)
%loads the parameter structure based on the tags and labview parameters in
%the experimental folder and parameters.csv settings file
    starting_column_index = 4;
    type_index = 3;
    parameters = load('EigenVectors.mat'); %load eigenvectors for eigenworms
    parameters.SampleRate = 30; %this will be overwritten from an experiment
    if nargin < 1
        tags = {};
    else
        %load the tags
        if exist([folder_name, filesep, 'tags.txt'],'file')
            %this is an image folder with tags
            tags = textread([folder_name, filesep, 'tags.txt'], '%s', 'delimiter', ' ');

            % load the labview parameters first, then replace with the latest csv
            param_table = readtable([folder_name, filesep, 'labview_parameters.csv'],'ReadVariableNames',false);
            for parameter_index = 2:size(param_table,1)
                value = param_table{parameter_index,starting_column_index}{1,1};
                value_type = param_table{parameter_index,type_index}{1,1};
                param_name = param_table{parameter_index,1}{1,1};
                if ~isempty(value)
                    if strcmp(value_type, 'String')
                        %string
                        parameters.(param_name) = value;
                    elseif strcmp(value_type, 'Numeric')
                        %numeric value, check if array
                        k = strfind(value,';');
                        if isempty(k)
                            parameters.(param_name) = str2double(value);
                        else
                            value_cells = strsplit(value, ';');
                            double_values = zeros(1,numel(value_cells));
                            for index = 1:numel(value_cells)
                                double_values(index) = str2double(value_cells{index});
                            end
                            parameters.(param_name) = double_values;
                        end
                    end
                end
            end
        else
            tags = {};
        end
    end

    param_table = readtable('parameters.csv','ReadVariableNames',false);

    %load the default parameters first
    for tag_index = starting_column_index:size(param_table,2)
        current_tags = param_table{1,tag_index}{1,1};
        current_tags = strsplit(current_tags,';');
        if tag_index == starting_column_index || all(ismember(current_tags, tags)) 
            %load the default and the correct tags
            for parameter_index = 2:size(param_table,1)
                value = param_table{parameter_index,tag_index}{1,1};
                value_type = param_table{parameter_index,type_index}{1,1};
                param_name = param_table{parameter_index,1}{1,1};
                if ~isempty(value)
                    if strcmp(value_type, 'String')
                        %string
                        parameters.(param_name) = value;
                    elseif strcmp(value_type, 'Numeric')
                        %numeric value
                        parameters.(param_name) = str2double(value);
                    end
                end
            end            
        end
    end

    if ischar(parameters.Mask) && exist(parameters.Mask, 'file')
       %get the mask
       parameters.Mask = imread(parameters.Mask); 
    end

%     if ischar(parameters.power500) && exist(parameters.power500, 'file')
%        %get the power distribution
%        load(parameters.power500);
%        parameters.power500 = power500; 
%     end
%     
   
end

