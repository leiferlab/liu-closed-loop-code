clear
clc
close all

cd('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/reviewer_comments/tracking_frame_rate/folder_name_list_for_fig_S2/all_data')

matfile_names = dir('*.mat');
folder_list=[];
for i=1:size(matfile_names,1)
    matfile_names(i).name
    load(matfile_names(i).name)
    folder_list=[folder_list folders];
    folders=[];
end

folders=[];
folders=folder_list;
