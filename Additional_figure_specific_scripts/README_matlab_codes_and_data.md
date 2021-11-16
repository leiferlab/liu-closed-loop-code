## Here are all the details to access MATLAB code and data related to the figure 4 of the publication "Closed-loop targeted optogenetic stimulation of C. elegans populations"

## 1)  List of MATLAB codes and explanations: ##

a)  `Liu_etal_closed_loop_and_open_loop_analysis_1.m` is the primary code which is used to generate the information presented in figure 4, table 1 (Closed-loop optogenetic row), and supplementary rable S1. Before running this code make sure to add the github repo to the path by selecting "Add to Path" and then selecting "Selected Folders and Subfolders". 

(i) In the closed loop experiments, this code generates a variable called `count_number_of_reversal_folder_turns`. These variables (saved as .mat) are then fed to the `bar_plot_AML67_1.m` or `bar_plot_AML470_1.m` code to generate bar plots (i.e. figure 4). The variables `count_number_of_reversal_folder_turns` has four columns. The first column represents the folder ID, the second column is the worm track id, the third column is the corresponding frame of the track and the fourth column has 0s and 1s. 0 means no reversal is detected while 1 means a reversal is detected. 
(ii) Note: In the open loop datasets, the variable `count_number_of_reversal_folder_rails` has five columns. The first column represents the folder ID, the second column is the worm track id, the third column is the corresponding frame of the track and the fourth column has 0s and 1s. 0 means no reversal is detected while 1 means a reversal is detected. The fifth column is the frame number of stim onset corresponding to the actual time in the experimental assay.

b)  `Liu_etal_counting_stim_on_turn_onset_in_open_loop_experiment_1.m` is the special case of the above code (`Liu_etal_closed_loop_and_open_loop_analysis_1.m`). It was used only to find the number of stim on turn onset events in open-loop assay for filling the information in table 1 (Open-loop optogenetic row for this work).

c)  `bar_plot_AML67_1.m` is the code to plot the bar plot for AML67 data i.e. figure 4a. It reads the .mat files for both open loop and closed loop datasets for AML67. These datasets are present in the folder `data_to_plot_probability_of_reversal/AML67`.

d)  `bar_plot_AML470_1.m` is the code to plot the bar plot for AML470 data i.e. figure 4b. It reads the .mat files for both open loop and closed loop datasets for AML470. These datasets are present in the folder `data_to_plot_probability_of_reversal/AML470`.

e)  `bootstrap_mean_and_ci.m` is the MATLAB function to generate 95% CI using bootstrap. These functions are called in the above codes `bar_plot_AML67_1.m` and `bar_plot_AML470_1.m`.

f) `find_stim_on_turn_onsets.m` is a matlab function which makes sure that we are only detecting stims which are delivered after the worm's ellipse ratio goes below the ellipse ratio threshold. Thus it makes sure that we are only selecting stim on turn onsets.

g) `GenerateFigures_S2.m` is the code to generate the supplementary figure S2 and the data corresponding to Stimulation durign turning vs forward enteries in supplementary table S1 (columns cumulative recording length, animals per frame).

h) `get_behavior_triggers.m` is a function used by the above codes.

i) `make_list_of_all_folders.m` is a code to create a cell of all the folders with specific criteria in a given folder (for e.g. all open loop AML67 experiments).

j) to_plot_habituation_3.m is the code to plot the habituation curve for both AML67 and AML470 as shown in supplementary figure S5.


## 2)  Details of the dataset ##
1) The final .mat files generated after running the `Liu_etal_closed_loop_and_open_loop_analysis_1.m` code are arranged into closed loop and open loop assay folder for both AML67 and AML470 strain in the folder `data_to_plot_probability_of_reversal/AML67`.

2) The same code is ran for the open loop dataset for 80uW/mm^2 and 3 sec conditions for AML67 and AML470. The output .mat files are stores in `data_to_plot_habituation_curve` fodler.

3) For the ease of users, a list of all the folders for specific experimental conditions (for e.g. all AML67 open loop experiments) are saved as .mat files in the folder `folder_name_list_for_fig_S2`. The user can load the corresponding .mat file in `GenerateFigures_S2.m` to get the numbers specific to a strain and experimental condition.
