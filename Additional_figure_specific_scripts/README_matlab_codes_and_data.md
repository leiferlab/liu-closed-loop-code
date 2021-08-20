## Here are all the details to access MATLAB code and data related to the figure 4 of the publication "Closed-loop targeted optogenetic stimulation of C. elegans populations"

## 1)  List of MATLAB codes and explanations: ##

a)  `Liu_etal_closed_loop_and_open_loop_analysis_1.m` is the primary code which is used to generate the information presented in figure 4 and table 1 (Closed-loop optogenetic row). Before running this code make sure to add the github repo to the path by selecting "Add to Path" and then selecting "Selected Folders and Subfolders". 

(i) This code generate variables called `count_number_of_reversal_folder_turns` for closed loop experiments and `count_number_of_reversal_folder_rails` for open loop experiments. These variables (saved as .mat) are then fed to the `bar_plot_AML67_1.m` or `bar_plot_AML470_1.m` code to generate bar plots (i.e. figure 4).These two variables `count_number_of_reversal_folder_turns` and `count_number_of_reversal_folder_rails` have four columns. The first column represents the folder ID, the second column is the worm track id, the third column is the corresponding frame id and the fourth column has 0s and 1s. 0 means no reversal is detected while 1 means a reversal is detected.

b)  `Liu_etal_counting_stim_on_turn_onset_in_open_loop_experiment_1.m` is the special case of the above code (`Liu_etal_closed_loop_and_open_loop_analysis_1.m`). It was used only to find the number of stim on turn onset events in open-loop assay for filling the information in table 1 (Open-loop optogenetic row for this work).

c)  `bar_plot_AML67_1.m` is the code to plot the bar plot for AML67 data i.e. figure 4a. It reads the .mat files for both open loop and closed loop datasets for AML67. These datasets are present in the folder `data_to_plot_probability_of_reversal/AML67`.

d)  `bar_plot_AML470_1.m` is the code to plot the bar plot for AML470 data i.e. figure 4b. It reads the .mat files for both open loop and closed loop datasets for AML470. These datasets are present in the folder `data_to_plot_probability_of_reversal/AML470`.

e)  `bootstrap_mean_and_ci.m` is the MATLAB function to generate 95% CI using bootstrap. These functions are called in the above codes `bar_plot_AML67_1.m` and `bar_plot_AML470_1.m`.

f) `find_stim_on_turn_onsets.m` is a matlab function which makes sure that we are only detecting stims which are delivered after the worm's ellipse ratio goes below the ellipse ratio threshold. Thus it makes sure that we are only selecting stim on turn onsets.

## 2)  Details of the dataset ##
The final .mat files generated after running the `Liu_etal_closed_loop_and_open_loop_analysis_1.m` code are arranged into closed loop and open loop assay folder for both AML67 and AML470 strain in the folder `data_to_plot_probability_of_reversal/AML67`.

