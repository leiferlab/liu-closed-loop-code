# This is the code used to analyze the data presented in: Liu, Kumar, Sharma and Leifer, "A high-throughput method to deliver targeted optogenetic stimulation to moving *C. elegans* population" available at [https://arxiv.org/abs/2109.05303](https://arxiv.org/abs/2109.05303) and forthcoming in *PLOS Biology*.
The raw data used in this publication is available on ieee-dataport at [https://ieee-dataport.org/open-access/data-high-throughput-method-deliver-targeted-optogenetic-stimulation-moving-c-elegans](https://ieee-dataport.org/open-access/data-high-throughput-method-deliver-targeted-optogenetic-stimulation-moving-c-elegans) 

Please see the supplement of the publication for the implementation details.
Additional insights can be found in the final chapter of ["C. ELEGANS BEHAVIORS AND THEIR MECHANOSENSORY DRIVERS"](https://dataspace.princeton.edu/handle/88435/dsp01tt44pq78z)

The code is roughly divided into 3 parts: Real-time recording and stimulation in LabVIEW, post-processing in MATLAB, and generating figures in MATLAB

## Section 1: LabVIEW Real-time

Requirements:
 - LabVIEW 2019 64-bit Windows
 - LabVIEW package manager
 - LabVIEW OpenG add-on
 - HDF5 v1.8.18+ (latest v1.8 but not v1.10)
 - h5labview

To run the real-time LabVIEW code, open the LabVIEW project at `\LabviewVIs\ProjectAPI.lvproj` and then select the appropriate experimental protocol vi. For example, the head/tail stimulation is done with RunHeadandTailRailswithDelays.vi.

Associated hardware in the correct configuration is required for the LabVIEW code to run properlly.

## Section 2: Post Processing 

Requirements:
MATLAB 2016a

To conduct the post-processing on the raw datasets, run `\ProcessDateDirectory.m` and select the dataset when prompted.

## Section 3: Generating figures

a) Please read the `instruction_to_generate_figures.csv` to get detailed explanation of which MATLAB script is used to generate a specific figure. It also mentions the combination of "Tags" which needs to be selected by the user to select the correct experimental folders. Additionally, please make sure to add the github repo to the path by selecting "Add to Path" and then selecting "Selected Folders and Subfolders". 

b) After all the dataset that has gone through post-processing, run `plot_figures_2_3_S1_S3_S4_S6.m` to generate figures 2, 3, S1, S3, S4, and S6 as shown in the publication by providing the corresponding input dataset folder using the "Tags".   

c) `plot_figures_4_S5.m` is the primary code which is used to generate the information presented in figure 4, figure S5, table 1 (Closed-loop optogenetic row), and supplementary table S1. 

d)`Liu_etal_counting_stim_on_turn_onset_in_open_loop_experiment_1.m` is the special case of the above code (`plot_figures_4_S5.m`). It was used only to find the number of stim on turn onset events in open-loop assay for filling the information in table 1 (Open-loop optogenetic row for this work).

e) `plot_figure_S2.m` is the code to generate the supplementary figure S2 and the data corresponding to Stimulation during turning vs forward enteries in supplementary table S1 (columns cumulative recording length, animals per frame).