# Code for "Closed-loop targeted optogenetic stimulation of C. elegans populations"

Please see the supplement of the publication for the implementation details.

Additional insights can be found in the final chapter of ["C. ELEGANS BEHAVIORS AND THEIR MECHANOSENSORY DRIVERS"](https://dataspace.princeton.edu/handle/88435/dsp01tt44pq78z)

The code is roughly divided into 3 parts: Real-time recording and stimulation in LabVIEW, post-processing in MATLAB, and generating figures in MATLAB

# LabVIEW Real-time

requirements:
 - LabVIEW 2019 64-bit Windows
 - LabVIEW package manager
 - LabVIEW OpenG add-on
 - HDF5 v1.8.18+ (latest v1.8 but not v1.10)
 - h5labview

To run the real-time LabVIEW code, open the LabVIEW project at \LabviewVIs\ProjectAPI.lvproj and then select the appropriate experimental protocol vi. For example, the head/tail stimulation is done with RunHeadandTailRailswithDelays.vi.

Associated hardware in the correct configuration is required for the LabVIEW code to run properlly.

# Post Processing 

requirements:
MATLAB 2016a

To conduct the post-processing on the raw datasets, run \ProcessDateDirectory.m and select the dataset when prompted.

# Generating figures

After all the dataset has gone through post-processing, run \GenerateFigures.m to generate the figures as shown in the publication.