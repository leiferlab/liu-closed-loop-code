# Code for "Closed-loop targeted optogenetic stimulation of C. elegans populations"

Please see the supplement of the publication for the implementation details.

Additional insights can be found in the final chapter of ["C. ELEGANS BEHAVIORS AND THEIR MECHANOSENSORY DRIVERS"](https://dataspace.princeton.edu/handle/88435/dsp01tt44pq78z)

LabVIEW requirements:
LabVIEW 2019 64-bit Windows
LabVIEW package manager
LabVIEW OpenG add-on
HDF5 v1.8.18+ (latest v1.8 but not v1.10)
h5labview

Post analysis requirements:
MATLAB 2016a

To run the real-time LabVIEW code, open the LabVIEW project at \LabviewVIs\ProjectAPI.lvproj
Associated hardware in the correct configuration is required for proper functionality.

To conduct the post-processing on the raw datasets, run \ProcessDateDirectory.m and select the dataset when prompted.

After all the dataset has gone through post-processing, run \GenerateFigures.m to generate the figures as shown in the publication.