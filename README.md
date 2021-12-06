## Abstract

Modern spatial transcriptomics techniques can reach ever-finer resolutions, and their resulting data
makes it possible to investigate the cellular interaction spaces between different tissues. 
In particular, the dynamics of these interaction spaces are well suited to the concept of RNA velocity. 
Here, we present a pipeline to map RNA velocities of these interaction spaces.

## Overview

![Flowchart of the interaction velocity pipeline](https://i.imgur.com/DALNgWj.png)

Scripts [`1_deconvolution.R`](1_deconvolution.R) and [`2_detect_markers.R`](2_detect_markers.R) run the `RCTD` and `Seurat` sections, and need to be run in that sequence.
After `velocyto` has also been run, [`3_velocity.py`](3_velocity.py) can be used to map velocities.

[`compare_markers.R`](compare_markers.R) and [`plot_overviews.R`](plot_overviews.R) are miscellaneous scripts to generate statistics and plots. 
