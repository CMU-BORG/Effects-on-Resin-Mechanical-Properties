# README for Dataset Described In:

**Impacts of sterilization method on material properties of 3D printable resins with prolonged exposure to cell culture environments**

AS Liao, MJ Bennington, K Dai, AB Irez, A Sun, S Schaffer, B Chopra, JM Seok, R Adams, YJ Zhang, VA Webster-Wood

2025

## Overview:

This work investigates the impact of common sterilization techniques on the mechanical properties of 3D printed resins after being exposed to mammalian cell culture conditions. The analysis in this paper is based on data reported in the data descriptor Liao <i>et al.</i>, Sci. Data, 2025 [doi](https://doi.org/10.1038/s41597-025-05738-7). This repository contains the code to reproduce the results reported in this paper, along with copies of the relevant data from Liao <i>et al.</i> 2025 and the generate figures presented in this paper. These files are contained within the directory [Mechanical Analysis](https://github.com/CMU-BORG/Effects-on-Resin-Mechanical-Properties/tree/main/Mechanical%20Analysis).

The [Mechanical Analysis](https://github.com/CMU-BORG/Effects-on-Resin-Mechanical-Properties/tree/main/Mechanical%20Analysis) directory contains two subdirectories:

 * [Liao et al. 2025 - Scientific Data Datasets](https://github.com/CMU-BORG/Effects-on-Resin-Mechanical-Properties/tree/main/Mechanical%20Analysis/Liao%20et%20al.%202025%20-%20Scientific%20Data%20Datasets): This subdirectory contains copies of the relevant mechanical datasets from Liao <i>et al.</i> 2025 used to conduct the analysis reported in this paper.

 * [Output Figures](https://github.com/CMU-BORG/Effects-on-Resin-Mechanical-Properties/tree/main/Mechanical%20Analysis/Output%20Figures): This subdirectory contains the figures generated during this analysis.

In addition to these subdirectory, there are three MATLAB scripts used to conduct the analysis in this paper:

 * [MechanicalAnalysis.m](https://github.com/CMU-BORG/Effects-on-Resin-Mechanical-Properties/blob/main/Mechanical%20Analysis/MechanicalAnalysis.m): This script generates the figures reported in this paper.

 * [RigidResinStressCalculations.m](https://github.com/CMU-BORG/Effects-on-Resin-Mechanical-Properties/blob/main/Mechanical%20Analysis/RigidResinStressCalculations.m): This script recalculates the stress-strain data for the rigid resin data from the raw force-displacement data reported in Liao <i>et al.</i> 2025.

 * [StatisticalAnalysis_SoftResins.m](https://github.com/CMU-BORG/Effects-on-Resin-Mechanical-Properties/blob/main/Mechanical%20Analysis/StatisticalAnalysis_SoftResins.m): This script conducts the bootstrapped statistical analysis of the soft resin Young's moduli.