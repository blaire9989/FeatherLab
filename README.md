## FeatherLab
Open-source code base of the SIGGRAPH Asia 2024 paper, "Appearance Modeling of Iridescent Feathers with Diverse Nanostructures." For more details on our wave simulations and appearance modeling pipeline, please refer to the paper.

This repository contains our procedural geometric models for feather barbules, our fast approximate and slow reference wave simulators for computing barbule scattering, as well as our feather BRDF generation pipeline.

### Code Base Overview
This repository is organized with four sub-directories, containing code written in $\texttt{C++}$ and $\texttt{MATLAB}$, used for different stages in iridescent feather appearance modeling. The numbered directories, namely $\texttt{1-Geometry}$, $\texttt{2-FastSim}$, and $\texttt{3-BRDF}$, contain code for procedural barbule geometry modeling (our first step for appearance modeling), efficient barbule scattering simulation (our second step), and iridescent feather BRDF generation (our third step). The directory $\texttt{FullWave}$ contains our full-wave simulation code for computing barbule scattering, which we used as our reference tool when developing our fast, approximate simulator, but did not use for appearance modeling due to its expensive nature.

### Tutorial and Support
The guidelines for using code from each stage of our appearance modeling are listed below, and can be found in the respective sub-directories:

Procedural geometric modeling: please see description [here](https://github.com/blaire9989/FeatherLab/blob/main/1-Geometry/README.md).

Fast wave simulation: please see tutorial for building [here](https://github.com/blaire9989/FeatherLab/blob/main/2-FastSim/README.md).

BRDF generation pipeline: please see description [here](https://github.com/blaire9989/FeatherLab/blob/main/3-BRDF/README.md).

Reference simulation: please see tutorial [here](https://github.com/blaire9989/FeatherLab/blob/main/FullWave/README.md), while keeping in mind that these simulations are slow.

For questions or collaborations, please contact the authors at yy735@cornell.edu.
