### Approximate Wave Simulations on Barbules
This directory contains our fast, approximate wave simulation code that computes scattering from different types of barbules. As verified with our full-wave simulations, our approximate simulation results are almost as accurate as the reference solutions, while the approximate methods gives rise to almost $1000 \times$ speedup.

#### Input Data
All the input data to our simulations are in the folders $\texttt{data1}-\texttt{data7}$. In our framework, the easiest way to generate these data folders is running our provided $\texttt{MATLAB}$ scripts, $\texttt{rockdoveProduction.m}$, $\texttt{starlingProduction.m}$, $\texttt{bronzewingProduction.m}$, $\texttt{hummingbirdProduction.m}$, $\texttt{mallardProduction.m}$, $\texttt{magpieProduction.m}$, and $\texttt{peacockProduction.m}$. 

By default, these scripts generate geometric models for multiple variants of each type of barbules; for instance, our $\texttt{data2}$ folder has 8 subfolders, representing European starling feathers of 8 different colors! Users should feel free to comment out different parts of the aforementioned $\texttt{MATLAB}$ scripts, to avoid generate barbules models they do not need. Users are also welcome to explore using different input geometric parameters to the $\texttt{MATLAB}$ functions.

Upon generating the data folders $\texttt{data1}-\texttt{data7}$, users can simply move these folders to our current $\texttt{2-FastSim}$ directory (as where they are for the time being). Inspecting one of the subfolders named after a barbule (e.g. $\texttt{rockdove1Forest}$) reveals three or four sub-directories: $\texttt{geometry}$, $\texttt{render}$, $\texttt{visual}$, and an additional $\texttt{coefs}$ sub-directory that contains photonic crystal reflectivity tables for photonic-type barbules. 

All the files containing the input geometric models are in the $\texttt{geometry}$ sub-directories, and these sub-directories in our repository already each contain geometric model files for 50 barbule instances (users are welcome to generate and use new ones). Moreover, the $\texttt{render}$ and $\texttt{visual}$ sub-directories are pregenerated to store the output files.

#### Command Line Arguments
Our wave simulation program takes the following command line arguments:

$\texttt{-a}$: Specifies the $\theta_i$ parameter that describes the incident direction (see our paper for the definition of $\theta_i$). In our framework, **this parameter must be an integer between 1 and 20**, and we have $\cos \theta_i = 0.05a$.

$\texttt{-b}$: Specifies the $\phi_i$ parameter that describes the incident direction (see our paper for the definition of $\phi_i$). In our framework, **this parameter must be an integer between 1 and 20**, and we have $\cos \phi_i = 0.05b$.
