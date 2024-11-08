### Approximate Wave Simulations on Barbules
This directory contains our fast, approximate wave simulation code that computes scattering from different types of barbules. As verified with our full-wave simulations, our approximate simulation results are almost as accurate as the reference solutions, while the approximate methods gives rise to almost $1000 \times$ speedup.

#### Input Data
All the input data to our simulations are in the folders $\texttt{data1}-\texttt{data7}$. In our framework, the easiest way to generate these data folders is running our provided $\texttt{MATLAB}$ scripts, $\texttt{rockdoveProduction.m}$, $\texttt{starlingProduction.m}$, $\texttt{bronzewingProduction.m}$, $\texttt{hummingbirdProduction.m}$, $\texttt{mallardProduction.m}$, $\texttt{magpieProduction.m}$, and $\texttt{peacockProduction.m}$. 

By default, these scripts generate geometric models for multiple variants of each type of barbules; for instance, our $\texttt{data2}$ folder has 8 subfolders, representing European starling feathers of 8 different colors! Users should feel free to comment out different parts of the aforementioned $\texttt{MATLAB}$ scripts, to avoid generating barbules models they do not need. Users are also welcome to explore using different input geometric parameters to the $\texttt{MATLAB}$ functions.

Upon generating the data folders $\texttt{data1}-\texttt{data7}$, users can simply move these folders to our current $\texttt{2-FastSim}$ directory (as where they are for the time being). Inspecting one of the subfolders named after a barbule (e.g. $\texttt{rockdove1Forest}$) reveals three or four sub-directories: $\texttt{geometry}$, $\texttt{render}$, $\texttt{visual}$, and an additional $\texttt{coefs}$ sub-directory that contains photonic crystal reflectivity tables for photonic-type barbules. 

All the files containing the input geometric models are in the $\texttt{geometry}$ sub-directories, and these sub-directories in our repository already each contain geometric model files for 50 barbule instances (users are welcome to generate and use new ones). Moreover, the $\texttt{render}$ and $\texttt{visual}$ sub-directories are pregenerated to store the output files.

#### Command Line Arguments
Our wave simulation program takes the following command line arguments:

-a: Specifies the $\theta_i$ parameter that describes the incident direction (see our paper for the definition of $\theta_i$). In our framework, **this parameter must be an integer between 1 and 20**, and we have $\cos \theta_i = 0.05a$.

-b: Specifies the $\phi_i$ parameter that describes the incident direction (see our paper for the definition of $\phi_i$). Here, **this parameter must be an integer between 1 and 20**, and we have $\cos \phi_i = 0.05b$.

-i: Specifies the type of barbule. 1: rock dove barbules; 2: European starling barbules; 3: bronzewing barbules; 4: hummingbird barbules; 5: mallard barbules; 6: magpie barbules; 7: peacock barbules.

-l: Selects the simulated wavelength. **This argument is passed as an integer, with 50 possible values: 400, 406, 412, ..., 688, 694**, representing wavelengths in nanometers. 

-n: The number of instances simulated for each type barbule. We usually simulate 50 instances at a time.

-v: A binary argument with the value of 0 or 1. When $v=0$, we are in the material characterizing mode and our output files will contain BRDF parameters (as described in our paper) used for rendering. When $v=1$, we are in the trial mode, and our output files can be processed to generate some example scattering patterns (see below).

-z: The name of the simulated barbule, which needs to match the name of a subfolder in some $\texttt{dataX}$ folder. One example barbule name is $\texttt{rockdove1Forest}$.

After building the wave simulation program using
```
make
```
running the following command
```
./main -a 20 -b 20 -i 1 -l 400 -n 50 -v 0 -z rockdove1Forest
```
simulates the scattering of 50 green-colored rock dove barbules under normal incidence of 400nm wavelength light, and computes relevant BRDF parameters.

#### Trial Mode
When trying out the simulation code for the first time, users are recommended to start from the trial mode: $v=1$. Running a small set of simulations in this mode can generate nice, visualizable scattering patterns from barbules. Users can try running the following commands to test the system:
```
for l in {400..694..6}
do
  ./main -a 20 -b 20 -i 1 -l ${l} -n 50 -v 1 -z rockdove1Forest
done
```
will generate 50 output files in the $\texttt{visual}$ sub-directory in the $\texttt{data1/rockdove1Forest}$ folder.

Users can further verify their results by moving the entire $\texttt{data1}$ folder into the $\texttt{3-BRDF}$ directory and run the $\texttt{MATLAB}$ script:
```
load('D65.mat');
a = 20;
b = 20;
colors = viewPattern(1, "rockdove1Forest", 20 * (a - 1) + b, D65);
imshow(colors);
```
Users should expect to get an image that looks like the following:
<p align="center">
  <img src="https://github.com/blaire9989/FeatherLab/blob/main/3-BRDF/example.jpg" alt="gray" style="width:300px;"/>
</p>

The five sub-images above represent: the average scattering pattern (in RGB) from 50 rock dove barbules, the average scattering pattern fitted to an analytical form, and three individual scattering patterns from different BRDF instances. Each pattern that comes from our 2.5D simulation is a 1D function of the outgoing azimuthal angle $\phi_o$, as described in our paper, and in our framework we simulate a total of 400 outgoing azimuthal angles, uniformly distributed between $-90\degree$ and $90\degree$.

#### Material Characterizing Mode
For fully characterzing the scattering properties of one kind of barbule, users can simulate our generated barbule instances under all 50 wavelengths and 400 incident directions required in our system (and use $v = 0$). Theoretically, this can be done by running our simulation commands in loops:
```
for a in {1..20..1}
do
  for b in {1..20..1}
  do
    for l in {400..694..6}
    do
      ./main -a ${a} -b ${b} -i 1 -l ${l} -n 50 -v 0 -z rockdove1Forest
    done
  done
done
```
**However, users should NEVER attempt to run these commands on a single machine, as the underlying computations will take multiple days.** Users are instead recommended to submit jobs on a distributed cluster, such that these $20 \times 20 \times 50$ simulations run concurrently on many machines.

After the entire set of simulations terminates, users should expect to find 20000 BRDF parameters files in the $\texttt{render}$ sub-directory of the simulated barbule folder, as well as a file named $\texttt{noise.binary}$. These output files will be processed and compressed in the next stage of our appearance modeling pipeline, as discussed in the next tutorial.
