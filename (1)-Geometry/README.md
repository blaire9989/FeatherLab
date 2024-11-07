### Procedural Geometric Modeling
We procedurally generate 2D geometric models for representing barbules cross sections from seven types of birds with iridescent feathers: 

(1) Rock dove

(2) European starling

(3) The common bronzewing

(4) Anna's hummingbird (the code can model many other hummingbird barbules as well)

(5) The common mallard

(6) Black-billed magpie

(7) Indian peafowl (the common peacock)

The code for generating 2D cross sections are written as $\texttt{MATLAB}$ functions. Each run of each $\texttt{MATLAB}$ function creates one specific instance of a barbule cross section. The characteristic features, such as average layer thicknesses and average melanosome diameters, in each type of barbule are referenced from the ornithology literature. These average layer thicknesses and average melanosome sizes are approximately invariant from barbule to barbule, while each generated barbule instance is unique. For instance, each specific melanosome in each barbule has a randomized diameter, and the layer structures are modeled as random rough surfaces.

The $\texttt{MATLAB}$ scripts with names in the form of $\texttt{xxxProduction.m}$ generates 50 barbule instances of a specific type, and stores the output files in appropriately created folders and subfolders. Users can enter each of the 7 numbered sub-directories (e.g. $\texttt{(1)-rockdove}$) from $\texttt{MATLAB}$, and run the production MATLAB scripts. As an example, the command
```
rockdoveProduction
```
creates a folder named $\texttt{data1}$ with 6 subfolders, each containing geometric models for a type of rock dove barbule with a specific color—different shades of green and purple in this case. As will be explained in the next tutorial, these $\texttt{dataX}$ folders ($X = 1, 2, ..., 7$) will be dragged into the directory containing our wave simulation code and used as simulation inputs.
