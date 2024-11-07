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

The $\texttt{MATLAB}$ scripts with names in the form of "xxxProduction.m"
