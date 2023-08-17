# RootWaterUptake
Matlab files for the concurrent solution of soil and root water transport within/towards  a root with length dependent axial and radial conductivities.
Tested with Matlab R2020a Update 6

Main file run.m calls either: 
RootStretch   = Iterated matrix method of solution, or                 
Bvpsuite2      =  DAE type reference solution using BVPSuite. 
Install BVPsuite from  https://github.com/NumODEsTUW/bvpsuite2.0.       

For more details see: 
Jan Graefe, Richard Pauwels, Michael Bitterlich
Water flow within and towards plant roots – a new concurrent solution
In silico Plants (in revision)
