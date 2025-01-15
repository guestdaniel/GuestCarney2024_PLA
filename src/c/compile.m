% Compile source code into MEX functions.  Requires C compiler.
% Run "mex -setup" first.
mex model_IHC.c complex.c
mex model_Synapse_2023.c complex.c 
mex model_Synapse_v2025a.c complex.c 