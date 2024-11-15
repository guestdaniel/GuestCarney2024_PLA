# Guest and Carney (2024, JASA)
This is the repository for code necessary to reproduce the analyses and figures reported in:

```
Guest, D. R. and Carney, L. H. (2024). "A fast and accurate approximation of power-law adaptation for auditory computational models (L)." The Journal of the Acoustical Society of America, XX(XX), XX—XX, doi:XXX.
```

Corresponding author: Daniel Guest (daniel_guest@urmc.rochester.edu, https://github.com/guestdaniel)

# Files and paths
```
.  
├── figs                          # Figure raster files (.png) and layout files (.svg)
├── src                           # Primary source directory
│   ├── c                         # Modified source code in C/Mex/MATLAB for the AN models
│   ├── matlab                    # MATLAB code to run AN model simulations
│   |   └── figures               # Folder containing MATLAB scripts to generate Figs 1e–1g
│   ├── genfigs.jl                # Julia script to generate Figs 1a–1d
│   └── PowerlawApproximation.jl  # Julia package definition file
├── README.md                     # This documentation file
├── cfg.R                         # Script to provide config shared across all R files
├── LICENSE                       # License file for the code contained in this repository
└── Project.toml                  # Julia environment management file
```

# Instructions
## Reproducing paper figures
If you would like to reproduce the figures in the paper, there are two sets of instructions.
The first set handles reproducing Figs 1a–1d, which are implemented in Julia and demonstrate the theory behind the approximation and the results of the heuristic and optimize weight-selection strategies.
The second set handles reproducing Figs 1e–1g, which are implemented in MATLAB and demonstrate the results of the approximation in terms of speed and accuracy at the model-output stage. 

### Figs 1a–1d (Julia)
1. This repository is a Julia environment.
Start a Julia REPL and then navigate to this folder and activate the environment (press `]` in the REPL and then type `activate .`).
2. Use the script `src/genfigs.jl` to run the different script components that generate the various subfigures of Fig 1.
The first run will probably take quite some time as dependencies are installed and code is precompiled.
3. The resulting `.png` files should be saved to the `figs` folder.
Note that the names of the image files correspond to their final names in the paper for Figs 1a–1c (see the `saveplot` commands in the script).
Fig 1d is composed of the individual `.png` files generated at the end of the script and named `fig3a_*.png`.
4. The Julia code is fairly short and only a few key functions are actually used in the figure-generation process, but the code is also sparsely documented; if you have any questions, please reach out to the corresponding author.

### Figs 1e–1g (MATLAB)
1. Add this entire folder to your MATLAB path.
2. Unzip `zbc2014_with_fast_pla.zip` into its own folder and add that folder to your path.
3. Run `compile.m`, which will compile the `.c` files that implement the model into `.mex*` files that MATLAB can use.
4. Run the scripts contained in `src/matlab/figures` to generate the relevant figures. 
Note that some of these figures can take quite a while to run because true power-law adaptation becomes very slow for long-duration simulations.
`fig_main_a_v2.m` generates Fig 1g, `fig_main_b.m` generates Fig 1f, and `fig_main_c.m` generates Fig 1e.

## Using the new approximation
If you are a user of the Zilany, Bruce, and Carney (2014) auditory-nerve model and would like to use the updated power-law adaptation approximation reported in the paper, follow these instructions:
1. Unzip `zbc2014_with_fast_pla.zip` into its own folder and add that folder to your path.
2. Run `compile.m`, which will compile the `.c` files that implement the model into `.mex*` files that MATLAB can use.
3. As in previous variants of the 2014 model, `model_IHC` accepts a sound-pressure waveform and parameter values as input and returns an IHC response. The IHC is unchanged from previous versions.
4. Instead of calling `model_Synapse` as in previous versions of the 2014 model to compute the AN model, you now call `model_Synapse_2023`. 
It has similar inputs and outputs to the original `model_Synapse`, except that it additionally accepts vectors for time constants and weights for the slow and fast power-law adaptation pathway approximations described in Guest and Carney (2024). 
5. We provide a wrapper, `sim_an_zbc2023.m`, that makes handling this change somewhat easier.
The wrapper can either be used to avoid the need to learn the new system, as it handles most of the work for you, or you can work from the wrapper code as an example to make appropriate changes in your own model code if you so wish.
In the wrapper, arguments that were previously passed to `model_Synapse` as positional arguments are now passed as name-value pairs (e.g., `fiber_type=2`).
The argument `implnt`, which controls the implementation of power-law adaptation, has several options in the wrapper.
These options are documented in `sim_an_zbc2023.m` itself, but also below for convenience:
- `implnt=0` uses the original approximation scheme from the earlier papers.
- `implnt=1` uses a direct implementation of power-law adaptation.
- `implnt=2` uses the numerically optimized weights as reported in Guest and Carney (2024).
- `implnt=3` uses the heuristic weights as reported in Guest and Carney (2024).

