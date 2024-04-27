# SpatialKcorePercolation
SpatialKcorePercolation is a Python project for generating spatial networks and performing k-core percolation simulations to analyze percolation behavior. It includes functions for generating spatial networks with various parameters, performing k-core percolation simulations, and plot the simulation results.

## SpatialKcorePercolation Project Overview
- **Objective:** The **SpatialKcorePercolation** project is designed to analyze the robustness of spatial k-core systems by generating spatial networks and performing k-core percolation simulations. 

- **Repository Content:** This repository contains code for network generation, k-core percolation simulations and visualization. All code is organized according to the figures presented in the associated paper
[Nucleation phenomena and extreme vulnerability of spatial k-core systems](https://doi.org/10.48550/arXiv.2311.13579). 

- **Organization:** Each directory contains the numerical simulation and plot code necessary to generate the corresponding figure, and the parameters included in the current code represent those used to generate the figures in the paper.

- **Additional Functionalities:** Additional commonly used functionalities are written in the `utils` directory, specifically in the `kcorePercolation.py` file.

- **Replicating Results:** If you want to replicate the results presented in the paper, simply download the pre-generated data from [Mendeley Data](10.17632/jkvk97nfjc.1) and directly execute the corresponding `plot_figurex.py`. The data and code paths match, so you only need to place the data in the same directory as the code to run it.

- **Note:** Generating the simulation data for each figure may require a significant amount of time due to the larger system $$L \times L = 10^6$$, which also depending on your machine's configuration. 

## Code
### Network Generation
**Usage:** 
1. Clone or download this repository to your local machine.
2. Navigate to the `/figure1/code/` directory.
3. Open the `NetworkGeneration.py` script in a text editor.
4. Locate the `rootpath` and `packagepath` variable at the beginning of the script.
5. Modify the `rootpath` and `packagepath` variable to specify the local directory where you want to save the generated networks.
6. Save the changes and close the text editor.
7. Run the `NetworkGeneration.py` script.
This script generates networks with various parameters, including different zeta values, average degrees, and lattice lengths. 

### K-Core Percolation Simulation
**Execute Numerical Simulations**:
- After generating the spatial networks, execute the k-core percolation simulations by running the code located in each directory, e.g. `\figure1\code\KCorePercolationSimulation.py`, `\figure2\code\binary_seachpc_kcore.py`. 

**Usage:**
1. Ensure that the required spatial networks are generated.
2. Navigate to the directory containing the code for the desired figure.
3. Verify that the path configurations in the code are correct.
4. Run the code to execute the k-core percolation simulations and generate the data.

### Plotting Figures
**Usage:**
1. Open the Python script containing the plotting functions (e.g., `\figure3\code\plot_figure3.py`) in a text editor.
2. Ensure that the `packagepath`, `networkpath`, `resultpath` are correctly configured.
3. Execute the plotpinfty function with the appropriate parameters to generate the plot.

## Prerequisites
Before running the code, make sure the following are installed:

- Python 3.x
- NumPy
- NetworkX
- Matplotlib
- Seaborn
- Pandas
- itertools (included in Python standard library)
- multiprocessing (included in Python standard library)

## Using Pre-generated Data
Alternatively, if you prefer not to generate the networks yourself and run k-core percolation simulation due to the long runtime, pre-generated network and k-core percolation data is available in the Mendeley dataset. 
Follow the steps below to download the data:
1. Access the [Mendeley dataset](10.17632/jkvk97nfjc.1).
2. Download the network data files corresponding to your desired parameters.
3. Extract the downloaded files to a directory of your choice.

### Loading Pre-generated Data
To load pre-generated network data and k-core percolation simulation data into your Python environment, you can use the load function from the kcorePercolation module. 
Here's an example of how to do it:
### Example: Loading pre-generated network data
```python 
import sys
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

network_data = kp.load('../figure1/network/NetID0_avgk10_zeta7_spatialNet.pkl')## Replace '../figure1/network/NetID0_avgk10_zeta7_spatialNet.pkl' with the actual path to the downloaded network data file.
```

## License
> This project is licensed under the MIT License. Feel free to use, modify, and distribute the code as needed.
> If you use this code in your research, please cite the following paper:
> **Nucleation phenomena and extreme vulnerability of spatial k-core systems**
> - [Link to the paper](https://doi.org/10.48550/arXiv.2311.13579)
> ```
> @misc{xue2023nucleation,
>       title={Nucleation phenomena and extreme vulnerability of spatial k-core systems}, 
>       author={Leyang Xue and Shengling Gao and Lazaros K. Gallos and Orr Levy and Bnaya Gross and Zengru Di and Shlomo Havlin},
>       year={2023},
>       eprint={2311.13579},
>       archivePrefix={arXiv},
>       primaryClass={physics.soc-ph}
> }
> ```



