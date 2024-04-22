# SpatialKcorePercolation
SpatialKcorePercolation is a Python project for generating spatial networks and performing k-core percolation simulations to analyze percolation behavior. It includes functions for generating spatial networks with various parameters, performing k-core percolation simulations, and analyzing simulation results.

## Introduction
This project is designed to generate spatial networks, perform k-core percolation simulations, analyze simulation results, and plot figures. It consists of multiple scripts organized into different folders.

## Network Generation
1. **Generate Networks**: 
Before running the simulation, generate spatial networks using appropriate parameters. Adjust the parameters in the script (`NetworkGeneration.py`) according to your requirements. If you would like to replicate the results or analyze the generated networks, follow the instructions below.
2. **Path Configuration**:
Ensure that the path to the `kcorePercolation` package is correctly specified in the script for importing custom utility functions.
3. **Usage**:
After setting the path and parameters, execute the script to generate spatial networks. This step is necessary before running the k-core percolation simulation.

### Generating Networks
To generate networks using this code:
1. Clone or download this repository to your local machine.
2. Navigate to the code directory.
3. Open the `NetworkGeneration.py` script in a text editor.
4. Locate the `rootpath` and `packagepath` variable at the beginning of the script.
5. Modify the `rootpath` and `packagepath` variable to specify the local directory where you want to save the generated networks.
6. Save the changes and close the text editor.
7. Run the `NetworkGeneration.py` script.
This script generates networks with various parameters, including different zeta values, average degrees, and lattice lengths. 
Note that the generation process may take a significant amount of time depending on your machine's specifications and the chosen parameters.

## K-Core Percolation Simulation
1. **Run K-Core Percolation Simulation**:
   After generating networks, run the k-core percolation simulation by executing the script (`KCorePercolation.py`). The simulation iterates over different parameters such as k-core value, zeta value.
2. **Analysis**:
   Once the simulation is complete, analyze the results to understand the percolation behavior. Results include the largest connected component size and the number of iterations.
3. **Downloading Pre-computed Results**:
   If running the simulation takes too much time, pre-computed percolation results are available for download in the Mendely dataset associated with this project. Use these results for analysis without running the code.

## Plotting Figures
Once you have completed the network generation and k-core percolation simulation, you can plot the figures to visualize the results.
1. Open the Python script containing the plotting functions (e.g., PlottingFunctions.py) in a text editor.
2. Locate the plotpinfty function within the script.
3. Ensure that the function parameters are correctly configured。
4. Execute the plotpinfty function with the appropriate parameters to generate the plot.

## Prerequisites
Before running the code, make sure you have the following installed:
·Python 3.x
·NumPy
·multiprocessing (standard library)

## Using Pre-generated Data
Alternatively, if you prefer not to generate the networks yourself and run k-core percolation simulation due to the long runtime, pre-generated network and k-core percolation data is available in the Mendeley dataset. 
Follow the steps below to download the data:
1. Access the Mendeley dataset at [10.17632/jkvk97nfjc.1].
2. Download the network data files corresponding to your desired parameters.
3. Extract the downloaded files to a directory of your choice.

## Loading Pre-generated Data
To load pre-generated network data and k-core percolation simulation data into your Python environment, you can use the load function from the kcorePercolation module. 
Here's an example of how to do it:
from utils import kcorePercolation as kp
# Example: Loading pre-generated network data
network_data = kp.load('../figure1/network/NetID0_avgk10_zeta7_spatialNet.pkl')
Replace '../figure1/network/NetID0_avgk10_zeta7_spatialNet.pkl' with the actual path to the downloaded network data file.

License
This project is licensed under the MIT License - feel free to use, modify, and distribute the code as needed.
If you use this code in your research, please cite the following papers:
Nucleation phenomena and extreme vulnerability of spatial k-core systems
https://doi.org/10.48550/arXiv.2311.13579
@misc{xue2023nucleation,
      title={Nucleation phenomena and extreme vulnerability of spatial k-core systems}, 
      author={Leyang Xue and Shengling Gao and Lazaros K. Gallos and Orr Levy and Bnaya Gross and Zengru Di and Shlomo Havlin},
      year={2023},
      eprint={2311.13579},
      archivePrefix={arXiv},
      primaryClass={physics.soc-ph}
}


