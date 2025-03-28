# Modeling resource consumption in the US air transportation system via minimum-cost percolation

This repository contains the code for
- Minsuk Kim, C. Tyler Diggans, and Filippo Radicchi, [Modeling resource consumption in the US air transportation system via minimum-cost percolation](TBD_URL), Journal XXXX.
- [Preprint (arXiv)](TBD_URL) (not there yet)

- BibTex entry: To be added
    <!-- ```
    @article{kim2024shortest,
    title={Shortest-Path Percolation on Random Networks},
    author={Kim, Minsuk and Radicchi, Filippo},
    journal={Physical Review Letters},
    volume={133},
    number={4},
    pages={047402},
    year={2024},
    publisher={APS}
    }
    ``` -->

# How to retrieve the publicly available data
Our work rely on publicly available data from various sources:
- [Airline On-Time Performance Data (flight schedule data)](https://www.transtats.bts.gov/Tables.asp?QO_VQ=EFD&QO_anzr=Nv4yv0r%FDb0-gvzr%FDcr4s14zn0pr%FDQn6n&QO_fu146_anzr=b0-gvzr)
- [DB1B Data (sold tickets data)](https://www.transtats.bts.gov/Tables.asp?QO_VQ=EFI&QO_anzr=Nv4yv0r%FDb4vtv0%FDn0q%FDQr56v0n6v10%FDf748rB%FD%FLQOEO%FM&QO_fu146_anzr=b4vtv0%FDn0q%FDQr56v0n6v10%FDf748rB)
- [FAA Aircraft Database](https://registry.faa.gov/aircraftinquiry)
- [OpenFlights Database](https://openflights.org/data.php)
- [Census Population Data](https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-apdens-wpp-2015-r11-4.11)

# Installation and Requirements

## Environment Setup
There are two ways to set up the environment:

### Option 1: Using environment.yml (Recommended)
```bash
# Create and activate conda environment with all dependencies
conda env create -f environment.yml
conda activate mcp-us-airline

# Install the project package
pip install -e .
```

### Option 2: Manual Installation
```bash
# Create and activate a new conda environment
conda create -n mcp-us-airline python=3.10
conda activate mcp-us-airline

# Install the project package with all dependencies
pip install -e .
```

<!-- ## Required Python Packages
The following main packages are required (all dependencies are handled by conda):
- Python 3.10
- numpy >= 1.26
- pandas >= 2.2
- matplotlib >= 3.8
- seaborn >= 0.13
- networkx >= 3.2
- geopandas >= 0.14
- folium >= 0.15
- jupyter >= 1.0
- ipython >= 8.21
- scikit-learn >= 1.4
- pulp >= 2.7
- scipy >= 1.12 -->

<!-- ## Install Project Library
```
# From the root directory of the project
pip install -e .
``` -->

<!-- I'm not sure if this is necessary though... -->
## System Requirements 
- Python 3.10 or higher
- C compiler (gcc recommended) for MCP model compilation
- At least 32GB RAM recommended for large dataset processing

# How to generate the input files for the minimum-cost-percolation (MCP) model
The raw data retrieved from public sources should be preprocessed before running the MCP model. We provide a Jupyter notebook file `generate-input-files.ipynb` which does the job. To properly run this notebook, you need to install this library `mcp_us_airline` from this repository.


# How to run the MCP model?
We provide a `makefile` that compiles the source code of the MCP model, which is written in `C language`. You can compile it by typing:

```bash
make
```

To execute the model:
```bash
./min_cost_percolation.out <DEMAND> <FLIGHTS_INFO> <FCN> <OUTPUT_FILE> <COST>
```

Required arguments:
- `<DEMAND>`: File containing origin-destination demand data
- `<FLIGHTS_INFO>`: File containing flight details for the flight-connection-network (FCN)
- `<FCN>`: Adjacency list representation of the FCN
- `<OUTPUT_FILE>`: Name of the file where results will be saved
- `<COST>`: Cost function to minimize (integer 1-3)
  - 1: Length
  - 2: Duration
  - 3: Seat-availability-cost


Example usage:
```bash
./min_cost_percolation.out demand.txt flights.txt network.txt results.txt 1
```


# Misc.
We also provide Jupyter notebook files which allow you to generate the figures in the main text.

<!-- # Future README Improvements

## Priority Improvements
- [x] Add direct links to data sources once available
- [ ] Include paper citation details after publication
- [ ] Add example input/output file formats
- [ ] Document expected runtime and memory requirements

## Additional Enhancements
- [ ] Data Section
    - Add data years used in the study
    - Specify required data format
    - Include sample data snippets
    - Add data size estimates

- [ ] Installation Section
    - List dependencies and versions
    - Add mcp_us_airline library installation steps
    - Specify minimum system requirements
    - Add troubleshooting tips

- [ ] Documentation
    - Add brief methodology description
    - Include example output visualization
    - Document all output file formats
    - Add comments about parallel processing capabilities

- [ ] Usage Examples
    - Add complete workflow example
    - Include sample scripts
    - Document common use cases
    - Add performance optimization tips

- [ ] Optional Additions
    - Contributing guidelines
    - License information
    - Contact information
    - Known limitations
    - Acknowledgments section
    - Related publications/projects


 -->
