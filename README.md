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
- [On-Time Performance Data (flight schedule data)](TBD_URL)
- [DB1B Data (sold tickets data)](TBD_URL)
- [SEDAC Data(population data)](TBD_URL)
- [FAA Aircraft database](TBD_URL)
- [Openflight database](TBD_URL)


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
```


# Misc.
We also provide Jupyter notebook files which allows you to generate the figures in the main text.


# Future README Improvements

## Priority Improvements
- [ ] Add direct links to data sources once available
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
