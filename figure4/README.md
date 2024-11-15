# Figure 4: HLA diversity analysis

This subdirectory contains the script for generating the plots for Figure 4 of the manuscript.

## How to Run

To generate the plots, you can either:
- Run the script from the command line:
  ```bash
  python plot_hla_diversity.py
  ```
- Or load the script into an IDE (e.g., Spyder) for step-by-step execution and exploration.

## Input Files

The script requires the following input file:

- **HLA Genotype and Frequency Data** (`../data/hla_calls.xlsx`): This file includes two sheets:
  - **Sheet "genotypes"**: Contains HLA genotypes from our study, from [Immel et al. (2021)](https://www.nature.com/articles/s42003-020-01627-4), and from the [1000 Genomes Project](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt).
  - **Sheet "allelefrequencies.net"**: Provides allele frequency data for modern Germans, available on [Allele Frequencies Net](https://www.allelefrequencies.net/pop6001c.asp?pop_id=3767).

## Package requirements

This script was tested with `scikit-bio` version 0.5.8. If you encounter any issues, please verify that you are using this specific version. 
