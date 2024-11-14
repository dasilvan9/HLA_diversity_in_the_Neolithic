# Figure 2: HLA frequencies

This subdirectory contains the script for generating the plots for Figure 2 of the manuscript.

## How to Run

To generate the plots, you can either:
- Run the script from the command line:
  ```bash
  python plot_hla_frequency.py
  ```
- Or load the script into an IDE (e.g., Spyder) for step-by-step execution and exploration.

## Input Files

The script requires the following input file:

- **HLA Genotype and Frequency Data** (`../data/hla_calls.xlsx`): This file includes two sheets:
  - **Sheet "genotypes"**: Contains HLA genotypes from our study, from Immel et al. (2021), and from the 1000 Genomes Project. The 1000 Genomes HLA data can be accessed [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt).
  - **Sheet "allelefrequencies.net"**: Provides allele frequency data for modern Germans, available on [Allele Frequencies Net](https://www.allelefrequencies.net/pop6001c.asp?pop_id=3767).
