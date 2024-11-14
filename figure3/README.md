# Figure 3: local ancestry inference in the MHC region

This subdirectory contains the script for generating the plots for Figure 3 of the manuscript.

## How to Run

To generate the plots, you can either:
- Run the script from the command line:
  ```bash
  python plot_local_ancestry.py
  ```
- Or load the script into an IDE (e.g., Spyder) for step-by-step execution and exploration.

## Input Files

The script requires the following input file:

- **RFmix2 results** (`../data/Results_LF.WHG.G22.tsv.zip`): A zipped file containing the merged RFmix2 output for all chromosomes. The RFmix2 results per chromosome (`*.fb.tsv`) were merged, and genome-wide Z-scores were calculated for analysis.

- **HLA genotypes** (`../data/AdditionalFile1.xlsx`): Supplementary data from the manuscript, with HLA genotypes per sample for early farmers (EF) and late farmers (LF) in the `"DataS1"` sheet.

## RFmix2 Workflow

To rerun RFmix2, a simple Snakemake workflow is available. Use the following command, specifying the number of cores you’d like to use:
  ```bash
  snakemake --cores <number_of_cores>
  ```

