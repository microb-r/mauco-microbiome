# Mauco Microbiome 

Repository for the Maucio microbiome project. The goal is to have a centralized repository for the code and analysis of the project, to allow reproducibility and easy share of protocols, analysis, results, etc.

The structure of this repository is:

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered (not on the repo, because of its size)
	| |- process/     # cleaned data, will not be altered once created;
	|
	|- code/          # any programmatic code
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
