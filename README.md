# Summary
Multivariate Data Split - Program to create subsets of a data sets that are likely to be not statistically significantly different.

# Description
Subsets are picked based on the psuedo gower distance, where the distance between subjects picked in each iteration is minimized.


## Input

	a) phenotype file (only numerical values, i.e. please convert information such as sex (M/F) into 0/1)

	b) phenotype list to be used for comparison (see example phenotype_list.csv file)

	c) Number of subsets


## Output

CSV file with added column including subset IDs.


## Execution

./subset_data.py INFILE OUTFILE NUM_GROUPS (PHENOTYPE_LIST)


## Example

./subset_data.py phenotypes.csv phenotypes_split.csv 5 phenotype_list.csv
 
# Requires
numpy, 
scipy

# Citations
Schirmer et al. "Rich-Club organization: an important determinant of functional outcome after acute ischemic stroke." BioRxiv (2019): 545897.

