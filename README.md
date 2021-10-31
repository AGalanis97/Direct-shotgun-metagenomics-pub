# Direct shotgun metagenomics captures species abundance of honey samples
This repository contains all files and code to reproduce the results in the publication [add link].
<p>When running the code from scripts, there is no need to change working directories - all paths are relative to the R project file located in this repository.</p>

## Folders
### 1. Direct-shotgun-metagenomics-pub/Figures/Mocks_analysis/Data_fig_2/
This folder contains the results from analysing simulated mock samples. The software used is in the filename and the file contains both the known abundance of the mock and the one reported by the software. For simplicity we have also offered the true abundance in standalone files:
mock_genus_plants.csv: reports the abundance of the plant-only mock at the genus level
mock_others.csv: reports the abundance of the mock that contains other Eukaryotic organisms than plants
mock_community.csv: reports the abundance of the community mock

### 2. Direct-shotgun-metagenomics-pub/Simulated_mocks/kraken2_threshold_evaluation/
This folder contains the results from analysed the additional mocks that were created to evaluate the confidence threshold for kraken2. For each condition 3 mocks were created (A,B,C). The filename contains information about the mocks:
110: These are the mocks with a fragment size of 110 bp.
200: These are the mocks with a fragment size of 200 bp.
23: These are the mocks with 23 species.
64: These are the mocks with 64 species.
1M: These are the mocks with 1 million reads.
8M: These are the mocks with 8 million reads.

Each folder contains the results from the kraken2 threshold evaluation, while in the top level there is the known abundance of each one of the samples.

### 3. Direct-shotgun-metagenomics-pub/Figures
This folder contains subfolders and scripts named after the figure they reproduce. Required data exists in each folder.
