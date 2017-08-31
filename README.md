# polypharmacology-db
This app is a wrapper for the Children's Tumor Foundation drug interaction database. 

The purpose of this app is to facilitate exploration of drug-target interaction databases. This app currently contains the Children's Tumor Foundation Drug-Target Database, licensed from Evotec, which summarizes activity data deposited in ChEMBL and inactivity data deposited in Pubchem.


How does PPDB work?

PPDB leverages structural information of molecules and the associated target annotations to build a drug-target map based on chemical similarity between molecules. PPDB includes drug-target interactions collated by Evotec, as well as a subset of those available in the DGIdb app. Examples of use-cases for this include:

- prediction of molecular targets for novel molecules based on structural similarity

- identification of off targets for molecules of interest

- facilitating polypharmacologic drug discovery
