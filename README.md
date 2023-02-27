# Drug-Target Explorer and Companion Database

This repo contains the code for an interactive web interface that was developed using the R Shiny platform and cheminformatics R packages. The app enables the end user with a specific query molecule to search a database of experimentally-derived drug-target interactions. The database can be queried using drug names or structures. A similarity parameter allows the user to expand their search to other structurally related molecules. The structure of the query molecule is compared to every database molecule, and structurally similar molecules and targets are presented in interactive tabular and network-based forms for in-depth exploration.

The app also performs enrichment analysis on the target lists, and allows the user to evaluate structure-activity relationships for drug response data using publicly-available in vitro datasets. If a user has a target of interest, they can search for molecules that bind that target and explore the resulting data interactively. 

Potential use cases include:
- prediction of molecular targets for novel molecules based on structural similarity

- identification of off targets for molecules of interest

- facilitating polypharmacologic drug discovery

To learn more about this project and access the underlying database, please go to www.synapse.org/dtexplorer. 


## Deployment to ShinyApps

- Enable workflows in the GitHub repository
- Under [secrets](https://github.com/Sage-Bionetworks/polypharmacology-db/settings/secrets/actions) click 'New repository secret'
- Enter secrets for `RSCONNECT_USER`, `RSCONNECT_TOKEN`, and `RSCONNECT_SECRET`, the values for which are saved in Sage's LastPass.
- Trigger the GitHub action.
- Check out the app here: https://sagebio.shinyapps.io/polypharmacology-db-staging.
- After verifying correctness, create a Git branch named release*, e.g., `release-1.0`.
- The application will become available at https://sagebio.shinyapps.io/polypharmacology-db
