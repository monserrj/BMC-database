# README.md `BMC database`

Space designed to create the Bacterial microcompartment database. This database is being designed with the objective of facilitating access to BMC protein information for research purposes. Proteins include both structural proteins (e.g.: hexamers) and proteins with other functions that are part of the bacterial microcompartment, such as a chaperone or signature enzyme

## ER Diagram

There is an ER diagram of the current design at the link below:

- [https://dbdiagram.io/d/DB_STRUCTURE_052025-682494635b2fc4582f9453e4]

## TABLE OF CONTENT

- .github : useful github instructions
- data : Data generated (raw) or downloaded for database development
- docs : Important documents such as training material
- scripts : code in development to create database
- .gitignore
- _config.yml
- LICENSE
- README.md

## INSTALLATION INSTRUCTIONS

Python environment (Miniforge / conda-forge)

This repository uses a Python 3.11 environment with dependencies from the community-maintained conda-forge channel.

Create the environment:

conda env create -f conda_envs/BMC_database_313_forge.yml

Activate the environment:

conda activate BMC_database_313_forge

Environment migration (October 2026)

The original Windows and macOS Anaconda environment specifications were replaced with a conda-forge-only configuration.

The replacement environment uses updated scientific and SQL-related dependencies and is not an exact reproduction of the original environments.

The historical specifications remain accessible through Git history.

## USAGE INSTRUCTIONS

To be added when database is further developed

## LICENCE (UNSURE YET, left template one for reference)

This repository is licensed under the Creative Commons Attribution-Sharealike 4.0 licence. Please see the [`LICENSE` file]([./LICENSE](https://github.com/sipbs-compbiol/template_bioinformatics_project/blob/master/LICENSE)) in this repository for further detail.
