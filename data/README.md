# README.md `data`

This subfolder contains _data_. This is an ambiguous term and might mean many things, including:

- data you generate
- data you have downloaded from elsewhere
- "raw" data: data as it was collected from a machine
- "cleaned" data: data that has been processed for analysis, for example by correcting mistakes, removing outliers and "bad" data, or normalisation

The way you organise your data should suit the way you work, but there are principles of _good practice_ that can help avoid problems later on, and make it easier to reanalyse the data or otherwise share your work.

## Good practice

- Store the raw data by itself, and never modify the raw data directly
  - You should always have a "chain of evidence" reproducibly establishing the steps of your analysis
  - You may need to reanalyse your data, or "clean" it, in a different way in future. You cannot do this if the original data has been modified
  - Keeping the raw data in its own folder keeps it separate from other files that might be modified
- Use plain text formats wherever possible
  - Plain text will always be readable in some way; binary and proprietary files probably will not
  - `git` works better with plain text than binary files - so it is easier to see if a file has been changed (accidentally or not) if it is plain text.
- Keep the "cleaned" and processed data separate from the raw data, in its own folder
  - Data processing should - wherever possible - be handled by a script or other repeatable code/software process. This helps ensure reproducibility and establish a "chain of evidence" for your work
  - Processed data will be the basis for your analyses
- Keep reference data separate from the data you generate yourself
  - You might, for example, need to download genome sequences from a remote database, like NCBI. As the remote database is updated you might want to repeat your analyses. This will be less prone to accidental problems if you do not mix your reference and self-generated datasets.
- Results should be kept completely separate from data - not in _this_ folder
  - The general principle is that anything you generate in your analysis can be regenerated, but raw data cannot
  - Processed/cleaned data is a starting point for your analysis, and you may make many attempts at your analysis. 

### Subfolder structure

One way to organise your data is suggested in this repository:

```bash
.
├── README.md
├── cleaned
│   └── README.md
├── raw
│   └── README.md
└── reference
    └── README.md
```

Separate subfolders are used for each of the `raw`, `cleaned` and `reference` datasets. A standard variation on this is if you try more than one data processing method (e.g. normalisation vs removing outliers). In that case you might want to have further subfolders in the `cleaned` data subfolder, and give them more meaningful names (e.g. `normalised`, `no_outliers`).