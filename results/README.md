# README.md `results`

The results folder might contain any kind of output that could be generated from your processed (or maybe raw) data, including:

- dataframes or tables
- statistical models
- figures

There are many ways to organise your work, and you should choose a way that makes sense for you and your project. A common principle is to have a different subfolder for each kind of analysis you do, or by the date the analysis was run, or on the sample set that was run - or some combination of these. 

As an example for this template, we provide a folder structure that reflects one way to organise data covering two kinds of analysis: "conservation" and "positive selection." These analyses are each performed more than once, on different dates, and the results are recorded separately, as shown below.

```bash
.
├── README.md
├── conservation
│   ├── 2022-05-14
│   │   └── README.md
│   ├── 2022-05-17
│   │   └── README.md
│   └── README.md
└── positive_selection
    ├── 2022-04-08
    │   └── README.md
    ├── 2022-05-01
    │   └── README.md
    └── README.md
```