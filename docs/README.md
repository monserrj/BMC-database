# README.md `docs`

This folder is intended to hold _documentation_ for my database repository and project. "Documentation" contents include or will include:

- A training folder from the sqlite carpentry course with Python examples
- A detailed account (or accounts) of your project's analysis workflow (to be developed)
  - essentially, this would be an online Methods section for your work

## Good practice information to always have into account

- It is useful to keep the different kinds of documentation in different subfolders, to help you navigate the project files
- Documentation file should be clearly named, so you don't have to open them to know the general contents; this will help you find relevant information again later in the project
- Plain text works better with `git` than binary files like `Word` (the version control aspect if `git` works with plain text, but not binary)
  - Use plain text files for data and document storage where possible
- There are a number of tools and resources that can make your documentation easier to write and navigate
  - All the `README.md` files in this template, including their formatting, are written in [Markdown](https://www.markdownguide.org). This is a plain text format.
  - Text written in Markdown can be converted simply into many other formats, including [the webpages for this repository](https://sipbs-compbiol.github.io/template_bioinformatics_project/), and Word documents, using tools such as [pandoc](https://pandoc.org).

### Subfolder structure

One way to organise your data is suggested in this repository:

```bash
.
├── README.md
├── methods
│   └── README.md
├── papers
│   └── README.md
└── thesis
    └── README.md
```

Separate folders are used for `methods` (a contemporary written account of the analysis), `papers`, published references that are relevant to the work, and `thesis` (the dissertation describing the work).
