# Snakemake 101

## What is a workflow manager?

A workflow manager is a tool to develop computational workflows that are reproducible by aiding in setting up, performing and monitoring defined sequences of commands. It is a structured approach to keep track of the commands used to analyze data, including parameters.

Bioinformatics field relies heavily on pipelines to analyze genomic data; usually you perform a sequence of analyses on several files. Hence, there are many workflow managers niche to the field. But workflow management is a general need in many other applications as well. 

### Useful definitions
- **Implicit framework:** the user defines rules or recipes for performing operations on files separately. It forces the users to think more carefully about filenames rather than about the process.
- **Explicit framework:** the user defines the rule rule topology, the order is fixed and tasks simply refer to each other rather than being a target naming scheme.
- Implicit wildcard rule: define available file transformations based on file suffixes[^leipzig2017].
- Configuration framework: pipeline framework that uses a configuration-based, rather than convention-based, approach to describe tasks. Requires a fixed XML file to describe individual job run instances and their dependencies.
- Class-based framework: class-based implementations closely bound to an existing code library, rather than executables. 

### Some examples of workflow managers
- Good old scripts
- Make
- Snakemake
- Nextflow
- Galaxy
- Jupyter notebooks

## Snakemake: a Python-esque make 
Snakemake essentially builds on the implicit wildcard rule approach of Make, and it extends its capabilities by allowing the use of Python in a pipeline.  It was developed to create scalable bioinformatics and genomics pipelines, although it can be generalized to other applications as well.

Just like Make, its goal is to produce a set of requested output files based on predefined rules and steps. You can also think of it as cooking or baking.

### Noteworthy features of Snakemake
- You can describe workflows using a human readable, Python based language.
- If some steps of your workflow have already been run, Snakemake can recognize that and avoid rerunning the same analyses.
- It can accommodate both serial and parallel jobs since each "work units" in a workflow can be run independently of one another.
- It makes debugging easier since it keeps track of all files generated, you can identify which steps in your workflow have failed.
- Integration with conda allows you to define conda environments for both the whole workflow and individual steps. 
- There are cloud integrations established to run Snakemake workflows seamlessly in different cloud platforms. 
- Workflows can be scaled to server, cluster, grid and cloud environments, without modifying the workflow itself.
- It has an active user and developer base.


## Useful links
- [Awesome pipeline repository](https://github.com/pditommaso/awesome-pipeline): a curated list of tools for creating pipelines
- [The official Snakemake documentation](https://snakemake.readthedocs.io/en/stable/): an extensive documentation that includes in depth tutorials and edge cases.
- [An Introduction to Snakemake with R for Economics](https://lachlandeer.github.io/snakemake-econ-r-tutorial/index.html)
- [Sustainable data analysis with Snakemake](https://doi.org/10.12688/f1000research.29032.2): A general overview article about Snakemake published in F1000Research
- [A review of bioinformatic pipeline frameworks.](https://doi.org/10.1093/bib/bbw020): A review article published in Briefings in Bionformatics, focused on Bionformatics pipelines


[^leipzig2017]: Leipzig, J. A review of bioinformatic pipeline frameworks. *Briefings in Bioinformatics*. 2017. 18 (3). https://doi.org/10.1093/bib/bbw020