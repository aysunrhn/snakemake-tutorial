# Snakemake 101

## What is a workflow manager?

A workflow manager is a tool to develop computational workflows that are reproducible by aiding in setting up, performing and monitoring defined sequences of commands. It is a structured approach to keep track of the commands used to analyze data, including parameters.

Bioinformatics field relies heavily on pipelines to analyze genomic data; usually you perform a sequence of analyses on several files. Hence, there are many workflow managers niche to the field. But workflow management is a general need in many other applications involving data to support reproducible analysis.

### Some examples of workflow managers
- Good old scripts
- Make
- Snakemake
- Nextflow
- Galaxy
- Airflow
- Jupyter notebooks

Galaxy, KNIME, and Watchdog have GUI where you can define and execute workflows, making them the easiest to pick up and use without prior coding experience.

Most of the workflow management systems are code-based; workflows in Andruil, Balsam, Hyperloom, Ruffus, SciPipe, SCOOP, COMPs and JUDI are defined using classes and functions borrowed from other programming languages, such as Python and Scala. Airflow, Nextflow, Snakemake, ClusterFlow and BigDataScript, on the other hand, have their own domain specific langauge (DSL) to specify workflows. They offer the advantage of being easier to read and understand.

## Snakemake: a Python-esque make 
Snakemake essentially builds on the implicit wildcard rule approach of Make, and it extends its capabilities by allowing the use of Python in a pipeline.  It was developed to create scalable bioinformatics and genomics pipelines, although it can be generalized to other applications as well.

Just like Make, its goal is to produce a set of requested output files based on predefined rules and steps. You can also think of it as cooking or baking.

### Noteworthy features of Snakemake
- If Python and Make were to have a baby.
- You can describe workflows using a human readable, Python based language.
- It has built-in caching: if some steps of your workflow have already been run, Snakemake can recognize that and avoid rerunning the same analyses.
- It can accommodate both serial and parallel jobs since each "work units" in a workflow can be run independently of one another.
- It makes debugging easier since it keeps track of all files generated, you can identify which steps in your workflow have failed.
- Integration with conda allows you to define conda environments for both the whole workflow and individual steps. 
- There are cloud integrations established to run Snakemake workflows seamlessly in different cloud platforms. 
- Workflows can be scaled to server, cluster, grid and cloud environments, without modifying the workflow itself.
- It has an active user and developer base.

### Why not Airflow?
Airflow has become the go-to tool for building data pipelines, it has started in Airbnb and was later adopted by the Apache project. Being backed by the Apache project, and a growing community of contributors, it's the most popular workflow manager in software engineering.

Pros:
- It is the most popular choice of workflow management system in software engineering.
- It has a bigger community support and hence more detailed tutorials and documentation.
- It is richer in features, especially enterprise integrations for industry.
- It is more production friendly.

Cons:
- It has a complex infrastructure and a steeper learning curve.
- It is more difficult to share pipelines between different environments.
- It lacks native support for conda.
- It has a more convoluted approach to parallel computation.
- It has limited support for basic cluster operations and cloud computing, except for Google Cloud and Astronomer.

> [!TIP] Reproducible research 
> Snakemake is a good option to introduce researchers to. Airflow produces well defined production level workflows that are meant to be run continually, and hence Snakemake is much better suited for the needs of researchers who want to run reproducible analyses for a particular project or an appplication. 

# Useful links
- [Awesome pipeline repository](https://github.com/pditommaso/awesome-pipeline): a curated list of tools for creating pipelines
- [The official Snakemake documentation](https://snakemake.readthedocs.io/en/stable/): an extensive documentation that includes in depth tutorials and edge cases.
- [An Introduction to Snakemake with R for Economics](https://lachlandeer.github.io/snakemake-econ-r-tutorial/index.html)
- [Sustainable data analysis with Snakemake](https://doi.org/10.12688/f1000research.29032.2): A general overview article about Snakemake published in F1000Research
- [A review of bioinformatic pipeline frameworks.](https://doi.org/10.1093/bib/bbw020): A review article published in Briefings in Bionformatics, focused on Bionformatics pipelines
- [Airflow vs Snakemake](https://learn.flowdeploy.com/blog/airflow-vs-snakemake): A direct comparison between Airflow and Snakemake for data pipelines, while promoting FlowDeploy.