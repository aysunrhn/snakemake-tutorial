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
Snakemake essentially builds on the implicit wildcard rule approach of Make, and it extends its capabilities by allowing the use of Python in a pipeline. Just like Make, its goal is to produce a set of requested output files based on predefined rules and steps. You can also think of it as cooking or baking.

Although it was originally developed to create scalable bioinformatics and genomics pipelines, it can be generalized to other applications as well. It has become a crucial tool in reproducible research; being cited more than 12 times per week in 2023, it has been used extensively in scientific publications in several different fields. Currently, it has over one million downloads on Conda.

### Noteworthy features of Snakemake
- If Python and Make were to have a baby.
- You can describe workflows using a human readable, Python based language.
- It has built-in caching: if some steps of your workflow have already been run, Snakemake can recognize that and avoid rerunning the same analyses.
- It can accommodate both serial and parallel jobs since each "work units" in a workflow can be run independently of one another.
- It makes debugging easier since it keeps track of all files generated, you can identify which steps in your workflow have failed.
- Integration with conda allows you to define conda environments for both the whole workflow and individual steps. 
- You can incorporate tools or methods written in different scripting languages.
- There are cloud integrations established to run Snakemake workflows seamlessly in different cloud platforms. 
- Workflows can be scaled to server, cluster, grid and cloud environments, without modifying the workflow itself.
- It has an active user and developer base in mainly bioinformatics, and other scientific research community. It's development is driven by the needs of scientists and their needs of reproducible research.

### What about Airflow?
Airflow has become the go-to tool for building data pipelines, it has started in Airbnb and was later adopted by the Apache project. Being backed by the Apache project, and a growing community of contributors, it's the most popular workflow manager in software engineering.

Pros:
- It is the most popular choice of workflow management system in software engineering.
- It has a bigger community support and hence more detailed tutorials and documentation.
- It is richer in features, especially enterprise integrations for industry.
- It integrates well with databases, APIs, cloud services and data warehouse like Google BigQuery, AWS S3 and Snowflake.
- It is more production friendly with built-in scheduling and monitoring features.

Cons:
- It has a complex infrastructure and a steeper learning curve.
- It is more difficult to share pipelines between different environments.
- It lacks native support for conda.
- It has a more convoluted approach to parallel computation.
- It has limited support for HPC.

### What is special about bioinformatics workflows?
In bioinformatics, a workflow is a collection of steps run in series to transform raw data input to processed results (figures, insights, decisions etc.). Each step can be made up of different tools, programs, parameters, databases and dependencies/requirements. 

There are 3 main reasons why bioinformatics workflows are different than those data engineering workflows in the industry:
1. Differences in data type, shape and scale
   - Bioinformatics datasets are typically very large and come from various sources (DNA sequences, RNA sequencing, proteomics data, imagining data)
   - Different data types have different structures and different preprocessing/analysis needs
2. Differences in programs and tooling
   - Bioinformatics pipelines often involve numerous steps with intricate dependencies
   - Different steps use different, specialized tools 
   - Open source, highly specialized tools that are not meant to be integrated natively, there are no software packages or standalone platforms to run analyses from start to finish 
   - New algorithms, tools and reference databases are updated frequently, researchers need flexible pipelines that can be adapted easily
   - Many bioinformatics tasks and tools are computationally expensive (genome assembly, alignment, sequence search) and require HPC
3. Community support behind bioinformatics workflow managers and open source software
   - Field of bioinformatics has a strong emphasis on open science, reproducible and transparent research. All of which are achieved using workflow managers 

### A summary of workflow management tools discussed

| **Feature**                  | **Snakemake**                  | **Nextflow**                 | **Airflow**                  |
|------------------------------|-------------------------------|-----------------------------|------------------------------|
| **Target Audience**           | Scientists, researchers        | Bioinformaticians, scientists | Data engineers, DevOps        |
| **Ease of Use**               | High                          | Medium                      | Low                          |
| **Reproducibility Focus**     | Strong                        | Medium                      | Low                          |
| **HPC Support**               | Excellent                     | Good                        | Minimal                      |
| **Cloud-Native Support**      | Moderate                      | Strong                      | Excellent                    |
| **Community**                 | Scientific community           | Bioinformatics              | Data engineering, cloud-native |



> [!TIP] Reproducible research 
> Snakemake is a good option to introduce researchers to. Airflow produces well defined production level workflows that are meant to be run continually, and hence Snakemake is much better suited for the needs of researchers who want to run reproducible analyses for a particular project or an appplication. 

### Other fields that already use or can benefit from Snakemake
- [Large, parallel deep learning experiments using Snakemake](https://waterdata.usgs.gov/blog/snakemake-for-ml-experiments): The USGS demonstrates how Snakemake can be used to perform deep learning experiments on environmental data
- [Workflow managers in high-energy physics - Enhancing analyses with Snakemake: Physics and high-energy particle research:](https://archive.fosdem.org/2024/events/attachments/fosdem-2024-3415-workflow-managers-in-high-energy-physics-enhancing-analyses-with-snakemake/slides/22450/FOSDEM-HEP-workflow-managers_hcQO9j3.pdf) Large scale HEP experiments with complex dependencies on cloud
- [Neuroscience:](https://www.nature.com/articles/s41598-024-77615-z) Using Snakemake to process and analyze MRI data- 
- [Ecology and evolutionary research:](https://doi.org/10.1111/2041-210X.14113) Article outlining the specific benefits of using Snakemake to manage data analysis workflows in ecology and evolutionary research, and demonstrating its use through case studies
- Earth and Climate Science: Preprocessing large datasets from climate models (NetCDF format) for regional climate projections
- Automating preprocessing EEG/MEG data, MRI and fMRI data analysis and predictive modeling for thousands of neuroimaging files
- Chemistry: Running large-scale molecular dynamics simulations with different parameter combinations and collecting results in a reproducible framework


# Useful links
- [Awesome pipeline repository](https://github.com/pditommaso/awesome-pipeline): a curated list of tools for creating pipelines
- [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/): an extensive documentation that includes in depth tutorials and edge cases
- [An Introduction to Snakemake with R for Economics](https://lachlandeer.github.io/snakemake-econ-r-tutorial/index.html)
- [Sustainable data analysis with Snakemake](https://doi.org/10.12688/f1000research.29032.2): A general overview article about Snakemake published in F1000Research
- [A review of bioinformatic pipeline frameworks.](https://doi.org/10.1093/bib/bbw020): A review article published in Briefings in Bionformatics, focused on Bionformatics pipelines
- [Airflow vs Snakemake](https://learn.flowdeploy.com/blog/airflow-vs-snakemake): A direct comparison between Airflow and Snakemake for data pipelines, while promoting FlowDeploy.
- [Why are bioinformatics workflows different?](https://www.bsiranosian.com/bioinformatics/why-are-bioinformatics-workflows-different/)