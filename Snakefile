configfile: "config.yaml"

rule all:
    input:
        "output/visualizations/price_distribution.png",
        "output/visualizations/availability_by_neighborhood.png",
        "output/visualizations/reviews_per_month.png"

rule preprocess:
    input:
        "data/AB_NYC_2019.csv"
    output:
        "output/cleaned_data.csv"
    script:
        "scripts/preprocess.py"

rule analyze:
    input:
        "output/cleaned_data.csv"
    output:
        "output/summary.csv"
    script:
        "scripts/analyze.py"

rule visualize:
    input:
        "output/summary.csv"
    output:
        "output/visualizations/price_distribution.png",
        "output/visualizations/availability_by_neighborhood.png",
        "output/visualizations/reviews_per_month.png"
    script:
        "scripts/visualize.py"