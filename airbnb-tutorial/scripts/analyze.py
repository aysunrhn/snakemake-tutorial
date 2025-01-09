import pandas as pd

# Load cleaned data
data = pd.read_csv(snakemake.input[0])

# Group data by neighborhood group
summary = data.groupby("neighbourhood_group").agg(
    avg_price=("price", "mean"),
    avg_availability=("availability_365", "mean"),
    total_reviews=("reviews_per_month", "sum")
).reset_index()

# Save summary statistics
summary.to_csv(snakemake.output[0], index=False)