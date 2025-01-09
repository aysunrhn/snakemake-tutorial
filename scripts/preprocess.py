import pandas as pd

# Load raw data
data = pd.read_csv(snakemake.input[0])

# Drop rows with missing values in critical columns
data = data.dropna(subset=["name", "host_name", "neighbourhood_group", "price"])

# Filter out listings with unrealistic prices (e.g., over $1,000)
data = data[data["price"] <= 1000]

# Normalize column names
data.columns = [col.strip().lower().replace(" ", "_") for col in data.columns]

# Save cleaned data
data.to_csv(snakemake.output[0], index=False)