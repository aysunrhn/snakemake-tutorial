import pandas as pd

# Load preprocessed training data
X_train = pd.read_csv(snakemake.input[0])

# Dynamically decide which models to train
selected_models = []
if X_train.shape[0] > 1000:  # Example condition: large dataset
    selected_models.extend(["Ridge", "Lasso"])
else:  # Small dataset
    selected_models.append("LinearRegression")

# Save selected models to file
with open(snakemake.output[0], "w") as f:
    for model in selected_models:
        f.write(model + "\n")