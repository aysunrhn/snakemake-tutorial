import pandas as pd
from sklearn.metrics import mean_squared_error
import os

# Load test data
X_test = pd.read_csv(snakemake.input[2])  # Features for testing
y_test = pd.read_csv(snakemake.input[3]).values.ravel()  # Target values for testing

# Evaluate models
results = {}
for model_file in snakemake.input.models:
    model_name = os.path.basename(model_file).split(".")[0]
    coefficients = pd.read_csv(model_file)  # Load model coefficients
    
    # Ensure coefficients are in the correct shape
    coef = coefficients["Coefficient"].values  # Extract as 1D array

    # Make predictions
    y_pred = X_test @ coef  # Perform the dot product between X_test and coefficients

    # Calculate Mean Squared Error (MSE)
    mse = mean_squared_error(y_test, y_pred)
    results[model_name] = mse

# Save evaluation results
pd.DataFrame.from_dict(results, orient="index", columns=["MSE"]).reset_index().rename(columns={"index": "Model"}).to_csv(snakemake.output[0], index=False)
