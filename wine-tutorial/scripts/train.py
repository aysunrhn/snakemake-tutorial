import pandas as pd
from sklearn.linear_model import LinearRegression, Ridge, Lasso
import os

# Load preprocessed data
X_train = pd.read_csv(snakemake.input[0])
y_train = pd.read_csv(snakemake.input[1])
# y_train = pd.read_csv(snakemake.input[1]).values.ravel()

# Determine model type from wildcard
model_type = snakemake.wildcards.model

if model_type == "LinearRegression":
    model = LinearRegression()
elif model_type == "Ridge":
    model = Ridge(alpha=1.0)
elif model_type == "Lasso":
    model = Lasso(alpha=0.1)
else:
    raise ValueError(f"Unknown model type: {model_type}")

# Train the model
model.fit(X_train, y_train)

# Save the model coefficients
coefficients = pd.DataFrame(model.coef_, columns=["Coefficient"])
coefficients.to_csv(snakemake.output[0], index=False)