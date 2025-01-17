import pandas as pd
import matplotlib.pyplot as plt

# Load model evaluation results
results_df = pd.read_csv(snakemake.input[0])

# Create a bar chart for model performance
plt.figure()
plt.bar(results_df["Model"], results_df["MAE"], color="skyblue")
plt.title("Model Performance (MAE)")
plt.xlabel("Model")
plt.ylabel("Mean Absolute Error")
plt.xticks(rotation=45)
plt.tight_layout()

# Save the visualization
plt.savefig(snakemake.output[0])
plt.close()
