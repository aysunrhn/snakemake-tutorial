import pandas as pd
import matplotlib.pyplot as plt

# Load model evaluation results
results_df = pd.read_csv(snakemake.input[0])

# Create a bar chart for model performance
plt.figure(figsize=(10, 6))
plt.bar(results_df["Model"], results_df["MSE"], color="skyblue")
plt.title("Model Performance (MSE)")
plt.xlabel("Model")
plt.ylabel("Mean Squared Error")
plt.xticks(rotation=45)
plt.tight_layout()

# Save the visualization
plt.savefig(snakemake.output[0])
plt.close()
