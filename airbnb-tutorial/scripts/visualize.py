import pandas as pd
import matplotlib.pyplot as plt

# Load summary statistics
summary = pd.read_csv(snakemake.input[0])

# 1. Price Distribution Visualization
plt.figure(figsize=(8, 6))
plt.bar(summary["neighbourhood_group"], summary["avg_price"], color="skyblue")
plt.title("Average price by neighborhood group")
plt.xlabel("Neighborhood group")
plt.ylabel("Average price (USD)")
plt.tight_layout()
plt.savefig("output/visualizations/price_distribution.png")
plt.close()

# 2. Availability Visualization
plt.figure(figsize=(8, 6))
plt.bar(summary["neighbourhood_group"], summary["avg_availability"], color="salmon")
plt.title("Average availability by neighborhood group")
plt.xlabel("Neighborhood group")
plt.ylabel("Average availability (days per year)")
plt.tight_layout()
plt.savefig("output/visualizations/availability_by_neighborhood.png")
plt.close()

# 3. Reviews Per Month Visualization
plt.figure(figsize=(8, 6))
plt.bar(summary["neighbourhood_group"], summary["total_reviews"], color="gray")
plt.title("Total reviews by neighborhood group")
plt.xlabel("Neighborhood group")
plt.ylabel("Total reviews")
plt.tight_layout()
plt.savefig("output/visualizations/reviews_per_month.png")
plt.close()
