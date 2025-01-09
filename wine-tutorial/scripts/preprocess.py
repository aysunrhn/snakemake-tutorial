import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# Load the raw dataset
data = pd.read_csv(snakemake.input[0], sep=";")

# Handle missing values (if any)
data = data.dropna()

# Separate features and target variable
X = data.iloc[:, :-1]  # Features (all columns except the last one)
y = data.iloc[:, -1]   # Target variable (last column)

# Normalize the feature data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Convert scaled data back to a DataFrame
X_scaled = pd.DataFrame(X_scaled, columns=X.columns)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# Save the processed data to CSV files
X_train.to_csv(snakemake.output[0], index=False)
X_test.to_csv(snakemake.output[1], index=False)
y_train.to_csv(snakemake.output[2], index=False)
y_test.to_csv(snakemake.output[3], index=False)