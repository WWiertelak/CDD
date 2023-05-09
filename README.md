# CDD: Computational Drug Discovery

CDD is a Python library for downloading compound data from the ChEMBL database, generating molecular fingerprints using the PadelPy library, and training a Random Forest model to predict pIC50 values.

## Requirements

- chembl_webresource_client
- pandas
- numpy
- padelpy
- scikit-learn
- tqdm

## Installation

Install the required libraries using pip:

```
pip install chembl_webresource_client pandas numpy padelpy scikit-learn tqdm
```

## Usage

```python
from cdd import CDD
```

### Initialize the CDD object

```python
cdd = CDD()
```

### Show the list of targets from the ChEMBL database

```python
targets = cdd.show_targets("gene_name")
```

### Select a target to download data for

```python
cdd.select_target(gene="gene_name", index=0)
```

Alternatively, use the ChEMBL ID:

```python
cdd.select_target(chembl_id="CHEMBLXXXX")
```

### Generate molecular fingerprints for the downloaded compounds

```python
cdd.describe()
```

### Train the Random Forest model using the generated fingerprints and pIC50 values

```python
score, size = cdd.fit()
```

### Predict the pIC50 of a molecule using the trained model

```python
predicted_pIC50 = cdd.predict("CHEMBLXXXX")
```

## Functions

### `pIC50(input: pd.DataFrame) -> pd.DataFrame`

Converts IC50 values to pIC50 (-log10(IC50)) and adds a new column to the DataFrame.

### `norm_value(input: pd.DataFrame) -> pd.DataFrame`

Normalizes IC50 values by capping them at 100,000,000 and adds a new column to the DataFrame.

## Class `CDD`

### Methods

- `__init__(self)`: Initializes the CDD object.

- `show_targets(self, gene: str) -> pd.DataFrame`: Returns a list of targets from the ChEMBL database.

- `select_target(self, gene: str=None, index: int=None, chembl_id: str=None)`: Downloads a list of molecules from the ChEMBL database for the selected target, cleans the data, and loads it into the model.

- `describe(self)`: Generates molecular fingerprints for the list of molecules using PadelPy and saves the result in a .CSV file.

- `fit(self) -> List[float, int]`: Trains the Random Forest model using the molecular fingerprints and pIC50 values, and returns the model score and the size of the database in number of molecules.

- `predict(self, chembl_id: str) -> float`: Predicts the pIC50 of a molecule using the trained model and returns the predicted pIC50 value.
