import pandas as pd

# Load the CSV file and skip the initial lines containing metadata
file_path = './comsol_output.csv'
column_names = ['X', 'Y', 'Z', 'BX', 'BY', 'BZ', 'BMOD/HMOD']
data = pd.read_csv(file_path, comment='%', skiprows=9, names=column_names)

print(data)
