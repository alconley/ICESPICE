import pandas as pd

# Load the CSV file and skip the initial lines containing metadata
file_path = '~/Downloads/comsol_output_1_2mm_grid.csv'
column_names = ['X', 'Y', 'Z', 'BX', 'BY', 'BZ', 'BMOD/HMOD']
data = pd.read_csv(file_path, comment='%', skiprows=9, names=column_names)

# convert x y and y from millimeters to meters
data['X'] = data['X'] / 1000
data['Y'] = data['Y'] / 1000
data['Z'] = data['Z'] / 1000

# First, create a copy of the dataframe with the specified column transformations
formatted_data = data.copy()
formatted_data = formatted_data.fillna(0) # had to add otherwise it would not work on ubunutu lol but mac it would smh 

# Convert coordinates from millimeters to meters and format them in scientific notation
formatted_data['X'] = formatted_data['X'].apply(lambda x: f"{x:.12E}")
formatted_data['Y'] = formatted_data['Y'].apply(lambda y: f"{y:.12E}")
formatted_data['Z'] = formatted_data['Z'].apply(lambda z: f"{z:.12E}")
formatted_data['BX'] = formatted_data['BX'].apply(lambda bx: f"{bx:.12E}")
formatted_data['BY'] = formatted_data['BY'].apply(lambda by: f"{by:.12E}")
formatted_data['BZ'] = formatted_data['BZ'].apply(lambda bz: f"{bz:.12E}")
formatted_data['BMOD/HMOD'] = formatted_data['BMOD/HMOD'].apply(lambda ratio: f"{ratio:.12E}")

unique_x_values = formatted_data['X'].nunique()
unique_y_values = formatted_data['Y'].nunique()
unique_z_values = formatted_data['Z'].nunique()

# look through formatted data and if there are any values that are NAN, replace them with 0

# Prepare header and footer
header_info = f"""
          {unique_x_values}          {unique_y_values}         {unique_z_values}
 1 X                                                                              
 2 Y                                                                              
 3 Z                                                                              
 4 BX                                                                             
 5 BY                                                                             
 6 BZ                                                                             
 7 BMOD/HMOD                                                                       
00 [METRE] """ 
 
# Write to file
output_file_path = './MiniOrange3D.TABLE'
with open(output_file_path, 'w') as file:
    file.write(header_info + '\n')
    formatted_data.to_csv(file, header=False, index=False, sep=' ', mode='a')



# Original file
# Load the CSV file and skip the initial lines containing metadata
# file_path = './MiniOrange3D_Backup.TABLE'
# column_names = ['X', 'Y', 'Z', 'BX', 'BY', 'BZ', 'BMOD/HMOD']
# data = pd.read_csv(file_path, comment='%', skiprows=10, names=column_names, delim_whitespace=True)
# print(data)

# # Calculate the number of unique values in the 'X' column
# unique_x_values = data['X'].nunique()
# unique_y_values = data['Y'].nunique()
# unique_z_values = data['Z'].nunique()

# # Display the number of unique 'X' values
# print(f"Number of unique X values: {unique_x_values}")
# print(f"Number of unique Y values: {unique_y_values}")
# print(f"Number of unique Z values: {unique_z_values}")