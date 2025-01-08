import pandas as pd

# Load the CSV file and skip the initial lines containing metadata
# file_path = './python_scripts/new_field.csv'

# file_path = '/home/alconley/Downloads/comsol_output_5N42_1x1x16in_x50_y50_z70_res1_2mm.csv'
# file_path = '/home/alconley/Downloads/comsol_output_3N42_1x1x16in_x50_y50_z70_res1_2mm.csv'
file_path = '/Users/alconley/OneDrive - Florida State University/ICESPICE/5N42_1x1x1_8/geant4/comsol_output_5N42_1x1x1_8in_withmounts_1_2mm_grid.csv'

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
output_file_path = './comsol_output_5N42_1x1x8in_with_mounts_x50_y50_z70_res1_2mm.mag'
with open(output_file_path, 'w') as file:
    file.write(header_info + '\n')
    formatted_data.to_csv(file, header=False, index=False, sep=' ', mode='a')