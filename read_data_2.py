#%%
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  

# Define the path to your test folder
folder_path = 'test'

# Read eigenvalues file and compute frequency in Hz
csv = pd.read_csv(os.path.join(folder_path, 'eigenvalues.txt'), sep='\s+', skiprows=1)[1:]
csv['FREQUENCY_HZ'] = np.float16(csv['RADIANS']) / (2 * np.pi)

# Get all mode files in the folder (excluding eigenvalues file)
files = [f for f in os.listdir(folder_path) 
         if os.path.isfile(os.path.join(folder_path, f)) 
         and not f.startswith('eigen')]

# Sort files to ensure mode1, mode2, ... order
files.sort()

# Only keep the first two mode files
files = files[:2]
output_folder = os.path.join(folder_path, "plots")
os.makedirs(output_folder, exist_ok=True)
# Increase font size globally
plt.rcParams.update({'font.size': plt.rcParams['font.size'] + 6})
frequency = [0.52,0.76]
# Iterate over each mode file (only first two)
for i, file in enumerate(files):
    file_path = os.path.join(folder_path, file)

    # Read the mode shape file into a DataFrame
    df = pd.read_csv(file_path, sep='\s+', skipinitialspace=True)

    # Extract the mode number from the file name
    mode_number = int(file.replace('mode', '').replace('.txt', ''))

    # Get the corresponding frequency in Hz
    frequency_hz = frequency[i]

    # Plot mode shape components
    plt.figure(figsize=(10, 6))
    plt.plot(df['T1'], label='T1')
    plt.plot(df['T2'], label='T2')
    plt.plot(df['T3'], label='T3')
    plt.plot(df['R1'], label='R1')
    plt.plot(df['R2'], label='R2')
    plt.plot(df['R3'], label='R3')

    plt.title(f'Mode Shape - Mode {mode_number} (Frequency: {frequency_hz:.2f} Hz)')
    plt.xlabel('GRID ID')
    plt.ylabel('Norm. Displacement')
    plt.legend()
    plt.grid(True)
    # Save the figure
    save_path = os.path.join(output_folder, f'mode_{mode_number}.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.tight_layout()
    plt.show()