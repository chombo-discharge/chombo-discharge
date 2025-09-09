import os
import re
import numpy as np

# Folder containing .in files
folder_path = '.'

# Regex patterns to match relevant lines
time_step_pattern = re.compile(r"Driver::Time step report -- Time step #(\d+)")
time_pattern = re.compile(r"Time\s+=\s+([0-9.eE+-]+)")
dt_pattern = re.compile(r"dt\s+=\s+([0-9.eE+-]+)")
cfl_ito_pattern = re.compile(r"CFL \(Ito\)\s+=\s+([0-9.eE+-]+)")

# Iterate over each .in file
for filename in os.listdir(folder_path):
    if filename.endswith('.in'):
        file_path = os.path.join(folder_path, filename)
        output_path = os.path.join(folder_path, filename.rsplit('.', 1)[0] + '.out')

        extracted_data = []

        with open(file_path, 'r') as f:
            time_step = None
            time = None
            dt = None
            cfl_ito = None

            for line in f:
                ts_match = time_step_pattern.search(line)
                time_match = time_pattern.search(line)
                dt_match = dt_pattern.search(line)
                cfl_ito_match = cfl_ito_pattern.search(line)

                if ts_match:
                    time_step = int(ts_match.group(1))
                    time = None
                    dt = None
                    cfl_ito = None
                elif time_step is not None and time is None and time_match:
                    time = float(time_match.group(1))
                elif time_step is not None and dt is None and dt_match:
                    dt = float(dt_match.group(1))
                elif time_step is not None and cfl_ito is None and cfl_ito_match:
                    cfl_ito = float(cfl_ito_match.group(1))

                # Collect the data when all values are found
                if (time_step is not None and time is not None and
                        dt is not None and cfl_ito is not None):
                    extracted_data.append([time_step, time, dt, cfl_ito])
                    time_step = None  # reset for the next block

        if extracted_data:
            arr = np.array(extracted_data)
            # Save as text: time_step (int), time, dt, cfl_ito (floats)
            np.savetxt(output_path, arr, fmt=['%d', '%.6e', '%.6e', '%.6e'],
                       header="time_step time dt cfl_ito", comments='')
            print(f"Extracted and saved raw data to: {output_path}")
        else:
            print(f"No data found in {filename}")
import os
import re
import numpy as np
from scipy.signal import savgol_filter

# Define the folder containing .in files
folder_path = '.'

# Regex patterns
time_step_pattern = re.compile(r"Driver::Time step report -- Time step #(\d+)")
time_pattern = re.compile(r"Time\s+=\s+([0-9.eE+-]+)")
dt_pattern = re.compile(r"dt\s+=\s+([0-9.eE+-]+)")
cfl_ito_pattern = re.compile(r"CFL \(Ito\)\s+=\s+([0-9.eE+-]+)")

# Smoothing parameters
window_length = 101  # Must be odd and <= number of samples
polyorder = 1      # Polynomial order for the filter

# Process each .in file
for filename in os.listdir(folder_path):
    if filename.endswith('.in'):
        file_path = os.path.join(folder_path, filename)
        output_path = os.path.join(folder_path, filename.rsplit('.', 1)[0] + '.out')

        extracted_data = []
        with open(file_path, 'r') as f:
            time_step = None
            time = None
            dt = None
            cfl_ito = None

            for line in f:
                ts_match = time_step_pattern.search(line)
                time_match = time_pattern.search(line)
                dt_match = dt_pattern.search(line)
                cfl_ito_match = cfl_ito_pattern.search(line)

                if ts_match:
                    time_step = int(ts_match.group(1))
                    time = None
                    dt = None
                    cfl_ito = None
                elif time_step is not None and time is None and time_match:
                    time = float(time_match.group(1))
                elif time_step is not None and dt is None and dt_match:
                    dt = float(dt_match.group(1))
                elif time_step is not None and cfl_ito is None and cfl_ito_match:
                    cfl_ito = float(cfl_ito_match.group(1))

                # Collect the data when all are found
                if (time_step is not None and time is not None and
                        dt is not None and cfl_ito is not None):
                    extracted_data.append([time_step, time, dt, cfl_ito])
                    time_step = None

        if extracted_data:
            arr = np.array(extracted_data)

            # Apply Savitzky-Golay smoothing to the third column (dt)
            if arr.shape[0] >= window_length:
                arr[:, 2] = savgol_filter(arr[:, 2], window_length=window_length, polyorder=polyorder)
                arr[:, 3] = savgol_filter(arr[:, 3], window_length=window_length, polyorder=polyorder)                
            else:
                print(f"Skipping smoothing for {filename} — not enough data points")

            # Save result: first column as integer, others as float
            np.savetxt(output_path, arr, fmt=['%d', '%.6e', '%.6e', '%.6e'],
                       header="time_step time dt cfl_ito", comments='')
            print(f"Smoothed and saved: {output_path}")
        else:
            print(f"No data found in {filename}")
