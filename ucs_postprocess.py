# UCS Post-Processing Script
# Reads wall force and position data from the simulation output and
# computes the stress-strain curve for the uniaxial compressive strength test.
#
# Output file: uniStressStrain.dat
#   Column 1: strain (%)
#   Column 2: stress (units of 1e-5)

from math import *

# ----------------------------- Configuration ---------------------------------
# Output file for stress-strain data
output_file = "uniStressStrain.dat"

# Input data files produced by the ESyS-Particle simulation
force_file    = "data/out_wallForce.dat"
position_file = "data/out_wallPosition.dat"

# Sample cross-sectional area (pixels): 65 x 65
area = 65 * 65

# Sampling interval (read every 'inter'-th line from the data files)
inter = 10
# ----------------------------- Read data files --------------------------------
try:
    with open(force_file, 'r') as fh1:
        fh1data = fh1.readlines()
except IOError:
    raise IOError("Could not open wall force file: {}".format(force_file))

try:
    with open(position_file, 'r') as fh2:
        fh2data = fh2.readlines()
except IOError:
    raise IOError("Could not open wall position file: {}".format(position_file))

# ----------------------------- Initialise geometry ----------------------------
# Determine the initial sample height from the first position record.
# Each line contains two wall positions: top wall (index 1) and bottom wall (index 4).
k2 = fh2data[0].split()
l1 = float(k2[1])   # top-wall Y position
l2 = float(k2[4])   # bottom-wall Y position
height = l1 - l2    # initial sample height (mm)

N = len(fh1data)    # total number of data records

# ----------------------------- Compute stress-strain --------------------------
# Open the output file in write mode
with open(output_file, 'w') as k0_file:
    for i in range(0, N, inter):
        k1 = fh1data[i].split()
        k2 = fh2data[i].split()

        # Wall forces: f1 = bottom wall, f2 = top wall (Y-components)
        f1 = float(k1[1])
        f2 = float(k1[4])

        # Stress: average axial force divided by cross-sectional area (units of 1e-5)
        stress = 0.5 * (f2 - f1) / area

        # Current wall positions
        l1 = float(k2[1])
        l2 = float(k2[4])

        # Strain: relative change in sample height (dimensionless)
        strain = (height - (l1 - l2)) / height

        # Write strain (%) and stress (1e-5 units) to match experimental data format
        info = str(strain * 100) + '\t' + str(stress / 1e5) + '\n'
        k0_file.write(info)

        # Stop once axial strain exceeds 1 %
        if strain * 100 > 1:
            break
