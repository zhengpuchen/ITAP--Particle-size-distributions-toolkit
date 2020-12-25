# psd.py
# Analyze the properties of particle size distributions.
# Usage:  python3 ./distribution.py input_file.csv

import os, sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import psd

def ReadCumulativeDistribution(filename):
    # Read an input file with a cumulative distribution.
    print("Reading the cumulative distribution file: ", filename)
    with open(filename) as csv_file:
        data = []
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            data.append([float(row[0]), float(row[1])])
        npdata = np.array(data)
    return npdata[:,0], npdata[:,1]

def ReadFrequencyDistribution(filename):
    # Read an input file with a frequency distribution.
    print("Reading the frequency distribution file: ", filename)
    with open(filename) as csv_file:
        data = []
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            data.append([float(row[0]), float(row[1])])
        npdata = np.array(data)
    return npdata[:,0], npdata[0:len(npdata[:,1])-1,1]

def ReadSieveData(filename):
    # Read an input file with sieve data
    print("Reading the sieve file: ", filename)
    with open(filename) as csv_file:
        data = []
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            data.append([float(row[0]), float(row[1]), float(row[2])])
            npdata = np.array(data)
    sieve_size = npdata[:,0]
    mass_empty = npdata[:,1]
    mass_with_material = npdata[:,2]
    mass_retained = mass_with_material - mass_empty
    mass_total = sum(mass_retained)

    mass_fraction_coarser = np.zeros([len(npdata[:,0])],dtype=float)
    mass_coarser = 0
    for i in range(len(mass_retained)):
        mass_coarser += mass_retained[i]
        mass_fraction_coarser[i] = mass_coarser/mass_total
    sieve_size = np.flip(sieve_size)
    mass_fraction_coarser = np.flip(mass_fraction_coarser)
    mass_fraction_finer = 1 - mass_fraction_coarser
    
    return sieve_size, mass_fraction_finer


# Read the command line to get the input file name.
if len(sys.argv) < 2:
    sys.exit('Usage: %s [stl file]' % sys.argv[0])
if not os.path.exists(sys.argv[1]):
    sys.exit('ERROR: File %s was not found.' % sys.argv[1])

# Read a frequency distribution input file.
x, q = ReadFrequencyDistribution(sys.argv[1])
f = psd.Distribution(x=x, q0=q)

# Read a cumulative distribution input file.
#x, Q = ReadCumulativeDistribution(sys.argv[1])
#f = psd.Distribution(x=x, Q2=Q)

# Read a sieve input file.
#x, Q = ReadSieveData(sys.argv[1])
#f = psd.Distribution(x=x, Q3=Q)

# Write output files with frequency and cumulative distribution data.
f.WriteFrequencyDistribution()
f.WriteCumulativeDistribution()

# Generate frequency distribution data.
print("f.type = ", f.type)
print("f.by = ", f.by)
print("f.x0 = ", f.x0, "units")
print("f.x1 = ", f.x1, "units")
print("f.x_avg = ", f.x_avg, "units")
print("f.q = ", f.q, "1/units")
print("f.x = ", f.x, "units")
print("f.Q = ", f.Q)
print("f.GetD(0) = ", f.GetD(0), "units")
print("f.GetD(5) = ", f.GetD(5), "units")
print("f.GetD(10) = ", f.GetD(10), "units")
print("f.GetD(50) = ", f.GetD(50), "units")
print("f.GetD(90) = ", f.GetD(90), "units")
print("f.GetD(95) = ", f.GetD(95), "units")
print("f.GetD(100) = ", f.GetD(100), "units")
print("f.GetD_Median() = ", f.GetD_Median(), "units")
print("f.GetSpan() = ", f.GetSpan())
print("f.GetD_Mode() = ", f.GetD_Mode(), "units")
print("f.GetD_ArithmeticMean() = ", f.GetD_ArithmeticMean(), "units")
print("f.GetD_QuadraticMean() = ", f.GetD_QuadraticMean(), "units")
print("f.GetD_CubicMean() = ", f.GetD_CubicMean(), "units")
print("f.GetD_HarmonicMean() = ", f.GetD_HarmonicMean(), "units")
print("f.GetD_GeometricMean() = ", f.GetD_GeometricMean(), "units")
print("f.GetSauterMeanDiameter() = ", f.GetSauterMeanDiameter(), "units")
print("f.GetDeBroukereMeanDiameter() = ", f.GetDeBrouckereMeanDiameter(), "units")
print("f.GetStandardDeviation() = ", f.GetStandardDeviation(), "units")
print("f.GetSkewness() = ", f.GetSkewness(), "units")
print("f.GetExcessKurtosis() = ", f.GetExcessKurtosis(), "units")
print("f.GetDistributionMoment(2) = ", f.GetDistributionMoment(2))
print("f.FractionInInterval(0, 200) = ", f.FractionInInterval(0, 200))

f0= psd.Distribution(x=f.x, q0=f.GetFrequencyDistributionByNumber())
f2= psd.Distribution(x=f.x, q2=f.GetFrequencyDistributionByArea())
f3= psd.Distribution(x=f.x, q3=f.GetFrequencyDistributionByVolume())

#Other ways to get the Sauter Mean Diameter and DeBroukere Mean Diameter
M40 = f0.GetDistributionMoment(4, 0)
M30 = f0.GetDistributionMoment(3, 0)
M20 = f0.GetDistributionMoment(2, 0)
SMD = M30/M20
VMD = M40/M30
print("M30/M20 = ", SMD)
print("M40/M20 = ", VMD)

Mm13 = f0.GetDistributionMoment(-1, 3)
M13 = f3.GetDistributionMoment(1, 3)
SMD = 1/Mm13
VMD = M13
print("Mm13 = ", SMD)
print("M13 = ", VMD)

print(" Done")


# Make plots
fig,axs = plt.subplots(3)
axs[0].plot(f.x, 1-f.Q, '-o', label="1-F")
axs[0].plot(f0.x, 1-f0.Q, '-x', label="1-F0")
axs[0].plot(f2.x, 1-f2.Q, '-x', label="1-F2")
axs[0].plot(f3.x, 1-f3.Q, '-*', label="1-F3")
axs[0].set_xlabel('size [units]')
axs[0].set_ylabel('mass fraction coarser [-]')
axs[0].legend()

axs[1].plot(f.x, f.Q, '-o', label="F")
axs[1].plot(f0.x, f0.Q, '-x', label="F0")
axs[1].plot(f2.x, f2.Q, '-*', label="F2")
axs[1].plot(f3.x, f3.Q, '-*', label="F3")
axs[1].set_xlabel('size [units]')
axs[1].set_ylabel('mass fraction finer [-]')
axs[1].legend()

axs[2].plot(f.x_avg[:], f.q, '-o', label="f")
axs[2].plot(f0.x_avg[:], f0.q, '-x', label="f0")
axs[2].plot(f2.x_avg[:], f2.q, '-x', label="f2")
axs[2].plot(f3.x_avg[:], f3.q, '-x', label="f3")
axs[2].set_xlabel('average size [units]')
axs[2].set_ylabel('frequency [1/units]')
axs[2].legend()
plt.show()

# To Do:
# -  Fit to common distributions:  normal, log-normal, Rosin-Rammler


print("Done.")






