# importing the required modules
import matplotlib.pyplot as plt
import numpy as np
import AgeCompartmentalModel as ageCompModel
import MatrixMath as matMath

# Macros
CSV_DATA_FILE_NAME = "us-contact-matrix.csv"

# COVID-19 Contact Matrix - PLOS Computational Biology Study
incubationRates = []
recoveryRates = []
infectionRates = []
contactMatrix = matMath.getSymmetricContactMatrix(CSV_DATA_FILE_NAME)

# Incubation Rate for COVID-19 - 5.2
# Recovery Rates for COVID-19 - 0.11

for ctrVariable in range(16):
    incubationRates.append(5.2)
    recoveryRates.append(0.11)
    infectionRates.append(0.5)

populations = [19736, 20212, 20827, 20849, 21254, 23277, 21932, 21443, 19584, 20345, 20355, 21163, 20592, 17356, 14131, 21300]
time = 150
dt = time / 1000
acm = ageCompModel.AgeCompartmentalModel(contactMatrix, time, populations, infectionRates, recoveryRates, incubationRates, dt)
acm.graphCumulativeModels()