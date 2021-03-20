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
    infectionRates.append(0.23)

populations = [19.736, 20.212, 20.827, 20.849, 21.254, 23.277, 21.932, 21.443, 19.584, 20.345, 20.355, 21.163, 20.592, 17.356, 14.131, 21.300]
time = 364
dt = time / 1000
acm = ageCompModel.AgeCompartmentalModel(contactMatrix, time, populations, infectionRates, recoveryRates, incubationRates, dt)
acm.graphCumulativeModels()