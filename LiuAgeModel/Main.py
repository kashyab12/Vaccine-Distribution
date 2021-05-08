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
vaccineDistributions = []
Dtrans = []
Dvulner = []
T = 50
contactMatrix = [[8.27, 1.395, 4.165, 1.51, 0.715],
                 [1.395, 5.65, 2.385, 1.83, 0.895],
                 [4.165, 2.385, 6.55, 3.425, 1.383],
                 [1.51, 1.83, 3.425, 4.2, 2.055],
                 [0.715, 0.895, 1.383, 2.055, 2.66]]

# Incubation Rate for COVID-19 - 5.2
# Recovery Rates for COVID-19 - 0.11

for ctrVariable in range(5):
    incubationRates.append(0.25)
    recoveryRates.append(0.334)
    vaccineDistributions.append([[T,200,1/150]])

for value in contactMatrix[0]:
    Dtrans.append([[T,200, value/sum(contactMatrix[0])/100]])
    print(value, " ", value/sum(contactMatrix[0]))

infectionRates = [0.434, 0.158, 0.118, 0.046, 0.046]

for value in infectionRates:
    Dvulner.append([[T, 200, value/sum(infectionRates)/100]])

populations = [0.94, 0.91, 2.3, 1.86, 0.85]
time = 364
dt = time / 1000
acm = ageCompModel.AgeCompartmentalModel(contactMatrix, time, populations, infectionRates, recoveryRates, incubationRates, dt, Dvulner)
acm.graphCumulativeModels()
# acm.graphCumulativeModels(True)
acm.graphIndividualModels()
