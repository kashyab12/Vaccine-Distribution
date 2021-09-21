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
Dtrans = [[], [], []]
Dvulner = []
startT = 50
endT=startT+1
Q = 0.3 - 0.115
time = 150
dt = time/1000
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

for value in contactMatrix[0]:
    Dtrans[0].append([[0,1, value/sum(contactMatrix[0]) *10* Q/(endT-startT)]])
    Dtrans[1].append([[50, 51, value / sum(contactMatrix[0]) * 10 * Q / (endT - startT)]])
    Dtrans[2].append([[100, 101, value / sum(contactMatrix[0]) * 10 * Q / (endT - startT)]])

infectionRates = [0.434, 0.158, 0.118, 0.046, 0.046]

for value in infectionRates:
    Dvulner.append([[startT, endT, value/sum(infectionRates)*Q*10/(endT-startT)]])

populations = [0.94, 0.91, 2.3, 1.86, 0.85]
#acmVuln = ageCompModel.AgeCompartmentalModel(contactMatrix, time, populations, infectionRates, recoveryRates, incubationRates, dt, Dvulner)
#acmTrans = ageCompModel.AgeCompartmentalModel(contactMatrix, time, populations, infectionRates, recoveryRates, incubationRates, dt, Dtrans)
#acmVuln.graphIndividualModels()
#acmTrans.graphCumulativeModels()
# acm.graphIndividualModels()

acm = []
for dTrans in Dtrans:
    acm.append(ageCompModel.AgeCompartmentalModel(contactMatrix, time, populations, infectionRates, recoveryRates, incubationRates, dt, dTrans))

def compareInfected(models):
    for model in models:
        plt.plot(np.arange(0, time, dt), model.cumulativeI, label="Infected", color="blue")
    plt.xlabel("Time(days)")
    plt.ylabel("Population (million)")
    plt.title("Comparing infected of different deployments")
    plt.show()


compareInfected(acm)