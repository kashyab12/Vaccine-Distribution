# importing the required modules
import matplotlib.pyplot as plt
import numpy as np


class AgeCompartmentalModel:
    def __init__(self, contactMatrix, time, populations, infectionRates, recoveryRates, dt):
        self.contactMatrix = contactMatrix
        self.time = time
        self.populations = populations
        self.infectionRates = infectionRates
        self.recoveryRates = recoveryRates
        self.dt = dt
        self.ageBrackets = []
        self.t = np.arange(0, time, dt)
        index = 0
        for contactVector in contactMatrix:
            self.ageBrackets.append(AgeBracket(self.ageBrackets, populations[index], infectionRates[index], incubationRates[index], recoveryRates[index], contactVector))
            index += 1

        for increment in self.t:
            if increment == 0:
                continue
            for ageBracket in self.ageBrackets:
                ageBracket.stepIncrement(self.dt, 0.0001)
            for ageBracket in self.ageBrackets:
                ageBracket.updateIncrement()

    # Graphs each age bracket's set of diff eqs in a separate window.
    def graphIndividualModels(self):
        for ageBracket in self.ageBrackets:
            plt.xlabel("Time(days)")
            plt.ylabel("Population (million)")
            plt.plot(self.t, ageBracket.pastS, label="Susceptible", color="blue")
            plt.plot(self.t, ageBracket.pastE, label="Exposed", color="orange")
            plt.plot(self.t, ageBracket.pastI, label="Infected", color="red")
            plt.plot(self.t, ageBracket.pastR, label="Recovered", color="gray")
            plt.plot(self.t, ageBracket.pastV, label="Vaccinated", color="green")

            plt.title("SEIRV Model")
            plt.show()

# TODO finish this method which totals each ageBracket's seperate equations into one for the total population
    def graphCumulativeModels(self):
        index = 0;

# TODO It may be better to just have S,E,I,R,V as a 2D array to make the code cleaner. This is not a high priority
# Each age bracket has a set of diff eqs that are used to model the actual trajectories of the disease
class AgeBracket:
    def __init__(self, brackets, population, infectionRate, incubationRate, recoveryRate, contactVector):
        self.brackets = brackets
        self.population = population
        self.infectionRate = infectionRate
        self.incubationRate = incubationRate
        self.recoveryRate = recoveryRate
        self.contactVector = contactVector
        self.newS, self.newE, self.newI, self.newR, self.newV = 0, 0, 0, 0, 0
        self.S, self.E, self.I, self.R, self.V = self.population, 0, 1 / self.population, 0, 0
        self.pastS, self.pastE, self.pastI, self.pastR, self.pastV = [], [], [], [], []
        self.pastS.append(self.S)
        self.pastE.append(self.E)
        self.pastI.append(self.I)
        self.pastR.append(self.R)
        self.pastV.append(self.V)

    # This attempts to apply the equation in figure 2 on page 1158
    def getInfectionRisk(self):
        sum, index = 0, 0
        for bracket in self.brackets:
            sum += self.contactVector[index] * bracket.I / bracket.population
        return (1 / len(self.brackets)) * sum * self.S / self.population * self.infectionRate

# TODO This is no incrementing correctly. Must find bug
    # Increments the set of diff eqs using Euler's method. Ref figure 1 on page 1157
    def stepIncrement(self, dt, dV):
        ir = self.getInfectionRisk()
        self.newS = self.S + ((-ir) * (self.S - dV) - dV) * dt
        self.newE = self.E + ((-self.incubationRate) * self.E + ir * (self.S - dV)) * dt
        self.newI = self.I + ((-self.recoveryRate) * self.I + ir * self.incubationRate) * dt
        self.newR = self.R + (self.recoveryRate * self.I) * dt
        self.newV = self.V + dV * dt

    # updateIncrement is needed since each set of diff eqs depends on the rest of the age brackets.
    # This means that the actual increment can be performed after all the diff eqs have their new values calculated
    def updateIncrement(self):
        self.pastS.append(self.S)
        self.pastE.append(self.E)
        self.pastI.append(self.I)
        self.pastR.append(self.R)
        self.pastV.append(self.V)
        self.S = self.newS
        self.E = self.newE
        self.I = self.newI
        self.R = self.newR
        self.V = self.newV


# TODO this is data from the h1n1 influenza breakout. We will need to replace this with COVID-19 data
contactMatrix = np.array([[8.27, 1.395, 4.165, 1.51, 0.715],
                          [1.395, 5.65, 2.83, 1.83, 0.895],
                          [4.165, 2.385, 6.55, 3.425, 1.383],
                          [1.51, 1.83, 3.425, 4.2, 2.055],
                          [0.715, 0.895, 1.383, 2.055, 2.66]], dtype=float)

incubationRates = [0.25, 0.25, 0.25, 0.25, 0.25]
recoveryRates = [0.334, 0.334, 0.334, 0.334, 0.334]
infectionRates = [0.434, 0.158, 0.118, 0.046, 0.046]
populations = [0.94, 0.91, 2.30, 1.86, 0.85]
time = 365
dt = time / 1000
acm = AgeCompartmentalModel(contactMatrix, time, populations, infectionRates, recoveryRates, dt)
acm.graphIndividualModels()