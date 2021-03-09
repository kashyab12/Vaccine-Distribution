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
        index = 0
        for contactVector in contactMatrix:
            self.ageBrackets.append(AgeBracket(self.ageBrackets, populations[index], infectionRates[index], incubationRates[index], recoveryRates[index], contactVector))
            index += 1

        day = 0
        while day < time:
            for ageBracket in self.ageBrackets:
                ageBracket.stepIncrement(self.dt, 0.0001)
            day += 1


class AgeBracket:
    def __init__(self, brackets, population, infectionRate, incubationRate, recoveryRate, contactVector):
        self.brackets = brackets
        self.population = population
        self.infectionRate = infectionRate
        self.incubationRate = incubationRate
        self.recoveryRate = recoveryRate
        self.contactVector = contactVector
        self.S, self.E, self.I, self.R, self.V = self.population, 0, 1 / self.population, 0, 0

    def getInfectionRisk(self):
        sum, index = 0, 0
        for bracket in self.brackets:
            sum += self.contactVector[index] * bracket.I / bracket.population
        return (1 / len(self.brackets)) * sum * self.S / self.population * self.infectionRate

    def stepIncrement(self, dt, dV):
        ir = self.getInfectionRisk()
        newS = self.S + ((-ir) * (self.S - dV) - dV) * dt
        newE = self.E + ((-self.incubationRate) * self.E + ir * (self.S - dV)) * dt
        newI = self.I + ((-self.recoveryRate) * self.I + ir * self.incubationRate) * dt
        newR = self.R + (self.recoveryRate * self.I) * dt
        newV = self.V + dV * dt
        self.S = newS
        self.E = newE
        self.I = newI
        self.R = newR
        self.V = newV


contactMatrix = np.array([[8.27, 1.395, 4.165, 1.51, 0.715],
                          [1.395, 5.65, 2.83, 1.83, 0.895],
                          [4.165, 2.385, 6.55, 3.425, 1.383],
                          [1.51, 1.83, 3.425, 4.2, 2.055],
                          [0.715, 0.895, 1.383, 2.055, 2.66]], dtype=float)

incubationRates = [0.25, 0.25, 0.25, 0.25, 0.25]
recoveryRates = [0.334, 0.334, 0.334, 0.334, 0.334]
infectionRates = [0.434, 0.158, 0.118, 0.046, 0.046]
populations = [0.94, 0.91, 2.30, 1.86, 0.85]
time = 100
dt = time / 1000
acm = AgeCompartmentalModel(contactMatrix, time, populations, infectionRates, recoveryRates, dt)
print(acm.ageBrackets[0].I)