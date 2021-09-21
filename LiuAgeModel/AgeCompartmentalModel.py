# importing the required modules
import matplotlib.pyplot as plt
import numpy as np


# TODO It may be better to just have S,E,I,R,V as a 2D array to make the code cleaner. This is not a high priority
# Each age bracket has a set of diff eqs that are used to model the actual trajectories of the disease
class AgeBracket:
    def __init__(self, brackets, population, infectionRate, incubationRate, recoveryRate, contactVector, dVDistribution):
        self.brackets = brackets
        self.population = population
        self.infectionRate = infectionRate
        self.incubationRate = incubationRate
        self.recoveryRate = recoveryRate
        self.contactVector = contactVector
        self.dVDistribution = dVDistribution
        self.newS, self.newE, self.newI, self.newR, self.newV = 0, 0, 0, 0, 0
        self.S, self.E, self.I, self.R, self.V = self.population, 0, 0.0000001, 0, 0
        self.Sv, self.Ev, self.Iv, self.Rv, self.Vv = 0, 0, 0, 0, 0
        self.newSv, self.newEv, self.newIv, self.newRv, self.newV = 0, 0, 0, 0, 0
        self.pastS, self.pastE, self.pastI, self.pastR, self.pastV = [], [], [], [], []
        self.pastS.append(self.S)
        self.pastE.append(self.E)
        self.pastI.append(self.I)
        self.pastR.append(self.R)
        self.pastV.append(self.V)

    def getdV(self, day, dt):
        index = 0
        dv = 0
        for vacCoord in self.dVDistribution:
            if vacCoord[0] < day < vacCoord[1]:
                dv += vacCoord[2] * dt
        return dv


    def getInfectionRisk(self):
        sum, index = 0, 0
        for bracket in self.brackets:
            sum += self.contactVector[index] * bracket.I / bracket.population
        return (1 / len(self.brackets)) * sum * self.S / self.population * self.infectionRate

    # Increments the set of diff eqs using Euler's method. Ref figure 1 on page 1157
    def stepIncrement(self, dt, day):
        ir = self.getInfectionRisk()
        self.newS = self.S + ((-ir) * (self.S - self.getdV(day, dt)) - self.getdV(day, dt)) * dt
        self.newE = self.E + ((-self.incubationRate) * self.E + ir * (self.S - self.getdV(day, dt))) * dt
        self.newI = self.I + ((-self.recoveryRate) * self.I + self.incubationRate * self.E) * dt
        self.newR = self.R + (self.recoveryRate * self.I) * dt
        self.newV = self.V + self.getdV(day, dt) * dt

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


class AgeCompartmentalModel:
    def __init__(self, contactMatrix, time, populations, infectionRates, recoveryRates, incubationRates, dt, dvDistributions):
        self.contactMatrix = contactMatrix
        self.time = time
        self.populations = populations
        self.infectionRates = infectionRates
        self.recoveryRates = recoveryRates
        self.incubationRates = incubationRates
        self.dt = dt
        self.ageBrackets = []
        self.dvDistributions = dvDistributions
        self.t = np.arange(0, time, dt)
        index = 0
        self.totalpopulation = 0;

        # Initialize Brackets
        for contactVector in contactMatrix:
            self.ageBrackets.append(AgeBracket(self.ageBrackets, self.populations[index], self.infectionRates[index], self.incubationRates[index], self.recoveryRates[index], contactVector, self.dvDistributions[index]))
            self.totalpopulation += populations[index]
            index += 1

        # Run Euler's Method on Model
        day = 0
        for increment in self.t:
            if increment == 0:
                continue
            for ageBracket in self.ageBrackets:
                ageBracket.stepIncrement(self.dt, day)
            for ageBracket in self.ageBrackets:
                ageBracket.updateIncrement()
            day += self.dt

        index = 0;
        self.cumulativeS, self.cumulativeE, self.cumulativeI, self.cumulativeR, self.cumulativeV = np.empty(
            self.t.__len__()), np.empty(self.t.__len__()), np.empty(self.t.__len__()), np.empty(
            self.t.__len__()), np.empty(self.t.__len__())
        while index < self.t.__len__():
            for ageBracket in self.ageBrackets:
                self.cumulativeS[index] += ageBracket.pastS[index]
                self.cumulativeE[index] += ageBracket.pastE[index]
                self.cumulativeI[index] += ageBracket.pastI[index]
                self.cumulativeR[index] += ageBracket.pastR[index]
                self.cumulativeV[index] += ageBracket.pastV[index]
            index += 1

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

    def graphCumulativeModels(self, inpercent=False):
        plt.xlabel("Time(days)")
        if inpercent:
            plt.ylabel("Population (%)")
            plt.plot(self.t, self.cumulativeS/self.totalpopulation, label="Susceptible", color="blue")
            plt.plot(self.t, self.cumulativeE/self.totalpopulation, label="Exposed", color="orange")
            plt.plot(self.t, self.cumulativeI/self.totalpopulation, label="Infected", color="red")
            plt.plot(self.t, self.cumulativeR/self.totalpopulation, label="Recovered", color="gray")
            plt.plot(self.t, self.cumulativeV/self.totalpopulation, label="Vaccinated", color="green")
        else:
            plt.ylabel("Population (millions)")
            plt.plot(self.t, self.cumulativeS, label="Susceptible", color="blue")
            plt.plot(self.t, self.cumulativeE, label="Exposed", color="orange")
            plt.plot(self.t, self.cumulativeI, label="Infected", color="red")
            plt.plot(self.t, self.cumulativeR, label="Recovered", color="gray")
            plt.plot(self.t, self.cumulativeV, label="Vaccinated", color="green")
        plt.legend(loc="best")
        plt.title("Cumulative")
        plt.savefig("result.jpg")
        plt.show()