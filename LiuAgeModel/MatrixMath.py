import matplotlib.pyplot as plt
import numpy as np

# This function helps Compute a Symmetric Matrix given a Square Matrix.
# csvDataFile - CSV Data File where the Original Matrix is parsed from.
def getSymmetricContactMatrix(csvDataFile):
        
    matrixDataFile = open(csvDataFile)
    numpy_matrix = np.loadtxt(matrixDataFile, delimiter = ",", skiprows = 1)
        
    transposeNumpyMatrix = np.transpose(numpy_matrix)

    matrixSum = np.add(numpy_matrix, transposeNumpyMatrix)
    symmetricMatrix = matrixSum / 2

    return symmetricMatrix;