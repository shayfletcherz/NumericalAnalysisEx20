# https://github.com/shayfletcherz/NumericalAnalysisEx20.git

# EX20

from datetime import datetime

def LUdecomposition(matrix, b):
    inversL, L = upperMatrix(matrix)  # calculate matrix and her inverse
    print("L = ")
    printMatrix(L)
    U = multMat(inversL, matrix)  #calculate upper matrix  (U)
    print("U = ")
    printMatrix(U)
    inversU = FromUpperToInvers(U, createIdentityMatrix(len(U)))  #calculate upper matrix inverse
    x = multMat(inversL, b)  #resulte vector
    x = multMat(inversU, x)  #resulte vector
    print("X = ")
    printMatrixSolution(x)  #showing solution X

# function for inverse matrix
def FromUpperToInvers(A, inverseMat):
    elemntarMat = createIdentityMatrix(len(A))
    for i in range(len(A) - 1, -1, -1):
        for j in range(i):  # running on the row upper the pivot
            elemntarMat[j][i] = -(A[j][i] / A[i][i])
            A = multMat(elemntarMat, A)
            inverseMat = multMat(elemntarMat, inverseMat)
            elemntarMat[j][i] = 0
        if A[i][i] != 1:  # turning pivot to 1
            elemntarMat[i][i] = 1 / A[i][i]
            A = multMat(elemntarMat, A)
            inverseMat = multMat(elemntarMat, inverseMat)
            elemntarMat[i][i] = 1
    return inverseMat


def upperMatrix(matrix):
    identityMat1 = createIdentityMatrix(len(matrix))
    identityMat2 = createIdentityMatrix(len(matrix))
    for i in range(len(matrix)):
        if matrix[i][i] == 0:  # check if the pivot is 0
            for j in range(i + 1, len(matrix)):
                if matrix[j][i] != 0:  # check if the numbers under the pivot are 0
                    identityMat1 = multMat((identityMat1, linesSwaps(matrix, i, j)))
                    identityMat2 = multMat((linesSwaps(matrix, i, j)), identityMat2)  # swap and multiply rows
                    matrix = multMat((linesSwaps(matrix, i, j)), matrix)
                    break
        if matrix[i][i] != 0:
            for j in range(i + 1, len(matrix)):
                identity = createIdentityMatrix(len(matrix))
                identity[j][i] = -(matrix[j][i] / matrix[i][i])
                identityMat2 = multMat(identity, identityMat2)  # inverse matrix
                matrix = multMat(identity, matrix)
                identity[j][i] *= -1  # Changing values for finding matrix
                identityMat1 = multMat(identityMat1, identity)
    return identityMat2, identityMat1

# the function will return an elementary matrix after swapping the rows
def linesSwaps(matrixA, line1, line2):
    idendityMax = createIdentityMatrix(len(matrixA))

    #changes in row 1
    temp = idendityMax[line1][line1]
    idendityMax[line1][line1] = idendityMax[line2][line1]
    idendityMax[line2][line1] = temp

    #changes in row 2
    temp = idendityMax[line2][line2]
    idendityMax[line2][line2] = idendityMax[line1][line2]
    idendityMax[line1][line2] = temp
    return idendityMax

# create Identity Matrix
def createIdentityMatrix(size):
    identityMat = newMat(size, size)
    for index in range(size):  # checking the diagonally
        identityMat[index][index] = 1
    return identityMat

# multiply matrix
def multMat(matrixA, matrixB):
    if len(matrixA[1]) == len(matrixB):
        C = newMat(len(matrixA), len(matrixB[0]))
        for i in range(len(C)):
            for j in range(len(C[1])):
                for k in range(len(matrixB)):
                    C[i][j] += matrixA[i][k] * matrixB[k][j]
        return C
    else:
        return None  # can not multiply

def iterativGuassSeidel(matrix, b, epsilon, checkDominence):
    print("Iterative Guass Seidel:\n")
    converge = False  # will return false if the matrix dose not convergence
    x = newMat(len(b), 1)
    print("The results are:\nThe first guess is: ")
    printMatrix(x)
    for k in range(99):  # work maximum until 99 iterations
        oldX1 = x[0][0]
        for i in range(len(x)):
            if matrix[i][i] == 0:  # for not dividing by 0
                return None
            temp = b[i][0] / matrix[i][i]
            for j in range(len(x)):  # calculating new values in new iteration
                if i != j:
                    temp -= (matrix[i][j] * x[j][0]) / matrix[i][i]
            x[i][0] = temp  # update value
        print("The result of the " + str(k + 1) + " iteration is: ")
        printMatrix(x)
        if abs(oldX1 - x[0][0]) < epsilon:  # check terms of epsilon
            converge = True
            break
    if converge is True:
        print("\n")
        if checkDominence is False:
            print("\n")
        print("The final result is: x = ")
        printMatrixSolution(x)
        print("There are : " + str(k + 1) + " iterations")
    else:
        print("Linear equations does not converge.")


# printing the solution matrix
def printMatrixSolution(matrix):
    day = str(datetime.today().day)
    hour = str(datetime.now().hour)
    minute = str(datetime.now().minute)

    print("\n")
    for i in range(len(matrix)):
        x = str(matrix[i])
        x = x.replace("]","")
        if i is len(matrix) - 1:
            print(" " + x + "00000" + day + hour + minute + "]" + "]")  # add the wanted addition to the solution
        elif i == 0:
            print("[" + x + "00000" + day + hour + minute + "]")
        else:
            print(" " + x + "00000" + day + hour + minute + "]")
    print("\n")


# printing matrix
def printMatrix(matrix):
    print("\n")
    for i in range(len(matrix)):
        if i is len(matrix) - 1:
            print(" " + str(matrix[i]) + "]")
        elif i == 0:
            print("[" + str(matrix[i]))
        else:
            print(" " + str(matrix[i]))
    print("\n")

# create zeros matrix
def newMat(numRow, numCol):
    mat = []
    for i in range(numRow):
        mat.append([])  # new row
        for j in range(numCol):
            mat[i].append(0)
    return mat

# function to check if the diadonal is dominant
def dominantDiagonal(matrix):
    for i in range(len(matrix)):
        lineSum = 0
        for j in range(len(matrix)):
            if i != j:
                lineSum += matrix[i][j]
        if matrix[i][i] <= lineSum:
            return False
    print("There is a dominant diagonal")
    return True

# function to copy matrix so the origin would not change
def copyMat(matrix):
    matrixB = newMat(len(matrix), len(matrix[0]))  # create a zero matrix of the same size as A
    for i in range(len(matrix)):  # copy A
        for j in range(len(matrix[0])):
            matrixB[i][j] = matrix[i][j]
    return matrixB

# function to change the matrix to be with dominant diagonal
def createDominantDiagonal(matrix, b):
    max = 0
    maxIndex = 0
    for i in range((len(matrix))):  # counting the highest value on the diagonal
        sum = 0
        for j in range(len(matrix)):  # running on row
            sum += abs(matrix[i][j])  # summing up all the values on the row
            if abs(matrix[i][j]) > max:  # checking for highest value
                max = abs(matrix[i][j])
                maxIndex = j
        if (sum - max) <= max:  # check if the highest value fit for the dominant diagonal
            matrix = manualSwapCol(matrix, maxIndex,  i)  # swap rows of the current diagonal
        else:  # if not search for the highest value on the column
            max = 0
            maxIndex = 0
            for j in range(len(matrix)):
                if abs(matrix[j][i]) > max:
                    max = abs(matrix[j][i])
                    maxIndex = j
            if rowSum(matrix[j]) - max <= max:  # check if the highest value fit for the dominant diagonal
                matrix, b = manualSwapRow(matrix, b, i, maxIndex)
            else:
                print("No dominant diagonal")
                return None, None
    return matrix, b

# function will return matrix after swapping rows
def manualSwapRow(matrix, b, row1, row2):
    if row2 < len(matrix) and row1 < len(matrix):
        temp = matrix[row1]
        matrix[row1] = matrix[row2]
        matrix[row2] = temp
        if b is not None:  # checking if the vector is not empty and then swapping
            temp = b[row1]
            b[row1] = b[row2]
            b[row2] = temp
    return matrix, b

# function will return matrix after swapping columns
def manualSwapCol(matrix, col1, col2):
    if col2 < len(matrix) and col1 < len(matrix):
        for i in range(len(matrix)):
            temp = matrix[i][col1]
            matrix[i][col1] = matrix[i][col2]
            matrix[i][col2] = temp
    return matrix

# returning all the values as Absolute values
def rowSum(line):
    lineSum = 0
    for index in range(len(line)):  # run over all the line`s members
        lineSum += abs(line[index])
    return lineSum


def driver():
    A = [[5, 1, 10],
         [10, 8, 1],
         [4, 10, -5]]

    b = [[1.5],
         [-7],
         [2]]

    epsilon = 0.0001

    domDiagonal = dominantDiagonal(A)
    if domDiagonal is False: # checking for dominant diagonal
        copyA = copyMat(A)
        copyB = copyMat(b)
        copyA, copyB = createDominantDiagonal(copyA, copyB)
        if (copyA is not None) and (copyB is not None): # checking that the matrix are not empty
            A = copyA
            b = copyB
            domDiagonal = True

    LUdecomposition(A, b)
    iterativGuassSeidel(A, b, epsilon, domDiagonal)


driver()