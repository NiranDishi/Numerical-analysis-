def MatUcreation(A):  # Create matrix U and L and send them back
    L = copy_matrix(A)
    size = len(A)
    global MatL
    for i in range(len(A)):
        for j in range(len(A)):
            if j == i:
                L[i][j] = 1
            else:
                L[i][j] = 0
    Yehida = copy_matrix(L)

    for pivot in range(size):
        for i in range(pivot + 1, size):
            if (A[pivot][pivot] == 0):
                print(A)
                exit()
            Yehida[i][pivot] = float(A[i][pivot] / A[pivot][pivot] * (-1))
            A = MatrixMultiply(Yehida, A)
            L = MatrixMultiply(L, Inversion(copy_matrix(Yehida), LU=1))
            Yehida[i][pivot] = 0

    MatL=L
    return (A)


def Inversion(A,
              LU=0):  # creates A-1 by a formula learned at class , we noticed to take care about Gauss Unstablize problem
    InverseA = copy_matrix(A)
    size = len(A)
    global MatB

    if LU != 1:
        for i in range(
                size):  # swap rows by MAX pivot for ignoring the gauss unstablize  החלפת שורות כדי לתקן את האי היציבות של גאוס עם התייחסות לוקטור הפתרון כלומר עמודה בי
            max = abs(A[i][i])
            for j in range(i, size):

                if abs(A[j][i]) > max:
                    temp = A[j]
                    A[j] = A[i]
                    A[i] = temp
                    temp = MatB[j]
                    MatB[j] = MatB[i]
                    MatB[i] = temp
                    max = abs(A[i][i])

    for i in range(len(A)):  # create a matrix with diagonal of values and other indexes are 0
        for j in range(len(A)):
            if j == i:
                InverseA[i][j] = 1
            else:
                InverseA[i][j] = 0
    Yehida = copy_matrix(InverseA)
    for pivot in range(size):
        for i in range(pivot):
            if A[pivot][pivot] == 0:
                return InverseA
            Yehida[i][pivot] = float(A[i][pivot] / A[pivot][pivot] * (-1))
            A = MatrixMultiply(Yehida, A)
            InverseA = MatrixMultiply(Yehida, InverseA)
            Yehida[i][pivot] = 0

        for i in range(pivot + 1, size):
            if A[pivot][pivot] == 0:
                return InverseA
            Yehida[i][pivot] = float(A[i][pivot] / A[pivot][pivot] * (-1))
            A = MatrixMultiply(Yehida, A)
            InverseA = MatrixMultiply(Yehida, InverseA)
            Yehida[i][pivot] = 0
    # now we got A with diagonal of not 1,lets make it 1

    ValueOfDividnes = 1
    for i in range(len(A)):
        if A[i][i] == 1:
            continue
        else:
            ValueOfDividnes = A[i][i]
            for j in range(len(A)):
                A[i][j] = A[i][j] / ValueOfDividnes
                InverseA[i][j] = InverseA[i][j] / ValueOfDividnes

    return InverseA


def MatrixMultiply(A, B): # function that calcualte matrix multiply
    if isinstance(B[0], list): # in case matrix is not a vector
        sizeA = len(A)
        newMat = []
        for i in range(sizeA):  # A for loop for row entries
            a = []
            for j in range(sizeA):  # A for loop for column entries
                a.append(0)
            newMat.append(a)

        for i in range(sizeA):
            for j in range(sizeA):
                for k in range(sizeA):
                    newMat[i][j] = newMat[i][j] + A[i][k] * B[k][j]
        return newMat
    else:# in case matrix is vector
        sizeA = len(A)
        newVect = [0] * sizeA
        for i in range(sizeA):
            for j in range(sizeA):
                newVect[i] = float(newVect[i]) + float(A[i][j]) *float(B[j])
        return newVect

def Determinant(A, total=0):
    #choose a row for calculate the determinant by each number in the row
    pivots = list(range(len(A)))

    #final condition of recursion , if 2X2 mat so calculate the determinant by the formula
    if len(A) == 2 and len(A[0]) == 2:
        detAmount = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        return detAmount

    # define submatrix for focus column 
    for fc in pivots:  # for each focus column
        # find the submatrix 
        As = copy_matrix(A)  # make a hard copy of original mat
        As = As[1:]  #  remove  first row
        height = len(As)

        for i in range(height):
            As[i] = As[i][0:fc] + As[i][fc + 1:]#     remove the pivot column elements

        sign = (-1) ** (fc % 2)  # F) 
        #  pass submatrix recursively
        sub_det = Determinant(As)
        # total all returns from recursion
        total += sign * A[0][fc] * sub_det

    return total


def copy_matrix(s): # creates a copy of s mat and return it
    MatA = []
    MatSize = len(s)
    for i in range(MatSize):  # A for loop for row entries
        a = []
        for j in range(MatSize):  # A for loop for column entries
            a.append(int(s[i][j]))
        MatA.append(a)
    return MatA



MatSize=3
MatA = [[2,4,6],[8,1,2],[2,4,8]]
MatL=0
MatB= [4,2,6]



if Determinant(MatA)!=0: #chekc if singular and in that case solution is by A^-1 regular formula
    InverseMat=Inversion(MatA)

    vetcX=MatrixMultiply(InverseMat,MatB)
    print("")
    print ("The X vector solution by X=A^-1*B  method is :  ")
    for i in range(len(vetcX)):
        print("[%.4f]"%vetcX[i])
else: # solution is by L and U matrixes

    MatU=MatUcreation(MatA)
    print("MATRIX U IS:",MatU)
    print("MATRIX L IS:",MatL)
    print("THE MULTIFICTAION OF MATRIX L AND MATRIX U IS (LU METHOD) :" ,MatrixMultiply(MatL,MatU))





    

