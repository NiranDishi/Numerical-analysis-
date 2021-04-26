def createD(A):
    D=copy_matrix(A)
    for i in range(MatSize):
        for j in range(MatSize):
            if j != i:
                D[i][j] = 0

    global MatD
    MatD=D

def createU(A):
    U = copy_matrix(A)
    for i in range(MatSize-1,-1,-1):
        for j in range(i,-1,-1):
            U[i][j] = 0
    global MatU
    MatU = U

def createL(A):
    L = copy_matrix(A)
    for i in range(MatSize):
        for j in range(i,MatSize):
            L[i][j]=0
    global MatL
    MatL = L




def checkG(A):
    x=1



def copy_matrix(s): # creates a copy of s mat and return it

    MatSize = len(s)
    MatA = [0]*MatSize
    if isinstance(s[0],list)==False:
        for i in range(MatSize) :
            MatA[i]=s[i]
        return MatA
    for i in range(MatSize):  # A for loop for row entries
        a = []
        for j in range(MatSize):  # A for loop for column entries
            a.append(int(s[i][j]))
        MatA.append(a)
    return MatA

def yakobi(A,B):
    isDominant=dominantDiagnosal(A)
    if(isDominant==0):
        print("the Matrix isnt Mitkaneset")
        exit()
    OldVariables=[0,0,0]
    NewVariables = [0, 0, 0]


    for i in range(MatSize):
        NewVariables[i]=B[i]

        for j in range(i+1,MatSize,1):
            NewVariables[i]+=(-1)*A[i][j]*OldVariables[j]

        for k in range(0,i,1):
            NewVariables[i] += (-1) * A[i][k]*OldVariables[k]

        NewVariables[i]=NewVariables[i]/A[i][i]

    COUNTER=1
    while(abs(OldVariables[0]-NewVariables[0])>0.00001 and abs(OldVariables[1]-NewVariables[1])>0.00001 and abs(OldVariables[2]-NewVariables[2])>0.00001):
        COUNTER+=1
        OldVariables=copy_matrix(NewVariables)
        for i in range(MatSize):
            NewVariables[i] = B[i]

            for j in range(i + 1, MatSize, 1):
                NewVariables[i] += (-1) * A[i][j] * OldVariables[j]

            for k in range(0, i, 1):
                NewVariables[i] += (-1) * A[i][k] * OldVariables[k]

            NewVariables[i] = NewVariables[i] / A[i][i]
    print("After  {}  Iterations , by Yakobi Method - answers are : ".format(COUNTER))
    print(NewVariables)


def zaidel(A,B):
    isDominant = dominantDiagnosal(A)
    if (isDominant == 0):
        print("the Matrix isnt Mitkaneset")
    else:
        OldVariables = [0, 0, 0]
        NewVariables = [0, 0, 0]

        for i in range(MatSize):
            NewVariables[i] = B[i]

            for j in range(i + 1, MatSize, 1):
                NewVariables[i] += (-1) * A[i][j] * OldVariables[j]

            for k in range(0, i, 1):
                NewVariables[i] += (-1) * A[i][k] * OldVariables[k]

            NewVariables[i] = NewVariables[i] / A[i][i]


        COUNTER = 1
        while (abs(OldVariables[0] - NewVariables[0]) > 0.00001 and abs(
                OldVariables[1] - NewVariables[1]) > 0.00001 and abs(OldVariables[2] - NewVariables[2]) > 0.00001):
            COUNTER += 1
            OldVariables= copy_matrix(NewVariables)
            for i in range(MatSize):
                NewVariables[i] = B[i]

                for j in range(i + 1, MatSize, 1):
                    NewVariables[i] += (-1) * A[i][j] * NewVariables[j]

                for k in range(0, i, 1):
                    NewVariables[i] += (-1) * A[i][k] * NewVariables[k]

                NewVariables[i] = NewVariables[i] / A[i][i]
    print("After  {}  Iterations , by Zaidel Method - answers are : ".format(COUNTER))
    print(NewVariables)


def pivoting(A,B):
    for i in range(MatSize):  # swap rows by MAX pivot for ignoring the gauss unstablize  החלפת שורות כדי לתקן את האי היציבות של גאוס עם התייחסות לוקטור הפתרון כלומר עמודה בי
        max = abs(A[i][i])
        for j in range(i, MatSize):

            if abs(A[j][i]) > max:
                temp = A[j]
                A[j] = A[i]
                A[i] = temp
                temp = B[j]
                B[j] = B[i]
                B[i] = temp
                max = abs(A[i][i])
    global MatB
    MatB=B
    return A

def dominantDiagnosal(A):
    pivot=0
    sum=0
    for i in range (MatSize):
        pivot =abs(A[i][i])
        sum=0
        for j in range(MatSize):
            if j!=i:
                sum+=abs(A[i][j])
        if pivot<sum:
            print("There is No Dominant Diagnosal\n")
            return 0
    print("There is Dominant Diagnosal\n")
    return 1

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

def Inversion(A):  # creates A-1 by a formula learned at class , we noticed to take care about Gauss Unstablize problem
    InverseA = copy_matrix(A)
    size = len(A)
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



#######MAIN#######
MatSize=3
MatA = [[4,2,0],[2,10,4],[0,4,5]]
MatB= [2,6,5]
MatU =0
MatD=0
MatL=0


print("Pivoting Procsessing...")
#MatA=pivoting(MatA,MatB)




print("Wich Method You Want To Use To Solve Matrix?\n1-Yakobi Method\n2-Zaidel Method")
choise=input()
if choise=="1":
    yakobi(MatA,MatB)
else:
    zaidel(MatA,MatB)
