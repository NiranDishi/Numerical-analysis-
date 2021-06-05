def copy_matrix(s): # creates a copy of s mat and return it
    MatA = []
    MatSize = len(s)
    for i in range(MatSize):  # A for loop for row entries
        a = []
        for j in range(MatSize):  # A for loop for column entries
            a.append(int(s[i][j]))
        MatA.append(a)
    return MatA



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









def Inversion(A,
              LU=0):  # creates A-1 by a formula learned at class , we noticed to take care about Gauss Unstablize problem
    InverseA = copy_matrix(A)
    size = len(A)
    global MatB



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











def LinearMethod(List,Number):

    start=0
    end=0
    index=0
    for i in range(len(List)-1):
        if Number>List[i][1]and Number<List[i+1][1]:
            start=List[i][1]
            end=List[i+1][1]
            index=i
            break
    return (((List[i][0]-List[i+1][0])/(start-end))*Number+((start*List[i+1][0]-end*List[i][0])/(start-end)))

def PolynomalMethod(List,Number):

     Mat = []
     MatSize = len(List)

     for i in range(MatSize):
         a = []
         for j in range(MatSize):
             if(j==0):
                 a.append(1)
             else:
                 a.append(List[i][1]**j)
         Mat.append(a)

     yVector=[0]*MatSize
     for i in range(MatSize):
         yVector[i]=List[i][0]

     matrixMultiplyResult=MatrixMultiply(Inversion(Mat,0),yVector)
     sum=0
     for i in range(len(matrixMultiplyResult)):
         if i==0:
             sum=sum+matrixMultiplyResult[i]
         else:
             sum=sum+matrixMultiplyResult[i]*Number**i

     return sum


def LaGrangeMethod(List,Number):
    SUM=0

    for i in range(len(List)):
        mult = 1
        for j in range(len(List)):
            if j==i:
                continue
            else:
                mult = mult*((Number-List[j][1])/(List[i][1]-List[j][1]))
        mult=mult*List[i][0]
        SUM=SUM+mult
    return SUM

def P(start,end,i,tempdict,iters,List):


        return tempdict ['p' + str(start) + ',' + str(end)]



def  NevilleMethod(List,Number):


        temp=0
        for i in range(len(List)):
            temp=List[i][0]
            List[i][0]=List[i][1]
            List[i][1]=temp

        dict = {}
        tempdict = {}
        n = len(List)
        iters = 1
        iterationNumer=0
        while len(dict) != 1:
            for i in range(n - iters):
                if iterationNumer==0:
                    temp = ((Number - List[i][0]) * List[i + iters][1] - (Number - List[i + iters][0]) * List[i][1]) / (
                            List[i + iters][0] - List[i][0])
                    tempdict['p' + str(i) + ',' + str(i + iters)] = temp
                else:
                    temp = ((Number - List[i][0]) * P(i+1,i+iters,i,dict,iters,List) - (Number - List[i + iters][0]) * P(i,i+iters-1,i,dict,iters,List) ) / (
                            List[i + iters][0] - List[i][0])
                    tempdict['p' + str(i) + ',' + str(i + iters)] = temp

            iterationNumer=1
            dict = tempdict
            tempdict = {}
            iters += 1
        temp = list(dict.items())
        return temp[0][1]




#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]
numbersTable = [[0, 1], [0.112463, 1.2], [0.167996, 1.3],[0.222709,1.4]] #####[F(X),X]
#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]#####[F(X),X]


numberThatNotInTheTable = 1.28######### X VALUE
print("Which method you want to use?")
print("1.Linear Method")
print("2.Polynomial Method")
print("3.La Grange Method")
print("4.Neville Method")

methodChose=input()

if(methodChose=="1"):
    print("The Point {} value, is : {}".format(numberThatNotInTheTable,LinearMethod(numbersTable,numberThatNotInTheTable)) )

if(methodChose=="2"):
   print("The Point {} value, is : {}".format(numberThatNotInTheTable,PolynomalMethod(numbersTable,numberThatNotInTheTable)) )

if(methodChose=="3"):
    print("The Point {} value, is : {}".format(numberThatNotInTheTable,LaGrangeMethod(numbersTable,numberThatNotInTheTable)) )

if(methodChose=="4"):
    print("The Point {} value, is : {}".format(numberThatNotInTheTable,NevilleMethod(numbersTable,numberThatNotInTheTable)) )

