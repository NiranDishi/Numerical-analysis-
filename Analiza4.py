import math
from sympy.utilities.lambdify import lambdify
import sympy as sp



def secant_method(polynom,startPoint,endPoint,epsilon):
    x = sp.symbols('x')
    xrMinusOne=startPoint
    xr=endPoint
    f=toFunc(polynom)
    f=lambdify(x,f)
    xrPlusOne=(xrMinusOne*f(xr)-xr*f(xrMinusOne))/(f(xr)-f(xrMinusOne))
    counter=1
    while(abs(xr-xrPlusOne)>epsilon):
        if counter>100:
            print("bad function")
            return 0
        xrMinusOne=xr
        xr=xrPlusOne
        xrPlusOne=(xrMinusOne*f(xr)-xr*f(xrMinusOne))/(f(xr)-f(xrMinusOne))
        counter+=1

    print("after {} iterations the root is {}".format(counter,xrPlusOne))





def Newton_Raphson(polynom,startPoint,endPoint,epsilon):
    x = sp.symbols('x')
    xr = 0
    xrPlusOne=0
    fx = 0
    fxTAG = 0
    xr=(startPoint+endPoint)/2
    fx=functionOF(polynom,xr)
    fxTAG=Difrencial(polynom)
    fxTAG = lambdify(x, fxTAG)
    fxTAG=fxTAG(xr)
    xrPlusOne=xr-(fx/fxTAG)
    eterationNumber=0
    while abs(xr-xrPlusOne)>epsilon:
        eterationNumber=eterationNumber+1
        xr=xrPlusOne
        fx = functionOF(polynom, xr)
        fxTAG = Difrencial(polynom)
        fxTAG = lambdify(x, fxTAG)
        fxTAG = fxTAG(xr)
        xrPlusOne = xr - (fx / fxTAG)
    print("after: {} iterations we found the root:  {} ".format(eterationNumber,xrPlusOne))


def toFunc(polynom):

    x = sp.symbols('x')


    my_f=0
    maala=len(polynom)-1

    for i in range(len(polynom)):
        my_f+=polynom[i]*x**(maala-i)


    return my_f

def Difrencial(polynom):

    x = sp.symbols('x')


    my_f=0
    maala=len(polynom)-1

    for i in range(len(polynom)):
        my_f+=polynom[i]*x**(maala-i)


    my_f1 = sp.diff(my_f, x)


    return my_f1



def Bisection_Method(polynom,startPoint,endPoint,epsilon):
    c=0
    a=startPoint
    b=endPoint
    global eterationNumber

    while (abs(b-a)>epsilon):
        c = (a + b) / 2
        #eterationNumber=eterationNumber+1

        if functionOF(polynom,a)*functionOF(polynom,c)>0:
            a=c
        else:
            b=c

    print("After {} iterations we Found that one of the roots is :".format(eterationNumber))
    print(a)



def functionOF(polynom,ABC):
    result = 0  # ITHUL
    polynomMaala = len(polynom)
    misparHofshi = polynom[polynomMaala - 1]

    for i in range(polynomMaala - 1):
        result = result + polynom[i] * math.pow(ABC, polynomMaala - 1)
        polynomMaala = polynomMaala - 1
    result = result + misparHofshi
    return result


######MAIN######

print("please enter the MAALA of your equivalent ")
numberOfRoots=int(input())

polynom=[0]*(numberOfRoots+1)
for i in range(numberOfRoots+1):
    if (i==numberOfRoots):
        print("please enter the Eivar Hofshii of your equivalent: ")
        polynom[i]= int(input())
    else:
        print("please enter the x^{}  of your equivalent: ".format(numberOfRoots-i))
        polynom[i]=int( input())

print(("enter the start point please :"))
startPoint=int(input())
print(("enter the end point please :"))
endPoint=int(input())
print(("enter the epsilon please : (should be 0.0001)"))
epsilon=float(input())


iterationNumber=0 #FIRST STATE
a=startPoint #HASAMA
b=endPoint #HASAMA


result=0 # ITHUL
polynomMaala=len(polynom)
misparHofshi=polynom[polynomMaala-1]

for i in range(polynomMaala-1):
    result=result+polynom[i]*math.pow(a,polynomMaala-1)
    polynomMaala=polynomMaala-1
result=result+misparHofshi

rootsArray=[0]*numberOfRoots
NextResault=0

a=a+0.1
eterationNumber=0
print("Which method do you want to use ?")
print("Click 1 For HALFING METHOD")
print("Click 2 For NUTHON RAFSON METHOD")
print("Click 3 For MEITAR METHOD")
choose=int(input())

flag=0
while a<=b:
    eterationNumber+=1
    polynomMaala = len(polynom)

    for i in range(polynomMaala - 1):
        NextResault = NextResault + polynom[i] * math.pow(a, polynomMaala - 1)
        polynomMaala = polynomMaala - 1
    polynomMaala = len(polynom)
    misparHofshi = polynom[polynomMaala - 1]
    NextResault = NextResault + misparHofshi

    if(a>-0.05 and a<0.05):
        print("After {} iterations we Found that one of the roots is :\n0".format(eterationNumber))
    if(result*NextResault<0):
        if(choose==1):
            Bisection_Method(polynom,a,a-0.1,epsilon)
        if(choose==2):
            Newton_Raphson(polynom, startPoint, endPoint, epsilon)
        if (choose == 3):
            secant_method(polynom, startPoint, endPoint, epsilon)


    if(iterationNumber+1>-1 * sp.log(epsilon / (b - a)) / sp.log(2)):
        print("Method Error")
    z=-1 * sp.log(epsilon / (b - a)) / sp.log(2)

    result=NextResault
    NextResault=0
    a=a+0.1











