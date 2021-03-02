import sympy
from sympy import* 
import numpy as np
from numpy import*
import math
from scipy.stats import*
from scipy.optimize import*
from scipy.optimize import fsolve

u = 78     #estar entre 50 a 250
CV = 0.44    #estar entre 0.1 a 0.6
A = 567     #estar entre 100 a 1000
h = 1.34       #estar entre 0.1 a 1.5
b1 = 20
L = 4
s=CV*u

alfa= 0.01 #entre 0 a 0.5
beta= 0.05 #entre 0 a 1
t=1
parada = 1e-3

R1= []
Q1= []

C_R1=[]
C_Q1=[]

C_R2=[]
C_Q2=[]

C2_RQ=[]
C2_QR=[]

def x (r,q):
    return ((r+q-u*L)/(((s**2)*L)**0.5))
def H (x):
    return ((0.5)*((x**2+1)*(1-norm.cdf(x))-x*norm.pdf(x)))
def G (x):
    return (norm.pdf(x) -x*(1-norm.cdf(x)))
#Punto Inicial
Q_0 = ((2*u*A)/h)**0.5
R_k = lambda R: ((h)+(h+b1)*((((s**2)*L)**0.5)/Q_0)*(G((R+Q_0-u*L)/(L*s**2)**0.5)-G((R-u*L)/(L*s**2)**0.5)))
R_0 = fsolve(R_k, Q_0)[0]
#Agregar los valores iniciales a las listas
R1.append(R_0)
Q1.append(Q_0)

#FUNCION C Y SUS DERIVADAS
def C (R,Q):
    return (h*(R+(Q/2)-u*L)+(h+b1)*(((s**2)*L)/Q)*(H(x(R,0))-H(x(R,Q)))+(A*u/Q))
def C_R (R,Q):
    return (h+(h+b1)*((((s**2)*L)**0.5)/Q)*(G(x(R,Q))-G(x(R,0))))
def C_Q (R,Q):
    return ((h/2)-(A*u/(Q**2))-(h+b1)*((((s**2)*L)/(Q**2)))*(H(x(R,0))-H(x(R,Q))-(Q/(((s**2)*L)**0.5))*G(x(R,Q))))
def C_RR (R,Q):
    return (((h+b1)/Q)*(norm.cdf(x(R,Q))-norm.cdf(x(R,0))))
def C_QQ (R,Q):
    return (((2*A*u)/(Q**3))-(h+b1)*(((s**2)*L)/(Q**3))*(H(x(R,0))-H(x(R,Q))+(((s**2)*L)**0.5)*(G(x(R,Q)))+((Q**2)/(((s**2)*L)**0.5))*(norm.cdf(x(R,Q))-1)))
def C_RQ (R,Q):
    return ((-1)*(h+b1)*((((s**2)*L)*0.5)/(Q**2))*(G(x(R,Q))-G(x(R,0))-(Q/(((s**2)*L)*0.5))*(norm.cdf(x(R,Q))-1)))
def C_QR (R,Q):
    return ((-1)*(h+b1)*((((s**2)*L)**0.5)/(Q**2))*(G(x(R,Q))-G(x(R,0))-(Q/(((s**2)*L)**0.5))*(norm.cdf(x(R,Q))-1)))

#Agregar los valores iniciales a las listas
C_R1.append(C_R(R1[0],Q1[0]))
C_Q1.append(C_Q(R1[0],Q1[0]))
C_R2.append(C_RR(R1[0],Q1[0]))
C_Q2.append(C_QQ(R1[0],Q1[0]))
C2_RQ.append(C_RQ(R1[0],Q1[0]))
C2_QR.append(C_QR(R1[0],Q1[0]))

l_deter2=[]
l_criterio=[]

deter2 = ((1)/((C_R2[-1])*(C_Q2[-1])-(C2_QR[-1])*(C2_RQ[-1])))
criterio = (deter2)*(C_R1[-1]*(C_R1[-1]*C_Q2[-1] - C_Q1[-1]*C2_RQ[-1]) + C_Q1[-1]*(C_Q1[-1]*C_R2[-1] - C_R1[-1]*C2_QR[-1]))

l_deter2.append(deter2)
l_criterio.append(criterio)
print (criterio)

n=0
j=0
while ((l_criterio[-1]**2)/2) > parada and n<6:
    t=1
    #Calculo de delta X
    determinante = ((1)/((C_R2[-1])*(C_Q2[-1])-(C2_QR[-1])*(C2_RQ[-1])))
    d_R = ((-1)*(determinante)*((C_Q2[-1])*(C_R1[-1]) - (C2_QR[-1])*(C_Q1[-1])))
    d_Q = ((-1)*(determinante)*(C_R2[-1]*C_Q1[-1] - C2_RQ[-1]*C_R1[-1]))
    #Calculo del nuevo punto
    R = R1[-1] + t*d_R
    Q = Q1[-1] + t*d_Q
    #Calculo de la funcion C con nuevo punto
    C1 = C(R,Q)
    #Calculo de la primera derivada
    Derivada_R = C_R(R1[-1],Q1[-1])
    Derivada_Q = C_Q(R1[-1],Q1[-1])
    #Calculo de la funcion a comparar
    F = (C(R1[-1],Q1[-1]) + alfa*t*((Derivada_R*d_R) + (Derivada_Q*d_Q)))
    #Contador
    n=n+1
    while C1 > F:
        t=t*beta
        R = R1[-1] + t*d_R
        Q = Q1[-1] + t*d_Q
        C1 = C(R,Q)
        F = (C(R1[-1],Q1[-1]) + alfa*t*((Derivada_R*d_R) + (Derivada_Q*d_Q)))
        j=j+1

    R = R1[-1] + t*d_R
    Q = Q1[-1] + t*d_Q
    R1.append(R)
    Q1.append(Q)
    C_R1.append(C_R(R,Q))
    C_Q1.append(C_Q(R,Q))
    C_R2.append(C_RR(R,Q))
    C_Q2.append(C_QQ(R,Q))
    C2_RQ.append(C_RQ(R,Q))
    C2_QR.append(C_QR(R,Q))

    deter2 = ((1)/((C_R2[-1])*(C_Q2[-1])-(C2_QR[-1])*(C2_RQ[-1])))
    criterio = (deter2)*(C_R1[-1]*(C_R1[-1]*C_Q2[-1] - C_Q1[-1]*C2_RQ[-1]) + C_Q1[-1]*(C_Q1[-1]*C_R2[-1] - C_R1[-1]*C2_QR[-1]))
    l_deter2.append(deter2)
    l_criterio.append(criterio)

#print ('valor de R:', R1[-1], 'valor de Q:', Q1[-1])
#print ('valor de C:', C(R1[-1],Q1[-1]))
#print ('Con una cantidad de iteraciones de:', n)
#print (j)
#print (l_criterio[-1])
#print ('----------------------------')
print ('R es', R1)
print ('----------------------------')
print ('Q es', Q1)
#print ('----------------------------')
#print ('C_R es', C_R1)
#print ('C_Q es', C_Q1)
#print ('----------------------------')
#print ('C_RR es', C_R2)
#print ('C_QQ es', C_Q2)
#print ('----------------------------')
#print ('C_RQ es', C2_RQ)
#print ('C_QR es', C2_QR)


