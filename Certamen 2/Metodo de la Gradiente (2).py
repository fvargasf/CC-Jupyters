#Método de minimización

    #Llamar las bibliotecas
import sympy
from sympy import*
import numpy as np
from scipy.stats import*
from scipy.optimize import*
import math

#Parámetros
u = 100     #estar entre 50 a 250
CV = 0.5    #estar entre 0.1 a 0.6
A = 200     #estar entre 100 a 1000
h = 1       #estar entre 0.1 a 1.5
b1 = 20
L = 4

#Variables
#R = sympy.symbols('R')
#Q = sympy.symbols('Q')
#R=1
#Q=1

    #Los X para la función H
    #x_1 = ((R - u*L)/((((CV*u)**2)*L)**(1/2)))
    #x_2 = (R + Q - u*L)/((((CV*u)**2)*L)**(1/2))

        #Funciones
    #H_1= ((1/2)*((((x_1)**2)+1)*(1-norm.cdf(float(x_1)))-float(x_1)*norm.pdf(float(x_1))))
    #H_2= ((1/2)*((((x_2)**2)+1)*(1-norm.cdf(float(x_2)))-float(x_2)*norm.pdf(float(x_2))))

    #G_1= norm.pdf(float(x_1))-(float(x_1))*(1-norm.cdf(float(x_1)))
    #G_2= norm.pdf(float(x_2))-(float(x_2))*(1-norm.cdf(float(x_2)))

##ITERACIONES
alfa= 0.01
beta= 0.05
t=1
#PUNTO INICIAL
#R= 1
#Q= 1

def gradiente (R_actual=1, Q_actual=1, err=0.00001):
    t=symbols('t', real=True)
    R=1
    Q=1
    x_1 = ((R - u*L)/((((CV*u)**2)*L)**(1/2)))
    x_2 = (R + Q - u*L)/((((CV*u)**2)*L)**(1/2))
    H_1= ((1/2)*((((x_1)**2)+1)*(1-norm.cdf(float(x_1)))-float(x_1)*norm.pdf(float(x_1))))
    H_2= ((1/2)*((((x_2)**2)+1)*(1-norm.cdf(float(x_2)))-float(x_2)*norm.pdf(float(x_2))))
    G_1= norm.pdf(float(x_1))-(float(x_1))*(1-norm.cdf(float(x_1)))
    G_2= norm.pdf(float(x_2))-(float(x_2))*(1-norm.cdf(float(x_2)))
    C= h*(R+(Q/2)-u*L)+(h+b1)*((((CV*u)**2)*L)/Q)*(H_1 - H_2)*(A*u/Q)
    C_R= h+(h+b1)*((((CV*u)**2)*L)/Q)*(G_2 - G_1)
    C_R= C_R.evalf(subs={float(R):float(R_actual),float(Q):float(Q_actual)})
    C_Q= ((h/2)-(A*u/(Q**2))-(h+b1)*((((CV*u)**2)*L)/(Q**2))*(H_1 - H_2 - (Q/((((CV*u)**2)*L)**(1/2))*G_2)))
    C_Q= C_Q.evalf(subs={float(R):float(R_actual),float(Q):float(Q_actual)})
    while abs(C_R) > err or abs(C_Q) > err:
        m= C.subs([(R,R_actual + t*C_R),(Q,Q_actual + t*C_Q)])
        n= diff(m,t)
        tn= solve(n,t) [0]
        R_actual = R_actual + tn*C_R
        Q_actual = Q_actual + tn*C_Q
        C_R= h+(h+b1)*((((CV*u)**2)*L)/Q)*(G_2 - G_1)
        C_R= C_R.evalf(subs={R:R_actual,Q:Q_actual})
        C_Q= ((h/2)-(A*u/(Q**2))-(h+b1)*((((CV*u)**2)*L)/(Q**2))*(H_1 - H_2 - (Q/((((CV*u)**2)*L)**(1/2))*G_2)))
        C_Q= C_Q.evalf(subs={R:R_actual,Q:Q_actual})
        print (R_actual, Q_actual)
    


                                                    
    








