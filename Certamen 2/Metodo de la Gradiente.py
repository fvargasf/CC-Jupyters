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
R0 = sympy.symbols('R0')
Q0 = sympy.symbols('Q0')
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
t=1
#PUNTO INICIAL
#R= 1
#Q= 1


Q0=sqrt(2*u*A/h)
R0="R0"
x_1 = int((R0 - u*L)/((((CV*u)**2)*L)**(1/2)))
x_2 = int((R0 + Q0 - u*L)/((((CV*u)**2)*L)**(1/2)))
H_1= ((1/2)*(((x_1**2)+1)*(1-norm.cdf(x_1))-(x_1)*norm.pdf(x_1)))
H_2= ((1/2)*(((x_2**2)+1)*(1-norm.cdf(x_2))-x_2*norm.pdf(x_2)))





print (H_1)
print (H_2)










                                                    
    








