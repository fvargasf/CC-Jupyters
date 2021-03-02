import numpy as np
from scipy.stats import norm
from scipy.optimize import fsolve

u = 100     #estar entre 50 a 250
CV = 0.5    #estar entre 0.1 a 0.6
s = CV*u
A = 200     #estar entre 100 a 1000
h = 1       #estar entre 0.1 a 1.5
b1 = 20
L = 4

n_iter = 0
H = lambda x: 0.5*((x**2+1)*(1-norm.cdf(x)) - x*norm.pdf(x))
G = lambda x: norm.pdf(x) - x*(1-norm.cdf(x))

Q_k = lambda R, Q: ( (2*A*u)/h + (((2*(h+b1)*L*s**2)**0.5)/h)*(H(((R-u*L)/(L*s**2)**0.5))- H(((R+Q-u*L)/(L*s**2)**0.5)) - (Q/(L*s**2)**0.5)*((R+Q-u*L)/(L*s**2)) ) )**0.5 
C = lambda R, Q: h*(R+0.5*Q-u*L)+(h+b1)*(((L*s**2)**0.5)/Q)*(H((R-u*L)/(L*s**2)**0.5) - H((R+Q-u*L)/(L*s**2)**0.5)) + A*u/Q

Q_0 = ((2*u*A)/h)**0.5

R_0 = lambda R: h+(h+b1)*(((L*s**2)**0.5)/Q_0)*(G((R+Q_0-u*L)/((L*s**2)**0.5))-G((R-u*L)/((L*s**2)**0.5)))
R_0 = fsolve(R_0, Q_0)[0]

Q_n = [Q_0]
R_n = [R_0]

Q_n.append(Q_k(R_0, Q_0))

R_k = lambda R: h+(h+b1)*(((L*s**2)**0.5)/Q_k(R_0, Q_0))*(G((R+Q_k(R_0, Q_0)-u*L)/((L*s**2)**0.5))-G((R-u*L)/((L*s**2)**0.5)))
R_n.append(fsolve(R_k, R_n[-1])[0])

n_iter = 0
while abs(C(R_n[len(R_n)-1], Q_n[len(Q_n)-1])- C(R_n[len(R_n)-2], Q_n[len(Q_n)-2])) > 1e-3 and n_iter < 20:
    Q_n.append(Q_k(R_n[len(R_n)-1], Q_n[len(Q_n)-1]))
    
    R_k = lambda R: h+(h+b1)*(((L*s**2)**0.5)/Q_k(R_n[len(R_n)-1], Q_n[len(Q_n)-1]))*(G((R+Q_k(R_n[len(R_n)-1], Q_n[len(Q_n)-1])-u*L)/((L*s**2)**0.5))-G((R-u*L)/((L*s**2)**0.5)))
    R_n.append(fsolve(R_k, R_n[-1])[0])
    
    n_iter = n_iter + 1
print("Q:", Q_n[len(Q_n)-1], "\nR:", R_n[len(R_n)-1], "\nC:", C(R_n[len(R_n)-1], Q_n[len(Q_n)-1]) ,"\nNumero iteraciones:", n_iter)
