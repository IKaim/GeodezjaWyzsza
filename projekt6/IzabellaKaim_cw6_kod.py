import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math as m
matplotlib.use('TkAgg')

def hirvonen(x,y,z):
    e2 = 0.00669438002290
    a = 6378137
    p = np.sqrt(x**2 + y**2)
    phi = np.arctan(z/(p*(1-e2)))
    
    while True:
        phi_stare = phi
        N = a / m.sqrt(1-(e2*(m.sin(phi)**2)))
        h = p/np.cos(phi) - N
        phi = np.arctan(z/(p * (1-(N*e2)/(N-h))))

        if abs(phi-phi_stare) < (0.000001/206265):
            break

    N = a / m.sqrt(1-(e2*(m.sin(phi)**2)))
    h = p/np.cos(phi) - N 
    lam = np.arctan2(y,x)

    return (phi, lam, h)

def NEU(f, l):
    n = np.array([
        -m.sin(f)*m.cos(l),
        -m.sin(f)*m.sin(l),
        m.cos(f)])
    e = np.array([
        -m.sin(f),
        m.cos(f),
        0 ])
    u = np.array([
        m.cos(f)*m.cos(l),
        m.cos(f)*m.sin(l),
        m.sin(f)])

    return [n,e,u]


def countRneu(f, l):
    Rneu = np.column_stack((NEU(f, l)))

    return Rneu


dane = np.genfromtxt('model.txt')
czas = dane[:,2]
xyz = dane[:,3:6]


dxyz = xyz - np.mean(xyz, axis=0)   

xm, ym ,zm = np.mean(xyz, axis = 0)
fi, lam, h = hirvonen(xm, ym ,zm)
R = countRneu(fi, lam)
neu_all = []
for xyz1 in dxyz:
    neu = R.T @ xyz1
    neu_all.append(neu)
neu_all = np.array(neu_all)


# różnice wartości względem średniej zarejestrowanych wartości
x = czas - czas[0]
y_model_all = []
y = dxyz[:,0]
for y in neu_all.T:
    A = np.column_stack((x, np.ones(len(x))))
    xx = np.linalg.inv(A.T@A)@(A.T@y)

    y_model = A@xx
    y_model_all.append(y_model)
y_model_all = np.array(y_model_all).T


fig, ax = plt.subplots(3,1,figsize=(20,10))
for i in range(3):
    ax[i].plot(x, neu_all[:,i])
    ax[i].plot(x, y_model_all[:,i], color = 'k')
#plt.xlabel('Time [years]')
#plt.ylabel('[m]')
plt.show()


