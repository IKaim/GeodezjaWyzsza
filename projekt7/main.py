import math as m
import numpy as np

a_kras = 6378245
e2_kras = 0.00669342

def blh2xyz(f, l, h):
    a = 6378245
    e2 = 0.00669342
    f = np.deg2rad(f)
    l = np.deg2rad(l)
    N = a / (m.sqrt(1 - e2*np.sin(f)*np.sin(f)))
    x = (N+h)*np.cos(f)*np.cos(l)
    y = (N+h)*np.cos(f)*np.sin(l)
    z = (N*(1-e2)+h)*np.sin(f)
    return [x, y, z]

punkty = []

with open('dataKRAS.txt', 'r') as plik:
    for linia in plik:
        linia = linia.strip().split()
        punkty.append(linia)

#print(punkty)

def kras2grs80(xk, yk, zk):
    k = -0.84078048 * 10**(-6)
    alfa = 1.73888389 * 10**(-6)
    beta = 0.25613864 * 10**(-6)
    gam = -4.08960007 * 10**(-6)

    B = np.array([[k, gam, -beta],
                  [-gam, k, alfa],
                  [beta, -alfa, k]])
    xyz = np.array([[xk], [yk], [zk]])
    #xyz = np.transpose(xyz)
    T = [33.4297, -1146.5746, -76.2865]
    #T = np.array([[-33.4297], [1146.5746], [76.2865]])
    xk = xk-T[0]
    yk = yk-T[1]
    zk = zk-T[2]
    xg = xk + k*xk + gam*yk + (-beta)*zk
    yg = yk + (-gam)*xk + k*yk + (alfa)*zk
    zg = zk + beta*xk + (-alfa)*yk + k*zk
##    w1 = np.multiply(B, xyz)
##    w2 = np.add(w1, xyz)
##    w3 = np.add(w2, T)
##    print(w3)
    return xg, yg, zg

for p in punkty:
##    x = np.deg2rad(float(p[1]))
##    y = np.deg2rad(float(p[2]))
##    z = np.deg2rad(float(p[3]))
    x = float(p[1])
    y = float(p[2])
    z = float(p[3])
    
    x,y,z = blh2xyz(x,y,z)
    kras = kras2grs80(x, y, z)
    print(kras)

#pkt1 = 3673743.94493585  1410085.19668459  5002888.19185362


def Np(f, a, e2):
    N = a / (m.sqrt(1 - e2*np.sin(f)*np.sin(f)))
    return N

def hirvonen(x,y,z,a,e2):
    p = np.sqrt(x**2 + y**2)
    phi = np.arctan(z/(p * (1-e2)))
    
    while True:
        phi_stare = phi
        N = Np(phi, a, e2) #funkcja Np
        h = p/np.cos(phi) - N
        phi = np.arctan(z/(p * (1 - (N*e2)/(N+h))))
        if abs(phi-phi_stare) < (0.000001/206265):
            break
    N = Np(phi, a, e2) #funkcja Np
    h = p/np.cos(phi) - N
    lam = np.arctan2(y,x)
    return phi, lam, h

x = 5000000
y = 1000000
z = 6000000
a = 6378137
e2 = 0.00669438002290
f, l, h = hirvonen(x,y,z,a,e2)
#print('Hirvonen: ' , np.rad2deg(f), np.rad2deg(l), h)
