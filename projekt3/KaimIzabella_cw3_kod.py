import matplotlib
import math as m
import numpy as np
import matplotlib.pyplot as plt
import geographiclib
matplotlib.use('TkAgg')


a = 6378137
e2 = 0.00669438002290


def kivioji(fi, li, Ai, s):
    n = round(s/1000)
    
    ds = s/n
    for i in range(n):
        Mi = (a*(1 - e2)) / (m.sqrt((1 - e2*m.sin(fi)*m.sin(fi))**3))
        Ni = a / (m.sqrt(1 - e2*np.sin(fi)*np.sin(fi)))

        dfi = (ds * m.cos(Ai)) / Mi
        dAi = (ds * m.sin(Ai) * m.tan(fi)) / Ni

        fm = fi + (1/2)* dfi
        Am = Ai + (1/2)* dAi

        Mm = (a*(1 - e2)) / (m.sqrt((1 - e2 * np.sin(fm)*np.sin(fm))**3))
        Nm = (a)/ (m.sqrt(1 - e2*np.sin(fm)*np.sin(fm)))

        dfm = (ds * m.cos(Am)) / Mm
        dlm = (ds * m.sin(Am)) / (Nm * m.cos(fm))
        dAm = (ds * m.sin(Am) * m.tan(fm)) / Nm

        fi = fi + dfm
        li = li + dlm
        Ai = Ai + dAm
        i += 1
        
    return fi, li, Ai

def vincenty(BA,LA,BB,LB):
    b = a * np.sqrt(1-e2)
    f = 1-b/a
    dL = LB - LA
    UA = np.arctan((1-f)*np.tan(BA))
    UB = np.arctan((1-f)*np.tan(BB))
    L = dL
    while True:
        sin_sig = np.sqrt((np.cos(UB)*np.sin(L))**2 +\
                          (np.cos(UA)*np.sin(UB) - np.sin(UA)*np.cos(UB)*np.cos(L))**2)
        cos_sig = np.sin(UA)*np.sin(UB) + np.cos(UA) * np.cos(UB) * np.cos(L)
        sig = np.arctan2(sin_sig,cos_sig)
        sin_al = (np.cos(UA)*np.cos(UB)*np.sin(L))/sin_sig
        cos2_al = 1 - sin_al**2
        cos2_sigm = cos_sig - (2 * np.sin(UA) * np.sin(UB))/cos2_al
        C = (f/16) * cos2_al * (4 + f*(4 - 3 * cos2_al))
        Lst = L
        L = dL + (1-C)*f*sin_al*(sig+C*sin_sig*(cos2_sigm+C*cos_sig*(-1 + 2*cos2_sigm**2)))
        if abs(L-Lst)<(0.000001/206265):
            break
    
    u2 = (a**2 - b**2)/(b**2) * cos2_al
    A = 1 + (u2/16384) * (4096 + u2*(-768 + u2 * (320 - 175 * u2)))
    B = u2/1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    d_sig = B*sin_sig * (cos2_sigm + 1/4*B*(cos_sig*(-1+2*cos2_sigm**2)\
            - 1/6 *B*cos2_sigm * (-3 + 4*sin_sig**2)*(-3+4*cos2_sigm**2)))
    sAB = b*A*(sig-d_sig)
    A_AB = np.arctan2((np.cos(UB) * np.sin(L)),(np.cos(UA)*np.sin(UB) - np.sin(UA)*np.cos(UB)*np.cos(L)))
    A_BA = np.arctan2((np.cos(UA) * np.sin(L)),(-np.sin(UA)*np.cos(UB) + np.cos(UA)*np.sin(UB)*np.cos(L))) + np.pi

    return sAB, np.rad2deg(A_AB), np.rad2deg(A_BA)



f1 = np.deg2rad(57) #54 + 12 * 0.25
l1 = np.deg2rad(22) #16 + 12 * 0.5

s12 = 20000 #m
s23 = 35000 #m
s34 = 20000 #m
s41 = 35000 #m

A12 = np.deg2rad(0)
A23 = np.deg2rad(90)
A34 = np.deg2rad(180)
A41 = np.deg2rad(270)

f2 = kivioji(f1, l1, A12, s12)[0]
l2 = kivioji(f1, l1, A12, s12)[1]

f3 = kivioji(f2, l2, A23, s23)[0]
l3 = kivioji(f2, l2, A23, s23)[1]

f4 = kivioji(f3, l3, A34, s34)[0]
l4 = kivioji(f3, l3, A34, s34)[1]

f1p = kivioji(f4, l4, A41, s41)[0]
l1p = kivioji(f4, l4, A41, s41)[1]

print('Wsp punktu 1: ', np.rad2deg(f1), np.rad2deg(l1))
print('Wsp punktu 2: ', np.rad2deg(f2), np.rad2deg(l2))
print('Wsp punktu 3: ', np.rad2deg(f3), np.rad2deg(l3))
print('Wsp punktu 4: ', np.rad2deg(f4), np.rad2deg(l4))
print('Wsp punktu 1*: ', np.rad2deg(f1p), np.rad2deg(l1p))


l1p = np.rad2deg(l1p)
f1p = np.rad2deg(f1p)
l2 = np.rad2deg(l2)
f2 = np.rad2deg(f2)
l3 = np.rad2deg(l3)
f3 = np.rad2deg(f3)
l4 = np.rad2deg(l4)
f4 = np.rad2deg(f4)
l1 = np.rad2deg(l1)
f1 = np.rad2deg(f1)

plt.scatter(l1,f1)
plt.text(l1,f1, s='1')
plt.scatter(l2,f2)
plt.text(l2,f2, s='2')
plt.scatter(l3,f3)
plt.text(l3,f3, s='3')
plt.scatter(l4,f4)
plt.text(l4,f4, s='4')
plt.scatter(l1p,f1p)
plt.text(l1p,f1p, s='1*')
plt.show()

from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84
p = geod.Polygon()

p.AddPoint(f4, l4)
p.AddPoint(f3, l3)
p.AddPoint(f2, l2)
p.AddPoint(f1, l1)


num, perim, area = p.Compute()

print('Pole:', area/1000000, '[km^2]')

