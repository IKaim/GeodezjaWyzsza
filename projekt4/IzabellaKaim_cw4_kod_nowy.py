import matplotlib
import math as m
import numpy as np
import matplotlib.pyplot as plt
import geographiclib
from pyproj import Transformer
matplotlib.use('TkAgg')

global a, e2, ep2
a = 6378137
e2 = 0.006694380022900
ep2 = 0.0067394909677548

def deg2dms(dd):
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    return dms

def fN(f, a, e2):
    N =  a / np.sqrt(1 - e2 * (np.sin(f) ** 2))
    return N

def fM(f, a, e2):
    M = a * (1 - e2) / np.sqrt((1 - e2 * (np.sin(f) ** 2)) ** 3)
    return M

def fm(y, r):
    m = 1 + (y**2)/(2*(r**2)) + (y**4)/(24*(r**4))
    return m

def gk2u1992(xgk, ygk):
    x = xgk*0.9993 - 5300000
    y = ygk*0.9993 + 500000
    return x, y

def gk2u2000(xgk, ygk, nr):
    x = 0.999923 * xgk
    y = 0.999923 * ygk + 500000 + nr*1000000
    return x, y

def R(f):
    a = 6378137
    e2 = 0.006694380022900
    M = a * (1 - e2) / np.sqrt((1 - e2 * (np.sin(f) ** 2)) ** 3)
    N =  a / np.sqrt(1 - e2 * (np.sin(f) ** 2))
    return m.sqrt(M*N)

def R2(lat1, lat2):
    a = 6378137
    e2 = 0.006694380022900
    f = lat1 + (lat2 - lat1)/2
    M = a * (1 - e2) / np.sqrt((1 - e2 * (np.sin(f) ** 2)) ** 3)
    N =  a / np.sqrt(1 - e2 * (np.sin(f) ** 2))
    return m.sqrt(M*N)

def geo2gk(f, l, lon_0):
    f = np.deg2rad(f)
    l = np.deg2rad(l)
    lon_0 = np.deg2rad(lon_0)

    a = 6378137
    e2 = 0.006694380022900
    b2 = a**2 * (1 - e2)
    ep2 = (a**2 - b2) / b2
    d_lon = l - lon_0
    A0 = 1 - e2/4 - 3*e2**2/64 - 5*e2**3/256
    A2 = (3/8)*(e2 + e2**2/4 + 15*e2**3/128)
    A4 = (15/256)*(e2**2 + 3*e2**3/4)
    A6 = 35*e2**3/3072

    N = fN(f, a, e2)
    t = np.tan(f)
    n2 = ep2 * m.cos(f)*m.cos(f)
    
    sigma = a*(A0*f - A2*m.sin(2*f) + A4*m.sin(4*f)-A6*m.sin(6*f))

    x = sigma + ((d_lon**2)/2) * N * m.sin(f) * m.cos(f) * ( 1 + ((d_lon**2)/12) * (m.cos(f)**2) * (5 - (t**2) + 9*n2 + 4*(n2**2)) + ((d_lon**4)/360) * (m.cos(f)**4) * (61 - 58*(t**2)+(t**4) + 270*n2 - 330*n2*(t**2)))
    y = d_lon * N * m.cos(f)
    y = y * (1 + ((d_lon**2)/6)*(m.cos(f)**2)*(1 - (t**2) + n2)+((d_lon**4)/120) * (m.cos(f)**4)*(5 - 18*(t**2) + (t**4)+ 14*n2 - 58*n2*(t**2)))

    return x, y

def distance(x1,y1, x2, y2):
    return m.sqrt((x2-x1)**2 + (y2-y1)**2)

def gk2geo(xk, yk, lon_0):
    a = 6378137
    e2 = 0.006694380022900
    b2 = a**2 * (1 - e2)
    ep2 = (a**2 - b2) / b2
    
    A0 = 1 - e2/4 - 3*e2**2/64 - 5*e2**3/256
    A2 = (3/8)*(e2 + e2**2/4 + 15*e2**3/128)
    A4 = (15/256)*(e2**2 + 3*e2**3/4)
    A6 = 35*e2**3/3072
    
    f_prev = xk/(a * A0)

    sigma = a*(A0*f_prev - A2*m.sin(2*f_prev) + A4*m.sin(4*f_prev)-A6*m.sin(6*f_prev))
    fi = f_prev + (xk - sigma)/(a * A0)

    min_dif = np.deg2rad(0.000001/3600)

    while True:
        if abs(fi - f_prev) < min_dif:
            break
        else:
            f_prev = fi
            sigma = a*(A0*f_prev - A2*m.sin(2*f_prev) + A4*m.sin(4*f_prev)-A6*m.sin(6*f_prev))
            fi = f_prev + (xk - sigma)/(a * A0)

    return fi

def fil2gamma(f, l, lon_0):
    f = np.deg2rad(f)
    l = np.deg2rad(l)
    lon_0 = np.deg2rad(lon_0)
    d_lon = l - lon_0
    a = 6378137
    e2 = 0.006694380022900
    b2 = a**2 * (1 - e2)
    ep2 = (a**2 - b2)/ b2
    n2 = ep2 * m.cos(f)*m.cos(f)
    t = np.tan(f)

    g = d_lon * np.sin(f) + ((d_lon**3)/3 * np.sin(f) * (np.cos(f)**2) * (1 + 3 * n2 + 2 * (n2**2))) + (((d_lon**5)/15 * np.sin(f) * (np.cos(f)**4) * (2 - t ** 2)))

    return g

def xy2gamma(x, y, lon_0):
    a = 6378137
    e2 = 0.006694380022900
    b2 = a**2 * (1 - e2)
    A0 = 1 - (e2/4) - (3 * e2**2 / 64) - (5 * e2**3 / 256)

    f = gk2geo(x, y, lon_0)
    
    n2 = ep2 * m.cos(f)*m.cos(f)
    t = np.tan(f)
    N =  a / np.sqrt(1 - e2 * (np.sin(f) ** 2))

    gamma = y / N*t*(1 - (y**2 / (3 * N**2) * (1 + t**2 - n2 - 2*n2**2)) + (y**4 / (15 * N**4) * (2 + 5 * t**2 + 3 * t**4)))

    return gamma
    
def redukcjaAB(xA, yA, xB, yB, f1, f2):
    a = 6378137
    e2 = 0.006694380022900

    fm = (f1 + f2) / 2
    Rm = R(fm)

    red = (xB - xA) * (2 * yA + yB) / (6 * Rm**2)

    return red

plik = open('wsp.txt', 'w')


punkty =[[57.0, 22.0],
        [57.179595310234625, 22.0],
        [57.17826148114223, 22.578698075251012],
        [56.998666132425434, 22.578698075251012]]

lat = [57.0, 57.179595310234625, 57.17826148114223, 56.998666132425434, 57.0]
lon = [22.0, 22.0, 22.578698075251012, 22.578698075251012, 22.0]



pkt_2000 = []
for i in range(len(lon)):
    lon_0 = 21
    nr = 7
    xgk,ygk = geo2gk(lat[i], lon[i], lon_0)
    x,y = gk2u2000(xgk, ygk, nr)
    pkt_2000.append([x,y])
    

pkt_1992 = []
for i in range(len(lon)):
    lon_0 = 19
    xgk,ygk = geo2gk(lat[i], lon[i], lon_0)
    x,y = gk2u1992(xgk, ygk)
    pkt_1992.append([x,y])

print("1992:\n", pkt_1992)
print()
print("2000:\n", pkt_2000)

##x = [el[0] for el in pkt_1992]
##y = [el[1] for el in pkt_1992]
##
##plt.scatter(x, y)
##plt.show()
##
##x = [el[0] for el in pkt_2000]
##y = [el[1] for el in pkt_2000]
##
##plt.scatter(x, y)
##plt.show()


ss2000 = np.array([m.dist(pkt_2000[i], pkt_2000[i-1]) for i in range(1,len(pkt_2000))])

ss1992 = np.array([m.dist(pkt_1992[i], pkt_1992[i-1]) for i in range(1,len(pkt_1992))])

m0_2000 = 0.999923
m0_1992 = 0.9993

gk_2000 = []
for i in range(len(lon)):
    lon_0 = 21
    nr = 7
    xgk,ygk = geo2gk(lat[i], lon[i], lon_0)
    gk_2000.append([xgk,ygk])

gk_1992 = []
for i in range(len(lon)):
    lon_0 = 19
    xgk,ygk = geo2gk(lat[i], lon[i], lon_0)
    gk_1992.append([xgk,ygk])

sgk2000 = ss2000/m0_2000
Rs = []
for i in range(len(lat) - 1):
    Rs.append(R2(lat[i], lat[i+1]))

rs2000 = []
for i in range(len(ss2000)):
    r = ss2000[i] * (gk_2000[i][1]**2 + gk_2000[i][1]*gk_2000[i + 1][1] + gk_2000[i + 1][1]**2) / (6 * Rs[i]**2)
    rs2000.append(r)
          
ss_elip2000 = []
for i in range(len(rs2000)):
          ss_elip2000.append(sgk2000[i] - rs2000[i])

print()
print("Długości odcinków - układ 2000 \n", ss_elip2000)
print()

sgk1992 = ss1992/m0_1992
Rs = []
for i in range(len(lat) - 1):
    Rs.append(R2(lat[i], lat[i+1]))

rs1992 = []
for i in range(len(ss1992)):
    r = ss1992[i] * (gk_1992[i][1]**2 + gk_1992[i][1]*gk_1992[i + 1][1] + gk_1992[i + 1][1]**2) / (6 * Rs[i]**2)
    rs1992.append(r)
          
ss_elip1992 = []
for i in range(len(rs1992)):
          ss_elip1992.append(sgk1992[i] - rs1992[i])

print("Długości odcinków - układ 1992 \n", ss_elip1992)
print()

katy_2000 = []
for i in range(len(gk_2000)-1):
    xA, yA = gk_2000[i]
    xB, yB = gk_2000[i+1]
    alfa = np.arctan2(yB-yA, xB-xA)
    #alfa = deg2dms(np.rad2deg(alfa))
    katy_2000.append(alfa)

print("Kąty 2000: \n", katy_2000)
    
katy_1992 = []
for i in range(len(gk_1992)-1):
    xA, yA = gk_1992[i]
    xB, yB = gk_1992[i+1]
    alfa = np.arctan2(yB-yA, xB-xA)
    #alfa = deg2dms(np.rad2deg(alfa))
    katy_1992.append(alfa)

print("Kąty 1992: \n", katy_1992)

print()

gamma_2000 = []
gamma_1992 = []
for i in range(len(lat) - 1):
    gamma_2000.append(fil2gamma(lat[i], lon[i], 21))
    gamma_1992.append(fil2gamma(lat[i], lon[i], 19))

print("Zbieżności południków 1992: \n", gamma_1992) #, [np.rad2deg(gamma_1992[i]) for i in range(len(gamma_1992))])
print()
print("Zbieżności południków 2000: \n", gamma_2000)
print()

red_1992 = []
red_2000 = []
for i in range(len(gk_1992) - 1):
    xA, yA = gk_1992[i]
    xB, yB = gk_1992[i+1]
    redAB = redukcjaAB(xA, yA, xB, yB, lat[i], lat[i+1])
    red_1992.append(redAB)

print("Reduckje kierunków 1992: \n", red_1992)
print()

for i in range(len(gk_2000) - 1):
    xA, yA = gk_2000[i]
    xB, yB = gk_2000[i+1]
    redAB = redukcjaAB(xA, yA, xB, yB, lat[i], lat[i+1])
    red_2000.append(redAB)

print("Reduckje kierunków 2000: \n", red_2000)
print()

az_elip_2000 = []
az_elip_1992 = []

for i in range(len(red_1992)):
    az = katy_1992[i] + gamma_1992[i] + red_1992[i]
    az = np.rad2deg(az)
    az = deg2dms(az)
    az_elip_1992.append(az)

print("Azymuty odcinków na elipsoidzie - 1992: \n", az_elip_1992)
print()

for i in range(len(red_2000)):
    az = katy_2000[i] + gamma_2000[i] + red_2000[i]
    az = np.rad2deg(az)
    az = deg2dms(az)
    az_elip_2000.append(az)

print("Azymuty odcinków na elipsoidzie - 2000: \n", az_elip_2000)
print()

punkty_92 = pkt_1992[0:4]
#pole PL-1992
PP = 0
for i in range(0, 4):
    if i == 0:
        PP += punkty_92[i][0] * (punkty_92[1][1] - punkty_92[3][1])
    elif i == 3:
        PP += punkty_92[i][0] * (punkty_92[0][1] - punkty_92[i-1][1])
    else:
        PP += punkty_92[i][0] * (punkty_92[i+1][1] - punkty_92[i-1][1])
Pole_1992 = PP/2

#kontrola
mPP = 0
for i in range(0, 4):
    if i == 0:
        mPP += punkty_92[i][1] * (punkty_92[1][0] - punkty_92[3][0])
    elif i == 3:
        mPP += punkty_92[i][1] * (punkty_92[0][0] - punkty_92[i-1][0])
    else:
        mPP += punkty_92[i][1] * (punkty_92[i+1][0] - punkty_92[i-1][0])
mPole_1992 = np.abs(mPP/2)

punkty_2000 = pkt_2000[0:4]
#pole PL-2000
PP = 0
for i in range(0, 4):
    if i == 0:
        PP += punkty_2000[i][0] * (punkty_2000[1][1] - punkty_2000[3][1])
    elif i == 3:
        PP += punkty_2000[i][0] * (punkty_2000[0][1] - punkty_2000[i-1][1])
    else:
        PP += punkty_2000[i][0] * (punkty_2000[i+1][1] - punkty_2000[i-1][1])
Pole_2000 = PP/2

#kontrola
mPP = 0
for i in range(0, 4):
    if i == 0:
        mPP += punkty_2000[i][1] * (punkty_2000[1][0] - punkty_2000[3][0])
    elif i == 3:
        mPP += punkty_2000[i][1] * (punkty_2000[0][0] - punkty_2000[i-1][0])
    else:
        mPP += punkty_2000[i][1] * (punkty_2000[i+1][0] - punkty_2000[i-1][0])
mPole_2000 = np.abs(mPP/2)

print('pole 1992:\n', Pole_1992, 'kon:', mPole_1992)
print('\npole 2000:\n', Pole_2000, 'kon:', mPole_2000)


transformer = Transformer.from_crs(4326, 3035)

pkt_LAEA = []
for i in range(len(lon)-1):
    x, y = transformer.transform(lat[i], lon[i])
    pkt_LAEA.append([x,y])
    print(x,y)

PP = 0
for i in range(0, 4):
    if i == 0:
        PP += pkt_LAEA[i][0] * (pkt_LAEA[1][1] - pkt_LAEA[3][1])
    elif i == 3:
        PP += pkt_LAEA[i][0] * (pkt_LAEA[0][1] - pkt_LAEA[i-1][1])
    else:
        PP += pkt_LAEA[i][0] * (pkt_LAEA[i+1][1] - pkt_LAEA[i-1][1])

Pole_LAEA = PP/2
 
print("Pl-LAEA:", Pole_LAEA)
