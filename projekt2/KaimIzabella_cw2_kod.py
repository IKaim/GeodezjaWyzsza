import matplotlib
import math as m
import numpy as np
import matplotlib.pyplot as plt
import pyproj as pp
from datetime import datetime
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
#import cartopy.crs as ccrs
#import cartopy.io.img_tiles as cimgt
matplotlib.use('TkAgg')

def deg2dms(dd):
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    return dms

def dms2rad(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg)
    return rad

def geo2xyz(lat, lon, h):
    a = 6378137
    e2 = 0.00669438002290
    N = a / (m.sqrt(1 - e2*np.sin(lat)*np.sin(lat)))
    x = (N+h)*np.cos(lat)*np.cos(lon)
    y = (N+h)*np.cos(lat)*np.sin(lon)
    z = (N*(1-e2)+h)*np.sin(lat)
    xyz = [x, y, z]
    return np.array(xyz)

def geo2neu(f1, l1, h1, f2, l2, h2 ):
    f1 = np.deg2rad(f1)
    l1 = np.deg2rad(l1)
    f2 = np.deg2rad(f2)
    l2 = np.deg2rad(l2)
    R = np.array([[-np.sin(f1)*np.cos(l1), -np.sin(l1), np.cos(f1)*np.cos(l1)],
                  [-np.sin(f1)*np.sin(l1), np.cos(l1), np.cos(f1)*np.sin(l1)],
                  [np.cos(f1), 0, np.sin(f1)]])
    x1, y1, z1 = geo2xyz(f1, l1, h1, a, e2)
    x2, y2, z2 = geo2xyz(f2, l2, h2, a, e2)
    x = np.array([[x2 - x1],
                  [y2 - y1],
                  [z2 - z1]])
    neu = np.dot(R.transpose(), x)
    return neu

def read_flightradar(file):
    with open(file, 'r') as f:
        i = 0
        size= []
        Timestamp = []; date = []; UTC = []; Latitude = []; Longitude = []; 
        Altitude = []; Speed = []; Direction = []
        for linia in f:
            if linia[0:1]!='T':
                splited_line = linia.split(',')
                size.append(len(splited_line))
                i+=1
                Timestamp.append(int(splited_line[0]))
                full_date = splited_line[1].split('T')
                date.append(list(map(int,full_date[0].split('-'))))
                UTC.append(list(map(int, full_date[1].split('Z')[0].split(':'))))
                Callsign = splited_line[2]
                Latitude.append(float(splited_line[3].split('"')[1]))
                Longitude.append(float(splited_line[4].split('"')[0]))
                Altitude.append(float(splited_line[5]))
                Speed.append(float(splited_line[6]))
                Direction.append(float(splited_line[7]))
    all_data = np.column_stack((np.array(Timestamp), np.array(date), np.array(UTC),
                                np.array(Latitude), np.array(Longitude), np.array(Altitude),
                                np.array(Speed), np.array(Direction)))            

    return all_data, i

def hms2deg(t):
    w = t[0] + t[1]/60 + t[2]/3600
    return w

def na_array(Warszawa, data):
    f_w = np.deg2rad(Warszawa[7])
    l_w = np.deg2rad(Warszawa[8])
    h_w = Warszawa[9] * 0.3048 + 135.4

    wspWarszawa = np.array(geo2xyz(f_w, l_w, h_w))
    
    wspSamolot = np.array([geo2xyz(dms2rad(deg2dms(data[i][7])),
                                   dms2rad(deg2dms(data[i][8])),
                                   data[i][9]*0.3048 + hel) for i in range(len(data))])

    return wspWarszawa, wspSamolot

def wektorsamolotu(wspSamolot, wspWarszawa):
    wektor = wspSamolot - wspWarszawa
    return wektor

def NEU(f_w, l_w):
    n = np.array([
        -m.sin(f_w)*m.cos(l_w),
        -m.sin(f_w)*m.sin(l_w),
        m.cos(f_w)])
    e = np.array([
        -m.sin(l_w),
        m.cos(l_w),
        0 ])
    u = np.array([
        m.cos(f_w)*m.cos(l_w),
        m.cos(f_w)*m.sin(l_w),
        m.sin(f_w)])

    return [n,e,u]

def RTneu(n,e,u):
    Rneu = np.column_stack((n,e,u))
    RTneu = Rneu.transpose()
    return RTneu

def xslneu(RTneu, wektor):
    xslneu = []
    for i in range(len(wektor)):
        xslneu.append(RTneu @ wektor[i])
    xn = []
    xe = []
    xu = []
    for i in range(len(xslneu)):
        xn.append(np.array(xslneu[i][0]))
        xe.append(np.array(xslneu[i][1]))
        xu.append(np.array(xslneu[i][2]))
    return [xn, xe, xu]

def dl_wektora(xn, xe, xu):
    s = np.sqrt(np.square(xn) + np.square(xe) + np.square(xu))
    return s

def azymut(xe, xn):
    Az = np.array(np.arctan2(xe,xn))
    #Az[Az<0] = Az[Az<0] + 2 * m.pi
    return Az

def wys(xn, xe, xu):
    h = np.array(np.arcsin(xu/dl_wektora(xn, xe, xu)))
    return h

data = read_flightradar('W61539_2dde28b4.csv')[0]
ile = read_flightradar('W61539_2dde28b4.csv')[1]

a = 6378137
e2 = 0.00669438002290
hel = 135.4 #metry

#Warszawa
Warszawa = data[0, :]
data = data[1:, :]

#print(data)
f_w = np.deg2rad(Warszawa[7])
l_w = np.deg2rad(Warszawa[8])
h_w = Warszawa[9] * 0.3048 + hel

startIndex = 0
endIndex = 0

for i, elem in enumerate(data):
    if elem[9] != 0:
        startIndex = i - 1
        break

for i in reversed(range(len(data))):
    if data[i][9] != 0:
        endIndex = i + 1
        break

data = data[startIndex:endIndex+1,:]

airportCoords, planeCoords = na_array(Warszawa, data)
wektor = wektorsamolotu(planeCoords, airportCoords)

n, e, u = NEU(np.deg2rad(Warszawa[7]), np.deg2rad(Warszawa[8]))
RTneu = RTneu(n, e, u)

xn, xe, xu = xslneu(RTneu, wektor)

s = dl_wektora(xn, xe, xu)
Az = azymut(xe, xn)
hs = wys(xn, xe, xu)

start = datetime(int(Warszawa[1]), int(Warszawa[2]), int(Warszawa[3]),int(Warszawa[4]), int(Warszawa[5]), int(Warszawa[6]))

times = [datetime(int(data[i][1]), int(data[i][2]), int(data[i][3]), int(data[i][4]), int(data[i][5]), int(data[i][6])) - start for i in range(len(data))]
times = np.array([time.seconds + (24*3600) * time.days for time in times])

times -= times[0]


underHorizon = []
for i, h in enumerate(hs):
    if h < 0 and data[i][7] > 53:
        print(h, deg2dms(data[i][7]), deg2dms(data[i][8]))
        print(times[i]/60)
        break
        
#print(hs[0])
wys = [ data[i][9]*0.3048 + hel for i in range(len(data))]
plt.plot(times/3600, wys)
plt.xlabel("godziny [h]")
plt.ylabel("wysokość [m]")
plt.show()

plt.plot(times/3600, s/1000)
plt.xlabel("godziny [h]")
plt.ylabel("odległość [km]")
plt.show()

dane = []
for i in range(len(data)):
    t = times[i]
    h = data[i][9]*0.3048 + hel
    lon = data[i][8]
    lat = data[i][7]
    dane.append({'czas': t, 'h': h, 'lat': lat, 'lon': lon})

dane = pd.DataFrame(dane)
##fig = px.line_mapbox(dane, lat='lat', lon='lon')
##fig.update_layout(mapbox_style="open-street-map")
##fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
##fig.show()

hsdeg = [np.rad2deg(hs[i]) for i in range(len(hs))]
plt.plot(times/3600, hsdeg)
plt.xlabel("godziny [h]")
plt.ylabel("wysokość horyzontalna [°]")
plt.show()


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(polar = True)
ax.set_theta_zero_location('N') # ustawienie kierunku północy na górze wykresu
ax.set_theta_direction(-1)

ax.set_yticks(range(0, 90+10, 10))

end = max(s)

yLabel = [str(round(int(end/9)/1000 * i, 3)) + "km" for i in range(10)]
ax.set_yticklabels(yLabel)
ax.set_rlim(0,90)

ax.scatter(Az, s/24000, s=1)
plt.show()
