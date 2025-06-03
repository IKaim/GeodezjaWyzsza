import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math as m
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
matplotlib.use('TkAgg')

def dms2deg(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    deg = d+m/60+s/3600
    return deg

def deg2dms(dd):
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    return dms

def hms2rad(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg*15)
    return rad

def dms2rad(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    
    deg = d+m/60+s/3600
    rad = np.deg2rad(deg)
    return rad

def hms2sec(hms):
    sec = hms[0]*3600 + hms[1] * 60 + hms[2]
    return sec

def sec2hms(s):
    hd = s/3600
    h = int(np.trunc(hd))
    m = int(np.trunc((hd-h) * 60))
    s = ((hd-h) * 60 - m) * 60
    hms = [h,abs(m),abs(s)]
    return hms

def rad2hms(rad):
    dd = np.rad2deg(rad)
    dd = dd/15
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    return dms

def rad2dms(rad):
    dd = np.rad2deg(rad)
    dd = dd
    deg = int(np.trunc(dd))
    mnt = int(np.trunc((dd-deg) * 60))
    sec = ((dd-deg) * 60 - mnt) * 60
    dms = [deg, abs(mnt), abs(sec)]
    return dms

def dms2hms(dms):
    sall = dms[0] * (4*60) + dms[1] * 4 + dms[2]/15    
    h = int(sall//3600)
    m = int((sall%3600)//60)
    s = sall%60
    return [h,m,s]

def julday(y,m,d,h):
    if m <= 2:
        y = y - 1
        m = m + 12
    jd = np.floor(365.25*(y+4716))+np.floor(30.6001*(m+1))+d+h/24-1537.5;
    return jd

def GMST(jd):
    T = (jd - 2451545) / 36525
    Tu = jd - 2451545
    g = 280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*T**2-T**3/38710000
    g = (g%360) / 15
    return g

def LST(y,m,d,h,lb,alfa):
    jd = julday(y,m,d,h)
    gm = GMST(jd)
    lst = gm*15 + lb
    if lst > 360:
        lst = lst - 360
    lst = np.deg2rad(lst)
    return lst

def A(fi, dek, t):
    h = m.asin((m.sin(fi)*m.sin(dek)) + m.cos(fi) * m.cos(dek) * m.cos(t))
    g = -np.cos(dek)*np.sin(t)
    d = np.cos(fi)*np.sin(dek)-np.sin(fi)*np.cos(dek)*np.cos(t)
    a = m.atan2(g,d)
    if a < 0:
        a = a + (2*m.pi)
    return a, h

def UTC2toUTC0(UTC2):
    return UTC2 - 2

def kh(lst, alfa):
    return lst - alfa

#numer gwiazdy 12 - RA FK5 456
alfa = hms2rad([12,16,31.815]) #alfa, rektascensja hms
dek = dms2rad([56,54,27.780]) #deklinacja dms

#warszawa
f_w = np.deg2rad(52)
l_w = np.deg2rad(21)

#rownik
f_r = np.deg2rad(0)
l_r = np.deg2rad(21)

hours = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
A_w = []
A_r = []
H_w = []
Z_w = []


dane_w = []
dane_r = []

for hour in hours: #warszawa
    utc = UTC2toUTC0(hour)
    lst = LST(2022,7,1,hour,l_w,alfa)
    if utc < 0:
        lst = LST(2022,6,30,hour,l_w,alfa)
    t = kh(lst,alfa)
    azymut, h = A(f_w, dek, t)
    r = 1
    z = m.acos(m.sin(f_w)*m.sin(dek) + m.cos(f_w) * m.cos(dek) * m.cos(t))
    gx = r * np.sin(z) * np.cos(azymut)
    gy = r * np.sin(z) * np.sin(azymut)
    gz = r * np.cos(z)
    A_w.append(azymut)
    H_w.append(90-np.rad2deg(h))
    #azymut = np.rad2deg(azymut)
     
    dane_w.append({'x': gx, 'y': gy, 'z': gz, 'azymut ': azymut, 't [hours]': hour, 'h [rad]': h})
    
dane_w = pd.DataFrame(dane_w)

#plt.scatter(A_w, H_w)
#plt.show()

#figury
fig = px.scatter_3d(dane_w, x='x', y = 'y', z='z', title = 'Ruch gwiazdy w Warszawie')
#fig.show()
fig_az_od_h = px.scatter(dane_w, x='azymut [rad]', y='h [rad]', title = 'Wykres zależności wysokości gwiazdy od azymutu')
#fig_az_od_h.show()
fig_az_od_t = px.scatter(dane_w, x='t [hours]', y = 'azymut [rad]', title = 'Wykres zależności azymutu od czasu dla Warszawy')
#fig_az_od_t.show()
fig_wys_od_t = px.scatter(dane_w, x='t [hours]', y = 'h [rad]', title = 'Wykres zależności wysokości od czasu dla Warszawy')
#fig_wys_od_t.show()

f_r = np.deg2rad(1)
A_r = []
H_r = []

for hour in hours: #rownik
    utc = UTC2toUTC0(hour)
    lst = LST(2022,7,1,hour,l_r,alfa)
##    if utc < 0:
##        lst = LST(2022,6,30,utc,l_r,alfa)
    t = kh(lst,alfa)
    azymut, h = A(f_r, dek, t)

    z = m.acos(m.sin(f_r)*m.sin(dek) + m.cos(f_r) * m.cos(dek) * m.cos(t))
    
    gx = r * np.sin(z) * np.cos(azymut)
    gy = r * np.sin(z) * np.sin(azymut)
    gz = r * np.cos(z)

    A_r.append(azymut)
    H_r.append(90-np.rad2deg(h))
    
    dane_r.append({'x': gx, 'y': gy, 'z': gz, 'azymut [rad]': azymut, 't [hours]': hour, 'h [rad]': h})
    
dane_r = pd.DataFrame(dane_r)

fig2 = px.scatter_3d(dane_r, x='x', y = 'y', z='z',  title = 'Ruch gwiazdy na równiku')
#fig2.show()
fig_az_od_h = px.scatter(dane_r, x='azymut [rad]', y='h [rad]', title = 'Wykres zależności wysokości gwiazdy od azymutu')
#fig_az_od_h.show()
fig_az_od_t = px.scatter(dane_r, x='t [hours]', y = 'azymut [rad]', title = 'Wykres zależności azymutu od czasu dla równika')
#fig_az_od_t.show()
fig_wys_od_t = px.scatter(dane_r, x='t [hours]', y = 'h [rad]', title = 'Wykres zależności wysokości od czasu dla równika')
#fig_wys_od_t.show()

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(polar = True)
ax.set_theta_zero_location('N') # ustawienie kierunku północy na górze wykresu
ax.set_theta_direction(-1)

ax.set_yticks(range(0, 90+10, 10)) # Define the yticks

yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
ax.set_yticklabels(yLabel)
ax.set_rlim(0,90)

xs = A_w
ys = H_w
ax.scatter(xs, ys, label='Warszawa')
xs = A_r
ys = H_r
ax.scatter(xs, ys, label='Równik')
ax.legend()
plt.show()



##fig = plt.figure(figsize = (10,10))
##ax = fig.add_subplot(projection = '3d')
### promien´ Ziemi
##r = 1
### siatka wspołrz˛ednych
##u, v = np.mgrid[0:(2 * np.pi+0.1):0.1, 0:np.pi:0.1]
##x = np.cos(u) * np.sin(v)
##y = np.sin(u) * np.sin(v)
##z = np.cos(v)
##z[z<0] = 0 # bez tego, narysowalibys´my cał ˛a kul˛e, a chcemy tylko półkul˛e
##ax.plot_surface(x,y,z, alpha = 0.1)
##
##
### narysowanie punktu na sferze
##gx = r * np.sin(A_w) * np.cos(H_w)
##gy = r * np.cos(A_w) * np.cos(H_w)
##gz = r * np.sin(H_w)
##ax.plot3D(gx,gy,gz)

##plt.show()
