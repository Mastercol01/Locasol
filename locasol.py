import numpy as np
import datetime as dat
import math
from matplotlib import pyplot as plt


def jd_diff(year, month, day, hour):
    delta=year-1949
    leap=int(delta/4)
    day2=dat.date(year,month,day).timetuple().tm_yday
    jd= 2432916.5 + delta*365 + leap + day2 + (hour/24)
    return jd

def ecliptic_coord(jd):
    n=jd-2451545.0
    L=280.460 + 0.9856474*n 
    L=L - (L/360).astype(int)*360        #(mean longitude in degs)
    g=357.528 + 0.9856003*n 
    g=g - (g/360).astype(int)*360        #(mean anomaly in degs)  
    l=L + 1.915*np.sin(np.deg2rad(g)) + 0.020*np.sin(2*np.deg2rad(g))
    l=l - (l/360).astype(int)*360        #(ecliptic longidute in degs) 
    ep=23.439-0.0000004*n                #(obliquity of the ecliptic in degs) 
    eclipcoord=np.array([n,L,g,l,ep])
    return eclipcoord

def equatorial_coord(eclipcoord):
    
    if (eclipcoord.ndim==2):
        epp=np.deg2rad(eclipcoord[4,:])
        ll=np.deg2rad(eclipcoord[3,:])
    else:
        epp=np.deg2rad(eclipcoord[4])
        ll=np.deg2rad(eclipcoord[3])
    
    ra=np.arctan2(np.cos(epp)*np.sin(ll),np.cos(ll))%(2*np.pi)
    dec=np.arcsin(np.sin(epp)*np.sin(ll))
    ra=np.rad2deg(ra)
    dec=np.rad2deg(dec)
    equacoord=np.array([ra, dec])
    return equacoord


def east_longitude(Longitude):
    if (Longitude[0]=='O'):
        eastlong=360-float(Longitude[1]) - float(Longitude[2])/60 - float(Longitude[3])/3600
    else:
        eastlong=float(Longitude[1]) + float(Longitude[2])/60 + float(Longitude[3])/3600
    return eastlong


def latitude2(Latitude):
    if (Latitude[0]=='S'):
        lat2=-float(Latitude[1]) - float(Latitude[2])/60 - float(Latitude[3])/3600
    else:
        lat2=float(Latitude[1]) + float(Latitude[2])/60 + float(Latitude[3])/3600
    return lat2
    


def local_coord_1(equacoord, eclipcoord, hour, eastlong):
    
    if (eclipcoord.ndim==2):
        gmst=6.697375+0.0657098242*eclipcoord[0,:]+hour
    else:
        gmst=6.697375+0.0657098242*eclipcoord[0]+hour
    
    gmst=gmst - (gmst/24).astype(int)*24
    
    lmst=gmst + eastlong/15
    lmst=lmst - (lmst/24).astype(int)*24
    
    ha=lmst*15-equacoord[0] 
    
    if (eclipcoord.ndim==2):
        
        for i in range(0,len(lmst)):
            
            if(i<=int(len(lmst)/2) and ha[i]>=0):
                ha[i]=ha[i]-360
            elif(i>int(len(lmst)/2) and ha[i]<=0): 
                ha[i]=ha[i]+360
    else:
        if(Hour<=12 and Hour>=0):
            ha=ha-360
        elif (Hour>12 and ha<=0):
            ha=ha+360
    
    loccoord1=np.array([ha,lmst,gmst])        
    
    return loccoord1        
    

def local_coord_2(loccord1, lat2, eclipcoord, equacoord):
    
    if (eclipcoord.ndim==2):
        dec=np.deg2rad(equacoord[1,:])
        ha=np.deg2rad(loccord1[0,:])
    else:
        dec=np.deg2rad(equacoord[1])
        ha=np.deg2rad(loccord1[0])
        
    lat=np.deg2rad(lat2)

    
    el=np.arcsin(np.sin(dec)*np.sin(lat) + np.cos(dec)*np.cos(lat)*np.cos(ha))
    az=np.arcsin(-np.cos(dec)*np.sin(ha)/np.cos(el))
    elc=np.arcsin(np.sin(lat)/np.sin(dec))
    
    el=np.rad2deg(el)
    az=np.rad2deg(az)
    elc=np.rad2deg(elc)
    
    
    if (eclipcoord.ndim==2):
        
        for i in range(0,len(el)):
            
            if(el[i]%360>=elc[i]):
                az[i]=180-az[i]
            elif(el[i]%360<=elc[i] and ha[i]>0): 
                az[i]=az[i]+360
                
    else:
        if(el%360>=elc):
            az=az-180
        elif (el%360<=elc and ha>0):
            az=az+360
    

    
    loccoord2=np.array([el,az,elc])
    
    return loccoord2
    
    
    
    
 #-------------------------------------------------------------------------------------------------------------------   
    
    
show=2
UTC=-5
Longitude=np.array(['O',75,33,48.92])
Latitude=np.array(['N',6,13,1])
Hi=6
Hf=18
Hour=np.linspace(Hi,Hf,(Hf-Hi)*1+1)-UTC
Date=np.array([2021,1,28])




jd=jd_diff(Date[0],Date[1],Date[2],Hour)
eclipcoord=ecliptic_coord(jd)
equacoord=equatorial_coord(eclipcoord)
eastlong=east_longitude(Longitude)
lat2=latitude2(Latitude)
loccord1=local_coord_1(equacoord, eclipcoord, Hour, eastlong)
loccoord2=local_coord_2(loccord1, lat2, eclipcoord, equacoord)


if show==1:
    
    plt.plot(Hour+UTC, loccoord2[0])
    plt.xlim(min(Hour+UTC),max(Hour+UTC))
    plt.ylim(min(loccoord2[0]),max(loccoord2[0])+1)
    plt.title('Hour vs Sun elevation')
    plt.xlabel('Hour (24h format)')
    plt.ylabel('Elevation (°)')

elif show==2:

    plt.plot(Hour+UTC, loccoord2[1])
    plt.xlim(min(Hour+UTC),max(Hour+UTC))
    plt.ylim(min(loccoord2[1]),max(loccoord2[1])+1)
    plt.title('Hour vs Sun azimuth')
    plt.xlabel('Hour (24h format)')
    plt.ylabel('Azimuth (°)')

