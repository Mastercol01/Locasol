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
    ep=23.439-0.0000004*n    #(obliquity of the ecliptic in degs) 
    eclipcoord=np.array([n,L,g,l,ep])
    return eclipcoord

z=np.linspace(0,24,25)
jd=jd_diff(2021,1,27,z)
eclipcoord=ecliptic_coord(jd)