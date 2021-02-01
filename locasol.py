import numpy as np
import datetime as dat
from matplotlib import pyplot as plt


def jd_calc(year, month, day, hour):                     #'jd_calc' calculates the current julian date (jd), given the current year, month, day and hour.
    delta=year-1949                                      #delta is the years passed since 1949 (including the current year)
    leap=int(delta/4)                                    #leap is the number of leap days that there have been since 1949
    day2=dat.date(year,month,day).timetuple().tm_yday    #'day2' represents the number of days thatb have passed since the beginning of the current year 
    jd= 2432916.5 + delta*365 + leap + day2 + (hour/24)  #'jd' is the current julian date  
    return jd





def ecliptic_coord(jd):                                                  #'ecliptic_coord' calculates the sun's position in the ecliptic coordinate system 
    n=jd-2451545.0                                                       #'n' is the difference between the current julian date and the julian date of January 1 2000 at 12:00pm UTC                                
    L=280.460 + 0.9856474*n                                              
    L=L - (L/360).astype(int)*360                                        #'L' is the mean longitude in degs         
    g=357.528 + 0.9856003*n 
    g=g - (g/360).astype(int)*360                                        #'g' is the mean anomaly in degs 
    l=L + 1.915*np.sin(np.deg2rad(g)) + 0.020*np.sin(2*np.deg2rad(g))
    l=l - (l/360).astype(int)*360                                        #'l' is the ecliptic longidute in degs 
    ep=23.439-0.0000004*n                                                #'ep' is the obliquity of the ecliptic in degs
    eclipcoord=np.array([n,L,g,l,ep])                                    #'eclipcoord' is an array containing all the important variables calculated in 'ecliptic_coord'
    return eclipcoord





def equatorial_coord(eclipcoord):                                  #'equatorial_coord' converts the sun's position from ecliptic coordinates to equatorial coordinates 
    
    if (eclipcoord.ndim==2):                                       #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
        epp=np.deg2rad(eclipcoord[4,:])
        ll=np.deg2rad(eclipcoord[3,:])
    else:
        epp=np.deg2rad(eclipcoord[4])
        ll=np.deg2rad(eclipcoord[3])
    
    ra=np.arctan2(np.cos(epp)*np.sin(ll),np.cos(ll))%(2*np.pi)
    dec=np.arcsin(np.sin(epp)*np.sin(ll))
    ra=np.rad2deg(ra)                                               #'ra' is the right ascension in degs from 0 to 360 (being eastward the positive direction)
    dec=np.rad2deg(dec)                                             #'dec' is the declination in degs from -90 to 90 (being nordwards the positive direction)
    equacoord=np.array([ra, dec])
    return equacoord                                                #'equacoord' is an array containing all the important variables calculated in 'equatorial_coord'





def east_longitude(Longitude):                                                                  # 'east_longitude' takes the longitude coordinate of a place on earth (in W/O, degs, min and secs) and converts it to the... 
    if (Longitude[0]=='O'):                                                                     #... to east longitude format (i.e from 0 to 360) in degrrees. 
        eastlong=360-float(Longitude[1]) - float(Longitude[2])/60 - float(Longitude[3])/3600
    else:
        eastlong=float(Longitude[1]) + float(Longitude[2])/60 + float(Longitude[3])/3600
    return eastlong






def latitude2(Latitude):                                                             #'Latitude2' takes the longitude coordinate of a place on earth (in N/S, degs, min and secs) and converts it to the...     
    if (Latitude[0]=='S'):                                                           #to degrees (i.e from -90 to 90)
        lat2=-float(Latitude[1]) - float(Latitude[2])/60 - float(Latitude[3])/3600
    else:
        lat2=float(Latitude[1]) + float(Latitude[2])/60 + float(Latitude[3])/3600
    return lat2
    





def local_coord_1(equacoord, eclipcoord, hour, eastlong):   #'local_coord_1' calculates the first round of necessary local coordinates from the equatorial coordinate system.
    
    if (eclipcoord.ndim==2):                                #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
        gmst=6.697375+0.0657098242*eclipcoord[0,:]+hour    
    else:
        gmst=6.697375+0.0657098242*eclipcoord[0]+hour
    
    if (eclipcoord.ndim==2):                                #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
        gmst=gmst - (gmst/24).astype(int)*24
        lmst=gmst + eastlong/15
        lmst=lmst - (lmst/24).astype(int)*24                #'gmst' is greenwhich mean siderial time in hours (1h= 15 degs)
    else:                                                   
        gmst=gmst - int(gmst/24)*24   
        lmst=gmst + eastlong/15
        lmst=lmst - int(lmst/24)*24                         #'lms' is the local siderial time in hours (1h= 15 degs)
    
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
            ha=ha+360                                                            #'ha' is the local hour angle in degs (from -90 to 90)
    
    loccoord1=np.array([ha,lmst,gmst])        
    
    return loccoord1                                                             #'loccoord1' is an array containing all the important variables calculated in 'local_coord_1'
    







def local_coord_2(loccord1, lat2, eclipcoord, equacoord):                        #'local_coord_2' calculates the second round of necessary local coordinates from the equatorial coordinate system.  
    
    if (eclipcoord.ndim==2):                                                     #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
        dec=np.deg2rad(equacoord[1,:])
        ha=np.deg2rad(loccord1[0,:])
    else:
        dec=np.deg2rad(equacoord[1])
        ha=np.deg2rad(loccord1[0])
        
    lat=np.deg2rad(lat2)                                                          #'lat' is just the latitude in radians 

    
    el=np.arcsin(np.sin(dec)*np.sin(lat) + np.cos(dec)*np.cos(lat)*np.cos(ha))    #'el' is the elevation degs (-90 to 90)  
    az=np.arcsin(-np.cos(dec)*np.sin(ha)/np.cos(el))                              
    elc=np.arcsin(np.sin(lat)/np.sin(dec))                                        #This formula was changed with respect to the Almanac. Otherwise it wouldn't have been right
    
    el=np.rad2deg(el)
    az=np.rad2deg(az)
    elc=np.rad2deg(elc)
    
    
    if (eclipcoord.ndim==2):                                                     #The first if-else statement is only here so the function can handle, both: single valued inputs and arrays. 
        
        for i in range(0,len(el)):
            
            if(el[i]%360>=elc[i]):                                               #Also, a little extra tampering with this part of the algorithim was necessary to make it work.
                az[i]=180-az[i]
            elif(el[i]%360<=elc[i] and ha[i]>0): 
                az[i]=az[i]+360
                
    else:
        if(el%360>=elc):
            az=180-az
        elif (el%360<=elc and ha>0):
            az=az+360
                                                                    
                                                                                  #'az' is the azimuth angle in degs (0 to 360)    
    loccoord2=np.array([el,az,elc])
    
    return loccoord2
    
    
    
def refract1(loccoord2):
    if (eclipcoord.ndim==2):
        el=np.deg2rad(loccoord2[0,:])
        correct=3.51561*(0.1594 + 0.0196*el + 0.00002*el**2)/(1 + 0.505*el + 0.0845*el**2)
        correct=correct/60
        el_correct=np.rad2deg(el)+correct
        
        return el_correct
        
    else:
        el=el=np.deg2rad(loccoord2[0])
        correct=3.51561*(0.1594 + 0.0196*el + 0.00002*el**2)/(1 + 0.505*el + 0.0845*el**2)
        correct=correct/60
        el_correct=np.rad2deg(el)+correct
        
        return el_correct
    
    
def refract2(loccoord2):
    if (eclipcoord.ndim==2):
        el=loccoord2[0,:]
        correct=1.02/np.tan(np.deg2rad(el + 10.3/(el+5.11)))
        correct=correct/60
        el_correct=el+correct
        
        return el_correct
    
    else:
        el=loccoord2[0]
        correct=1.02/np.tan(np.deg2rad(el + 10.3/(el+5.11)))
        correct=correct/60
        el_correct=el+correct
        
        return el_correct        


def refract3(loccoord2):
    if (eclipcoord.ndim==2):
        el=loccoord2[0,:]
        h=np.deg2rad(el)
        correct=np.zeros(len(h))
        
        for i in range(0, len(loccoord2[0,:])):
            if(el[i]<-0.575):
                correct[i]=(1/3600)*(-20.774/np.tan(h[i]))
            elif(el[i]>=-0.575 and el[i]<5):
                correct[i]=(1/3600)*(1735 - 518.2*el[i] + 103.4*el[i]**2 - 12.79*el[i]**3 + 0.711*el[i]**4)
            elif(el[i]>=5 and el[i]<85):
                correct[i]=(1/3600)*(58.1/np.tan(h[i]) - 0.07/np.tan(h[i])**3 + 0.000086/np.tan(h[i])**5)
            elif(el[i]>=85 and el[i]<=90):
                correct[i]=0
        
        el_correct=el+correct    
        
        return el_correct
            
    else:
        el=loccoord2[0]
        h=np.deg2rad(el)
        
        if(el<-0.575):
            correct=(1/3600)*(-20.774/np.tan(h))
        elif(el>=-0.575 and el<5):
            correct=(1/3600)*(1735 - 518.2*el + 103.4*el**2 - 12.79*el**3 + 0.711*el**4)
        elif(el>=5 and el<85):
            correct=(1/3600)*(58.1/np.tan(h) - 0.07/np.tan(h)**3 + 0.000086/np.tan(h)**5)
        elif(el>=85 and el<=90):
            correct=0
            
        el_correct=el+correct           
    
        return el_correct        
        
        
    
    
    
 #-------------------------------------------------------------------------------------------------------------------   
    
    
show=3
UTC=-5                                      #Timezone according to universal time
Longitude=np.array(['O',75,33,48.92])       #Longitude of the geographic point in question
Latitude=np.array(['N',6,13,1])             #Latitude of the geographic point in question  
Hi=6                                      #Intial hour
Hf=8                                     #Final hour
Hour=np.linspace(Hi,Hf,int(Hf-Hi)*10+1)-UTC   #Hour escalar/vector
#Hour=6+19/60-UTC
Date=np.array([2021,2,1])                  #Date wanted




jd=jd_calc(Date[0],Date[1],Date[2],Hour)
eclipcoord=ecliptic_coord(jd)
equacoord=equatorial_coord(eclipcoord)
eastlong=east_longitude(Longitude)
lat2=latitude2(Latitude)
loccord1=local_coord_1(equacoord, eclipcoord, Hour, eastlong)
loccoord2=local_coord_2(loccord1, lat2, eclipcoord, equacoord)
el_correct=refract3(loccoord2)
print(loccoord2[0],el_correct,loccoord2[1])


if eclipcoord.ndim==2:

    if show==1:                                           #Shows graph of elevation vs hour
        
        plt.plot(Hour+UTC, loccoord2[0])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(loccoord2[0]),max(loccoord2[0])+1)
        plt.title('Hour vs true Sun elevation')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('Elevation (°)')
    
    elif show==2:                                          #Shows graph of azimuth vs hour
    
        plt.plot(Hour+UTC, loccoord2[1])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(loccoord2[1]),max(loccoord2[1])+1)
        plt.title('Hour vs Sun azimuth')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('Azimuth (°)')
        
    elif show==3: 
        
        plt.plot(loccoord2[0], el_correct)
        plt.plot(loccoord2[0], loccoord2[0])
        plt.xlim(min(loccoord2[0]),max(loccoord2[0]))
        plt.ylim(min(el_correct),max(el_correct))
        plt.title('el vs el_correct')
        plt.xlabel('el (°)')
        plt.ylabel('el_correct (°)')
        
    elif show==4: 
        Refract=el_correct-loccoord2[0]
        plt.plot(loccoord2[0], Refract)
        plt.xlim(min(loccoord2[0]),max(loccoord2[0]))
        plt.ylim(min(Refract),max(Refract))
        plt.title('el vs Refract')
        plt.xlabel('el (°)')
        plt.ylabel('Refract (°)')
            

