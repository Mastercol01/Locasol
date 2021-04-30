import numpy as np
import datetime as dat
from matplotlib import pyplot as plt


#----------------------------------------CODE DESCRIPTION----------------------------------------------------------------
#This code allows you to calculate various parameters regarding the sun's position on the sky, as the day elapses.
#Given the geographical coordinates of the location for which the sun's trajectory in the sky is wished to be calculated
#and the date, along with the hour or range of hours, for which it is wished to be calculated, this code returns the approximate value of parameters such as:
#the sun's Hour Right Ascension, declination, Azimuth, Elevation/Altitude and Angle of incidence on a solar array, among other things.
#This code was primarly based on the Astronomical Almanac's Algorithm for approximate Solar position (1950-2050). However some parts of the algorithm had to be changed, since
#they were wrong or did not work. Having said that, enjoy.

#---------------------------------------PARAMETER DEFINITION-------------------------------------------------------------


#ENTER THE LOCATION'S NAME
LocaName='Medellín, Colombia'                                                  #LocalName: 'Descriptive name of the place for which calculations are wished to be performed'

#ENTER THE LOCATION'S GEOGRAPHIC COORDINATES
Longitude=np.array(['W',75,34,1])                                              #Longitude: [West/East, degrees, arcminutes, arcseconds]
Latitude=np.array(['N',6,13,1])                                                #Latitude: [North/South, degrees, arcminutes, arcseconds]

#ENTER THE LOCATIONS'S UTC TIMEZONE
UTC=-5                                                                         #UTC: timezone according to universal time

#ENTER THE LOCAL DATE 
Date=np.array([2021,4,30])                                                     #Date: [year, No. of the month, day of the month]  

#ENTER THE LOCAL TIME IN 24H FORMAT  
#(It can be either a single value or a 2D array)
Hour1D= 5 + (49/60) + (0/3600)                                                  #Hour1D: hour + (minutes/60) + (seconds/3600)
Hour2D=np.linspace(6,18,60)                                                    #Hour2D: np.linspace(initial hour, final hour, number of partitions )
HourDim=1                                                                      #HourDim=1 -> Hour1D becomes the input
                                                                               #HourDim=2 -> Hour2D becomes the input
                                                                                         
#ENTER THE SOLAR ARRAY'S RELATIVE POSITIONING  (Optional)
#(It can be either a single value or a 2D array, but, in the latter case, it must have the same dimensions as Hour2D)
slope1D=10                                                                     #Slope: Is the angle that the surface is tilted from the horizontal (deg)
saz1D=10                                                                       #Saz: Surface Azimuth of the solar panel. It is defined similarly as the azimuth for the sun (deg)          
slope2D=10*np.ones(len(Hour2D))                                        
saz2D=10*np.ones(len(Hour2D))                                                  #The 1D and 2D variants function exactly the same as explained in the section above for 'Hour'                                                  
slopesazDim=1
      
                                                
#CHOOSE WHICH GRAPH YOU WHISH TO SEE  (Applies only for Hour2D)  

show=5                  # show=-1 -> Nothing
                        # show=0 -> Hour vs true elevation
                        # show=1 -> Hour vs azimuth 
                        # show=2 -> true elevation vs refration-corrected elevation
                        # show=3 -> true elevation vs refration-correction term
                        # show=4 -> Hour vs corrected elevation 
                        # show=5 -> Angle of incidence on the solar panel 
                        # show=6 -> Hour angle vs Hour
                        # show=7 -> LMST vs Hour
                        # show=8 -> Angle of GMST vs Hour
                        # show=9 -> Right Ascension vs Hour
                        # show=10 -> Declination vs Hour
                                             
                        
#----------------------------------------DEFINITION OF MAIN FUNCTIONS-------------------------------------------------------------------


def jd_calc(year, month, day, hour):                         #'jd_calc' calculates the current julian date (jd), given the current year, month, day and hour.
    delta=year-1949                                          #delta is the years passed since 1949 (including the current year)
    leap=int(delta/4)                                        #leap is the number of leap days that there have been since 1949
    day2=dat.date(year,month,day).timetuple().tm_yday        #'day2' represents the number of days that have passed since the beginning of the current year 
    jd= 2432916.5 + delta*365 + leap + day2 + (hour/24)      #'jd' is the current julian date  
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
    if (Longitude[0]=='W'):                                                                     #... to east longitude format (i.e from 0 to 360) in degrrees. 
        eastlong=360-float(Longitude[1]) - float(Longitude[2])/60 - float(Longitude[3])/3600
    else:
        eastlong=float(Longitude[1]) + float(Longitude[2])/60 + float(Longitude[3])/3600
    return eastlong






def latitude2(Latitude):                                                                #'Latitude2' takes the longitude coordinate of a place on earth (in N/S, degs, min and secs) and converts it to the...     
    if (Latitude[0]=='S'):                                                              #to degrees (i.e from -90 to 90)
        lat2=-float(Latitude[1]) - float(Latitude[2])/60 - float(Latitude[3])/3600
    else:
        lat2=float(Latitude[1]) + float(Latitude[2])/60 + float(Latitude[3])/3600
    return lat2
    





def local_coord_1(equacoord, eclipcoord, hour, eastlong):                        #'local_coord_1' calculates the first round of necessary local coordinates from the equatorial coordinate system.
    
    if (eclipcoord.ndim==2):                                                     #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
        gmst=6.697375+0.0657098242*eclipcoord[0,:]+hour    
    else:
        gmst=6.697375+0.0657098242*eclipcoord[0]+hour
    
    if (eclipcoord.ndim==2):                                                      #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
        gmst=gmst - (gmst/24).astype(int)*24
        lmst=gmst + eastlong/15
        lmst=lmst - (lmst/24).astype(int)*24                                      #'gmst' is greenwhich mean siderial time in hours (1h= 15 degs)
    else:                                                   
        gmst=gmst - int(gmst/24)*24   
        lmst=gmst + eastlong/15
        lmst=lmst - int(lmst/24)*24                                               #'lmst' is the local siderial time in hours (1h= 15 degs)
    
    ha=gmst*15+eastlong-equacoord[0] 
        
    if (eclipcoord.ndim==2):
        
        ha=ha-(ha/360).astype(int)*360                                             #This part for the ha had to be included since it wasn't ther explicitly in the original algrithm
        for i in range(0, len(ha)):
                       if(ha[i]>=180):
                           ha[i]=ha[i]-360
    else: 
        ha=ha-int(ha/360)*360
        if(ha>=180):
            ha=ha-360
                                                                                   #'ha' is the local hour angle in degs (from -180 to 180)
    
    loccoord1=np.array([ha,lmst,gmst])        
    
    return loccoord1                                                               #'loccoord1' is an array containing all the important variables calculated in 'local_coord_1'
    







def local_coord_2(loccord1, lat2, eclipcoord, equacoord):                          #'local_coord_2' calculates the second round of necessary local coordinates from the equatorial coordinate system.  
    
    if (eclipcoord.ndim==2):                                                       #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
        dec=np.deg2rad(equacoord[1,:])
        ha=np.deg2rad(loccord1[0,:])
    else:
        dec=np.deg2rad(equacoord[1])
        ha=np.deg2rad(loccord1[0])
        
    lat=np.deg2rad(lat2)                                                            #'lat' is just the latitude in radians 

    
    el=np.arcsin(np.sin(dec)*np.sin(lat) + np.cos(dec)*np.cos(lat)*np.cos(ha))      #'el' is the elevation degs (-90 to 90)  
    elc=np.arcsin(np.sin(lat)/np.sin(dec))                                          #The formula for calculating az was completely changed, since the algorithm's formula was just wrong: pveducation.org/pvcdrom/properties-of-sunlight/azimuth-angle
    cos_az=(np.sin(dec)*np.cos(lat)-np.cos(dec)*np.sin(lat)*np.cos(ha))/np.cos(el)
    az=np.arccos(cos_az)    
    
    
    el=np.rad2deg(el)                                                                       
    az=np.rad2deg(az)
    elc=np.rad2deg(elc)
    ha=np.rad2deg(ha)
    
    
    if (eclipcoord.ndim==2):                                                          #The first if-else statement is only here so the function can handle, both: single valued inputs and arrays. 
        
        for i in range(0,len(el)):
            
            if(ha[i]<0):                                              
                az[i]=az[i]
            if(ha[i]>0): 
                az[i]=360-az[i]
                
    else:
        if(ha<0):
            az=az
        if (ha>0):
            az=360-az
                                                                    
                                                                                       #'az' is the azimuth angle in degs (0 to 360)    
    loccoord2=np.array([el,az,elc])
    
    return loccoord2
    
    
    
def refract1(loccoord2):                                                               #refract1 calculates a correction term for atmospheric refraction, based on the formula given by the The Astronomical Almanac.    
    if (eclipcoord.ndim==2):                                                           #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
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
    
    
def refract2(loccoord2):                                                               #refract2 calculates a correction term for atmospheric refraction, based on the formula given by Saemundsson. 
    if (eclipcoord.ndim==2):                                                           #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.
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


def refract3(loccoord2):                                                                  #refract3 calculates a correction term for atmospheric refraction, based on the formulas used by NOAA. 
    if (eclipcoord.ndim==2):                                                              #Empirically speaking, this is the best correction function so far.
        el=loccoord2[0,:]                                                                 #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.    
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
        
    
    
def AOI (loccoord2, el_correct, slope, saz):                       #AOI calculates the angle of incidence of solar rays onto the solar panel                            
                                                                   #This if-else statement is only here so the function can handle, both: single valued inputs and arrays.    
    if (loccoord2.ndim==2):     
        az=loccoord2[1,:]  
    else:
        az=loccoord2[1]  
              
    el_correct=np.deg2rad(el_correct)
    az=np.deg2rad(az)
    slope=np.deg2rad(slope)
    saz=np.deg2rad(saz)
    
    cos_aoi= np.sin(el_correct)*np.cos(slope) + np.cos(el_correct)*np.sin(slope)*np.cos(az-saz)
    
    aoi=np.rad2deg(np.arccos(cos_aoi))
    
    return aoi


    
    
    
    
 #------------------------------------CODE EXECUTION AND CALCULATION OF VARIABLES---------------------------------------------------   
    
if (HourDim==1):
    Hour=Hour1D
else:
    Hour=Hour2D
Hour=Hour-UTC

if (slopesazDim==1):
    slope=slope1D
    saz=saz1D
else:
    slope=slope2D
    saz=saz2D




jd=jd_calc(Date[0],Date[1],Date[2],Hour)
eclipcoord=ecliptic_coord(jd)
equacoord=equatorial_coord(eclipcoord)
eastlong=east_longitude(Longitude)
lat2=latitude2(Latitude)
loccord1=local_coord_1(equacoord, eclipcoord, Hour, eastlong)
loccoord2=local_coord_2(loccord1, lat2, eclipcoord, equacoord)
el_correct=refract3(loccoord2)
aoi=AOI(loccoord2, el_correct, slope, saz)


#-------------------------------------------RESULTS DISPLAY TEXT-----------------------------------------

print('----------------------------------------------------------------------------')
print('                           RESULTS DISPLAY')
print('----------------------------------------------------------------------------')

print('LOCATION:', LocaName )
print('LONGITUDE:', Longitude[0] , Longitude[1],'°', Longitude[2],'\'', Longitude[3],'\"', '             |', 'EAST LONGITUDE:', eastlong, '°'  )
print('LATITUDE:', Latitude[0] , Latitude[1],'°', Latitude[2],'\'', Latitude[3],'\"', '               |', 'LATITUDE:', lat2, '°'   )
print('TIMEZONE:', 'UTC', UTC)
print('DATE:', Date[2],'/',Date[1],'/',Date[0] )

if(HourDim==1):
    print('JULIAN DATE:',jd)
    print('LOCAL TIME:', int(Hour+UTC),'h' , int(Hour+UTC-int(Hour+UTC))*60,'min', ((Hour+UTC-int(Hour+UTC))- int(Hour+UTC-int(Hour+UTC)))*60,'sec')  
    print('----------------------------------------------------------------------------')
    print('HA:', loccord1[0], '°', '                |','HA:', int(loccord1[0]/15), 'h', int((loccord1[0]/15-int(loccord1[0]/15))*60), '\'', ((loccord1[0]/15-int(loccord1[0]/15))*60-int((loccord1[0]/15-int(loccord1[0]/15))*60))*60, '\"'   ) 
    print('LMST:', loccord1[1]*15, '°', '               |','LMST:', int(loccord1[1]), 'h', int((loccord1[1]-int(loccord1[1]))*60), '\'', ((loccord1[1]-int(loccord1[1]))*60-int((loccord1[1]-int(loccord1[1]))*60))*60, '\"'   ) 
    print('GMST', loccord1[2]*15, '°', '               |','GMST:', int(loccord1[2]), 'h', int((loccord1[2]-int(loccord1[2]))*60), '\'', ((loccord1[2]-int(loccord1[2]))*60-int((loccord1[2]-int(loccord1[2]))*60))*60, '\"'    ) 
    print('----------------------------------------------------------------------------')
    print('RIGHT ASCENSION:', equacoord[0], '°', '   |','RIGHT ASCENSION:', int(equacoord[0]/15), 'h', int((equacoord[0]/15-int(equacoord[0]/15))*60), '\'', ((equacoord[0]/15-int(equacoord[0]/15))*60-int((equacoord[0]/15-int(equacoord[0]/15))*60))*60, '\"'    ) 
    print('DECLINATION:', equacoord[1], '°' ) 
    print('AZIMUTH:', loccoord2[1], '°' ) 
    print('TRUE ELEVATION:', loccoord2[0], '°' ) 
    print('REFRACTION-CORRECTED ELEVATION:', el_correct, '°' ) 
    print('----------------------------------------------------------------------------')
    print('SLOPE:', slope, '°')  
    print('SURFACE AZIMUTH:', saz, '°' )  
    print('ANGLE OF INCIDENCE ON THE SOLAR PANEL:', aoi, '°' )
print('----------------------------------------------------------------------------')

#-------------------------------------------RESULTS DISPLAY GRAPHS-----------------------------------------

if eclipcoord.ndim==2:

    if show==0:                                           #Shows graph of elevation vs hour
        
        plt.plot(Hour+UTC, loccoord2[0])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(loccoord2[0]),max(loccoord2[0])+1)
        plt.title('True elevation vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('True elevation (°)')
    
    elif show==1:                                          #Shows graph of azimuth vs hour
    
        plt.plot(Hour+UTC, loccoord2[1])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(loccoord2[1]),max(loccoord2[1])+1)
        plt.title('Azimuth vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('Azimuth (°)')
        
    elif show==2:                                          #Refraction-corrected elevation vs True Elevation
        
        plt.plot(loccoord2[0], el_correct)
        plt.plot(loccoord2[0], loccoord2[0])
        plt.xlim(min(loccoord2[0]),max(loccoord2[0]))
        plt.ylim(min(el_correct),max(el_correct))
        plt.title('Refraction-corrected elevation vs True Elevation')
        plt.xlabel('el (°)')
        plt.ylabel('el_correct (°)')
        
    elif show==3:                                          #Shows graph of Refration-correction term vs True Elevation
        Refract=el_correct-loccoord2[0]
        plt.plot(loccoord2[0], Refract)
        plt.xlim(min(loccoord2[0]),max(loccoord2[0]))
        plt.ylim(min(Refract),max(Refract))
        plt.title(' Refration-correction term vs True Elevation ')
        plt.xlabel('True Elevation (°)')
        plt.ylabel('Refration-correction term (°)')
        
    elif show==4:                                           #Shows graph of Corrected sun elevation vs Hour
        plt.plot(Hour+UTC, el_correct)
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(el_correct),max(el_correct)+1)
        plt.title('Corrected sun elevation vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('Corrected Elevation (°)')
        
    elif show==5:                                           #Shows graph of Angle of incidence on the solar panel
        plt.plot(Hour+UTC, aoi)
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(aoi),max(aoi)+1)
        plt.title('Angle of incidence on the solar panel vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('Angle of incidence on the solar panel (°)')
        
    elif show==6:                                           #Shows graph of Hour angle vs Hour
        plt.plot(Hour+UTC, loccord1[0])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(loccord1[0]),max(loccord1[0])+1)
        plt.title('Hour angle vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('Hour angle (°)')
        
    elif show==7:                                           #Shows graph of LMST vs Hour
        plt.plot(Hour+UTC, loccord1[1])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(loccord1[1]),max(loccord1[1])+0.5)
        plt.title('LMST vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('LMST (h)')
        
    elif show==8:                                           #Shows graph of GMST vs Hour
        plt.plot(Hour+UTC, loccord1[2])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(loccord1[2]),max(loccord1[2])+0.5)
        plt.title('GMST vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('GMST (h)')
        
    elif show==9:                                           #Shows graph of Right Ascension vs Hour
        plt.plot(Hour+UTC, equacoord[0])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(equacoord[0]),max(equacoord[0]))
        plt.title('Right Ascension vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('Right Ascencion (°)')
        
    elif show==10:                                           #Shows graph of Declination vs Hour
        plt.plot(Hour+UTC, equacoord[1])
        plt.xlim(min(Hour+UTC),max(Hour+UTC))
        plt.ylim(min(equacoord[1]),max(equacoord[1]))
        plt.title('Declination vs Hour')
        plt.xlabel('Hour (24h format)')
        plt.ylabel('Declination (°)')
        


#-----------------------SOME PARAMETERS DEFINITIONS--------------------------------------------- 

#aoi= angle of incidence
#az= sun's azimuth
#ra= sun's right ascension
#dec= sun's declination
#el= sun's elevation
#elc= critical elevation angle
#el_correct= sun's corrected-elevation
#ha= local hour angle
#gmst= Greenwich mean siderial time
#lmst= Local mean siderial time
#n= jd-2451545.0
#L= mean longitude
#g= mean anomaly
#l= ecliptic longitude
#ep= obliquity of the ecliptic

#FOR 1D INPUTS:

#el=loccoord2[0]
#az=loccoord2[1]
#elc=loccoord2[2]
#ra=equacoord[0]
#dec=equacoord[1]
#ha=loccord1[0]
#lmst=loccord1[1]
#gmst=loccord1[2]
#n=eclipcoord[0]
#L=eclipcoord[1]
#g=eclipcoord[2]
#l=eclipcoord[3]
#ep=eclipcoord[4]


#FOR 2D INPUTS:

#el=loccoord2[0,:]
#az=loccoord2[1,:]
#elc=loccoord2[2,:]
#ra=equacoord[0,:]
#dec=equacoord[1,:]
#ha=loccord1[0,:]
#lmst=loccord1[1,:]
#gmst=loccord1[2,:]
#n=eclipcoord[0,:]
#L=eclipcoord[1,:]
#g=eclipcoord[2,:]
#l=eclipcoord[3,:]
#ep=eclipcoord[4,:]

print('NOTE: This code is a bit less accurate for the range of hours of dawn and dusk.')
 