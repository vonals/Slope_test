# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 09:14:36 2020

@author: Aizhong Ye
"""
import math
#1.read data
f = open('曹坪.asc','r')
lines = f.readlines()
f.close()
dem={}
for i in range(6):
    a,b=lines[i].split()
    dem[a]=float(b) 

# elevation=[[0]*int(dem['ncols'])]*int(dem['nrows'])
elevation=[]
for i in range(6,int(dem['nrows'])+6):
    a=lines[i].split()
    elevation1=[]
    for j in a:
        elevation1.append(float(j)) 
    elevation.append(elevation1)

#2.calculate the slope and aspect
slope=[[dem['NODATA_value']]*int(dem['ncols'])]*int(dem['nrows'])
aspect=[[dem['NODATA_value']]*int(dem['ncols'])]*int(dem['nrows'])
for i in range(1,int(dem['nrows'])-1):
    for j in range(1,int(dem['ncols'])-1): 
        flag1=[0]*9
        flag1[0]=1 if abs(elevation[i][j]-dem['NODATA_value']) < 1 else 0
        flag1[1]=1 if abs(elevation[i][j-1]-dem['NODATA_value'])<1 else 0
        flag1[2]=1 if abs(elevation[i][j+1]-dem['NODATA_value'])<1 else 0
        flag1[3]=1 if abs(elevation[i-1][j]-dem['NODATA_value'])<1 else 0
        flag1[4]=1 if abs(elevation[i-1][j-1]-dem['NODATA_value'])<1 else 0
        flag1[5]=1 if abs(elevation[i-1][j+1]-dem['NODATA_value'])<1 else 0
        flag1[6]=1 if abs(elevation[i+1][j]-dem['NODATA_value'])<1 else 0
        flag1[7]=1 if abs(elevation[i+1][j-1]-dem['NODATA_value'])<1 else 0
        flag1[8]=1 if abs(elevation[i+1][j+1]-dem['NODATA_value'])<1 else 0
        if sum(flag1)<1 :
            we=((elevation[i-1][j-1]+2*elevation[i][j-1]+elevation[i+1][j-1])-\
                (elevation[i-1][j+1]+2*elevation[i][j+1]+elevation[i+1][j+1]))/\
                8/dem['cellsize']            
            sn=((elevation[i+1][j+1]+2*elevation[i+1][j]+elevation[i+1][j-1])-\
                (elevation[i-1][j-1]+2*elevation[i-1][j]+elevation[i-1][j+1]))/\
                8/dem['cellsize']
            slope[i][j]=math.sqrt(we**2+sn**2)
            if abs(we)>0.001:   aspect=sn/we
            else :  aspect=sn/0.001
            #we=((e5+2e1+e8)-(e6+2e3+e7))/8cell
            #sn=((e7+2e4+e8)-(e5+2e2+e6))/8cell
            #slope=sqrt(we**2+sn**2)
            #aspect=sn/we
            '''
            5 2 6    -1-1  -1 0 -1+1
            1 0 3     0-1   0 0  0+1
            8 4 7     +1-1 +1 0 +1+1
            
            '''
        
        
        
        
        
        
        
        
        
        