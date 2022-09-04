#!/usr/bin/env python
# coding: utf-8

# In[1]:
from xml.dom.minidom import parse
import os
import pandas as pd
import numpy as np
import math
class Scenario:
    
    # Farm Parameters
    CT=0.8
    farmRad=0
    PRated=1500.0
    R=38.5
    eta=-500.0
    k=0.0750
    lam=140.86
    vCin=3.5
    vCout=20
    vRated=14
    width=0
    height=0
    wfe=0
    nturbines=0
    obstacles=0
    obstacles_df=0
    
    # Wind Resources
    c=np.empty((24))
    ks=np.empty((24))
    omega=np.empty((24))
    theta=np.empty((24,2))

    # Optimisation Parameters
    fac=np.pi/180
    coSinMidThetas=np.empty((0,0))
    rkRatio=0
    krRatio=0
    vints=np.empty((0))
    wblcdfValues=np.empty((0))
    wblcdfAccuracy=0
    cMax=0
    cMin=0
    atan_k=0
    trans_CT=0
    minDist=0

    def getScenario(self,sc):
        dir_path = "C:\\Users\\ggnar\\OneDrive\\Desktop\\MSc Data Science\\Dissertation\\Testing\\WindFlo"
        data_records = []
        mysc=sc

        with parse(f'{dir_path}/{mysc}') as xml_doc:
            root = xml_doc.documentElement


            angles = root.getElementsByTagName('angle')
            
            x=0
            for r in angles:
                self.c[x]=float(r.getAttributeNode("c").value)
                self.ks[x]=float(r.getAttributeNode("k").value)
                self.omega[x]=float(r.getAttributeNode("omega").value)
                self.theta[x][0]=x*15.0
                self.theta[x][1]=(x+1)*15.0
                x+=1

            data=[]
            data.append(self.c)
            data.append(self.ks)
            data.append(self.omega)
            data.append(self.theta)
            angles_df=pd.DataFrame(data).transpose()
            angles_df.columns=["c","k","omega","theta"]
            print("Wind Angles: ")
            print(angles_df)

            self.obstacles = root.getElementsByTagName('obstacle')
            xmin=[]
            xmax=[]
            ymin=[]
            ymax=[]

            if len(self.obstacles)>0:

                for r in self.obstacles:
                    xmin.append(int(r.getAttributeNode("xmin").value))
                    xmax.append(int(r.getAttributeNode("xmax").value))
                    ymin.append(int(r.getAttributeNode("ymin").value))
                    ymax.append(int(r.getAttributeNode("ymax").value))

                data=[]
                data.append(xmin)
                data.append(xmax)
                data.append(ymin)
                data.append(ymax)
                self.obstacles_df=pd.DataFrame(data).transpose()
                self.obstacles_df.columns=["xmin","xmax","ymin","ymax"]
                print("\nObstacles:")
                print(self.obstacles_df)

            print("\nOther Parameters:")
            self.width = int(root.getElementsByTagName('Width')[0].firstChild.data)
            print(self.width)
            self.height = int(root.getElementsByTagName('Height')[0].firstChild.data)
            print(self.height)
            self.nturbines = int(root.getElementsByTagName('NTurbines')[0].firstChild.data)
            print(self.nturbines)
            self.wfe = float(root.getElementsByTagName('WakeFreeEnergy')[0].firstChild.data)
            print(self.wfe)
            return angles_df, self.obstacles_df, self.width, self.height, self.nturbines, self.wfe

    def initOptpar(self):
        self.coSinMidThetas=np.empty((len(self.theta),2))
        for thets in range(0,len(self.theta)):
            t=(self.theta[thets][0]+self.theta[thets][1])/2.0*self.fac
            self.coSinMidThetas[thets][0]=math.cos(t)
            self.coSinMidThetas[thets][1]=math.sin(t)
        self.rkRatio=self.R/self.k
        self.krRatio=self.k/self.R
        self.vints=np.empty((int(2.0*self.vRated-7.0+1.0)))
        for i in range(0,len(self.vints)):
            self.vints[i]=3.5+i*0.5
        self.atan_k=math.atan(self.k)
        self.trans_CT=1.0-math.sqrt(1.0-self.CT)
        self.minDist=64.0*self.R*self.R

    def getCosMidThetas(self,idx):
        return self.coSinMidThetas[idx][0]
    def getSinMidThetas(self,idx):
        return self.coSinMidThetas[idx][1]

    def getWblcdfVints(self,c, vintIndex, ksIndex):
        return self.wblcdfValues[int((c-cMin)/self.wblcdfAccuracy)*(len(self.vints)+1)*len(self.ks)+self.vintIndex*len(self.ks)+self.ksIndex]
    def getWblcdfVrated(self,c, ksIndex):
        return self.wblcdfValues[int((c-cMin)/self.wblcdfAccuracy)*(len(self.vints)+1)*len(self.ks)+len(self.vints)*len(self.ks)+self.ksIndex]
    
    @staticmethod
    def wblcdf(x, sc, sh):
        return(1.0-math.exp(-Scenario.fastPow(x/sc,sh)))
    
    @staticmethod
    def fastPow(a,b):
        if(abs(b-2)<0.0001):
            return a*a
        elif(abs(b-1)<0.0001):
            return a
        elif(abs(b)<0.0001):
            return 1
        else:
            return pow(abs(a),abs(b))

