#!/usr/bin/env python
# coding: utf-8

# In[3]:
from Scenarios import Scenario
import numpy as np
import math

class Evaluator:
    energyCapture=0
    wakeFreeRatio=0
    energyCost=999999
    sc=0
    tspe=np.empty((0,0))
    tpos=np.empty((0,0))
    def initialise(self,scenario):
        self.sc=scenario

    
    
    # Compute energy output (calculating the wake)
    def evaluateLayout(self,layout):

        #nEvals?

        #copying the layout
        self.tpos=np.empty((len(layout),len(layout[0])))
        for i in range(0,len(layout)):
            for j in range(len(layout[i])):
                self.tpos[i][j]=layout[i][j]
             
        if self.checkConstraints(layout):
            #wind resource per turbine stored temporarily in tspe
            self.tspe=np.empty((len(self.sc.theta),len(self.tpos)))
            for turb in range(len(self.tpos)): # For each Turbine
                for thets in range(len(self.sc.theta)): # For each direction
                    #calculate wake
                    totalVdef=self.calWakeTurb(turb,thets)
                    cTurb=self.sc.c[thets]*(1.0-totalVdef)
                    # Annual power output per turbine per direction
                    tint=self.sc.theta[thets][1]-self.sc.theta[thets][0]
                    w=self.sc.omega[thets]
                    ki=self.sc.ks[thets]
                    totPow=0
                    for ghh in range(1,len(self.sc.vints)):
                        v=(self.sc.vints[ghh]+self.sc.vints[ghh-1])/2
                        P=self.powOutput(v)
                        prv=Scenario.wblcdf(self.sc.vints[ghh],cTurb,ki)-Scenario.wblcdf(self.sc.vints[ghh-1],cTurb,ki)
                        totPow+=prv*P
                    totPow+=self.sc.PRated*(1.0-Scenario.wblcdf(self.sc.vRated, cTurb, ki))
                    totPow*=tint*w
                    self.tspe[thets][turb]=totPow
                    self.energyCapture+=totPow
            self.wakeFreeRatio=self.energyCapture/(self.sc.wfe*len(self.tpos))
            return self.wakeFreeRatio
        else:
            self.energyCapture=0
            self.wakeFreeRatio=0
            self.tspe=0
            return 0

    def evaluate(self,layout):
        ct = 750000
        cs = 8000000
        m = 30
        r = 0.03
        y = 20
        com = 20000
        print("Reached Evaluate")

        wfr=self.evaluateLayout(layout)
        if (wfr<=0):
            return 99999
        n=len(layout)

        self.energyCost = (((ct*n+cs*np.floor(n/m))*(0.666667+0.333333*np.exp(-0.00174*n*n))+com*n)/
                      ((1.0-np.power(1.0+r,-y))/r)/(8760.0*self.sc.wfe*wfr*n))+0.1/n
        return self.energyCost

    def getEnergyOutputs():
        return self.tspe

    def getTurbineFitness():
        res=np.empty((len(self.tspe[0])))
        for i in range(len(res)):
            res[i]=0
            for j in range(len(self.tspe)):
                res[i]+=self.tspe[j][i]
            res[i]=res[i]/self.sc.wfe
        return res

    # Check constraints (avoid obstacles)
    def checkConstraints(self,layout):
        for i in range(len(layout)):
            if (layout[i][0]!=layout[i][0] or layout[i][1]!=layout[i][1] or layout[i][0]<0.0 or layout[i][1]<0.0 or layout[i][0]>self.sc.width  or layout[i][1]>self.sc.height):
                print("Layout Constraint")
                print("Turbine ",i,"(",layout[i][0],",",layout[i][1],") is invalid.")
                return False
            for j in range(len(self.sc.obstacles)):
                if (layout[i][0] > self.sc.obstacles_df.loc[j]["xmin"] and
                    layout[i][0] < self.sc.obstacles_df.loc[j]["xmax"] and
                    layout[i][1] > self.sc.obstacles_df.loc[j]["ymin"] and
                    layout[i][1] < self.sc.obstacles_df.loc[j]["ymax"]):
                    print("Obstacle Constraint")
                    print("Turbine ",i,"(",layout[i][0],",",layout[i][1],") is invalid.")
                    return False
            for j in range(len(layout)):
                if(i!=j):
                    dist=(layout[i][0]-layout[j][0])*(layout[i][0]-layout[j][0])+(layout[i][1]-layout[j][1])*(layout[i][1]-layout[j][1])
                    if(dist<self.sc.minDist):
                        print("Distance Constraint: ",dist)
                        print("Turbine ",i,"(",layout[i][0],",",layout[i][1],") is invalid.")
                        return False

        print("Valid Layout")
        return True

    def calWakeTurb(self,turb,thetidx):
        x=self.tpos[turb][0]
        y=self.tpos[turb][1]
        velDef=0
        for oturb in range(len(self.tpos)):
            if (oturb!= turb):
                xo=self.tpos[oturb][0]
                yo=self.tpos[oturb][1]
                beta=self.calculateBeta(x,y,xo,yo,thetidx)
                if(beta<self.sc.atan_k):
                    dij=self.calculateProjectedDist(x,y,xo,yo,thetidx)
                    curDef=self.calculateVelocityDeficit(dij)
                    velDef+=curDef*curDef
        return math.sqrt(velDef)
    
    def calculateBeta(self,xi,yi,xj,yj,thetidx):
        num=((xi-xj)*self.sc.getCosMidThetas(thetidx)+(yi-yj)*self.sc.getSinMidThetas(thetidx)+self.sc.rkRatio)
        a=xi-xj+self.sc.rkRatio*self.sc.getCosMidThetas(thetidx)
        b=yi-yj+self.sc.rkRatio*self.sc.getSinMidThetas(thetidx)
        denom=math.sqrt(a*a+b*b)
        return math.acos(num/denom)

    def powOutput(self,v):
        if (v<self.sc.vCin):
            return 0
        elif(v>=self.sc.vCin and v<=self.sc.vRated):
            return self.sc.lam*v+self.sc.eta
        elif(v<self.sc.vCout and v>self.sc.vRated):
            return self.sc.PRated
        else:
            return 0

    def calculateProjectedDist(self,xi,yi,xj,yj,thetidx):
        return abs((xi-xj)*self.sc.getCosMidThetas(thetidx)+(yi-yj)*self.sc.getSinMidThetas(thetidx))
    
    def calculateVelocityDeficit(self,dij):
        return self.sc.trans_CT/((1.0+self.sc.krRatio*dij)*(1.0+self.sc.krRatio*dij))

