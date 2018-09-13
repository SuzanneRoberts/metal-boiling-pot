# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 15:21:33 2017

@author: Roberts
Simplified, low, computational arc-heat model described by Fathi et al. (2015)
WARNING: arc length is specified in cm (not in the SI-unit of m)!

"""

import numpy as np, pylab as pl, unittest
from pylab import rcParams
rcParams['figure.figsize'] = 6, 4


class arc(object):
    # arc object containing values of arc properties that are assumed to be constant
    def __init__(self, cathodeSpotCurrentDensity, arcResistivity, averageArcTemperature, 
                 averageArcPressure, magneticPermeability, airMolarMass, 
                 anodeWorkFunction, electronCharge, anodeVoltageDrop, arcAirEnthalpyChange):
        
        self.j_k    = cathodeSpotCurrentDensity
        self.rho_a  = arcResistivity
        self.T      = averageArcTemperature
        self.P      = averageArcPressure
        self.mu_0   = magneticPermeability
        self.M      = airMolarMass
        self.O_an   = anodeWorkFunction
        self.e      = electronCharge
        self.U_an   = anodeVoltageDrop
        self.dh     = arcAirEnthalpyChange
        self.k_B    = 1.38064852e-23            # J/K, Boltzmann konstante
        self.R      = 8.314459848               # J/mol.K, universal gas constant
        
def arcShape(theI, theL, theArc):
    #========================#
    # Constants (Calculated) #
    #========================#

    r_k   = np.sqrt(theI/(np.pi*theArc.j_k))    # cm, cathode spot radius

    #========================================================#
    # Arc shape calculation (Section 2.1, Logar et al. 2015) #
    #========================================================#
    
    r = lambda z: (3.2*r_k) - (2.2*r_k*np.exp(-z/(5*r_k)))  # cm

    #=========================================================#
    # CAM radius calculation (Section 2.2, Logar et al. 2015) #
    #=========================================================#

    Vol = np.pi*(theL/3)*( (r_k)**2 + r(theL)**2 + (r(theL)*(r_k)) )   # cm^3
    R_arc = np.sqrt(Vol/(np.pi*theL))                          # cm
    
    return r_k, R_arc, r


def arcVoltage(theArc, r_k, theI, theL):
    
    #================================================================================================#
    # Arc voltage calculation (Section 2.3, Logar et al. 2015 | Section 2, Jones and Reynolds 2002 ) #
    #================================================================================================#
    
    a = 3.2*r_k 
    b = -2.2*r_k
    m = -1/(5*r_k)
    
    term1 = -1/(a**2 + a*b)
    term2 = 1/( a**2 + (a*b*np.exp(m*theL)) )
    term3 = np.log(a+b)/( a**2 )
    term4 = m*theL/(a**2)
    term5 = -np.log(a + (b*np.exp(m*theL)) )/(a**2)
    
    V_a = (theI*theArc.rho_a/(m*np.pi))*(term1 + term2 + term3 + term4 + term5)
    
    return V_a


def arcVelocity(theArc, r_k, R_arc, theI):
    # m not cm from here!!!

    #==========================================================#
    # MHD results calculation (Section 2.4, Logar et al. 2015) #
    #==========================================================#
    rho_k = theArc.P*theArc.M/(theArc.R*theArc.T)                 # kg/m^3, arc density
    v_arc = 0.5*theI/(np.pi*(R_arc/100)) * np.sqrt( (5*theArc.mu_0/3) * ((R_arc/100)**2/(r_k/100)**2 - 1) / (2*rho_k) )
    
    return v_arc, rho_k


def arcHeatDistribution(theI1, theL1, theArc1):
    r_k1, R_arc1, r1   = arcShape(theI1, theL1, theArc1)
    V_a1               = arcVoltage(theArc1, r_k1, theI1, theL1)
    v_arc1, rho_k1     = arcVelocity(theArc1, r_k1, R_arc1, theI1)
    
    #=====================================================================#
    # Arc energy dissipation calculation (Section 2.5, Logar et al. 2015) #
    #=====================================================================#
    
    P_a = V_a1*theI1
    
    #==================================================================================#
    # Arc electron flow energy dissipation calculation (Equation 3, Logar et al. 2015) #
    #==================================================================================#
    
    P_e = theI1*(theArc1.O_an + (5*theArc1.k_B*theArc1.T/(2*theArc1.e*np.pi*R_arc1**2)) + theArc1.U_an)
    
    #===============================================================================#
    # Arc convection energy dissipation calculation (Equation 4, Logar et al. 2015) #
    #===============================================================================#
    
    P_conv = np.pi*((R_arc1/100)**2)*rho_k1*v_arc1*(theArc1.dh)
    
    #===============================================================================#
    # Arc radiation energy dissipation calculation (Equation 12, Logar et al. 2015) #
    #===============================================================================#
    
    P_r = P_a - (P_conv + P_e)
    
    #===============================================#
    # Share of heat dissipated by arc electron flow #
    #===============================================#
    
    Q_P_e = P_e / P_a * 100
    
    #============================================#
    # Share of heat dissipated by arc convection #
    #============================================#
    
    Q_P_conv = P_conv / P_a * 100
    
    #============================================#
    # Share of heat dissipated by arc radiontion #
    #============================================#
    
    Q_P_r = P_r / P_a * 100
    
    return r_k1, R_arc1, V_a1, v_arc1, rho_k1, P_a, P_e, P_conv, P_r, Q_P_e, Q_P_conv, Q_P_r


def plotPowerAgainstArcLength(theI, theL, theArc, markerStr = '-', varStr = 'varStr', var1 = 'var1', var2 = 'var2', var3 = 'var3'):
    (cathodeSpotRadius, R_eff, avgArcVoltage, avgArcVelocity, avgArcDensity, 
     Ptotal, Pelec, Pconv, Prad, Pelec100, Pconv1str00, Prad100) = arcHeatDistribution(
             theI, theL, theArc)
    
    # plotting shape
    pl.figure(1)
    #pl.plot(L, cathodeSpotRadius, 'b')  # only a function of arc current, not of arc length
    pl.plot(theL, R_eff, 'y')
    pl.xlabel('Arc length [cm]')
    pl.ylabel('Effective arc radius [cm]')
    
    # plotting arc properties: arc voltage
    pl.figure(2)
    pl.plot(theL, avgArcVoltage, 'b' + markerStr)
    pl.xlabel('Arc length [cm]')
    pl.ylabel('Arc voltage [V]')
    
    # plotting arc properties: arc velocity
    pl.figure(3)
    pl.plot(theL, avgArcVelocity, 'y' + markerStr)
    pl.xlabel('Arc length [cm]')
    pl.ylabel('Arc velocity [m/s]')
    
    # plotting arc properties: arc density
    pl.figure(4)
    pl.plot(theL, avgArcDensity*np.ones(len(theL)), 'g' + markerStr) # not a function of arc current or arc length
    pl.xlabel('Arc length [cm]')
    pl.ylabel('Arc density [kg/m$^3$]')
    
    # plotting power percentages   
    pl.figure(5)
    pl.plot(theL, Pconv1str00, 'b' + markerStr)
    pl.plot(theL, Prad100, 'y' + markerStr)
    pl.plot(theL, Pelec100, 'g' + markerStr)
    pl.plot(theL, Prad100+Pconv1str00+Pelec100, 'r')
    pl.xlabel('Arc length [cm]')
    pl.ylabel('Power contribution [%]')
    pl.plot([], [], 'k')
    pl.plot([], [], 'k-.')
    pl.plot([], [], 'k:')
    pl.legend(['Conv','Rad','Elec','Total', var1, var2, var3])
    pl.savefig('autoPlot/'+varStr+'_pct.jpg', bbox_inches='tight')
    pl.savefig('autoPlot/'+varStr+'_pct.eps', bbox_inches='tight')
    
    # plotting power
    pl.figure(6)
    pl.plot(theL, Pconv/1e6, 'b' + markerStr)
    pl.plot(theL, Prad/1e6, 'y' + markerStr)
    pl.plot(theL, Pelec/1e6, 'g' + markerStr)
    pl.plot(theL, Ptotal/1e6, 'r' + markerStr)
    pl.plot([], [], 'k')
    pl.plot([], [], 'k-.')
    pl.plot([], [], 'k:')
    pl.legend(['Conv','Rad','Elec','Total', var1, var2, var3])
    pl.xlabel('Arc length [cm]')
    pl.ylabel('Power released from arc [MW]')
    pl.savefig('autoPlot/'+varStr+'_power.jpg', bbox_inches='tight')
    pl.savefig('autoPlot/'+varStr+'_power.eps', bbox_inches='tight')
    
    return Prad100, Prad
    

def plotPowerAgainstCurrent(theI, theL, theArc, markerStr = '-', varStr = 'varStr', var1 = 'var1', var2 = 'var2', var3 = 'var3'):
    (cathodeSpotRadius, R_eff, avgArcVoltage, avgArcVelocity, avgArcDensity, 
     Ptotal, Pelec, Pconv, Prad, Pelec100, Pconvstr100, Prad100) = arcHeatDistribution(
             theI, theL, theArc)
    
    # plotting shape
    pl.figure(7)
    pl.plot(theI/1e3, cathodeSpotRadius, 'b')  # only a function of arc current, not of arc length
    pl.plot(theI/1e3, R_eff, 'y')
    pl.legend(['Cath spot rad','R_eff'])
    pl.xlabel('Arc current [kA]')
    pl.ylabel('Radius [cm]')
    
    # plotting arc properties: arc voltage
    pl.figure(8)
    pl.plot(theI/1e3, avgArcVoltage, 'b' + markerStr)
    pl.xlabel('Arc current [kA]')
    pl.ylabel('Arc voltage [V]')
    
    # plotting arc properties: arc velocity
    pl.figure(9)
    pl.plot(theI/1e3, avgArcVelocity, 'y' + markerStr)
    pl.xlabel('Arc current [kA]')
    pl.ylabel('Arc velocity [m/s]')
    
    # plotting arc properties: arc density
    pl.figure(10)
    pl.plot(theI/1e3, avgArcDensity*np.ones(len(theI)), 'g' + markerStr) # not a function of arc current or arc length
    pl.xlabel('Arc current [kA]')
    pl.ylabel('Arc density [kg/m$^3$]')
    
    # plotting power percentages
    pl.figure(11)
    pl.plot(theI/1e3, Pconvstr100, 'b' + markerStr)
    pl.plot(theI/1e3, Prad100, 'y' + markerStr)
    pl.plot(theI/1e3, Pelec100, 'g' + markerStr)
    pl.plot(theI/1e3, Prad100+Pconvstr100+Pelec100, 'r')
    pl.xlabel('Arc current [kA]')
    pl.ylabel('Power contribution [%]')
    pl.plot([], [], 'k')
    pl.plot([], [], 'k-.')
    pl.plot([], [], 'k:')
    pl.legend(['Conv','Rad','Elec','Total', var1, var2, var3])
    pl.savefig('autoPlot/'+varStr+'_pct.jpg', bbox_inches='tight')
    pl.savefig('autoPlot/'+varStr+'_pct.eps', bbox_inches='tight')
    
    # plotting power
    pl.figure(12)
    pl.plot(theI/1e3, Pconv/1e6, 'b' + markerStr)
    pl.plot(theI/1e3, Prad/1e6, 'y' + markerStr)
    pl.plot(theI/1e3, Pelec/1e6, 'g' + markerStr)
    pl.plot(theI/1e3, Ptotal/1e6, 'r' + markerStr)
    pl.plot([], [], 'k')
    pl.plot([], [], 'k-.')
    pl.plot([], [], 'k:')
    pl.legend(['Conv','Rad','Elec','Total', var1, var2, var3])
    pl.xlabel('Arc current [kA]')
    pl.ylabel('Power released from arc [MW]')
    pl.savefig('autoPlot/'+varStr+'_power.jpg', bbox_inches='tight')
    pl.savefig('autoPlot/'+varStr+'_power.eps', bbox_inches='tight')
    
    return Prad100, Prad
    

#=============================================================================#
# Main program                                                                #
#=============================================================================#

# best guess arc        
bga = arc(3500,           # A/cm^2, Stromdichte des Brennflecks, Logar et al. zitiert Reynolds 2012
          0.0175,         # Ohm.cm, arc resistivity, Logar et al. zitiert Reynolds and Jones 2004
          16136,          # K, estimated average arc temperature, Logar et al.
          1200e3,         # Pa, estimated average arc pressure, Logar et al.
          4e-7 * np.pi,   # H/m, magnetic permeability of free space
          0.028971,       # kg/mol, molar mass of air
          4.2,            # V, work function for the anode, Bowman & Krüger, p12
          1.602e-19,      # J, electron charge (1 eV)
          6.6,             # V, anode voltage drop, Echterhof Folien: 10 V bis 30 V
          130e5 - 27.4e5)    # J/kg, specific enthalpy of the arc, using air at average arc temperature, Fig 8, Sanchez et al. 2009
                            # J/kg, specific enthalpy of the air surrounding the arc @ average air temperature, 1200 K (Gruber, 2015), Fig 8, Sanchez et al. 2009

# sensitivity study set-up
indVar   = 'current'
vStr     = 'temperature_'
vUnitStr = ' K'
v1 = 8000
v2 = 16136 
v3 = 24272 
v1str   = str(v1) + vUnitStr
v2str   = str(v2) + vUnitStr
v3str   = str(v3) + vUnitStr


if indVar == 'length':                            
   
    I = 15000                 #A, arc current
    L = np.linspace(5, 100, 75) # np.array([5, 20, 80]) #    #cm, arc length
    
    bga.T = v1
    Prad100v1, Pradv1 = plotPowerAgainstArcLength(I, L, bga, varStr=vStr+indVar, var1=v1str, var2=v2str, var3=v3str)
    bga.T = v2
    Prad100v2, Pradv2 = plotPowerAgainstArcLength(I, L, bga, '-.', varStr=vStr+indVar, var1=v1str, var2=v2str, var3=v3str)
    bga.T = v3
    Prad100v3, Pradv3 = plotPowerAgainstArcLength(I, L, bga, ':', varStr=vStr+indVar, var1=v1str, var2=v2str, var3=v3str)
    
    print('minDiffDown', np.min(np.abs(Pradv2 - Pradv1)) / 1e6 , 'MW')
    print('minDiffUp', np.min(np.abs(Pradv2 - Pradv3)) / 1e6 , 'MW')
    
    print('maxDiffDown', np.max(np.abs(Pradv2 - Pradv1)) / 1e6 , 'MW')
    print('maxDiffUp', np.max(np.abs(Pradv2 - Pradv3)) / 1e6 , 'MW')
    
    print('minDiffDown100', np.min(np.abs(Prad100v2 - Prad100v1)))
    print('minDiffUp100', np.min(np.abs(Prad100v2 - Prad100v3)))
    
    print('maxDiffDown100', np.max(np.abs(Prad100v2 - Prad100v1)))
    print('maxDiffUp100', np.max(np.abs(Prad100v2 - Prad100v3)))
    
    
if indVar == 'current':

    I2 = np.linspace(1000, 90000, 75) #A, arc current
    L2 = 30  # cm, arc length
    
    bga.T = v1
    Prad100v1, Pradv1 = plotPowerAgainstCurrent(I2, L2, bga, varStr=vStr+indVar, var1=v1str, var2=v2str, var3=v3str)
    bga.T = v2
    Prad100v2, Pradv2 = plotPowerAgainstCurrent(I2, L2, bga, '-.', varStr=vStr+indVar, var1=v1str, var2=v2str, var3=v3str)
    bga.T = v3
    Prad100v3, Pradv3 = plotPowerAgainstCurrent(I2, L2, bga, ':', varStr=vStr+indVar, var1=v1str, var2=v2str, var3=v3str)
    
    print('minDiffDown', np.min(np.abs(Pradv2 - Pradv1)) / 1e6 , 'MW')
    print('minDiffUp', np.min(np.abs(Pradv2 - Pradv3)) / 1e6 , 'MW')
    
    print('maxDiffDown', np.max(np.abs(Pradv2 - Pradv1)) / 1e6 , 'MW')
    print('maxDiffUp', np.max(np.abs(Pradv2 - Pradv3)) / 1e6 , 'MW')
    
    print('minDiffDown100', np.min(np.abs(Prad100v2 - Prad100v1)))
    print('minDiffUp100', np.min(np.abs(Prad100v2 - Prad100v3)))
    
    print('maxDiffDown100', np.max(np.abs(Prad100v2 - Prad100v1)))
    print('maxDiffUp100', np.max(np.abs(Prad100v2 - Prad100v3)))

#### classes - lees maandag verder
#### http://greenteapress.com/ModSimPy/ModSimPy.pdf
    
# unittests
class testArc(unittest.TestCase):
    def test_arcShape(self):
        r_k2, R_arc2, r2   = arcShape(15000, 30, bga)
        z = np.linspace(0,100,100)
        
        # assert that the arc shape is defined in the right direction
        # (smaller diameter at the kathode; larger diameter at the anode)
        # (z = 0: kathode)
        self.assertTrue(r2(100) > r_k2)
        
        pl.figure()
        pl.plot(z,r2(z))
        pl.show()
    
if __name__ == '__main__':
    unittest.main()