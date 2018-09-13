# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:02:44 2018

@author: Roberts
"""

import numpy as np, pylab as pl, matplotlib as mpl

mpl.rcParams['font.size'] = 12

maxPctCurrent = [25.86, 55.93, 26.65, 55.55, 35.59, 7.53, 11.15, 24.08]
maxPctLength = [26.13, 56.63, 28.83, 56.25, 36.34, 13.07, 19.35, 24.38]
x = 3*np.arange(len(maxPctCurrent)) + 0.5
x2 = 3*np.arange(len(maxPctCurrent)) - 0.5

parm = ('Temperature *16136 K \n [8000, 24272] K', 
        'Pressure *1200 kPa \n [10, 2390] kPa', 
        'Arc resistivity *0.0175 $\Omega$.cm \n [0.013, 0.03] $\Omega$.cm', 
        'Enthalpy *10260 kJ/kg \n [1000, 15000] kJ/kg', 
        'Cathode current density \n *3500 A/cm$^2$, [2300, 5000] A/cm$^2$', 
        'Anode work function \n *4.2 V, [1, 20] V', 
        'Anode voltage drop \n *6.6 V, [5, 30] V', 
        'Molar mass *0.0289 kg/mol \n [0.0289, 0.0718] kg/mol')

pl.figure() # current range
pl.bar(x,maxPctCurrent)
pl.bar(x2,maxPctLength)
pl.xticks(x, parm, rotation=75, horizontalalignment = 'right')
pl.subplots_adjust(bottom=0.15)

pl.ylabel('Maximum difference in the \n percentage of arc heat radiated')
pl.legend(['Arc current = [1 kA, 90 kA] \n arc length = 30 cm',
           'Arc current = 15 kA \n Arc length [5 cm, 100 cm]'])