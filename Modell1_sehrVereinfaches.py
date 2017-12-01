# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 09:49:11 2017

@author: Roberts
"""

# Erstes ganz vereinfachtes Lichtbogen Modell

#========#
# Inputs #
#========#

# Elektrische Spannung
U_bogen = 450        # V, Spannungsabfall über den Lichtbogen = Ofenspannung
U_kath  = 10         # V, Kathoden-Spannungsabfall 
U_an    = 30         # V, Anoden-Spannungsabfall
E_bogen = 0.8*1000   # V/m, Feldstärke des Bogens, AC

# Bogenmaterialeigenschaften
sigma_bogen = 50000 # S/m, ELektrische Leitfähigkeit, eigentlich sigma(T,p)


#==============#
# Berechnungen #
#==============#

# Lichtbogenlänge Berechnung
# U_bogen = U_kath + E_bogen*l + U_an
l_bogen = (U_bogen - U_kath - U_an)/E_bogen

# Lichtbogenstrom Berechnung 
# r_bogen =
# R_bogen =
# I_bogen = (E_bogen*l_bogen)/(l_bogen/R-bogen)


#===============================#
# Druck Outputs nach Bildschirm #
#===============================#

print('############################################')
print('Input values')
print('U_bogen = ', U_bogen, ' V')
print('U_kath = ', U_kath, ' V')
print('U_an = ', U_an, ' V')
print('E_bogen = ', E_bogen, ' V/m = ', E_bogen/1000, ' V/mm' )

print('############################################')
print('Lichtbogenlänge, l = ', l_bogen, ' m = ', l_bogen*1000, ' mm')

print('############################################')