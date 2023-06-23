# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 09:12:04 2023

@author: eleobar
"""

import numpy as np
from astropy.coordinates import solar_system_ephemeris, get_body#, SkyCoord

'''
This script computes the signal of a point source for the calibration of MISTRAL and
assign the coordinates to each calibrator.
Coordinates: 
	- obj = SkyCoord.from_name('A301') for object outside the solar system
	- create a SkyCoord array with the (ra,dec) in deg of the object
	- with solar_system_ephemeris.set('builtin'):
	    uranus = get_body('uranus', time, SRT)
'''
#AGGIUNGERE:
	#Titan, Pallas, Vesta
	#Calcolare il diametro angolare dei calibratori a seconda del giorno!!!


D_SRT = 60 #m
band = 30e9 #Hz
wl = 0.0033 #m 
nu = 90e9 #GHz
eff = 0.3
c = 3e8 #m/s
k_B = 1.38e-23 #J/K
beam = (wl/D_SRT)/(2*np.log(2)) #beam gaussiano
R_opt = 2e11 #rad/W
A_tele= np.pi*((D_SRT/2)**2) #m2

class Calibrator():
	def __init__(self, calibrator="", time="", loc=""):
		self.calibrator = calibrator
		if calibrator == "Uranus":
			self.T = 125 #K at 90 GHz (125+-9)
			self.ang_diameter = 3.7 #arcsec (3.4-3.7 arcsec)
			self.D = self.ang_diameter / 3600 #m

			self.A = np.pi * (self.D/2)**2 #m2
			with solar_system_ephemeris.set('builtin'):
			    self.coord = get_body('uranus', time, loc) 
			
		elif calibrator == "Neptune":
			self.T = 126 #K at 90 GHz (126+-9)
			self.ang_diameter = 2.4 #arcsec (2.2-2.4 arcsec)
			self.D = self.ang_diameter / 3600 #m
			self.A = np.pi * (self.D/2)**2 #m2
			with solar_system_ephemeris.set('builtin'):
			    self.coord = get_body('neptune', time, loc) 
			
		elif calibrator == "Mars":
			self.T = 190 #K at 90 GHz
			self.ang_diameter = 3.5 #arcsec (3.5-25.1 arcsec)
			self.D = self.ang_diameter / 3600 #m
			self.A = np.pi * (self.D/2)**2 #m2
			with solar_system_ephemeris.set('builtin'):
			    self.coord = get_body('mars', time, loc) 
			
		elif calibrator == "Ceres":
			self.T = 137 #K at 90 GHz (137 +- 25)
			self.ang_diameter = 0.854 #arcsec (0.854″ to 0.339″)
			self.D = self.ang_diameter / 3600 #m
			self.A = np.pi * (self.D/2)**2 #m2
			#with solar_system_ephemeris.set('builtin'):
			#    self.coord = get_body('ceres', time, loc) 
			
		elif calibrator == "Callisto":
			self.T = 95 #K at 90 GHz (95 +- 17)
			self.ang_diameter = (2 * np.arctan (4820.6 / 628.3e6))*3600 #arcsec 
			self.D = self.ang_diameter / 3600 #m
			self.A = np.pi * (self.D/2)**2 #m2
			with solar_system_ephemeris.set('builtin'):
			    self.coord = get_body('callisto', time, loc) 
			
		elif calibrator == "Ganymede":
			self.T = 136 #K at 90 GHz (136 +- 21)
			self.ang_diameter = 1.8 #arcsec (1.2″ to 1.8″)
			self.D = self.ang_diameter / 3600 #m
			self.A = np.pi * (self.D/2)**2 #m2
			with solar_system_ephemeris.set('builtin'):
			    self.coord = get_body('ganymede', time, loc) 
		'''	
		else:
        self.T = 0 #K at 90 GHz 
        self.ang_diameter = 0 #arcsec 
        self.D = self.ang_diameter / 3600 #m 
        self.d_ang_rad = ang_diameter*np.pi/(180*3600)# 
        self.AO = np.pi*A_tele*(np.tan(d_ang_rad/2))**2 
        self.A = np.pi * (self.D/2)**2 #m2 
        self.coord = None
		'''
		#cal = Calibrator(calibrator = "ganymede")
calibrator='neptune'        
T= 136 #K
ang_diam_arcsec = 1.8 #arcsec
ang_diam_rad = ang_diam_arcsec*np.pi/(180*3600) #rad
AO = np.pi*A_tele*(np.tan(ang_diam_rad/2))**2 #m2 sr

B_nu = 2* nu**2 * k_B * T / (c**2) #W/m2/Hz/sr
B = B_nu * band #W/m2/sr
P = B * AO #W
		
#in order to compute the phase (signal) read by the kids, we have to consider the optical responsivity
S_ph = P * R_opt * eff #rad
		
#conversione densità di flusso in Jansky: 1 Jy = 10^-26 W/m2/Hz
# -> approx 1 Jy = 0.5 rad
beam= wl/D_SRT
A_beam = np.pi*beam**2/(4*np.log(2))
P_Jy = B_nu*1e26*A_beam
		
print("-----------"+calibrator+"-----------------")
print("Angular diameter = "+str(ang_diam_arcsec)+" arcsec")
print(f"Phase = {S_ph:.5f} rad")
print("Flux density= "+str(P_Jy)+" Jy")
print("-------------------------------------")





