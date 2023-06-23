# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 11:41:58 2023

@author: eleobar
"""

'''
Suppose you are planning to visit the Sardinia Radio Telescope (Lat. 39°29′34″N - Long. 9°14′42″E - Alt. 600 m). 
You want to plan the observation at 9:00 pm local time, and you want 
to know if the calibrator/cluster/cluster pair will be up.

Minimum observational altitude = 30°
'''

import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style, quantity_support
import astropy.units as u
from astropy.coordinates import solar_system_ephemeris, AltAz, EarthLocation, Angle, get_body, get_sun, get_moon
from astropy.time import Time
from astroplan import Observer
from astroplan.plots import plot_sky

#Fot matplotlib, use a nicer set of plot parameters and set up support for plotting/converting quantities
plt.style.use(astropy_mpl_style)
quantity_support()


SRT_lat = Angle('39°29′34″')
SRT_lon = Angle('9°14′42″E')
SRT_alt = 600
SRT = EarthLocation(lat=SRT_lat, lon=SRT_lon, height=SRT_alt*u.m)
utcoffset = +2*u.hour  # CEST
time = Time('2024-2-1 21:00:00') - utcoffset

with solar_system_ephemeris.set('builtin'):
    uranus = get_body('uranus', time, SRT)

#example of observation at a certain date at a certain time
uranuscoord = uranus.transform_to(AltAz(obstime=time,location=SRT))
print('----------------------------------------------------------')
print('Date and Time of observation:', time+utcoffset)
print(f"Uranus's Altitude = {uranuscoord.alt:.2}")
print(f"Uranus's Azimuth = {uranuscoord.az:.2}")
print('----------------------------------------------------------')

#Find the airmass at 100 times evenly spaced between 10pm and 7am CEST:
midnight = Time('2024-2-1 00:00:00') - utcoffset
delta_midnight = np.linspace(-2, 10, 100)*u.hour
frame_night = AltAz(obstime=midnight+delta_midnight,
                          location=SRT)
uranuscoords_night = uranus.transform_to(frame_night)

uranusairmasss_night = uranuscoords_night.secz
'''
plt.plot(delta_midnight, uranusairmasss_night)
plt.xlim(-2, 10)
plt.ylim(1, 4)
plt.xlabel('Hours from CEST Midnight')
plt.ylabel('Airmass [Sec(z)]')
plt.show()
'''

#Find the visibility of a calibrator at 100 times evenly spaced between 12am of a day and 12am of the following day CEST:
delta_midnight = np.linspace(-12, 12, 1000)*u.hour
time_obs_night = midnight + delta_midnight
frame_obs_night = AltAz(obstime=time_obs_night, location=SRT)

sunaltazs_obs_night = get_sun(time_obs_night).transform_to(frame_obs_night)
mooncoord_obs_night = get_moon(time_obs_night).transform_to(frame_obs_night) 
uranuscoords_obs_night = uranus.transform_to(frame_obs_night)

plt.plot(delta_midnight, sunaltazs_obs_night.alt, color='gold', label='Sun')
plt.scatter(delta_midnight, mooncoord_obs_night.alt, c=mooncoord_obs_night.az, lw=0, s=8, label='Moon', cmap='gist_gray')# color=[0.75]*3,, ls='--'
plt.colorbar().set_label('Azimuth [deg]')
plt.axhline(30)
plt.axhline(80)
plt.axhspan(0, 30, facecolor='grey', alpha=0.2)
plt.axhspan(80, 90, facecolor='grey', alpha=0.2)
plt.scatter(delta_midnight, uranuscoords_obs_night.alt, c=uranuscoords_obs_night.az, label='Uranus', lw=0, s=8, cmap='viridis')
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, sunaltazs_obs_night.alt < -0*u.deg, color='0.5', zorder=0)
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, sunaltazs_obs_night.alt < -18*u.deg, color='k', zorder=0)
plt.colorbar().set_label('Azimuth [deg]')
plt.legend(loc='upper left')
plt.xlim(-12*u.hour, 12*u.hour)
plt.xticks((np.arange(13)*2-12)*u.hour)
plt.ylim(0*u.deg, 90*u.deg)
plt.xlabel('Hours from CEST Midnight')
plt.ylabel('Altitude [deg]')
plt.title(f'Visibily of Uranus\n({time_obs_night[0]} to {time_obs_night[-1]})', fontsize=12)
plt.show()

style_uranus = {"label" : 'Uranus', 'c' : 'paleturquoise',  'lw' : 0, 's' : 8}
loc= Observer(location=SRT)
moon = get_moon(time_obs_night)
style_moon = {"label" : 'Moon', 'c' : 'silver',  'lw' : 0, 's' : 8}
sun = get_sun(time_obs_night)
style_sun = {"label" : 'Sun', 'c' : 'gold',  'lw' : 0, 's' : 8}
plot_sky(uranus, loc, time_obs_night, style_kwargs=style_uranus)
plot_sky(moon, loc, time_obs_night, style_kwargs=style_moon)
plot_sky(sun, loc, time_obs_night, style_kwargs=style_sun) 
plt.legend(loc='center left', bbox_to_anchor=(1.15, 0.5))
plt.title(f'Target positions in the sky with respect to the observer’s location\n({time_obs_night[0]} to {time_obs_night[-1]})', fontsize=12)
plt.show()




