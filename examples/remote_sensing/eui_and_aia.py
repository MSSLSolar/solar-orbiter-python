"""
Carrington map reprojection
===========================
"""
###############################################################################
# Importing required modules.
import aiapy.calibrate
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
import sunpy_soar
import sunpy.map
import sunpy.sun.constants
from sunpy.net import Fido

from sunpy.net.attrs import Instrument, Level, Time, Wavelength
from sunpy_soar.attrs import Identifier

###############################################################################
# Create search attributes.
instrument = Instrument('EUI')
time = Time('2021-02-01', '2021-02-02')
level = Level(2)
identifier = Identifier('EUI-FSI174-IMAGE')

###############################################################################
# Do search.
result = Fido.search(instrument, time, level, identifier)
print(result)

###############################################################################
# Download first file.
files = Fido.fetch(result[0, 0])
print(files)

eui_map = sunpy.map.Map(files[0])


###############################################################################
# Get AIA map
instrument = Instrument('AIA')
time = Time('2021-02-01', '2021-02-02', eui_map.date)
wavelength = Wavelength(171*u.Angstrom)

result = Fido.search(instrument, time, wavelength)
print(result)

files = Fido.fetch(result[0, 0])
print(files)

aia_map = sunpy.map.Map(files[0])
aia_map = aiapy.calibrate.normalize_exposure(aia_map)
# This is an empirical correction factor to make the two maps have equal
# histograms in quiet Sun regions
aia_map = sunpy.map.Map(aia_map.data * 126 / 98, aia_map.meta)
aia_map.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)

###############################################################################
# Plot maps side by side
fig = plt.figure()
ax = fig.add_subplot(121, projection=eui_map)
eui_map.plot(axes=ax, norm=aia_map.plot_settings['norm'], vmin=0)
ax = fig.add_subplot(122, projection=aia_map)
aia_map.plot(axes=ax, vmin=0)

###############################################################################
# Cereate output header

# This is set deliberately low to reduce memory consumption
shape_out = (360, 720)

header = sunpy.map.make_fitswcs_header(shape_out,
                                       SkyCoord(0, 0, unit=u.deg,
                                                frame="heliographic_stonyhurst",
                                                obstime=eui_map.date),
                                       scale=[180 / shape_out[0],
                                              360 / shape_out[1]] * u.deg / u.pix,
                                       projection_code="CAR")
out_wcs = WCS(header)

###############################################################################
# Reproject
array, footprint = reproject_and_coadd([eui_map, aia_map], out_wcs, shape_out,
                                       reproject_function=reproject_interp)

outmap = sunpy.map.Map((array, header))
outmap.plot_settings = aia_map.plot_settings

###############################################################################
# Plot the reprojected map
fig = plt.figure()
ax = fig.add_subplot(projection=outmap)
outmap.plot(axes=ax)
plt.show()
