"""
Combining EUI and AIA maps
==========================
How to combine EUI and AIA maps to improve coverage of the Sun at any one time.

Since Solar Orbiter is not positioned at the Earth, it is possible to combine
remote sensing images from Solar Orbiter and Earth-based assets
(in this example, SDO/AIA) to improve the fraction of the Sun we can see at
any one time.

In this example we take an EUI 174A full sun image, and combine it with
an SDO/AIA 171A image to create an almost-complete map of the Sun. This
is done by reprojecting each map into a Carrington frame of reference,
and adding them together.
"""
###############################################################################
# First import the required modules.
import aiapy.calibrate
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
import sunpy_soar
from sunpy.coordinates import HeliographicCarrington, HeliographicStonyhurst
import sunpy.map
import sunpy.sun.constants
from sunpy.net import Fido

from sunpy.net.attrs import Instrument, Level, Time, Wavelength
from sunpy_soar.attrs import Identifier

###############################################################################
# We'll start by searching for an EUI image. This is done using the sunpy_soar
# package. We define the search attributes, do the search, and print the
# result.
instrument = Instrument('EUI')
time = Time('2021-02-01', '2021-02-02')
level = Level(2)
identifier = Identifier('EUI-FSI174-IMAGE')

result = Fido.search(instrument, time, level, identifier)
print(result)

###############################################################################
# The search returned a number of results. We'll just download the first one
# here, print the local path that it is downloaded to, and load it into a
# sunpy map.
files = Fido.fetch(result[0, 0])
print(files)
eui_map = sunpy.map.Map(files[0])

###############################################################################
# Next we'll search for an AIA map. Because `sunpy.net.Fido` provides a single
# interface for searching different data sources, the code is very similar to
# the code used earlier to get an EUI map.
instrument = Instrument('AIA')
# The third argument here tells sunpy to find a single map closest to that date
time = Time('2021-02-01', '2021-02-02', eui_map.date)
wavelength = Wavelength(171*u.Angstrom)

result = Fido.search(instrument, time, wavelength)
print(result)

files = Fido.fetch(result[0, 0])
print(files)
###############################################################################
# Now we'll load the downloaded AIA file into a sunpy map, normalise the data
# to the exposure time, apply an emperical correction factor
aia_map = sunpy.map.Map(files[0])
aia_map = aiapy.calibrate.normalize_exposure(aia_map)
# This is an empirical correction factor to make the two maps have equal
# histograms in quiet Sun regions
aia_map = sunpy.map.Map(aia_map.data * 126 / 98, aia_map.meta)
aia_map.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)

###############################################################################
# Plot the maps side by side.
fig = plt.figure()
ax = fig.add_subplot(121, projection=eui_map)
eui_map.plot(axes=ax, norm=aia_map.plot_settings['norm'], vmin=0)
ax = fig.add_subplot(122, projection=aia_map)
aia_map.plot(axes=ax, vmin=0)

###############################################################################
# Next we will reproject both maps on to a Carrington frame of reference.
# To do this we start by creating a FITS world coordinat system (WCS) header
# corresponding to the output coordinate frame.

# This is set deliberately low, and (at the cost of memory and processing time)
# can be increased to increase the output resolution
shape_out = (360, 720)
ref_coord_observer = HeliographicStonyhurst(0*u.deg, 0*u.deg, sunpy.sun.constants.radius,
                                            obstime=eui_map.date)
ref_coord = HeliographicCarrington(0*u.deg, 0*u.deg, sunpy.sun.constants.radius,
                                   obstime=eui_map.date,
                                   observer=ref_coord_observer)

header = sunpy.map.make_fitswcs_header(shape_out, ref_coord,
                                       scale=[180 / shape_out[0],
                                              360 / shape_out[1]] * u.deg / u.pix,
                                       projection_code="CAR")
out_wcs = WCS(header)

###############################################################################
# Next we reproject and add together the two maps
array, footprint = reproject_and_coadd([eui_map, aia_map], out_wcs, shape_out,
                                       reproject_function=reproject_interp)

outmap = sunpy.map.Map((array, header))
outmap.plot_settings = aia_map.plot_settings

###############################################################################
# Finally, we'll plot the reprojected map.
fig = plt.figure()
ax = fig.add_subplot(projection=outmap)
outmap.plot(axes=ax)
plt.show()
