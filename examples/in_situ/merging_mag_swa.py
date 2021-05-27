"""
Merging magnetic field and plasma data
======================================

How to merge magnetic field and plasma data on to a set of common time stamps.
"""
###############################################################################
# First, import the required packages.
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from heliopy.data.solo import download
import plasmapy.formulary
from sunpy.timeseries import GenericTimeSeries

###############################################################################
# Set some nice default plotting settings
from astropy.visualization import quantity_support
quantity_support()
matplotlib.rcParams['date.converter'] = 'concise'

###############################################################################
# Download one day of MAG data in normal mode.
start_time = "2020-07-07"
end_time = "2020-07-08"
mag_data = download(start_time, end_time, 'MAG-RTN-NORMAL', 'L2')
print(mag_data.columns)

###############################################################################
# Download one day of SWA/PAS ground moments.
swa_data = download(start_time, end_time, 'SWA-PAS-GRND-MOM', 'L2')
print(swa_data.columns)

###############################################################################
# Re-index the magnetic field data on to the plasma data timestamps.
# To do this the magnetic field `~sunpy.timeseries.GenericTimeSeries` is first
# converted to a `~pandas.DataFrame`, so the powerful time-series maniuplation
# tools present in `pandas` can be used. Then the data is re-indexed using
# a nearest neighbour method, and re-inserted into a
# `~sunpy.timeseries.GenericTimeSeries` to preserve the unit information
mag_df = mag_data.to_dataframe()
mag_reindexed = mag_data.to_dataframe().reindex(swa_data.index,
                                                method='nearest')
mag_reindexed = GenericTimeSeries(mag_reindexed, units=mag_data.units)

###############################################################################
# Plot the re-indexed data. Because the SWA data is a lower cadence than the
# MAG data, the re-indexed data has fewer data points.
fig, ax = plt.subplots()
ax.plot(mag_data.index, mag_data.quantity('B_RTN_0'),
        label='Original')
ax.plot(mag_reindexed.index, mag_reindexed.quantity('B_RTN_0'),
        label='Re-indexed')
ax.legend()

###############################################################################
# Now the plasma and magnetic field data are on similar timestamps, they can
# be combined into derived parameters, such as the plasma beta.
#
# First the magnetic field strength is calculated,
# and then the `plasmapy.formularly.beta` function is used to calculate the
# beta from the temperature, density, and magnetic field strength.
mod_B = np.sqrt(mag_reindexed.quantity('B_RTN_0')**2 +
                mag_reindexed.quantity('B_RTN_1')**2 +
                mag_reindexed.quantity('B_RTN_2')**2)
plasma_beta = plasmapy.formulary.beta(swa_data.quantity('T'),
                                      swa_data.quantity('N'),
                                      mod_B)

###############################################################################
# Finally, the plasma beta along with the
fig, axs = plt.subplots(nrows=4, sharex=True)
axs[0].plot(swa_data.index, swa_data.quantity('T'), label='Temperature')
axs[1].plot(swa_data.index, swa_data.quantity('N'), label='Density')
axs[2].plot(swa_data.index, mod_B, label='Magnetic field strength')

axs[3].plot(swa_data.index, plasma_beta, label=r'Plasma $\beta$')
axs[3].set_yscale('log')
axs[3].axhline(1, color='black', linewidth=1)

[ax.legend() for ax in axs]
plt.show()
