"""
Merging magnetic field and plasma data
======================================
"""
import matplotlib.pyplot as plt
from matplotlib import dates as mdates
import numpy as np

from heliopy.data.solo import download
import plasmapy.formulary
from sunpy.timeseries import GenericTimeSeries

###############################################################################
# Download 1 day of MAG data in normal mode.
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
mag_df = mag_data.to_dataframe()
mag_reindexed = mag_data.to_dataframe().reindex(swa_data.index,
                                                method='nearest')
mag_reindexed = GenericTimeSeries(mag_reindexed, units=mag_data.units)

###############################################################################
# Plot the re-indexed data.
fig, ax = plt.subplots()
ax.plot(mag_data.index, mag_data.quantity('B_RTN_0'),
        label='Original')
ax.plot(mag_reindexed.index, mag_reindexed.quantity('B_RTN_0'),
        label='Re-indexed')

###############################################################################
# Calculate plasma beta.
mod_B = np.sqrt(mag_reindexed.quantity('B_RTN_0')**2 +
                mag_reindexed.quantity('B_RTN_1')**2 +
                mag_reindexed.quantity('B_RTN_2')**2)
plasma_beta = plasmapy.formulary.beta(swa_data.quantity('T'),
                                      swa_data.quantity('N'),
                                      mod_B)

###############################################################################
# Plot the plasma beta
fig, ax = plt.subplots()
ax.plot(swa_data.index, plasma_beta, label='$\beta$')
ax.set_yscale('log')
ax.axhline(1, color='black', linewidth=1)
plt.show()
