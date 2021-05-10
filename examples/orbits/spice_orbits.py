"""
Orbits from SPICE kernels
=========================
Generating orbits from SPICE kernels.

In this example we'll download the Solar Orbiter SPICE kernel, and plot
its orbit for the first year.
"""


from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
import heliopy.data.spice as spicedata
import heliopy.spice as spice

###############################################################################
# Load the Solar Orbiter spice kernel. HelioPy will automatically fetch and
# load the latest kernel. We can then inspect the kernel to check what the
# date coverage is for Solar Orbiter.
kernels = spicedata.get_kernel('solo')
solo = spice.Trajectory('Solar Orbiter')
coverage = kernels[0].coverage(spice.Body('Solar Orbiter'))
print(coverage)

###############################################################################
# Next we define a set of times at which to sample the orbit.
starttime = coverage[0]
times = [starttime + timedelta(days=i) for i in range(1, 5 * 365)]

###############################################################################
# Generate positions. "IAU_SUN" is a Carrington frame of reference.
solo.generate_positions(times, 'Sun', 'IAU_SUN')
coords = solo.coords

###############################################################################
# Plot radial distance and elevation as a function of time
fig, axs = plt.subplots(3, 1, sharex=True)
for frame in ('heliographic_carrington', 'heliographic_stonyhurst'):
    coords = coords.transform_to(frame)
    axs[0].plot(times, coords.radius.to(u.au))
    axs[0].set_ylim(0, 1.1)
    axs[0].set_ylabel('r (AU)')

    axs[1].plot(times, coords.lat.to(u.deg))
    axs[1].set_ylabel('Latitude (deg)')

    axs[2].plot(times, coords.lon.to(u.deg), label=frame)
    axs[2].set_ylabel('Longitude (deg)')

axs[2].set_xlabel('Year')
axs[2].legend()
axs[0].set_title('Solar Orbiter orbit')
fig.align_ylabels()
plt.show()
