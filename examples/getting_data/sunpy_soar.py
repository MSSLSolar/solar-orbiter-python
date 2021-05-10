"""
Searching the Solar Orbiter archive (SOAR)
==========================================
"""
###############################################################################
# Importing required modules.
import sunpy_soar
from sunpy.net import Fido
from sunpy.net.attrs import Instrument, Level, Time
from sunpy_soar.attrs import Identifier

###############################################################################
# Create search attributes.
instrument = Instrument('EUI')
time = Time('2021-02-01', '2021-02-02')
level = Level(1)
identifier = Identifier('EUI-FSI174-IMAGE')

###############################################################################
# Do search.
result = Fido.search(instrument, time, level, identifier)
print(result)

###############################################################################
# Download first file.
files = Fido.fetch(result[0, 0])
print(files)
