from scimath.units import *

import properties as prop
reload(prop)

air = prop.ideal_gas()
air.T = 300. * temperature.K
air.P = 101. * pressure.kPa
air.set_TempPres_dependents()
