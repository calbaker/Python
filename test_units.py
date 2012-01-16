from scimath.units import *
from scimath.units.api import *

import properties as prop
reload(prop)

air = prop.ideal_gas()
air.T = UnitScalar(300., units=temperature.K)
air.P = UnitScalar(101., units=pressure.kPa)
air.c_p_air = air.get_c_p_air(air.T)
air.set_TempPres_dependents()
