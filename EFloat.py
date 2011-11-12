import math

# representation when printed
# Convert a string or a number to a floating point number, if
# possible.
# When asked to render itself for printing (via str() or print)
# it is normalized
# to engineering style notation at powers of 10 in multiples of 3
# (for micro, milli, kilo, mega, giga, etc.)

class EFloat(float):
    def __init__(self, value=0.0, prec=12):
        super(EFloat, self).__init__(value)
        self.precision = prec

    def _get_precision(self):
        return self._precision

    def _set_precision(self, p):
        self._precision = p
        self.format_string = "%3." + ("%d" % self._precision) + "fe% +d"
        return
        precision = property(_get_precision, _set_precision, doc="The number of decimal places printed")

    def _exponent(self):
        if self == 0.0:
            ret = 0
        else:
            ret = math.floor(math.log10(abs(self)))
        return ret

    def _mantissa(self):
        return self/math.pow(10, self._exponent())
    
    def _asEng(self):
        shift = self._exponent() % 3
        retval = self.format_string % (self._mantissa()*math.pow(10,shift), self._exponent() - shift)
        return retval

    def __str__(self):
        return self._asEng()

    def __repr__(self):
        return str(self)

    def __add__(self, x):
        return EFloat(float.__add__(self, float(x)))

    def __radd__(self, x):
        return EFloat(float.__add__(self, float(x)))

    def __mul__(self, x):
        return EFloat(float.__mul__(self, float(x)))

    def __rmul__(self, x):
        return EFloat(float.__mul__(self, float(x)))

    def __sub__(self, x):
        return EFloat(float.__sub__(self, float(x)))
    
    def __rsub__(self, x):
        return EFloat(float.__rsub__(self, float(x)))
    
    def __div__(self, x):
        return EFloat(float.__div__(self, float(x)))
