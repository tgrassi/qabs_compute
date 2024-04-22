

# *******************
def constant(constant_name):
    constants = {"clight": 2.99792458e10,
                 "hplanck_eV*s": 4.135667662e-15}
    return constants[constant_name]


# ****************
# convert micron to Hz
def micron_to_hz(xdata):
    import numpy as np
    return constant("clight") / (np.array(xdata) * 1e-4)


# *******************
def hz_to_ev(xdata):
    import numpy as np
    return constant("hplanck_eV*s") * np.array(xdata)


# *******************
def ev_to_micron(xdata):
    import numpy as np
    return constant("hplanck_eV*s") * constant("clight") / np.array(xdata) * 1e4


# **********************
def wavenumber_to_micron(xdata):
    import numpy as np
    return 1e4 / np.array(xdata)


# **********************
def thz_to_micron(xdata):
    import numpy as np
    xdata = np.array(xdata) * 1e12  # THz to Hz
    return constant("clight") / xdata * 1e4  # Hz to micron