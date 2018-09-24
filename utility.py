

# *******************
def constant(arg):
    constants = {"clight": 2.99792458e10,
                 "hplanck_eV*s": 4.135667662e-15}
    return constants[arg]


# ****************
# convert micron to Hz
def micron_to_hz(xdata):
    import numpy as np
    return constant("clight") / (np.array(xdata) * 1e-4)


# *******************
def hz_to_ev(xdata):
    return constant("hplanck_eV*s") * xdata
