from material import Material
from optical import Optical


class QabsManager:

    # ******************
    def __init__(self):
        self.materials = dict()
        self.opticals = dict()

    # ********************
    # load material and return optical object
    def load_material(self, fname, labels=None, name=None, force=False, units="micron"):
        import sys

        # default labs input
        default_labels = ["wlen", "real_eps1", "im_eps", "real_m1", "im_m"]

        # use default is labs is not set
        if labels is None:
            labels = default_labels

        # generate name automatically
        if name is None:
            name = "material_" + str(len(self.materials))

        if name in self.materials and not force:
            print "ERROR: material " + name + " already present!"
            print "Change name or use force=True"
            sys.exit()

        mat = Material(fname, labels, name, units)
        opt = Optical(mat)

        # add materials to dictionary
        self.materials[name] = mat
        self.opticals[name] = opt

        return opt

    # ******************
    def constant_material(self, wmin, wmax, val, ngrid=100, what="eps"):
        import numpy as np
        import sys

        mat = Material()
        mat.data["wlen"] = np.linspace(wmin, wmax, ngrid)

        if what == "eps":
            mat.data["real_eps"] = np.array([val.real] * ngrid)
            mat.data["im_eps"] = np.array([val.imag] * ngrid)
            e1 = mat.data["real_eps"]
            e2 = mat.data["im_eps"]
            mat.data["im_m"] = ((-e1 + np.sqrt(e1 ** 2 + e2 ** 2)) / 2.) ** 0.5
            mat.data["real_m"] = e2 / 2. / mat.data["im_m"]
        elif what == "m":
            mat.data["real_m"] = np.array([val.real] * ngrid)
            mat.data["im_m"] = np.array([val.imag] * ngrid)
            mat.data["real_eps"] = mat.data["real_m"]**2 - mat.data["im_m"]**2
            mat.data["im_eps"] = 2e0 * mat.data["real_m"] * mat.data["im_m"]
        else:
            sys.exit("ERROR: unknown optical property in constant material!")

        opt = Optical(mat)

        # add materials to dictionary
        name = "constant material"
        self.materials[name] = mat
        self.opticals[name] = opt

        return opt

    # *******************
    @staticmethod
    def vacuum_as(opt):
        from copy import deepcopy as copy
        opt_vac = copy(opt)
        opt_vac.materials = [opt_vac.materials[0]]
        opt_vac.materials[0].data["real_m"][:] = 1e0
        opt_vac.materials[0].data["im_m"][:] = 0e0
        opt_vac.materials[0].data["real_eps"][:] = 1e0
        opt_vac.materials[0].data["im_eps"][:] = 0e0
        return opt_vac

    # *****************
    # merge fname1 and fname 2 to fname_out, labs1 and labs2 are
    # the column labels (that must contain wlen assumed to be the
    # independent variable)
    @staticmethod
    def merge_data(fname1, fname2, labs1, labs2, fname_out):
        import numpy as np
        import sys

        # check to avid overwriting
        if fname1 == fname_out or fname2 == fname_out:
            sys.exit("ERROR: input and output file in file merger are suspiciously the same!")

        # check if independent variables in in both sets
        if "wlen" not in labs1 and "wlen" not in labs2:
            sys.exit("ERROR: when merging you must have wlen label!")

        # load data from file
        def load_data(fname, labs):
            data = {ll: [] for ll in labs}
            for row in open(fname):
                srow = row.strip()
                srow = srow.replace("\t", " ").replace(",", " ")
                if srow == "" or srow.startswith("#"):
                    continue
                arow = [float(x) for x in srow.split(" ") if x != ""]
                for ii, ll in enumerate(labs):
                    data[ll].append(arow[ii])
            return {k: np.array(v) for k, v in data.iteritems()}

        # create interpolator from data
        def create_interpolator(data):
            from scipy.interpolate import interp1d
            return {k: interp1d(data["wlen"], v) for k, v in data.iteritems() if k != "wlen"}

        # load data
        data1 = load_data(fname1, labs1)
        data2 = load_data(fname2, labs2)

        # create interpoltors
        interp1 = create_interpolator(data1)
        interp2 = create_interpolator(data2)

        data_merged = []
        wlen_merged = []
        count_total_points = count_interpolated_points = 0
        # interpolate data on the merged independent variables
        for wlen in np.append(data1["wlen"], data2["wlen"]):
            vrow = []
            interp_ok = True

            count_total_points += 1

            # intepolate data1 if wlen in range
            for lab in labs1:
                if lab != "wlen":
                    try:
                        vrow.append(interp1[lab](wlen))
                    except ValueError:
                        interp_ok = False

            # intepolate data2 if wlen in range
            for lab in labs2:
                if lab != "wlen":
                    try:
                        vrow.append(interp2[lab](wlen))
                    except ValueError:
                        interp_ok = False

            # skip non-interpolable data
            if not interp_ok:
                continue

            count_interpolated_points += 1

            # store data to sort and write after
            wlen_merged.append(wlen)
            data_merged.append(vrow)

        # write sorted data to file
        fout = open(fname_out, "w")
        labs_merged = [x for x in labs1 + labs2 if x != "wlen"]
        fout.write("# wlen " + " ".join(labs_merged) + "\n")
        for jj in np.argsort(wlen_merged):
            vrow = [wlen_merged[jj]] + data_merged[jj]
            fout.write(" ".join([str(x) for x in vrow]) + "\n")
        fout.close()

        print fname1 + " merged to " + fname2 + " in " + fname_out
        print "total points %d, interpolated points %d" % (count_total_points,
                                                           count_interpolated_points)

    # *****************
    # create new optical material from others
    def make_optical(self, opticals, name=None):
        import sys

        # check input
        if len(opticals) != 2:
            sys.exit("ERROR: you need 2 materials to create a composite one!")

        # generate a random name
        if name is None:
            name = "optical_" + str(len(self.opticals))

        # loop on opticals to concatenate materials
        materials = []
        for opt in opticals:
            materials += opt.materials

        # store optical to attribute dictionary
        opt = Optical(materials)
        self.opticals[name] = opt
        return opt

    # ********************
    @staticmethod
    def merge_opticals(opticals, method, frac):
        from scipy.interpolate import interp1d
        import numpy as np
        import sys

        if method != "Maxwell-Garnett":
            sys.exit("ERROR: only Maxwell-Garnett method available!")

        if len(opticals) != 2:
            sys.exit("ERROR: you can merge only two opticals!")

        wmin = 0e0
        wmax = 9e99
        wlens = []
        feps1 = []
        feps2 = []
        for optical in opticals:
            if len(optical.materials) != 1:
                sys.exit("ERROR: you can merge only opticals with a single material!")
            wlen = optical.materials[0].data["wlen"]
            wmin = max(wmin, min(wlen))
            wmax = min(wmax, max(wlen))
            wlens.append(wlen)

            eps1 = optical.materials[0].data["real_eps"]
            eps2 = optical.materials[0].data["im_eps"]
            feps1.append(interp1d(wlen, eps1))
            feps2.append(interp1d(wlen, eps2))

        wlen_all = sorted(set(np.concatenate(wlens)))
        wlen_all = [x for x in wlen_all if x >= wmin]
        wlen_all = [x for x in wlen_all if x <= wmax]
        wlen_all = np.array(wlen_all)

        eps1_m = []
        eps2_m = []
        for wlen in wlen_all:
            epsa = feps1[0](wlen) + 1j * feps2[0](wlen)
            epsb = feps1[1](wlen) + 1j * feps2[1](wlen)
            eterm = (epsb - epsa) / (epsb + 2e0 * epsa)
            # eps = epsa * (1e0 + 3e0 * frac * eterm / (1e0 - frac * eterm))
            eps = epsb * (epsa + 2e0 * epsb + 2e0 * frac*(epsa-epsb)) \
                / (epsa + 2e0 * epsb - frac * (epsa - epsb))
            eps1_m.append(eps.real)
            eps2_m.append(eps.imag)

        eps1_m = np.array(eps1_m)
        eps2_m = np.array(eps2_m)

        material = Material()
        material.data["wlen"] = wlen_all
        material.data["real_eps"] = eps1_m
        material.data["im_eps"] = eps2_m
        cconj = np.sqrt(eps1_m ** 2 + eps2_m ** 2)
        m1 = np.sqrt((cconj + eps1_m) / 2e0)
        m2 = np.sqrt((cconj - eps1_m) / 2e0)

        material.data["real_m"] = m1
        material.data["im_m"] = m2
        opt = Optical(material)

        return opt

    # *****************
    @staticmethod
    def merge_kappa(opticals, fractions):
        import sys
        import numpy as np
        from scipy.interpolate import interp1d

        # check if fractions and opticals are the same
        if len(fractions) != len(opticals):
            sys.exit("ERROR: when combining kappa, opticals and fractions should be the "
                     "same number!")

        # merge all the wlen
        wlen_all = np.concatenate([x.dust.data["wlen"] for x in opticals])

        # create interpolators
        interps = [interp1d(x.dust.data["wlen"], x.dust.data["kappa"]) for x in opticals]

        kappa_interp = []
        wlen_interp = []
        # loop to interpolate
        for wlen in wlen_all:
            interp_ok = True
            # loop on interpolators
            interpolated_kappa = 0e0
            for ii, interp in enumerate(interps):
                # try to interpolate otherwise skip
                try:
                    # interpolate at given wavelength
                    interpolated_kappa += interp(wlen) * fractions[ii]
                except ValueError:
                    interp_ok = False

            # store data
            if interp_ok:
                kappa_interp.append(interpolated_kappa)
                wlen_interp.append(wlen)

        # convert to np arrays
        wlen_interp = np.array(wlen_interp)
        kappa_interp = np.array(kappa_interp)

        # sort by wavelength
        idxs = np.argsort(wlen_interp)
        kappa_interp = kappa_interp[idxs]
        wlen_interp = wlen_interp[idxs]

        # create a new optical object and store kappa in its dust object
        # Note: this optical has no material
        combined = Optical(None)
        combined.dust.data["wlen"] = wlen_interp
        combined.dust.data["kappa"] = kappa_interp  # / mass_normalization

        return combined

    # ******************
    # write a database report
    def report(self):
        # loop on standard data and coating
        for name, material in self.materials.iteritems():
            data = material.data
            # if missing data skip
            if "wlen" not in data:
                continue
            # small table with data
            print "*********"
            print "Name:", name
            print "File:", data["fname"]
            print "Number of points:", len(data["wlen"])
            print "Range wavelength (micron): %e, %e" % (min(data["wlen"]), max(data["wlen"]))
            print "Range frequency (Hz): %e, %e" % (min(data["freq"]), max(data["freq"]))

    # ************************
    # perform a benchmark
    def benchmark(self, fname="benchmark_q.png"):

        # load material
        opt = self.load_material("data/eps_Sil.dat", labs=["wlen", "real_eps1", "im_eps",
                                                           "real_m1", "im_m"])

        # compute Qabs for given size, cm
        opt.compute_q(1e-7)
        # plot computed
        opt.plot_q(fname)

        # load Qabs from file to compare
        opt.load_q("data/Sil_21_1e3.dat")
        # over-plot loaded
        opt.add_plot_q(fname, linestyle="--")

    # *******************
    @staticmethod
    def clear_plots():
        import matplotlib.pyplot as plt
        plt.clf()
