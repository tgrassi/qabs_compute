from material import Material
from optical import Optical


class QabsManager:

    # ******************
    def __init__(self):
        self.materials = dict()
        self.opticals = dict()

    # ********************
    def load_material(self, fname, labs=None, name=None, force=False, units="micron"):
        import sys

        # default labs input
        default_labs = ["wlen", "real_eps1", "im_eps", "real_m1", "im_m"]

        # use default is labs is not set
        if labs is None:
            labs = default_labs

        # generate name automatically
        if name is None:
            name = "material_" + str(len(self.materials))

        if name in self.materials and not force:
            print "ERROR: material " + name + " already present!"
            print "Change name or use force=True"
            sys.exit()

        mat = Material(fname, labs, name, units)
        opt = Optical(mat)

        # add materials to dictionary
        self.materials[name] = mat
        self.opticals[name] = opt

        return opt

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