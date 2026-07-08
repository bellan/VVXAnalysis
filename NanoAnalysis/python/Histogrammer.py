###
# Class to bookeep histrograms
###

from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import ROOT
import array


class Histogrammer:
    def __init__(self):
        self.plots = {}
        self._profile = False
        self._category = None

        ROOT.TH1.SetDefaultSumw2(True)
        ROOT.TH1.AddDirectory(False)

    # =========================
    # Utilities
    # =========================
    def get(self, name):
        if name not in self.plots:
            raise KeyError(f"Histogram '{name}' not found")
        return self.plots[name]

    def clone(self, name, newname):
        if newname in self.plots:
            print(f"{newname} already exists")
            return
        self.plots[newname] = self.get(name).Clone(str(newname))

    def erase(self, name):
        return self.plots.pop(name, None) is None

    # =========================
    # BOOK 1D
    # =========================
    def book1D(self, name, title, nbins, xmin, xmax, hist_class=ROOT.TH1D):
        if name in self.plots:
            return self.plots[name]

        h = hist_class(
            str(name), str(title),
            int(nbins), float(xmin), float(xmax)
        )

        self.plots[name] = h
        return h

    def book1D_var(self, name, title, xbins, hist_class=ROOT.TH1D):
        if name in self.plots:
            return self.plots[name]

        xarr = array.array('d', xbins)

        h = hist_class(
            str(name), str(title),
            len(xbins) - 1, xarr
        )

        self.plots[name] = h
        return h

    def book1D_labels(self, name, title, labels, hist_class=ROOT.TH1D):
        if name in self.plots:
            return self.plots[name]

        h = hist_class(
            str(name), str(title),
            len(labels), 0, len(labels)
        )

        for i, lab in enumerate(labels):
            h.GetXaxis().SetBinLabel(i + 1, str(lab))

        self.plots[name] = h
        return h

    # =========================
    # FILL 1D
    # =========================
    def fill1D(self, name, title, nbins, xmin, xmax, value, weight=1.0):
        if self.profile_enabled():
            h = self.book2D(
                name, title,
                64, 0, 64,
                nbins, xmin, xmax
            )
            h.Fill(float(self._category), float(value), float(weight))
        else:
            h = self.book1D(name, title, nbins, xmin, xmax)
            h.Fill(float(value), float(weight))

    def fill1D_var(self, name, title, xbins, value, weight=1.0):
        h = self.book1D_var(name, title, xbins)
        h.Fill(float(value), float(weight))

    def fill1D_label(self, name, title, labels, value, weight=1.0):
        h = self.book1D_labels(name, title, labels)
        h.Fill(str(value), float(weight))

    # =========================
    # BOOK 2D
    # =========================
    def book2D(self, name, title,
               nbinsx, xmin, xmax,
               nbinsy, ymin, ymax,
               hist_class=ROOT.TH2D):

        if name in self.plots:
            return self.plots[name]

        h = hist_class(
            str(name), str(title),
            int(nbinsx), float(xmin), float(xmax),
            int(nbinsy), float(ymin), float(ymax)
        )

        self.plots[name] = h
        return h

    def book2D_var(self, name, title, xbins, ybins, hist_class=ROOT.TH2D):
        if name in self.plots:
            return self.plots[name]

        xarr = array.array('d', xbins)
        yarr = array.array('d', ybins)

        h = hist_class(
            str(name), str(title),
            len(xbins) - 1, xarr,
            len(ybins) - 1, yarr
        )

        self.plots[name] = h
        return h

    # =========================
    # FILL 2D
    # =========================
    def fill2D(self, name, title,
               nbinsx, xmin, xmax,
               nbinsy, ymin, ymax,
               xvalue, yvalue,
               weight=1.0):

        h = self.book2D(
            name, title,
            nbinsx, xmin, xmax,
            nbinsy, ymin, ymax
        )

        h.Fill(float(xvalue), float(yvalue), float(weight))

    def fill2D_var(self, name, title,
                   xbins, ybins,
                   xvalue, yvalue,
                   weight=1.0):

        h = self.book2D_var(name, title, xbins, ybins)
        h.Fill(float(xvalue), float(yvalue), float(weight))

    # =========================
    # PROFILE
    # =========================
    def set_profile(self, category):
        self._category = float(category)
        self._profile = True

    def profile_enabled(self):
        return (
            self._profile
            and self._category is not None
            and self._category >= 0
        )

    # =========================
    # WRITE
    # =========================
    def write(self, fout):
        fout.cd()
        for h in self.plots.values():
            h.Write()


class Histogrammers:
    def __init__(self, regions):
        self.histogrammers = {Regions[region]: Histogrammer() for region in regions}
        
    def checkRegions(self, region_word):
        active = [h for region, h in self.histogrammers.items()
                  if region_word & region == region]
        return _ActiveHistogrammers(active)

    def write(self, base_odir, analysis_name, sample_name):
        for region, h in self.histogrammers.items():
            odir = f"{base_odir}/{analysis_name}_{region.name}"
            os.makedirs(odir, exist_ok=True)
            outFile = ROOT.TFile.Open(f"{odir}/{sample_name}.root", "recreate")
            h.write(outFile)
            outFile.Close()


class _ActiveHistogrammers:
    def __init__(self, active):
        self._active = active

    def __getattr__(self, method_name):
        def dispatch(*args, **kwargs):
            for h in self._active:
                getattr(h, method_name)(*args, **kwargs)
        return dispatch
