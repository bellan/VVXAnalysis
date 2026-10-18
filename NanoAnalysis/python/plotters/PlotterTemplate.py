#!/usr/bin/env python3

import os
import ROOT

from array import array


# ============================================================
# CONFIGURATION
# ============================================================

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


# Input directory
input_dir = "../pathToRoot/"

# Output directory
output_dir = "../plots/"

# Input ROOT files
processes = {
    "signal": {
        "file": "sig.root",
        "label": "VZ#gamma",
        "color": ROOT.kRed,
        "is_signal": True,
    },

    "DY": {
        "file": "bkg1.root",
        "label": "Drell-Yan",
        "color": ROOT.kAzure + 1,
        "is_signal": False,
    },

    "ZG": {
        "file": "bkg2.root",
        "label": "Z#gamma",
        "color": ROOT.kOrange,
        "is_signal": False,
    },
}


# Variables to plot
variables = {
    "mllG": {
        "label": "m_{ll#gamma} [GeV]",
        "rebin": 1,
        "logy": False,
        "ymin": 0.0,
        "ymax": 10000,
        "sensitivity_max": 1.0,
    },

    "ZMass": {
        "label": "m_{ll} [GeV]",
        "rebin": 1,
        "logy": False,
        "ymin": 0.0,
        "ymax": 10000,
        "sensitivity_max": 1.0,
    },
}


# CMS plot labels
cms_label = "Work in progress"
lumi_label = "110 fb^{-1} (13.6 TeV)"

# Sensitivity plot configuration
sensitivity_min = 0.0

# Output formats
output_formats = ["png", "pdf"]


# ============================================================
# HISTOGRAM HANDLING
# ============================================================

def get_histogram(root_file, hist_name, new_name):
    """
    Retrieve an histogram from a ROOT file and return a detached clone.

    The clone is detached from the input file so that it remains
    available after the file is closed.
    """

    hist = root_file.Get(hist_name)

    if not hist:
        raise RuntimeError(
            f"Histogram '{hist_name}' not found "
            f"in file '{root_file.GetName()}'"
        )

    hist = hist.Clone(new_name)
    hist.SetDirectory(0)

    return hist


def configure_histogram(hist, color):
    """
    Configure the appearance of a histogram.
    """

    hist.SetFillColor(color)
    hist.SetLineColor(color + 2)
    hist.SetLineWidth(1)

    return hist


# ============================================================
# PLOT PRODUCTION
# ============================================================

def make_plot(var_name, var_config, histograms):

    print(f"Producing plot: {var_name}")

    var_label = var_config["label"]
    rebin = var_config["rebin"]
    logy = var_config["logy"]

    # --------------------------------------------------------
    # Rebin histograms
    # --------------------------------------------------------

    for hist in histograms.values():
        if rebin > 1:
            hist.Rebin(rebin)

    # --------------------------------------------------------
    # Build signal and background histograms
    # --------------------------------------------------------

    signal_hist = histograms["signal"]

    background_hists = [
        histograms[name]
        for name, config in processes.items()
        if not config["is_signal"]
    ]

    total_background = background_hists[0].Clone(
        f"{var_name}_total_background"
    )
    total_background.SetDirectory(0)

    for hist in background_hists[1:]:
        total_background.Add(hist)

    total_mc = total_background.Clone(f"{var_name}_total_mc")
    total_mc.SetDirectory(0)
    total_mc.Add(signal_hist)

    # --------------------------------------------------------
    # Create uncertainty band
    # --------------------------------------------------------

    mc_uncertainty = total_mc.Clone(f"{var_name}_mc_uncertainty")
    mc_uncertainty.SetDirectory(0)

    mc_uncertainty.SetFillColor(ROOT.kBlack)
    mc_uncertainty.SetFillStyle(3005)
    mc_uncertainty.SetMarkerStyle(1)
    mc_uncertainty.SetLineColor(ROOT.kBlack)
    mc_uncertainty.SetLineWidth(1)

    # --------------------------------------------------------
    # Canvas and pads
    # --------------------------------------------------------

    canvas = ROOT.TCanvas(
        f"canvas_{var_name}",
        var_name,
        800,
        1000,
    )

    main_pad = ROOT.TPad(
        f"main_pad_{var_name}",
        "Main pad",
        0.0, 0.30, 0.95, 1.0,
    )

    ratio_pad = ROOT.TPad(
        f"ratio_pad_{var_name}",
        "Sensitivity pad",
        0.0, 0.0, 0.95, 0.30,
    )

    main_pad.SetBottomMargin(0.10)
    ratio_pad.SetTopMargin(0.05)
    ratio_pad.SetBottomMargin(0.30)

    if logy:
        main_pad.SetLogy()

    main_pad.Draw()
    ratio_pad.Draw()

    # --------------------------------------------------------
    # Main pad: stacked histograms
    # --------------------------------------------------------

    main_pad.cd()

    stack = ROOT.THStack(f"stack_{var_name}", "")

    # Backgrounds first, signal last
    for name, config in processes.items():
        hist = histograms[name]

        stack.Add(hist)

    stack.Draw("HIST")

    stack.GetYaxis().SetTitle("Events")
    stack.GetXaxis().SetTitle(var_label)

    if logy:
        stack.SetMinimum(max(var_config["ymin"], 0.1))
    else:
        stack.SetMinimum(var_config["ymin"])

    stack.SetMaximum(var_config["ymax"])

    mc_uncertainty.Draw("E2 SAME")

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    legend = ROOT.TLegend(0.65, 0.67, 0.94, 0.89)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    for name, config in processes.items():
        legend.AddEntry(
            histograms[name],
            config["label"],
            "f",
        )

    legend.AddEntry(
        mc_uncertainty,
        "MC stat. unc.",
        "f",
    )

    legend.Draw()

    # --------------------------------------------------------
    # CMS labels
    # --------------------------------------------------------

    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextFont(62)
    cms_text.SetTextSize(0.05)

    cms_text.DrawLatex(
        0.10, 0.94,
        "CMS",
    )

    cms_extra = ROOT.TLatex()
    cms_extra.SetNDC()
    cms_extra.SetTextFont(52)
    cms_extra.SetTextSize(0.04)

    cms_extra.DrawLatex(
        0.22, 0.94,
        cms_label,
    )

    lumi_text = ROOT.TLatex()
    lumi_text.SetNDC()
    lumi_text.SetTextFont(42)
    lumi_text.SetTextSize(0.04)

    lumi_text.SetTextAlign(32)

    lumi_text.DrawLatex(
        0.94, 0.94,
        lumi_label,
    )

    # --------------------------------------------------------
    # Lower pad: signal sensitivity
    # --------------------------------------------------------

    ratio_pad.cd()

    # Create histograms
    sensitivity = signal_hist.Clone(f"{var_name}_sensitivity")
    sensitivity.SetDirectory(0)
    sensitivity.Reset()

    sensitivity_fwd = signal_hist.Clone(f"{var_name}_sensitivity_fwd")
    sensitivity_fwd.SetDirectory(0)
    sensitivity_fwd.Reset()

    sensitivity_bwd = signal_hist.Clone(f"{var_name}_sensitivity_bwd")
    sensitivity_bwd.SetDirectory(0)
    sensitivity_bwd.Reset()

    # Cumulative signal and background yields
    signal_cumulative_fwd = signal_hist.GetCumulative()
    signal_cumulative_fwd.SetDirectory(0)

    background_cumulative_fwd = total_background.GetCumulative()
    background_cumulative_fwd.SetDirectory(0)

    signal_cumulative_bwd = signal_hist.GetCumulative(False)
    signal_cumulative_bwd.SetDirectory(0)

    background_cumulative_bwd = total_background.GetCumulative(False)
    background_cumulative_bwd.SetDirectory(0)

    # Compute sensitivities
    for i_bin in range(1, signal_hist.GetNbinsX() + 1):

        # Bin-by-bin yields
        s = signal_hist.GetBinContent(i_bin)
        b = total_background.GetBinContent(i_bin)

        if s + b > 0:
            sensitivity.SetBinContent(i_bin, s / (s + b) ** 0.5)

        # Forward cumulative yields: bins 1 -> i_bin
        s_fwd = signal_cumulative_fwd.GetBinContent(i_bin)
        b_fwd = background_cumulative_fwd.GetBinContent(i_bin)

        if s_fwd + b_fwd > 0:
            sensitivity_fwd.SetBinContent(
                i_bin,
                s_fwd / (s_fwd + b_fwd) ** 0.5,
            )

        # Backward cumulative yields: bins i_bin -> N
        s_bwd = signal_cumulative_bwd.GetBinContent(i_bin)
        b_bwd = background_cumulative_bwd.GetBinContent(i_bin)

        if s_bwd + b_bwd > 0:
            sensitivity_bwd.SetBinContent(
                i_bin,
                s_bwd / (s_bwd + b_bwd) ** 0.5,
            )

    # --------------------------------------------------------
    # Style
    # --------------------------------------------------------

    sensitivity.SetLineColor(ROOT.kBlack)
    sensitivity.SetMarkerColor(ROOT.kBlack)
    sensitivity.SetMarkerStyle(20)
    sensitivity.SetMarkerSize(0.8)

    sensitivity_fwd.SetLineColor(ROOT.kBlue + 1)
    sensitivity_fwd.SetMarkerColor(ROOT.kBlue + 1)
    sensitivity_fwd.SetMarkerStyle(21)
    sensitivity_fwd.SetMarkerSize(0.8)

    sensitivity_bwd.SetLineColor(ROOT.kRed + 1)
    sensitivity_bwd.SetMarkerColor(ROOT.kRed + 1)
    sensitivity_bwd.SetMarkerStyle(22)
    sensitivity_bwd.SetMarkerSize(0.8)

    # --------------------------------------------------------
    # Axes
    # --------------------------------------------------------

    sensitivity.SetTitle("")
    sensitivity.GetYaxis().SetTitle("S/#sqrt{S+B}")
    sensitivity.GetXaxis().SetTitle(var_label)

    sensitivity.GetYaxis().SetRangeUser(
        sensitivity_min,
        var_config["sensitivity_max"],
    )

    sensitivity.GetYaxis().SetTitleSize(0.08)
    sensitivity.GetXaxis().SetTitleSize(0.08)

    sensitivity.GetYaxis().SetLabelSize(0.06)
    sensitivity.GetXaxis().SetLabelSize(0.06)

    sensitivity.GetYaxis().SetTitleOffset(0.5)

    # --------------------------------------------------------
    # Draw
    # --------------------------------------------------------

    sensitivity.Draw("E1")
    sensitivity_fwd.Draw("E1 SAME")
    sensitivity_bwd.Draw("E1 SAME")

    # Legend
    sensitivity_legend = ROOT.TLegend(0.55, 0.72, 0.94, 0.92)
    sensitivity_legend.SetBorderSize(0)
    sensitivity_legend.SetFillStyle(0)

    sensitivity_legend.AddEntry(
        sensitivity,
        "Bin by bin",
        "lep",
    )

    sensitivity_legend.AddEntry(
        sensitivity_fwd,
        "Cumulative forward",
        "lep",
    )

    sensitivity_legend.AddEntry(
        sensitivity_bwd,
        "Cumulative backward",
        "lep",
    )

    sensitivity_legend.Draw()
    
    # --------------------------------------------------------
    # Save plots
    # --------------------------------------------------------

    canvas.cd()
    canvas.Update()

    os.makedirs(output_dir, exist_ok=True)

    for extension in output_formats:

        output_file = os.path.join(
            output_dir,
            f"{var_name}.{extension}",
        )

        canvas.SaveAs(output_file)

        print(f"Saved: {output_file}")

    # Keep references alive until the canvas has been saved.
    return canvas


# ============================================================
# MAIN
# ============================================================

def main():

    # Open input files
    root_files = {}

    for name, config in processes.items():

        file_path = os.path.join(
            input_dir,
            config["file"],
        )

        root_file = ROOT.TFile.Open(file_path, "READ")

        if not root_file or root_file.IsZombie():
            raise RuntimeError(
                f"Cannot open ROOT file: {file_path}"
            )

        root_files[name] = root_file

    canvases = []

    try:

        for var_name, var_config in variables.items():

            histograms = {}

            for name, config in processes.items():

                histograms[name] = get_histogram(
                    root_files[name],
                    var_name,
                    f"{var_name}_{name}",
                )

                configure_histogram(
                    histograms[name],
                    config["color"],
                )

            canvas = make_plot(
                var_name,
                var_config,
                histograms,
            )

            canvases.append(canvas)

    finally:

        for root_file in root_files.values():
            root_file.Close()


if __name__ == "__main__":
    main()
