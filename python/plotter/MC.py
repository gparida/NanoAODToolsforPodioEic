from numbers import Integral
import ROOT
import argparse
import os
import glob
from re import search
from variables_topCR import *


# 8 visually distinct colors for ALP mass signals
SIGNAL_COLORS = [
    ROOT.kRed,
    ROOT.kBlue + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
    ROOT.kCyan + 2,
    ROOT.kViolet + 1,
    ROOT.kPink + 6,
]


class MakeHistograms(object):
    def __init__(self, RootFilePath, RootFileName, userWeight="1.0"):
        self.RootFileName = ROOT.TFile(RootFilePath + RootFileName + '.root')
        self.HistogramName = None
        self.userWeight = userWeight

    def CreateCutString(self, standardCutString, otherCuts, weighting):
        cuts = []
        if standardCutString:
            cuts.append(standardCutString)
        if otherCuts:
            cuts.extend(otherCuts)
        if cuts:
            cut_expr = ' && '.join('(' + c + ')' for c in cuts)
            return weighting + '*(' + cut_expr + ')'
        return weighting

    def StandardDraw(self, theFile, variable, standardCutString, additionalSelections,
                     histogramName, theWeight='1'):
        theTree = theFile.Get('events')
        print("\n\n\nHistogram for... =", histogramName)
        print(variable + '>>' + histogramName + '(' + variableSettingDictionary[variable] + ')',
              self.CreateCutString(standardCutString, additionalSelections, theWeight))
        theTree.Draw(variable + '>>' + histogramName + '(' + variableSettingDictionary[variable] + ')',
                     self.CreateCutString(standardCutString, additionalSelections, theWeight))
        try:
            theHisto = ROOT.gDirectory.Get(histogramName).Clone()
        except ReferenceError:
            theHisto = None
        self.HistogramName = theHisto


def clubHistograms(lst, histObjects):
    clubHist = None
    for name in lst:
        if histObjects[name].HistogramName is not None:
            if clubHist is None:
                clubHist = histObjects[name].HistogramName.Clone()
                continue
            clubHist.Add(histObjects[name].HistogramName)
    return clubHist


def MakeStackErrors(theStack):
    denominatorHistos = theStack.GetHists().At(0).Clone()
    denominatorHistos.Reset()
    for i in range(0, theStack.GetNhists()):
        denominatorHistos.Add(theStack.GetHists().At(i))
    theErrorHisto = denominatorHistos.Clone()
    theErrorHisto.Reset()
    for i in range(0, denominatorHistos.GetNbinsX() + 1):
        theErrorHisto.SetBinContent(i, denominatorHistos.GetBinContent(i))
        theErrorHisto.SetBinError(i, denominatorHistos.GetBinError(i))
    theErrorHisto.SetLineColor(0)
    theErrorHisto.SetLineWidth(0)
    theErrorHisto.SetMarkerStyle(0)
    theErrorHisto.SetFillStyle(3001)
    theErrorHisto.SetFillColor(15)
    return theErrorHisto


def alp_mass_label(filename):
    """Extract ALP mass from filename like 'ma_0.1_hadded' → 'm_{a}=0.1 GeV'."""
    m = search(r'ma_([\d.]+)', filename)
    return "m_{{a}}={} GeV".format(m.group(1)) if m else filename


def main():
    parser = argparse.ArgumentParser(description='Generate EIC good-electron control plots.')
    parser.add_argument('--batchMode', help='run in batch mode', action='store_true')
    parser.add_argument('--variables', nargs='+',
                        default=["gElectron_pt[0]", "gElectron_eta[0]", "gElectron_phi[0]",
                                 "gElectron_E[0]", "ngElectron",
                                 "ggenElectron_pt[0]", "ggenElectron_eta[0]"],
                        help='Variables to draw control plots for')
    parser.add_argument('--additionalSelections', '-C2', nargs='+',
                        default=["ngElectron>0"],
                        help='Additional selection cuts applied to all samples')
    parser.add_argument('--pause', help='pause after each plot', action='store_true')
    parser.add_argument('--standardCutString', '-C1', nargs='?',
                        default="", help='Base cut string (C1)')
    parser.add_argument('--changeHistogramBounds', nargs='?',
                        help='Override binning for all histograms (nbin,lo,hi)')
    parser.add_argument('--logScale', help='log scale Y axis', action='store_true')
    parser.add_argument('--Sub', required=True, help='Output subdirectory label')
    parser.add_argument('--Path', required=True, help='Path to input ROOT files (with trailing /)')
    parser.add_argument('--Weight', default='1',
                        help='Event weight expression (default: 1)')
    parser.add_argument('--lumi', default='ePIC Simulation, 10#times110 GeV',
                        help='Label printed at top right of canvas')
    parser.add_argument('--normalize', help='normalize all histograms to unit area',
                        action='store_true')
    args = parser.parse_args()

    ROOT.gStyle.SetOptStat(0)
    if args.batchMode:
        ROOT.gROOT.SetBatch(ROOT.kTRUE)

    fnames = glob.glob(args.Path + "*.root")
    print("Found files:", fnames)

    BkgNameList = []
    SignalNameList = []

    for file in fnames:
        filename = os.path.splitext(os.path.basename(file))[0]
        if search(r'ma_', filename):
            SignalNameList.append(filename)
        else:
            BkgNameList.append(filename)

    print("Background files:", BkgNameList)
    print("Signal files:    ", SignalNameList)
    print()

    for variable in args.variables:
        if variable not in variableSettingDictionary:
            print("No histogram settings for variable:", variable)
            continue
        if variable not in variableAxisTitleDictionary:
            print("No axis title for variable:", variable)
            continue

        if args.changeHistogramBounds is not None:
            variableSettingDictionary[variable] = args.changeHistogramBounds

        # --- Draw background histograms ---
        BkgObjects = {}
        for name in BkgNameList:
            BkgObjects[name] = MakeHistograms(args.Path, name)
        for name in BkgNameList:
            BkgObjects[name].StandardDraw(BkgObjects[name].RootFileName, variable,
                                          args.standardCutString, args.additionalSelections,
                                          name, theWeight=args.Weight)

        # --- Draw signal histograms ---
        SignalObjects = {}
        for name in SignalNameList:
            SignalObjects[name] = MakeHistograms(args.Path, name)
        for name in SignalNameList:
            SignalObjects[name].StandardDraw(SignalObjects[name].RootFileName, variable,
                                             args.standardCutString, args.additionalSelections,
                                             name, theWeight=args.Weight)

        # --- Combine all backgrounds into QED Compton stack ---
        QEDCompton_Histo = clubHistograms(BkgNameList, BkgObjects)

        # --- Collect valid signal histos in sorted mass order ---
        def _mass_key(name):
            m = search(r'ma_([\d.]+)', name)
            return float(m.group(1)) if m else 0.0

        Signal_Histos = []
        for name in sorted(SignalNameList, key=_mass_key):
            h = SignalObjects[name].HistogramName
            if h is not None:
                Signal_Histos.append((name, h.Clone()))

        # --- Print yields ---
        if QEDCompton_Histo is not None:
            print("QED Compton events (incl overflow):",
                  QEDCompton_Histo.Integral(0, QEDCompton_Histo.GetNbinsX() + 1))
        for name, h in Signal_Histos:
            print("Signal {:30s}: {:.1f}".format(name, h.Integral(0, h.GetNbinsX() + 1)))

        # --- Normalize to unit area ---
        if args.normalize:
            if QEDCompton_Histo is not None and QEDCompton_Histo.Integral() > 0:
                QEDCompton_Histo.Scale(1.0 / QEDCompton_Histo.Integral())
            for name, h in Signal_Histos:
                if h.Integral() > 0:
                    h.Scale(1.0 / h.Integral())

        # --- Skip if nothing to draw ---
        if QEDCompton_Histo is None and not Signal_Histos:
            print("WARNING: nothing to draw for", variable)
            del BkgObjects
            del SignalObjects
            continue

        # --- Style: background ---
        color_qed = "#4c9be8"
        if QEDCompton_Histo is not None:
            QEDCompton_Histo.SetFillColor(ROOT.TColor.GetColor(color_qed))
            QEDCompton_Histo.SetLineWidth(0)

        # --- Style: signals ---
        for i, (name, h) in enumerate(Signal_Histos):
            h.SetLineColor(SIGNAL_COLORS[i % len(SIGNAL_COLORS)])
            h.SetLineWidth(2)
            h.SetFillStyle(0)

        # --- Build background stack ---
        backgroundStack = ROOT.THStack('backgroundStack', 'backgroundstack')
        if QEDCompton_Histo is not None:
            backgroundStack.Add(QEDCompton_Histo, 'HIST')

        if backgroundStack.GetNhists() == 0:
            print("WARNING: empty background stack for", variable)
            del BkgObjects
            del SignalObjects
            continue

        backgroundStack_Errors = MakeStackErrors(backgroundStack)

        # --- Canvas ---
        theCanvas = ROOT.TCanvas("theCanvas", "theCanvas")
        theCanvas.Divide()
        plotPad = ROOT.gPad.GetPrimitive('theCanvas_1')
        plotPad.SetPad("pad1", "plot", 0, 0, 1, 1)
        plotPad.SetFillColor(0)
        plotPad.SetBorderMode(0)
        plotPad.SetBorderSize(1)
        plotPad.SetTickx(1)
        plotPad.SetTicky(1)
        plotPad.SetTopMargin(0.10)
        plotPad.SetBottomMargin(0.12)
        plotPad.SetLeftMargin(0.15)
        plotPad.SetRightMargin(0.20)
        plotPad.SetFrameFillStyle(0)
        plotPad.SetFrameLineStyle(0)
        plotPad.SetFrameLineWidth(1)
        plotPad.SetFrameBorderMode(0)
        plotPad.SetFrameBorderSize(1)
        if args.logScale:
            plotPad.SetLogy(1)
        plotPad.cd()
        plotPad.SetFrameLineWidth(1)
        plotPad.SetTickx()
        plotPad.SetTicky()

        sig_max = max((h.GetMaximum() for _, h in Signal_Histos), default=0.0)
        maxi = max(backgroundStack.GetMaximum(), sig_max)
        backgroundStack.SetMaximum(maxi + 0.3 * maxi)
        backgroundStack.SetMinimum(0.1)
        backgroundStack.Draw()
        backgroundStack_Errors.Draw('SAME e2')
        backgroundStack.SetTitle("")

        for _, h in Signal_Histos:
            h.Draw('SAME HIST')

        backgroundStack.GetYaxis().SetTitle("A.U." if args.normalize else "Events")
        backgroundStack.GetYaxis().SetTitleSize(0.04)
        backgroundStack.GetYaxis().SetLabelSize(0.03)
        backgroundStack.GetYaxis().SetTitleOffset(1.05)
        backgroundStack.GetXaxis().SetLabelSize(0.03)
        backgroundStack.GetXaxis().SetTitleSize(0.04)
        backgroundStack.GetXaxis().SetTitle(variableAxisTitleDictionary[variable])

        # --- Legend ---
        theLegend = ROOT.TLegend(0.80, 0.35, 1.0, 0.90, "", "brNDC")
        theLegend.SetTextSize(0.025)
        theLegend.SetLineWidth(0)
        theLegend.SetFillStyle(1001)
        theLegend.SetFillColor(0)
        theLegend.SetBorderSize(0)
        theLegend.SetTextFont(42)
        if QEDCompton_Histo is not None:
            theLegend.AddEntry(QEDCompton_Histo, 'QED Compton', 'f')
        for name, h in Signal_Histos:
            theLegend.AddEntry(h, alp_mass_label(name), 'l')
        theLegend.Draw('SAME')

        # --- Experiment and lumi labels ---
        cmsLatex = ROOT.TLatex()
        cmsLatex.SetNDC(True)
        cmsLatex.SetTextFont(61)
        cmsLatex.SetTextAlign(11)
        cmsLatex.SetTextSize(0.05)
        cmsLatex.DrawLatex(0.15, 0.91, "ePIC")
        cmsLatex.SetTextFont(52)
        cmsLatex.DrawLatex(0.15 + 0.07, 0.91, "Simulation")
        cmsLatex.SetTextFont(42)
        cmsLatex.SetTextAlign(31)
        cmsLatex.SetTextSize(0.035)
        cmsLatex.DrawLatex(0.80, 0.91, args.lumi)

        # --- Save ---
        output_dir = 'EICPlots/' + args.Sub + '/'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        fname_key = variableFileNameDictionary.get(variable,
                    variable.replace('[', '_').replace(']', '').replace('/', '_'))
        theCanvas.SaveAs(output_dir + fname_key + '.png')
        theCanvas.SaveAs(output_dir + fname_key + '.pdf')

        if args.pause:
            input("Press Enter to Continue...")

        del theCanvas
        del BkgObjects
        del SignalObjects


if __name__ == '__main__':
    main()
