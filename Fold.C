#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TApplication.h"
#include <iostream>

/*
This script can be read using 
root -x 'Fold.C("input_file.root", 2.0)'
... or load the script in ROOT and run the function with the parameters

For this example...
"input_file.root" = input ROOT file
2.0 = FWHM of the folding gaussian in keV

only the FWHM is significant in this script for caluclation, the other variables are for naming purposes

A canvas will appear but the top plots will not be there. You must comment out "Clean Up" lines and the "Exit ROOT" line to see the plots
*/

void Fold(const char* inputFileName, double fwhm);

void Fold(const char* inputFileName, double fwhm) {

    // Open the input ROOT file
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
        return;
    }

    // Get the histogram named "Esil" from the original file
    TH1F* srcHist = static_cast<TH1F*>(inputFile->Get("Esil"));
    if (!srcHist) {
        std::cerr << "Error: Could not find histogram Esil" << std::endl;
        return;
    }

    // Function for FWHM in MeV to σ conversion
    auto fwhmToSigma = [](double fwhm) {
        return (fwhm / 1000.0) / 2.3548;
    };

    // Define Gaussian smearing parameters
    // double fwhm = 2.0; // FWHM in keV
    TRandom3 rand;


    // Create a new histogram for the convoluted data
    TH1F* outHist = static_cast<TH1F*>(srcHist->Clone("folded"));
    outHist->SetTitle(Form("Folded with Gaussian (FWHM=%.2f keV)", fwhm));
    outHist->GetXaxis()->SetTitle("Energy [MeV]");
    outHist->GetYaxis()->SetTitle("Counts/keV");
    outHist->Reset(); // Clear the histogram
    outHist->SetOption("HIST"); // Set options to draw as bar charts

    srcHist->SetTitle("Energy Deposition in PIPS");
    srcHist->GetXaxis()->SetTitle("Energy [MeV]");
    srcHist->GetYaxis()->SetTitle("Counts/keV");

    // Loop over all bins in the original histogram
    int maxFilledBin = srcHist->FindLastBinAbove();
    for (int i = 1; i <= maxFilledBin; i++) {
        double counts = srcHist->GetBinContent(i);
        double binCenter = srcHist->GetBinCenter(i);

        // Calculate the standard deviation σ at this bin energy
        double sigma = fwhmToSigma(fwhm) * binCenter;

        // Smear the content of each bin
        for (int j = 0; j < counts; j++) {
            double smearedValue = rand.Gaus(binCenter, sigma);
            outHist->Fill(smearedValue);
        }
    }

    // Calculate the transmission probability: sum counts above the first bin
    double totalCounts = 0.0;
    double totalInteractionCounts = 0.0;
    for (int i = 1; i <= srcHist->GetNbinsX(); i++) {
        double binContent = srcHist->GetBinContent(i);
        totalCounts += binContent;
        if (i > 1) { // Exclude the first bin
            totalInteractionCounts += binContent;
        }
    }

    // Create a canvas with six pads (3x2 layout)
    TCanvas* c1 = new TCanvas("c1", "Plots", 1200, 900);
    c1->Divide(2, 1);  // Divide the canvas into 2 columns and 3 rows

    // Draw the original histogram in the first pad
    c1->cd(1);  // Select the first pad
    gPad->SetLogy();  // Set the y-axis to a logarithmic scale
    srcHist->SetLineColor(kBlue);
    srcHist->Draw();

    // Draw the convoluted histogram in the second pad
    c1->cd(2);  // Select the second pad
    gPad->SetLogy();  // Set the y-axis to a logarithmic scale
    outHist->SetLineColor(kRed);
    outHist->Draw("HIST");

    // Update and save the canvas
    c1->Update();
    c1->Modified();
}
