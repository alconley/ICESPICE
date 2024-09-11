#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TApplication.h"
#include <iostream>

void Fold(const char* inputFileName, double fwhm, bool saveToFile = false);

void Fold(const char* inputFileName, double fwhm, bool saveToFile) {

    // Open the input ROOT file
    TFile* inputFile = TFile::Open(inputFileName, "UPDATE");  // Open in UPDATE mode to allow writing
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
        return (fwhm / 1000.0) / 2.3548;  // FWHM in keV, convert to sigma
    };

    // Create a new histogram for the convoluted data
    TH1F* outHist = static_cast<TH1F*>(srcHist->Clone("folded"));
    outHist->SetTitle(Form("Folded with Gaussian (FWHM=%.2f keV)", fwhm));
    outHist->GetXaxis()->SetTitle("Energy [MeV]");
    outHist->GetYaxis()->SetTitle("Counts/keV");
    outHist->Reset();  // Clear the histogram content

    // Random number generator for smearing
    TRandom3 rand;

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

    // Create a canvas for plotting
    TCanvas* c1 = new TCanvas("c1", "Plots", 1200, 900);
    c1->Divide(2, 1);  // Divide canvas into 2 columns

    // Draw the original histogram in the first pad
    c1->cd(1);  // Select the first pad
    gPad->SetLogy();  // Set the y-axis to logarithmic scale
    srcHist->SetLineColor(kBlue);
    srcHist->Draw();

    // Draw the folded histogram in the second pad
    c1->cd(2);  // Select the second pad
    gPad->SetLogy();  // Set the y-axis to logarithmic scale
    outHist->SetLineColor(kRed);
    outHist->Draw("HIST");

    // Update and save the canvas
    c1->Update();
    c1->Modified();

    // Save the folded histogram to the same file if saveToFile is true
    if (saveToFile) {
        inputFile->cd();  // Ensure we're in the correct file
        outHist->Write();  // Save the folded histogram to the file
        std::cout << "Folded histogram saved to " << inputFileName << std::endl;

        // Clean up
        delete c1;
        inputFile->Close();
        delete inputFile;
    }


}
