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
root -x 'Plots.C("1000.0", "2.0", "70.0", "35.0", "1000")'
... or load the script in ROOT and run the function with the parameters

For this example...
1000.0 = electron energy in keV
2.0 = FWHM of the folding gaussian in keV
70.0 = f in mm
35.0 = g in mm
1000 = thickness in µm

only the FWHM is significant in this script for caluclation, the other variables are for naming purposes

A canvas will appear but the top plots will not be there. You must comment out "Clean Up" lines and the "Exit ROOT" line to see the plots
*/

void Plots(double energy, double fwhm, double f, double g, int thickness);

void Plots(double energy, double fwhm, double f, double g, int thickness) {

    // std::string custom_identifier = "";

    // Create a unique name for the convoluted histogram
    // std::string name = custom_identifier + Form("PIPS%i_f%.0fmm_g%.0fmm_e%.0fkeV", thickness, f, g, energy);
    std::string name = Form("PIPS%i_f%.0fmm_g%.0fmm_e%.0fkeV", thickness, f, g, energy);

    // Open the original ROOT file to read the histogram
    TFile* inputFile = TFile::Open("./ICESPICE.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file build/ICESPICE.root" << std::endl;
        return;
    }

    // Get the histogram named "Esil" from the original file
    TH1F* srcHist = static_cast<TH1F*>(inputFile->Get("Esil"));
    if (!srcHist) {
        std::cerr << "Error: Could not find histogram Esil" << std::endl;
        return;
    }

    // Open the new ROOT file with recreate option to save the output
    // output the file to the name of the histogram + ICEPSPICE
    std::string outputFileName = "ICESPICE_" + name + ".root";

    TFile* outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not open file build/ICESPICE_new.root" << std::endl;
        return;
    }

    std::cout << "Output file: " << outputFileName << std::endl;

    // Print the parameters
    std::cout << "Parameters:" << std::endl;
    std::cout << "\tElectron Energy: " << energy << " keV" << std::endl;
    std::cout << "\tFolded with FWHM: " << fwhm << " keV" << std::endl;
    std::cout << "\tf: " << f << " mm"<< std::endl;
    std::cout << "\tg: " << g << " mm"<< std::endl;
    std::cout << "\tDetector Thickness: " << thickness << " µm\n" << std::endl;


    // rename the original histogram
    srcHist->SetName(name.c_str());

    // Function for FWHM in MeV to σ conversion
    auto fwhmToSigma = [](double fwhm) {
        return (fwhm / 1000.0) / 2.3548;
    };

    // Define Gaussian smearing parameters
    // double fwhm = 2.0; // FWHM in keV
    TRandom3 rand;

    // Create a unique name for the convoluted histogram
    std::string histName = name + Form("_folded%.0fkeV", fwhm);

    // Create a new histogram for the convoluted data
    TH1F* outHist = static_cast<TH1F*>(srcHist->Clone(histName.c_str()));
    outHist->SetTitle(Form("Folded with Gaussian (FWHM=%.2f keV)", fwhm));
    outHist->GetXaxis()->SetTitle("Energy [MeV]");
    outHist->GetYaxis()->SetTitle("Counts/keV");
    outHist->Reset(); // Clear the histogram
    outHist->SetOption("HIST"); // Set options to draw as bar charts

    srcHist->SetTitle(Form("Energy Deposition in PIPS%i", thickness));
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
    double fullEnergyCounts = srcHist->GetBinContent(energy);
    double transmissionProbability = totalInteractionCounts / totalCounts * 100.0 / 2.0;
    double transmissionProbabilityFullEnergy = fullEnergyCounts / totalCounts * 100.0 / 2.0;
    double detectorEfficiency = fullEnergyCounts / totalInteractionCounts * 100.0;

    // Print out the statistics
    std::cout << "\tTotal Counts (2π): " << totalCounts << std::endl;
    std::cout << "\tCounts in Detector: " << totalInteractionCounts << std::endl;
    std::cout << "\tCounts at Full Energy: " << fullEnergyCounts << std::endl;
    std::cout << "\tTransmission Probability in 4π [%]: " << transmissionProbability << "%" << std::endl;
    std::cout << "\tFull Energy Deposited Transmission Probability in 4π [%]: " << transmissionProbabilityFullEnergy << "%" << std::endl;
    std::cout << "\tDetector Full Energy Deposited Efficiency [%]: " << detectorEfficiency << "%" << std::endl;

    // Create a canvas with six pads (3x2 layout)
    TCanvas* c1 = new TCanvas("c1", "Plots", 1200, 900);
    c1->Divide(2, 3);  // Divide the canvas into 2 columns and 3 rows

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

    // Create and draw the transmission probability on the third pad
    c1->cd(3);  // Select the third pad
    TGraph* transProbGraph = new TGraph();
    transProbGraph->SetPoint(0, energy / 1000.0, transmissionProbability); // Convert keV to MeV for x-axis
    transProbGraph->SetMarkerStyle(21);
    transProbGraph->SetMarkerColor(kBlue);
    transProbGraph->GetXaxis()->SetLimits(0, 2); // Set x-axis range from 0 to 2 MeV
    transProbGraph->GetYaxis()->SetRangeUser(0, 10); // Set y-axis range to include 0
    transProbGraph->GetXaxis()->SetTitle("Energy [MeV]");
    transProbGraph->GetYaxis()->SetTitle("Probability [%]");
    transProbGraph->SetTitle("4pi- Transmission Probability");
    transProbGraph->SetDrawOption("AP");  // Set draw option to "P" for points on
    transProbGraph->Draw("AP");

    // Create and draw the full energy transmission probability on the fourth pad
    c1->cd(4);  // Select the fourth pad
    TGraph* fullEnergyProbGraph = new TGraph();
    fullEnergyProbGraph->SetPoint(0, energy / 1000.0, transmissionProbabilityFullEnergy); // Convert keV to MeV for x-axis
    fullEnergyProbGraph->SetMarkerStyle(22);
    fullEnergyProbGraph->SetMarkerColor(kRed);
    fullEnergyProbGraph->GetXaxis()->SetLimits(0, 2); // Set x-axis range from 0 to 2 MeV
    fullEnergyProbGraph->GetYaxis()->SetRangeUser(0, 10); // Set y-axis range to include 0
    fullEnergyProbGraph->GetXaxis()->SetTitle("Energy [MeV]");
    fullEnergyProbGraph->GetYaxis()->SetTitle("Probability [%]");
    fullEnergyProbGraph->SetTitle("4pi- Full Energy Transmission Probability");
    fullEnergyProbGraph->SetDrawOption("AP");  // Set draw option to "P" for points on
    fullEnergyProbGraph->Draw("AP");

    // Create and draw the detector efficiency in the fifth pad
    c1->cd(5);  // Select the fifth pad
    TGraph* efficiencyGraph = new TGraph();
    efficiencyGraph->SetPoint(0, energy / 1000.0, detectorEfficiency); // Convert keV to MeV for x-axis
    efficiencyGraph->SetMarkerStyle(21);
    efficiencyGraph->SetTitle("Detector Efficiency");
    efficiencyGraph->GetXaxis()->SetLimits(0, 2); // Set x-axis range from 0 to 2 MeV
    efficiencyGraph->GetYaxis()->SetRangeUser(0, 100); // Set y-axis range to include 0
    efficiencyGraph->GetXaxis()->SetTitle("Energy [MeV]");
    efficiencyGraph->GetYaxis()->SetTitle("Detector Efficiency [%]");
    efficiencyGraph->SetDrawOption("AP");  // Set draw option to "P" for points on
    efficiencyGraph->Draw("AP");

    // Update and save the canvas
    c1->Update();
    c1->Modified();
    
    // Save the canvas to a file
    // c1->SaveAs("Plots.png");

    // Write the convoluted histogram and graphs into the same ROOT file
    outputFile->cd();           // Ensure we are in the right directory in the file
    outHist->Write();     // Write the convoluted histogram to the file
    // write the original histogram to the file
    srcHist->Write(name.c_str()); // Write the original histogram to the file
    std::string transProbName = "TransmissionProbability_" + name; 
    std::string fullEnergyProbName = "FEDTransmissionProbability_" + name; 
    std::string efficiencyName = "DetectorEfficiency_" + name; 
    transProbGraph->Write(transProbName.c_str()); // Write the transmission probability graph with the name of the histogram
    fullEnergyProbGraph->Write(fullEnergyProbName.c_str()); // Write the full energy transmission probability graph with the name of the histogram
    efficiencyGraph->Write(efficiencyName.c_str()); // Write the detector efficiency graph with the name of the histogram

    // These will be used to when hadding the files
    transProbGraph->Write(Form("A_Total_TransmissionProbability_PIPS%i_f%.0fmm_g%.0fmm", thickness, f, g)); // Write the transmission probability graph
    fullEnergyProbGraph->Write(Form("A_Total_FEDTransmissionProbability_PIPS%i_f%.0fmm_g%.0fmm", thickness, f, g)); // Write the full energy transmission probability graph
    efficiencyGraph->Write(Form("A_Total_DetectorEfficiency_PIPS%i_f%.0fmm_g%.0fmm", thickness, f, g)); // Write the detector efficiency graph

    // Clean up
    outputFile->Close();
    inputFile->Close();
    delete outputFile;
    delete inputFile;

    // Exit ROOT
    gApplication->Terminate(); // Alternatively, you can use gROOT->ProcessLine(".q"); if running from the command line
}
