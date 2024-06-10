#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <regex>
#include <tuple>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;
namespace fs = std::filesystem;


void sum_vector(vector<double>& v_destination, vector<double>& v_toSum)
{
    for (int i = 0; i < v_destination.size(); i++)
    {
        v_destination[i] += v_toSum[i];
    }
}

void divide_vector(vector<double>& v_destination, double divisor)
{
    for (int i = 0; i < v_destination.size(); i++)
    {
        v_destination[i] /= divisor;
    }
}

void baseline_correction(vector<double>& v, int n_start, int n_end)
{
    double sum = 0;
    for (int i = n_start; i < n_end; i++)
    {
        sum += v[i];
    }
    sum /= (n_end - n_start);
    for (int i = 0; i < v.size(); i++)
    {
        v[i] -= sum;
    }
}

void scale_vector(vector<double>& v, double scale)
{
    for (int i = 0; i < v.size(); i++)
    {
        v[i] *= scale;
    }
}

void detrending(vector<double>& v, int n_start, int n_end)
{
    double sum_yi = 0;
    double sum_xi = 0;
    double sum_yixi = 0;
    double sum_xi2 = 0;

    for (int i = n_start; i < n_end; i++)
    {
        sum_yi += v[i];
        sum_xi += i;
        sum_yixi += v[i] * i;
        sum_xi2 += i * i;
    }
    double avg_yi = sum_yi / (n_end - n_start);
    double avg_xi = sum_xi / (n_end - n_start);
    double avg_yixi = sum_yixi / (n_end - n_start);
    double avg_xi2 = sum_xi2 / (n_end - n_start);

    double cov_xy = avg_yixi - avg_yi * avg_xi;
    double var_x = avg_xi2 - avg_xi * avg_xi;

    double a = cov_xy / var_x;
    double b = avg_yi - a * avg_xi;

    for (int i = 0; i < v.size(); i++)
    {
        v[i] -= a * i + b;
    }
}



void ExtractTemplate(TTree* waves, vector<double> & templ)
{
    int recordLength;
    vector<double> waveform;
    vector<double> waveform_avg;



    waves -> SetBranchAddress("recordLength", &recordLength);
    waves -> GetEntry(0);
    waveform.resize(recordLength);
    waveform_avg.resize(recordLength);
    waves -> SetBranchAddress("waveform", &waveform[0]);

    for (int i = 0; i < waves -> GetEntries(); i++)
    {
        waves -> GetEntry(i);
        sum_vector(waveform_avg, waveform);
    }

    divide_vector(waveform_avg, waves -> GetEntries());
    baseline_correction(waveform_avg, 0, 1500);
    detrending(waveform_avg, 0, 1500);

    double max = waveform_avg[0];
    double min = waveform_avg[0];
    for (int i = 0; i < waveform_avg.size(); i++)
    {
        if (waveform_avg[i] > max)
        {
            max = waveform_avg[i];
        }
        if (waveform_avg[i] < min)
        {
            min = waveform_avg[i];
        }
    }

    // Determine the polarity of the waveform
    // If the maximum is the most distant from 0, the polarity is positive
    // If the minimum is the most distant from 0, the polarity is negative

    if (max > -min)
    {
        scale_vector(waveform_avg, 1/ max);
    }
    else
    {
        scale_vector(waveform_avg, -1/ min);
    }

    templ = waveform_avg;

}



void template_fitting(vector<double>& templ, vector<double>& wave, vector<double>& fit, int n_start, int n_end)
{
    // Fit is linear
    // The wave is Y[k]
    // The template is X[k]
    // The fit is Y[k] = aX[k] + b (between the samples n_start and n_end)

    double sum_yi = 0;
    double sum_xi = 0;
    double sum_yixi = 0;
    double sum_xi2 = 0;

    for (int i = n_start; i < n_end; i++)
    {
        sum_yi += wave[i];
        sum_xi += templ[i];
        sum_yixi += wave[i] * templ[i];
        sum_xi2 += templ[i] * templ[i];
    }

    double avg_yi = sum_yi / (n_end - n_start);
    double avg_xi = sum_xi / (n_end - n_start);
    double avg_yixi = sum_yixi / (n_end - n_start);
    double avg_xi2 = sum_xi2 / (n_end - n_start);

    double cov_xy = avg_yixi - avg_yi * avg_xi;
    double var_x = avg_xi2 - avg_xi * avg_xi;

    double a = cov_xy / var_x;
    double b = avg_yi - a * avg_xi;

    fit[0] = a;
    fit[1] = b;
}


double ADC_to_E_CH0(double ADC)
{
    double E0 = 0;
    double ADC0 = 0.;

    double E1 = 59.54e-3;
    double ADC1 = 472.0;
    
    double m = (E1 - E0) / (ADC1 - ADC0);
    double q = E0 - m * ADC0;
    
    return m * ADC + q;
}

double ADC_to_E_CH1(double ADC)
{
    double E0 = 0;
    double ADC0 = 0.;

    double E1 = 59.54e-3;
    double ADC1 = 700.2;
    
    double m = (E1 - E0) / (ADC1 - ADC0);
    double q = E0 - m * ADC0;
    
    return m * ADC + q;
}




int Analysis_DeltaE_E
()
{
    vector<string> filenames;
    filenames.push_back("/home/riccardo/Documenti/NUSES/DeltaE_E/Sr90/wave0.root");
    filenames.push_back("/home/riccardo/Documenti/NUSES/DeltaE_E/Sr90/wave1.root");

    vector<TFile*> files;
    vector<TTree*> waves;

    for (int i = 0; i < filenames.size(); i++)
    {
        TFile* file = new TFile(filenames[i].c_str(), "READ");
        TTree* wave = (TTree*)file -> Get("waves");
        files.push_back(file);
        waves.push_back(wave);
    }

    vector<double> templ0;
    vector<double> templ1;

    ExtractTemplate(waves[0], templ0);
    ExtractTemplate(waves[1], templ1);


    // Open a txt file to save the templates
    ofstream file0("template0.txt");
    ofstream file1("template1.txt");

    for (int i = 0; i < templ0.size(); i++)
    {
        file0 << templ0[i] << endl;
        file1 << templ1[i] << endl;
    }

    file0.close();
    file1.close();

    vector<double> templ0_2 = templ0;
    vector<double> templ1_2 = templ1;

    scale_vector(templ0_2, -1);


    TCanvas* c = new TCanvas("c", "c", 800, 800);
    TGraph* g0 = new TGraph(templ0_2.size());
    TGraph* g1 = new TGraph(templ1.size());

    for (int i = 0; i < templ0.size(); i++)
    {
        g0 -> SetPoint(i, i, templ0_2[i]);
        g1 -> SetPoint(i, i, templ1[i]);
    }

    TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg -> AddEntry(g0, "Template 0", "l");
    leg -> AddEntry(g1, "Template 1", "l");
    

    g0 -> SetLineColor(kRed);
    g1 -> SetLineColor(kBlue);

    g0 -> SetLineWidth(2.5);
    g1 -> SetLineWidth(2.5);

    g0 -> Draw("AL");
    g1 -> Draw("L SAME");
    leg -> Draw();


    int recordLength_ch0;
    int recordLength_ch1;

    int n_fit_start_ch0 = 1900;
    int n_fit_end_ch0 = 2870;
    int n_fit_start_ch1 = 1915;
    int n_fit_end_ch1 = 2130;


    vector<double> wave0;
    vector<double> wave1;
    vector<double> fit0(2);
    vector<double> fit1(2);

    waves[0] -> SetBranchAddress("recordLength", &recordLength_ch0);
    waves[1] -> SetBranchAddress("recordLength", &recordLength_ch1);

    waves[0] -> GetEntry(0);
    wave0.resize(recordLength_ch0);
    wave1.resize(recordLength_ch0);

    waves[0] -> SetBranchAddress("waveform", &wave0[0]);
    waves[1] -> SetBranchAddress("waveform", &wave1[0]);

    TFile *f = new TFile("Energies.root", "RECREATE");
    TTree *t = new TTree("Energies", "Energies");


    double E_ch0;
    double E_ch1;
    t -> Branch("E_thin", &E_ch0);
    t -> Branch("E_thick", &E_ch1);

    for (int i = 0; i < waves[0] -> GetEntries(); i++)
    {
        cout << "Processing event " << i << endl;
    
        waves[0] -> GetEntry(i);
        waves[1] -> GetEntry(i);


        baseline_correction(wave0, 0, 1500);
        baseline_correction(wave1, 0, 1500);

        detrending(wave0, 0, 1500);
        detrending(wave1, 0, 1500);

        template_fitting(templ0, wave0, fit0, n_fit_start_ch0, n_fit_end_ch0);
        template_fitting(templ1, wave1, fit1, n_fit_start_ch1, n_fit_end_ch1);

        E_ch0 = ADC_to_E_CH0(fit0[0]);
        E_ch1 = ADC_to_E_CH1(fit1[0]);

        cout << "E_ch0 = " << E_ch0 << endl;
        cout << "E_ch1 = " << E_ch1 << endl;
        cout << endl;

        t -> Fill();
    }

    t -> Write();

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
    t -> Draw("E_thin:E_thick", "", "colz");

    double xMin = TMath::MinElement(t -> GetEntries(), t -> GetV1());
    double xMax = TMath::MaxElement(t -> GetEntries(), t -> GetV1());
    double yMin = TMath::MinElement(t -> GetEntries(), t -> GetV2());
    double yMax = TMath::MaxElement(t -> GetEntries(), t -> GetV2());

    TH2D* h = new TH2D("h", "h", 200, xMin, xMax, 200, yMin, yMax);
    t -> Draw("E_thin:E_thick >> h", "", "colz");
    h -> GetXaxis() -> SetTitle("E_{thick} [MeV]");
    h -> GetYaxis() -> SetTitle("E_{thin} [MeV]");


    TCanvas* c3 = new TCanvas("c3", "c3", 800, 800);
    t -> Draw("log10(E_thin*(E_thin+E_thick)):(E_thin+E_thick)", "", "colz");

    double xMin2 = TMath::MinElement(t -> GetEntries(), t -> GetV2());
    double xMax2 = TMath::MaxElement(t -> GetEntries(), t -> GetV2());
    double yMin2 = TMath::MinElement(t -> GetEntries(), t -> GetV1());
    double yMax2 = TMath::MaxElement(t -> GetEntries(), t -> GetV1());

    TH2D* h2 = new TH2D("h2", "h2", 200, xMin2, xMax2, 200, yMin2, yMax2);
    t -> Draw("log10(E_thin*(E_thin+E_thick)):(E_thin+E_thick) >> h2", "", "colz");
    h2 -> GetXaxis() -> SetTitle("E_{total} [MeV]");
    h2 -> GetYaxis() -> SetTitle("log_{10}(E_{thin} * E_{total})");







    return 0;
}


