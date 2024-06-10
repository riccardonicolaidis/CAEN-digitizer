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
#include "TH1D.h"
#include "TF1.h"

using namespace std;

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

int CH0_AM()
{
    string fname_template_signal = "/home/riccardo/Documenti/NUSES/DeltaE_E/template0.txt";
    string fname_template_noise = "/home/riccardo/Documenti/NUSES/DeltaE_E/template1.txt";
    string file_signal = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH0/wave0.root";
    string file_noise = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH0/wave1.root";
    string file_output = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH0/Am_CH0.root";
    string file_results = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH0/Am_CH0_result.txt";


    int Am_source_on_ch = 0;

    if(Am_source_on_ch == 0)
    {
        fname_template_signal   = "/home/riccardo/Documenti/NUSES/DeltaE_E/template0.txt";
        fname_template_noise    = "/home/riccardo/Documenti/NUSES/DeltaE_E/template1.txt";
        file_signal             = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH0/wave0.root";
        file_noise              = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH0/wave1.root";
        file_output             = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH0/Am_CH0.root";
        file_results            = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH0/Am_CH0_result.txt";
    }

    if(Am_source_on_ch == 1)
    {
        fname_template_signal   = "/home/riccardo/Documenti/NUSES/DeltaE_E/template1.txt";
        fname_template_noise    = "/home/riccardo/Documenti/NUSES/DeltaE_E/template0.txt";
        file_signal             = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH1/wave1.root";
        file_noise              = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH1/wave0.root";
        file_output             = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH1/Am_CH1.root";
        file_results            = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am_CH1/Am_CH1_result.txt";

    }




    int n_start_signal = 1900;
    int n_end_signal = 2870;

    int n_start_noise = 1915;
    int n_end_noise = 2130;

    int n_baseline = 1200;
    
    ifstream file_template_signal(fname_template_signal);
    ifstream file_template_noise(fname_template_noise);

    vector<double> v_template_signal;
    vector<double> v_template_noise;

    double value;
    while (file_template_signal >> value)
    {
        v_template_signal.push_back(value);
    }

    while (file_template_noise >> value)
    {
        v_template_noise.push_back(value);
    }

    file_template_signal.close();
    file_template_noise.close();

    TGraph* gr_template_signal = new TGraph();
    TGraph* gr_template_noise = new TGraph();

    for (int i = 0; i < v_template_signal.size(); i++)
    {
        gr_template_signal -> SetPoint(i, i, v_template_signal[i]);
        gr_template_noise -> SetPoint(i, i, v_template_noise[i]);
    }

    TFile* f_signal = new TFile(file_signal.c_str(), "READ");
    TFile* f_noise = new TFile(file_noise.c_str(), "READ");

    TTree* waves_signal = (TTree*)f_signal -> Get("waves");
    TTree* waves_noise = (TTree*)f_noise -> Get("waves");

    waves_signal -> Print();
    waves_noise -> Print();


    // Now I want to create a TTree where to save the results of the analysis

    int recordLength_signal;
    int recordLength_noise;

    waves_signal -> SetBranchAddress("recordLength", &recordLength_signal);
    waves_noise -> SetBranchAddress("recordLength", &recordLength_noise);

    waves_signal -> GetEntry(0);
    waves_noise -> GetEntry(0);

    cout << "Record length signal: " << recordLength_signal << endl;
    cout << "Record length noise: " << recordLength_noise << endl;

    vector<double> v_wave_signal(recordLength_signal);
    vector<double> v_wave_noise(recordLength_noise);

    waves_signal -> SetBranchAddress("waveform", &v_wave_signal[0]);
    waves_noise -> SetBranchAddress("waveform", &v_wave_noise[0]);

    double E_signal = 0;
    double E_noise = 0;

    TFile* f_output = new TFile(file_output.c_str(), "RECREATE");
    TTree* Energies = new TTree("Energies", "Energies");

    Energies -> Branch("E_signal", &E_signal, "E_signal/D");
    Energies -> Branch("E_noise", &E_noise, "E_noise/D");

    for (int i = 0; i < waves_signal -> GetEntries(); i++)
    {
        waves_signal -> GetEntry(i);
        waves_noise -> GetEntry(i);

        vector<double> v_wave_signal_copy = v_wave_signal;
        vector<double> v_wave_noise_copy = v_wave_noise;

        baseline_correction(v_wave_signal_copy, 0, n_baseline);
        baseline_correction(v_wave_noise_copy, 0, n_baseline);

        detrending(v_wave_signal_copy, 0, n_baseline);
        detrending(v_wave_noise_copy, 0, n_baseline);

        vector<double> v_fit_signal(2);
        vector<double> v_fit_noise(2);

        template_fitting(v_template_signal, v_wave_signal_copy, v_fit_signal, n_start_signal, n_end_signal);
        template_fitting(v_template_noise, v_wave_noise_copy, v_fit_noise,  n_start_noise, n_end_noise);

        double a_signal = v_fit_signal[0];
        double b_signal = v_fit_signal[1];

        double a_noise = v_fit_noise[0];
        double b_noise = v_fit_noise[1];

        E_signal = a_signal;
        E_noise = a_noise;

        Energies -> Fill();
    }

    Energies -> Write();
    
    // Signal Histogram

    double E_signal_min = Energies -> GetMinimum("E_signal");
    double E_signal_max = Energies -> GetMaximum("E_signal");

    double E_noise_min = Energies -> GetMinimum("E_noise");
    double E_noise_max = Energies -> GetMaximum("E_noise");
    



    TH1D* h_E_signal = new TH1D("h_E_signal", "Signal Energy", 800, E_signal_min, E_signal_max);
    // Draw preliminary histogram
    // Find percentiles 10 and 90
    Energies -> Draw("E_signal >> h_E_signal", "", "goff");
    double prob[2] = {0.001, 0.995};
    double quantiles[2];
    h_E_signal -> GetQuantiles(2, quantiles, prob);
    E_signal_min = quantiles[0];
    E_signal_max = quantiles[1];
    h_E_signal = new TH1D("h_E_signal", "Signal Energy", 100, E_signal_min, E_signal_max);

    TH1D* h_E_noise = new TH1D("h_E_noise", "Noise Energy", 800, E_noise_min, E_noise_max);
    // Draw preliminary histogram
    // Find percentiles 10 and 90
    Energies -> Draw("E_noise >> h_E_noise", "", "goff");
    h_E_noise -> GetQuantiles(2, quantiles, prob);
    E_noise_min = quantiles[0];
    E_noise_max = quantiles[1];
    h_E_noise = new TH1D("h_E_noise", "Noise Energy", 100, E_noise_min, E_noise_max);


    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 800);
    Energies -> Draw("E_signal >> h_E_signal", "", "");
    h_E_signal -> Draw();
    
    if(Am_source_on_ch == 0)
    {
        h_E_signal -> SetTitle("Signal Energy - Am source on 300 um Si");
    }
    if(Am_source_on_ch == 1)
    {
        h_E_signal -> SetTitle("Signal Energy - Am source on 500 um Si");
    }




    TCanvas* c2 = new TCanvas("c2", "c2", 1200, 800);
    Energies -> Draw("E_noise >> h_E_noise", "", "");
    h_E_noise -> Draw();

    
    h_E_noise -> Fit("f2", "L", "", E_noise_min, E_noise_max);



  





    return 0;
}