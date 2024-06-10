// C/C++ script for the digitizer simulation
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>
#include <regex>
#include <tuple>


#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLine.h"


#include "DigitizerCAEN.C"

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


double template_fitting_2(vector<double>& templ, vector<double>& wave, int n_start, int n_end, int *offset_out)
{
    // Fit is linear
    // The wave is Y[k]
    // The template is X[k]
    // The fit is Y[k] = aX[k] + b (between the samples n_start and n_end)

    int offset = 1000;

    double max_outfilter = -999999;


    double sum_xi = 0;
    double sum_xi2 = 0;
        
    for(int i = n_start; i < n_end; i++)
    {
        sum_xi += templ[i];
        sum_xi2 += templ[i] * templ[i];
    }

    double avg_xi = sum_xi / (n_end - n_start);
    double avg_xi2 = sum_xi2 / (n_end - n_start);
    double var_x = avg_xi2 - avg_xi * avg_xi;


    for(int j = -offset; j < offset; j++)
    {



        double sum_yi = 0;
        double sum_yixi = 0;
    

        for (int i = n_start; i < n_end; i++)
        {
            sum_yi += wave[i-j];
            sum_yixi += wave[i-j] * templ[i];
        }

        double avg_yi = sum_yi / (n_end - n_start);
        double avg_yixi = sum_yixi / (n_end - n_start);
        double cov_xy = avg_yixi - avg_yi * avg_xi;

        double a = cov_xy / var_x;
        double b = avg_yi - a * avg_xi;
    
        if(a > max_outfilter)
        {
            max_outfilter = a;
            *offset_out = j;
        }
    }
    cout << max_outfilter << endl;
    return max_outfilter;
}




int TemplateFit()
{
    string fname_template = "/home/riccardo/Documenti/NUSES/DeltaE_E/Template_Gain1/template0.txt";
    string fname_toConvert = "/home/riccardo/Documenti/NUSES/DeltaE_E/Am/wave0.root";
    string fname_fit = fname_toConvert.substr(0, fname_toConvert.size() - 5) + "_fit.root";

    int n_baseline = 4000;

    int n_start = 4500;
    int n_end = 6000;

    double min_t = 0;
    double max_t = 0;

    vector<double> templ;
    vector<double> fit(2);


    ifstream file_template(fname_template);
    double value;
    while (file_template >> value)
    {
        templ.push_back(value);
    }
    file_template.close();



    for (int i = 0; i < templ.size(); i++)
    {
        if(templ[i] > max_t)
        {
            max_t = templ[i];
        }
        if(templ[i] < min_t)
        {
            min_t = templ[i];
        }
    }



    TFile* file = new TFile(fname_toConvert.c_str(), "READ");
    TTree* waves = (TTree*) file -> Get("waves");

    int recordLength;
    waves -> SetBranchAddress("recordLength", &recordLength);
    waves -> GetEntry(0);
    vector<double> wave(recordLength);
    waves -> SetBranchAddress("waveform", &wave[0]);

    
    double E = 0;
    TFile *file_fit = new TFile(fname_fit.c_str(), "RECREATE");
    TTree *tree_fit = new TTree("Energies", "Energies");
    tree_fit -> Branch("E", &E, "E/D");
    int offset = 0;
    tree_fit -> Branch("offset", &offset, "offset/I");



    int N = -1;
    if (N == -1)
    {
        N = waves -> GetEntries() -1;
    }

    for (int i = 0; i < N; i++)
    {
        cout << "Processing entry " << i << " of " << waves -> GetEntries() << endl;
        waves -> GetEntry(i);

        vector<double> wave_copy = wave;

        baseline_correction(wave_copy, 0, n_baseline);

        detrending(wave_copy, 0, n_baseline);

        vector<double> fit_values(2);

        //template_fitting(templ, wave_copy, fit_values, n_start, n_end);

        E =  template_fitting_2(templ, wave_copy, n_start, n_end, &offset);

        tree_fit -> Fill();
    }

    tree_fit -> Write();
    file_fit -> Close();


    // Plot Template and the samples used for the fit
    TCanvas* c = new TCanvas("c", "c", 1200, 600);

    TGraph* gr_temp = new TGraph();

    for (int i = 0; i < templ.size(); i++)
    {
        gr_temp -> SetPoint(i, i, templ[i]);
    }

    gr_temp -> SetLineColor(kBlue);
    gr_temp -> SetLineWidth(2);
    gr_temp -> Draw("AL");

    // Draw Vertical line at n_start and n_end and shaded area
    TLine* l_start = new TLine(n_start, min_t, n_start, max_t);
    TLine* l_end = new TLine(n_end, min_t, n_end, max_t);
    l_start -> SetLineColor(kRed);
    l_end -> SetLineColor(kRed);
    l_start -> SetLineWidth(2);
    l_end -> SetLineWidth(2);
    l_start -> Draw();
    l_end -> Draw();

    // Line at n_baseline
    TLine* l_baseline = new TLine(n_baseline, min_t, n_baseline, max_t);
    l_baseline -> SetLineColor(kGreen);
    l_baseline -> SetLineWidth(2);
    l_baseline -> Draw();

    string figname = fname_toConvert.substr(0, fname_toConvert.size() - 5) + "_fit.png";

    c -> SaveAs(figname.c_str());


    cout << n_baseline << endl;
    cout << n_start << endl;
    cout << n_end << endl;


    return 0;


}