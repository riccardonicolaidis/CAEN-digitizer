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
#include "TMath.h"
#include "TRandom3.h"


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





int MatchedFilter_Es()
{

    string fname_template = "/home/riccardo/Documenti/NUSES/DeltaE_E/Template_Gain1/template0.txt";

    vector<double> templ;


    ifstream file_template(fname_template);
    double value;
    while (file_template >> value)
    {
        value = -value;
        templ.push_back(value);
    }
    file_template.close();



    double max_t = TMath::MaxElement(templ.size(), &templ[0]);
    double min_t = TMath::MinElement(templ.size(), &templ[0]);


    int n_trigger = 0;
    double threshold = 0.1;
    bool trigger = false;

    for(int i = 0; i < templ.size(); i++)
    {
        if(templ[i] > threshold * max_t)
        {
            if(trigger == false)
            {
                n_trigger = i;
                trigger = true;
                break;
            }
        }
    }

    cout << "Trigger at sample: " << n_trigger << endl;

    vector<double> template_clean(templ.size() - n_trigger);

    for(int i = 0; i < template_clean.size(); i++)
    {
        template_clean[i] = templ[i + n_trigger];
    }


    TRandom3* rnd = new TRandom3(0);
    int N = template_clean.size() + 5000;
    
    vector<double> templated_extended(N);
    for(int i = 0; i < template_clean.size(); i++)
    {
        templated_extended[i] = template_clean[i];
    }
    for(int i = template_clean.size() - 1000; i < N; i++)
    {
        templated_extended[i] = templated_extended[i-1]* exp(-1/1000.);
    }
    
    vector<double> noisy_signal(N);
    vector<double> clean_signal(N);

    for(int i = 0; i < N; i++)
    {
        noisy_signal[i] = rnd -> Gaus(0, 1);
    }

    int n_start = 4000;



    for(int i = 0; i < templated_extended.size(); i++)
    {
        noisy_signal[(i + n_start)%N] += templated_extended[i];
        clean_signal[(i + n_start)%N] = templated_extended[i];
    }

    TGraph* gr_noisy_signal = new TGraph(N);
    for (int i = 0; i < N; i++)
    {
        gr_noisy_signal -> SetPoint(i, i, noisy_signal[i]);
    }

    TGraph* gr_clean_signal = new TGraph(N);
    for (int i = 0; i < N; i++)
    {
        gr_clean_signal -> SetPoint(i, i, clean_signal[i]);
    }




    vector<double> matched_filter(N);

    for(int i = 0; i < N; i++)
    {
        double sum = 0;
        for(int j = 0; j < templated_extended.size(); j++)
        {
            sum += noisy_signal[(i + j)%N] * templated_extended[j];
        }
        matched_filter[i] = sum;
    }

    TGraph* gr_matched_filter = new TGraph(N);
    for (int i = 0; i < N; i++)
    {
        gr_matched_filter -> SetPoint(i, i, matched_filter[i]);
    }



    TCanvas* c4 = new TCanvas("c4", "c4", 800, 1000);
    c4 -> Divide(1, 3);
    
    c4 -> cd(1);
    gr_clean_signal -> Draw("AL");
    TLine *l_trigger = new TLine(n_start, min_t, n_start, max_t);
    l_trigger -> SetLineColor(kRed);
    l_trigger -> SetLineWidth(2);
    l_trigger -> Draw();

    c4 -> cd(2);
    gr_noisy_signal -> Draw("AL");
    TLine *l_trigger2 = new TLine(n_start, TMath::MinElement(N, &noisy_signal[0]), n_start, TMath::MaxElement(N, &noisy_signal[0]));
    l_trigger2 -> SetLineColor(kRed);
    l_trigger2 -> SetLineWidth(2);
    l_trigger2 -> Draw();


    c4 -> cd(3);
    gr_matched_filter -> Draw("AL");
    TLine *l_trigger3 = new TLine(n_start, TMath::MinElement(N, &matched_filter[0]), n_start, TMath::MaxElement(N, &matched_filter[0]));
    l_trigger3 -> SetLineColor(kRed);
    l_trigger3 -> SetLineWidth(2);
    l_trigger3 -> Draw();








    return 0;
}