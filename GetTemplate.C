#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

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


int GetTemplate()
{
    string fname_root = "/home/riccardo/Documenti/NUSES/DeltaE_E/Sr90/wave0.root";
    string fname_template_out = "/home/riccardo/Documenti/NUSES/DeltaE_E/Sr90/template0.txt";


    TFile* file = new TFile(fname_root.c_str(), "READ");
    TTree* waves = (TTree*) file -> Get("waves");


    vector<double> templ;
    ExtractTemplate(waves, templ);

    ofstream file_out(fname_template_out.c_str());
    for (int i = 0; i < templ.size(); i++)
    {
        file_out << templ[i] << endl;
        cout << templ[i] << endl;
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    TGraph* gr = new TGraph(templ.size());
    for (int i = 0; i < templ.size(); i++)
    {
        gr -> SetPoint(i, i, templ[i]);
    }

    gr -> Draw("AL");
    string fname_graph_out = fname_template_out.substr(0, fname_template_out.size() - 4) + ".png";
    c1 -> SaveAs(fname_graph_out.c_str());


    return 0;

}

