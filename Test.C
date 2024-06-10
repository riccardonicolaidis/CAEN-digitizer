#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"


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

int Test()
{

    string fname = "/home/riccardo/Documenti/NUSES/DeltaE_E/wave0.root";

    TFile* file = new TFile(fname.c_str(), "READ");
    TTree* waves = (TTree*)file -> Get("waves");

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

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    TGraph* gr = new TGraph(recordLength);
    for (int i = 0; i < recordLength; i++)
    {
        gr -> SetPoint(i, i, waveform_avg[i]);
    }
    gr -> Draw("ALP");



    return 0;
}
