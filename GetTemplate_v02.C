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
#include "TString.h"
#include "TMath.h"

using namespace std;

void sum_vector(vector<double>& v_destination, vector<double>& v_toSum)
{
    for (int i = 0; i < v_destination.size(); i++)
    {
        v_destination[i] += v_toSum[i];
    }
}


void sum_vector_shifted(vector<double>& v_destination, vector<double>& v_toSum, int shift)
{
    for (int i = 0; i < v_destination.size(); i++)
    {
        if ((i + shift) >= 0 && (i + shift) < v_toSum.size())
        {
            v_destination[i] += v_toSum[i + shift];
        }
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


int CFD_detection(vector<double>& v, double CFD)
{
    double max = v[0];
    double min = v[0];

    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] > max)
        {
            max = v[i];
        }
        if (v[i] < min)
        {
            min = v[i];
        }
    }


    double threshold = min + CFD * (max - min);

    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] > threshold)
        {
            return i;
        }
    }

    return -1;
}





int GetTemplate_v02
(
    TString fName_TTree = "/home/riccardo/Documenti/NUSES/DeltaE_E/Cremat/Output/Template_wave0.root",
    TString fName_Template = "/home/riccardo/Documenti/NUSES/DeltaE_E/Cremat/Output/Template_wave0.txt",
    double postTrigger = 0.5,
    double CFD_fraction = 0.4,
    int discard_pts = 200
)
{
    TFile *f = new TFile(fName_TTree, "READ");
    TTree *waves = (TTree*)f -> Get("waves");





    int recordLength;
    int N_sum = 0;
    vector<double> waveform;
    vector<double> waveform_avg;

    waves -> SetBranchAddress("recordLength", &recordLength);
    waves -> GetEntry(0);
    waveform.resize(recordLength);
    waveform_avg.resize(recordLength);

    for(int i = 0; i < waveform_avg.size(); i++)
    {
        waveform_avg[i] = 0;
    }

    waves -> SetBranchAddress("waveform", &waveform[0]);

    int index_Trigger_theory = (int) TMath::Floor(recordLength * postTrigger);
    int index_Trigger = 0;


    for(int i = 0; i <waves -> GetEntries(); i++)
    {
        waves -> GetEntry(i);
        vector<double> waveform_tmp = waveform;
        //TCanvas *c = new TCanvas("c", "c", 800, 600);
        //TGraph *gr = new TGraph(waveform_tmp.size());

        // for (int i = 0; i < waveform_tmp.size(); i++)
        // {
        //     gr -> SetPoint(i, i, waveform_tmp[i]);
        // }

        // gr -> Draw("AL");

        // c -> SaveAs("waveform.pdf");

        //baseline_correction(waveform_tmp, 0, 10000);
        //detrending(waveform_tmp, 0, 10000);

        index_Trigger = CFD_detection(waveform_tmp, CFD_fraction);
        if(index_Trigger < index_Trigger_theory*0.4)
        {
            cout << "Event " << i << " - Trigger not found" << endl;
            continue;
        }

        double max = TMath::MaxElement(waveform_tmp.size(), &waveform_tmp[0]);
        if(max > 14000)
        {
            cout << "Saturation" << endl;
            continue;
        }

        if(max < 7000)
        {
            cout << "Low signal" << endl;
            continue;
        }




        sum_vector_shifted(waveform_avg, waveform, -index_Trigger_theory + index_Trigger);
        N_sum++;

        cout << "Event " << i << " - Trigger: " << index_Trigger << endl;

    }


    divide_vector(waveform_avg, N_sum);

    baseline_correction(waveform_avg, discard_pts+1, 1000);
    detrending(waveform_avg, discard_pts+1, 1000);



    for(int i = 0; i < discard_pts; i++)
    {
        waveform_avg[i] = 0;
    }





    int max = waveform_avg[0];

    for (int i = 0; i < waveform_avg.size(); i++)
    {
        if (waveform_avg[i] > max)
        {
            max = waveform_avg[i];
        }
    }

    scale_vector(waveform_avg, 1.0 / max);


    int index_negative = 0;

    for(int i = index_Trigger_theory; i < waveform_avg.size(); i++)
    {
        if(waveform_avg[i] < 0)
        {
            index_negative = i;
            break;
        }
    }

    for(int i = index_negative; i < waveform_avg.size(); i++)
    {
        waveform_avg[i] = 0;
    }



    ofstream file(fName_Template);
    for (int i = 0; i < waveform_avg.size(); i++)
    {
        file << waveform_avg[i] << endl;
    }
    file.close();

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    TGraph *gr = new TGraph(waveform_avg.size());
    for (int i = 0; i < waveform_avg.size(); i++)
    {
        gr -> SetPoint(i, i, waveform_avg[i]);
    }
    gr -> Draw("AL");
    c -> SaveAs(fName_Template + ".pdf");



    return 0;
}