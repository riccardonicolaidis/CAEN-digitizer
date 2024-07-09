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


int template_fitting_roll(vector<double>& templ, vector<double>& wave, vector<double>&fit, int n_roll_start, int n_roll_end, int increment, int index_id = 0, TString path_output_diagnostic = "")
{
    double a_max = 0;
    int N = templ.size();

    vector<double> coeff_a;

    int max_i = n_roll_start;
    for(int i = n_roll_start; i < n_roll_end; i += increment)
    {
        double sum_yi = 0;
        double sum_xi = 0;
        double sum_yixi = 0;
        double sum_xi2 = 0;
        int N_sum = 0;

        for (int j = 0; j < templ.size(); j++)
        {
            if(templ[(j + i + N)%N] < 0)
            {
                continue;
            }
            N_sum++;
            sum_yi += wave[(j + i + N)%N];
            sum_xi += templ[j];
            sum_yixi += wave[(j + i + N)%N] * templ[j];
            sum_xi2 += templ[j] * templ[j];
        }

        double avg_yi = sum_yi / N_sum;
        double avg_xi = sum_xi / N_sum;
        double avg_yixi = sum_yixi / N_sum;
        double avg_xi2 = sum_xi2 / N_sum;

        double cov_xy = avg_yixi - avg_yi * avg_xi;
        double var_x = avg_xi2 - avg_xi * avg_xi;

        double a = cov_xy / var_x;
        double b = avg_yi - a * avg_xi;

        coeff_a.push_back(a);

        //cout << "a: " << a << " b: " << b << endl;
        //cout << "a: " << a << " b: " << b << endl;
        //cout << "sum_yi: " << sum_yi << " sum_xi: " << sum_xi << " sum_yixi: " << sum_yixi << " sum_xi2: " << sum_xi2 << endl;
        if(a > a_max)
        {
            a_max = a;
            fit[0] = a;
            fit[1] = b;
            max_i = i;
        }
    }

    

    if(path_output_diagnostic != "")
    {
        TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
        TGraph *g1 = new TGraph(coeff_a.size());
        for(int i = 0; i < coeff_a.size(); i++)
        {
            g1 -> SetPoint(i, i, coeff_a[i]);
        }
        g1 -> Draw("ALP");
        c1 -> SaveAs(path_output_diagnostic + Form("Diagnostic_%d.png", index_id));
    }
    
    return max_i;
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




int MatchedFilter
(
    TString path_template = "/home/riccardo/Documenti/NUSES/DeltaE_E/Cremat/Output/Template_wave0.txt",
    TString path_waveforms = "/home/riccardo/Documenti/NUSES/DeltaE_E/Cremat/Output/Am_wave0.root",
    TString path_output = "/home/riccardo/Documenti/NUSES/DeltaE_E/Cremat/Output/MatchedFilter.root",
    double CFD = 0.6
)
{
    ifstream file_template(path_template);

    vector<double> templ;

    if (file_template.is_open())
    {
        double value;
        while (file_template >> value)
        {
            templ.push_back(value);
            cout << value << endl;
        }
    }
    else
    {
        cout << "Error: the template file could not be opened" << endl;
        return 1;
    }


    TFile* file_waveforms = new TFile(path_waveforms, "READ");
    TTree* waves = (TTree*) file_waveforms -> Get("waves");

    vector<double> waveform;
    int recordLength;
    waves -> SetBranchAddress("recordLength", &recordLength);
    waves -> GetEntry(0);
    waveform.resize(recordLength);
    waves -> SetBranchAddress("waveform", &waveform[0]);

    // Extend the template to the length of the waveform

    vector<double> templ_extended;
    templ_extended.resize(waveform.size());
    for (int i = 0; i < templ.size(); i++)
    {
        templ_extended[i] = templ[i];
    }


    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    for(int i = 0; i < 40;++i)
    {
        waves -> GetEntry(i);
        TGraph *gr = new TGraph();
        for(int j = 0; j < waveform.size(); j++)
        {
            gr -> SetPoint(j, j, waveform[j]);
        }

        double max_wvf = TMath::MaxElement(waveform.size(), &waveform[0]);



        if(i == 0)
        gr -> Draw("ALP");
        else if(max_wvf > 3000)
        gr -> Draw("ALPsame");
    }
    c1 -> SaveAs("waveforms.png");






    TFile *file_output = new TFile(path_output, "RECREATE");
    TTree *matchedFilter = new TTree("matchedFilter", "Matched Filter");

    
    int index_Trigger_template = CFD_detection(templ, CFD);


    vector<double> fit;
    fit.resize(2);

    matchedFilter -> Branch("a", &fit[0], "a/D");
    matchedFilter -> Branch("b", &fit[1], "b/D");
    int max_i;
    matchedFilter -> Branch("max_i", &max_i, "max_i/I");

    for (int i = 0; i < waves -> GetEntries(); i++)
    {
        waves -> GetEntry(i);
        int index_Trigger = CFD_detection(waveform, CFD);
        int shift = index_Trigger - index_Trigger_template;

        vector<double> waveform_tmp = waveform;


        bool diagnostic = false;

        TString diagnostic_path = "/home/riccardo/Documenti/NUSES/DeltaE_E/Cremat/Output/Diagnostics/";
        if(diagnostic)
        {
            diagnostic_path = "/home/riccardo/Documenti/NUSES/DeltaE_E/Cremat/Output/Diagnostics/";
        }
        else
        {
            diagnostic_path = "";
        }



        int template_fit_roll = 2500;
        int increment = 100;
        max_i = template_fitting_roll(templ_extended, waveform_tmp, fit, shift - template_fit_roll, shift + template_fit_roll, increment, i,
        diagnostic_path
        ); // 

        max_i = template_fitting_roll(templ_extended, waveform_tmp, fit, shift - template_fit_roll, shift + template_fit_roll, increment, i,
        diagnostic_path
        ); //

        //cout << "Fitting at sample: " << i << " a: " << fit[0] << " b: " << fit[1] << " max_i = "<< max_i << endl;
        template_fit_roll = 200;
        increment = 10;

        max_i = template_fitting_roll(templ_extended, waveform_tmp, fit, max_i - template_fit_roll, max_i + template_fit_roll, increment, i*2323,
        diagnostic_path
        );

        
        template_fit_roll = 20;
        increment = 1;

        int max_i_new = template_fitting_roll(templ_extended, waveform_tmp, fit, max_i - template_fit_roll, max_i + template_fit_roll, increment, i*2323,
        diagnostic_path
        );

        while(max_i_new != max_i)
        {
            max_i = max_i_new;
            max_i_new = template_fitting_roll(templ_extended, waveform_tmp, fit, max_i - template_fit_roll, max_i + template_fit_roll, increment, i*2323,
            diagnostic_path
            );
        }


        //cout << "Fitting at sample: " << i << " a: " << fit[0] << " b: " << fit[1] << " max_i = "<< max_i << endl;

        cout << "Fitting at sample: " << i << " a: " << fit[0] << " b: " << fit[1] << " max_i = "<< max_i << endl;
        matchedFilter -> Fill();
    }

    matchedFilter -> Write();
    


    return 0;
}