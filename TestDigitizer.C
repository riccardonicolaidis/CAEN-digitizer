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

#include "DigitizerCAEN.C"

int TestDigitizer()
{
    DigitizerCAEN* digitizer = new DigitizerCAEN();

    digitizer -> setNToProcess(-1);
    digitizer -> setPathDigitizerFileFolder("/media/riccardo/DATA/Sr90_300um_500um/RUN_0");
    digitizer -> setPathDestination("/home/riccardo/Documenti/NUSES/DeltaE_E/");
    digitizer -> setVerbosity(3);
    digitizer -> setProgressBar(true);

    digitizer -> startProcessing();
    return 0;
}
