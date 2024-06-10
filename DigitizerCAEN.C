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



// Filesystem namespace alias
namespace fs = std::filesystem;

/* ********************************************************************************************** */
/*                                         DATA STRUCTURE                                         */
/* ********************************************************************************************** */

struct Wave{
    int recordLength;           // Number of samples in the waveform
    int boardID;                // Board ID
    int channel;                // Channel number
    int eventNumber;            // Event number
    std::string pattern;        // Pattern
    long triggerTimeStamp;      // Trigger timestamp
    int dcOffset;               // DC offset
    std::vector<double> waveform;  // Waveform data
};

/* ********************************************************************************************** */
/*                                     EXAMPLE OF FILE EVENTS                                     */
/* ********************************************************************************************** */
/*
Record Length: 5000
BoardID: 31
Channel: 0
Event Number: 14322
Pattern: 0x0000
Trigger Time Stamp: 1595551493
DC offset (DAC): 0x3333
xxxxxx
xxxxxx
xxxxxx
xxxxxx
...
xxxxxx
Record Length: 5000
BoardID: 31
Channel: 0
Event Number: 14322
Pattern: 0x0000
Trigger Time Stamp: 1595551493
DC offset (DAC): 0x3333
...
*/





class DigitizerCAEN
{
public:
    /* ****************************************** VARIABLES ***************************************** */

    /* ******************************************* METHODS ****************************************** */
    DigitizerCAEN();
    ~DigitizerCAEN();

    // DBG PRINT
    void dbg_print(std::string TO_PRINT, int PRIORITY);
    void dbg_print(int TO_PRINT, int PRIORITY);
    void dbg_print(double TO_PRINT, int PRIORITY);
    void dbg_print(Wave wave, int PRIORITY);
    void dbg_print(std::vector<double> TO_PRINT, int PRIORITY);

    // FILE MANAGEMENT
    std::vector<std::string> getWaveFiles();
    std::vector<std::string> getFilesInDirectory();

    // Wave management
    Wave readSingleWave(std::ifstream &input);
    std::string readWaves(const std::string& filename);
    std::vector<std::string> processAllFiles();

    int quickScan();

    void startProcessing();

    /* ******************************************* GETTERS ****************************************** */
    int getVerbosity(){return VERBOSITY;};
    int getNToProcess(){return N_TO_PROCESS;};
    int getNWaveFiles(){return N_WAVE_FILES;};

    std::string getPathDigitizerFileFolder(){return PATH_DIGITIZER_FILE_FOLDER;};
    std::string getPathDestination(){return PATH_DESTINATION;};

    std::vector<std::string> getRootFiles(){return rootFiles;};

    /* ******************************************* SETTERS ****************************************** */
    void setVerbosity(int verbosity){VERBOSITY = verbosity;};
    void setPathDigitizerFileFolder(std::string path){PATH_DIGITIZER_FILE_FOLDER = path;};
    void setNToProcess(int n){N_TO_PROCESS = n;};
    void setPathDestination(std::string path){PATH_DESTINATION = path;};
    void setProgressBar(bool progress){PROGRESS_BAR = progress;};
private:
    /* ****************************************** VARIABLES ***************************************** */
    int VERBOSITY                           = 6;
    bool PROGRESS_BAR                       = true;
    int N_TO_PROCESS                        = 200;
    int N_WAVE_FILES                        = 0;
    int N_EVENTS                            = 0;
    std::string PATH_DIGITIZER_FILE_FOLDER  = "/media/riccardo/DATA/Sr90_300um_500um/RUN_0";
    std::string PATH_DESTINATION            = "/home/riccardo/Documenti/NUSES/DeltaE_E";

    std::vector<std::string> file_list;
    std::vector<std::string> countAndGetWaveFiles();
    std::vector<std::string> waveFiles;
    std::vector<std::string> rootFiles;
    std::vector<TTree*> trees;

};


/* ********************************************************************************************** */
/*                                      FUNCTION DESCRIPTION                                      */
/* ********************************************************************************************** */

DigitizerCAEN::DigitizerCAEN(){};
DigitizerCAEN::~DigitizerCAEN(){};

std::vector<std::string> DigitizerCAEN::getFilesInDirectory(){

    for (const auto& entry : fs::directory_iterator(PATH_DIGITIZER_FILE_FOLDER)) {
        if (fs::is_regular_file(entry.path())) {
            file_list.push_back(entry.path().string());
            dbg_print(entry.path().string(), 3);
        }
    }
    return file_list;
}

std::vector<std::string> DigitizerCAEN::countAndGetWaveFiles() {
    std::regex pattern(".*/wave(\\d+)\\.txt");
    N_WAVE_FILES = 0;
    for (const auto& file : file_list) {
        std::smatch match;
        if (std::regex_match(file, match, pattern)) {
            N_WAVE_FILES++;
            waveFiles.push_back(match[0]);            
            dbg_print(match[0], 3);
        }
    }
    return waveFiles;
}

std::vector<std::string> DigitizerCAEN::getWaveFiles()
{
    getFilesInDirectory();
    countAndGetWaveFiles();

    dbg_print("Counter Wave = ", 2);
    dbg_print(N_WAVE_FILES, 2);

    return waveFiles;
}

/* *************************************** Print functions ************************************** */

void DigitizerCAEN::dbg_print(std::string TO_PRINT, int PRIORITY)
{
    if(VERBOSITY > PRIORITY)
    {
        std::cout << TO_PRINT << std::endl;
    }
    return;
}


void DigitizerCAEN::dbg_print(int TO_PRINT, int PRIORITY)
{
    if(VERBOSITY > PRIORITY)
    {
        std::cout << TO_PRINT << std::endl;
    }
    return;
}

void DigitizerCAEN::dbg_print(double TO_PRINT, int PRIORITY)
{
    if(VERBOSITY > PRIORITY)
    {
        std::cout << TO_PRINT << std::endl;
    }
    return;
}

void DigitizerCAEN::dbg_print(std::vector<double> TO_PRINT, int PRIORITY)
{
    if(VERBOSITY > PRIORITY)
    {
        for(int i = 0; i < TO_PRINT.size(); i++)
        {
            std::cout << TO_PRINT[i] << "\t";
        }
        std::cout << std::endl;
    }
    return;
}

void DigitizerCAEN::dbg_print(Wave wave, int PRIORITY)
{
    dbg_print("Record Length: " + std::to_string(wave.recordLength), PRIORITY);
    dbg_print("Board ID: " + std::to_string(wave.boardID), PRIORITY);
    dbg_print("Channel: " + std::to_string(wave.channel), PRIORITY);
    dbg_print("Event Number: " + std::to_string(wave.eventNumber), PRIORITY);
    dbg_print("Pattern: " + wave.pattern, PRIORITY);
    dbg_print("Trigger Time Stamp: " + std::to_string(wave.triggerTimeStamp), PRIORITY);
    dbg_print("DC Offset: " + std::to_string(wave.dcOffset), PRIORITY);
    dbg_print("Waveform size: " + std::to_string(wave.waveform.size()), PRIORITY);


    for(int i = 0; i < wave.waveform.size(); i++)
    {
        dbg_print(wave.waveform[i], PRIORITY*2);
    }
    return;

}


/* ********************************************************************************************** */
/*                                   WAVE READING FROM TXT FILE                                   */
/* ********************************************************************************************** */

Wave DigitizerCAEN::readSingleWave(std::ifstream& input) {
    Wave wave;
    std::string line;

    bool HeaderProcessed = false;
    int HeaderLines = 0;
    int WaveformLines = 0;
    bool EventProcessed = false;
    
    while (!EventProcessed) {
        // If input is not valid, return the wave
        if (!input) {
            wave.recordLength =  -1;
            return wave;
        }
        std::getline(input, line);
        dbg_print("Reading line: " + line, 9);
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, ':')) {
            std::string value;
            std::getline(iss, value);
            if (key == "Record Length") {
                wave.recordLength = std::stoi(value);
                ++HeaderLines;
            } else if (key == "BoardID") {
                wave.boardID = std::stoi(value);
                ++HeaderLines;
            } else if (key == "Channel") {
                wave.channel = std::stoi(value);
                ++HeaderLines;
            } else if (key == "Event Number") {
                wave.eventNumber = std::stoi(value);
                ++HeaderLines;
            } else if (key == "Pattern") {
                wave.pattern = value;
                ++HeaderLines;
            } else if (key == "Trigger Time Stamp") {
                wave.triggerTimeStamp = std::stol(value);
                ++HeaderLines;
            } else if (key == "DC offset (DAC)") {
                wave.dcOffset = std::stoi(value);
                ++HeaderLines;
            }
            if (HeaderLines == 7) {
                HeaderProcessed = true;
                dbg_print("Header processed", 3);
                dbg_print(wave, 3);
            }
        } 
        
        if (HeaderProcessed && !EventProcessed) {
            // Leggi i campioni della waveform
            dbg_print("Reading waveform", 3);
            
            
            for (int i = 0; i < wave.recordLength; ++i) {
                // get the line with the integer
                std::getline(input, line);
                double sample;
                sample = std::stof(line);
                wave.waveform.push_back(sample);
            }
            EventProcessed = true;
            dbg_print("Waveform processed", 3);
        }
    }
    return wave;
}




std::string DigitizerCAEN::readWaves(const std::string& filename)
{
    // N_TO_PROCESS:
    // > 0: process N_TO_PROCESS waves
    // < 0: process all the waves in the file

    // Get the basename of the file
    std::string basename = filename.substr(filename.find_last_of('/') + 1);
    dbg_print("Reading waves from file: " + basename, 2);

    std::string basename_noext = basename.substr(0, basename.find_last_of('.'));

    // Replace the extension of the file with .root
    std::string rootFilename = PATH_DESTINATION + "/" + basename_noext + ".root";
    std::string txt_filtered = PATH_DESTINATION + "/" + basename_noext + "_filtered.txt";
    TFile *file = new TFile(rootFilename.c_str(), "RECREATE");
    TTree *waves = new TTree("waves", "Waveform data");

    dbg_print("Opening file: " + filename, 2);
    dbg_print("Creating ROOT file: " + rootFilename, 2);
    dbg_print("Creating TTree: waves", 2);
    


    int recordLength_evt;
    int boardID_evt;
    int channel_evt;
    int eventNumber_evt;
    std::string pattern_evt;
    int triggerTimeStamp_evt;
    int dcOffset_evt;
    std::vector<double> waveform_evt;

    waves->Branch("recordLength", &recordLength_evt, "recordLength/I");
    waves->Branch("boardID", &boardID_evt, "boardID/I");
    waves->Branch("channel", &channel_evt, "channel/I");
    waves->Branch("eventNumber", &eventNumber_evt, "eventNumber/I");
    waves->Branch("pattern", &pattern_evt, "pattern/C");
    waves->Branch("triggerTimeStamp", &triggerTimeStamp_evt, "triggerTimeStamp/L");
    waves->Branch("dcOffset", &dcOffset_evt, "dcOffset/I");
    waves->Branch("waveform", &waveform_evt[0], "waveform[recordLength]/D");

    dbg_print("Setting up branches", 2);



    std::ifstream
    input(filename);
    int counter = 0;
    
    dbg_print("Reading waves: begin while", 2);
    while (input) {
        dbg_print("Reading wave: begin", 3);
        Wave wave = readSingleWave(input);
        dbg_print("Reading wave: end", 3);
        if (input && wave.recordLength > 0) {
            recordLength_evt = wave.recordLength;
            boardID_evt = wave.boardID;
            channel_evt = wave.channel;
            eventNumber_evt = wave.eventNumber;
            pattern_evt = wave.pattern;
            triggerTimeStamp_evt = wave.triggerTimeStamp;
            dcOffset_evt = wave.dcOffset;
            waveform_evt = wave.waveform;
            waves -> SetBranchAddress("waveform", &waveform_evt[0]);
            dbg_print(waveform_evt, 12);

            waves->Fill();

            counter++;
            dbg_print("Wave " + std::to_string(counter) + " processed.", 3);
        }
        if(PROGRESS_BAR)
        {
            if(counter%100 == 0)
            {
                dbg_print("Events processed: " + std::to_string(counter), -1);
            }
        }

        if(counter >= N_TO_PROCESS && N_TO_PROCESS > 0)
        {
            break;
        }
    }
    file->Write();
    file->Close();
    return rootFilename;
}

std::vector<std::string> DigitizerCAEN::processAllFiles() {
    for (const auto& file : waveFiles) {
        dbg_print("Processing file: " + file, 2);
        rootFiles.push_back(readWaves(file));
    }
    return rootFiles;
}

int DigitizerCAEN::quickScan()
{
    std::ifstream input(waveFiles[0]);
    // Loop over the entire files and scan the number of events
    std::string line;
    while (std::getline(input, line)) {
        // If in the line there is "Record Length" then we have a new event
        if (line.find("Record Length") != std::string::npos) {
            N_EVENTS++;
            if(N_EVENTS%1000 == 0)
            {
                dbg_print("Events found: " + std::to_string(N_EVENTS), -1);
            }
        }
    }
    dbg_print("Events found: " + std::to_string(N_EVENTS), -1);
    return N_EVENTS;
}


void DigitizerCAEN::startProcessing()
{
    getWaveFiles();
    processAllFiles();
    return;
}