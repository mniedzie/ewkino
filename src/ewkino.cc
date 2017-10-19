//include ROOT classes
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TROOT.h"

//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <tuple>

//include other parts of the code
#include "../interface/treeReader.h"
#include "../interface/analysisTools.h"
#include "../interface/ewkinoTools.h"
#include "../plotting/plotCode.h"
#include "../plotting/tdrStyle.h"

//Temporary function to set _flavors correctly as it is buggy
void treeReader::setFlavors(){
    for(unsigned l = 0; l < _nL; ++l){
        if(l < _nMu) _lFlavor[l] = 1;
        else if(l < _nEle) _lFlavor[l] = 0;
        else _lFlavor[l] = 2;
    }
}

void treeReader::Analyze(){
    //Set CMS plotting style
    setTDRStyle();
    gROOT->SetBatch(kTRUE);
    //read samples and cross sections from txt file
    std::ifstream file("samples.txt");
    std::string line;
    std::vector<std::tuple<std::string, std::string, double> > samples;
    while(std::getline(file, line)){	
        samples.push_back(tools::readSampleLine(line));
    }
    //Always close file after usage
    file.close();
    for( auto it = samples.cbegin(); it != samples.cend(); ++it){
        std::cout << std::get<0>(*it) << "     " << std::get<1>(*it) << "      " << std::get<2>(*it) << std::endl;
    }	
    //info on kinematic distributions to plot
    std::vector< std::tuple < std::string, std::string, unsigned, double , double > > histInfo;
    //name      xlabel    nBins,  min, max
    histInfo = {
        std::make_tuple("sip3d", "SIP_{3D}", 100, 0, 8),
        std::make_tuple("dxy", "|d_{xy}| (cm)", 100, 0, 0.05),
        std::make_tuple("dz", "|d_{z}| (cm)", 100, 0, 0.1),
        std::make_tuple("miniIso", "miniIso", 100, 0, 0.4),
        std::make_tuple("leptonMva", "lepton MVA value", 100, -1, 1),
        std::make_tuple("ptRel", "P_{T}^{rel} (GeV)", 100, 0, 200),
        std::make_tuple("ptRatio", "P_{T}^{ratio}", 100, 0, 2),
        std::make_tuple("closestJetCsv", "closest jet CSV", 100, 0, 1),
        std::make_tuple("chargedTrackMult", "closest jet track multiplicity", 20, 0, 20),
        std::make_tuple("met", "E_{T}^{miss} (GeV)", 100, 0, 300),
        std::make_tuple("mll", "M_{ll} (GeV)", 200, 0, 200),
        std::make_tuple("leadPt", "P_{T}^{leading} (GeV)", 100, 25, 200),
        std::make_tuple("trailPt", "P_{T}^{trailing} (GeV)", 100, 15, 150),
        std::make_tuple("nVertex", "number of vertices", 100, 0, 100),
        std::make_tuple("nJets", "number of jets", 10, 0, 10),
        std::make_tuple("nBJets_CSVv2", "number of b-jets (CSVv2)", 5, 0, 8),
        std::make_tuple("nBJets_DeepCSV", "number of b-jets (Deep CSV)", 5, 0, 8)
    };

    //split histograms in run periods
    const unsigned nRuns = 7;
    const std::string runNames[nRuns] = {"all2017", "RunA", "RunB", "RunC", "RunD", "RunE", "RunF"};
    //split histograms in flavor combinations
    const unsigned nFlav = 4;
    const std::string flavNames[nFlav] = {"inclusive", "ee", "em", "mm"};

    const unsigned nDist = histInfo.size(); //number of distributions to plot
    //initialize vector holding all histograms
    std::vector< std::vector < std::vector< std::vector< TH1D* > > > > hists(nRuns);
    for(unsigned run = 0; run < nRuns; ++run){
        for(unsigned flav = 0; flav < nFlav; ++flav){
            hists[run].push_back(std::vector< std::vector < TH1D* > >() );
            for(unsigned dist = 0; dist < nDist; ++dist){
                hists[run][flav].push_back(std::vector < TH1D* >() );
                for(size_t sam = 0; sam < samples.size(); ++sam){
                    hists[run][flav][dist].push_back(nullptr);
                    hists[run][flav][dist][sam] = new TH1D( (const TString&) (std::get<1>(samples[sam]) + std::get<0>(histInfo[dist]) + flavNames[flav] + runNames[run] ), (const TString&) (std::get<1>(samples[sam]) + std::get<0>(histInfo[dist]) + flavNames[flav] + runNames[run]  + ";" + std::get<1>(histInfo[dist]) + ";Events" ),  std::get<2>(histInfo[dist]), std::get<3>(histInfo[dist]), std::get<4>(histInfo[dist]) );
                }
            }
        }
    }

    //tweakable options
    const double DataLuminosity = 18.90; //21.15; //in units of 1/fb
    const TString extra = ""; //for plot names

    //loop over all samples 
    for(size_t sam = 0; sam < samples.size(); ++sam){
        //read info from tree
        std::shared_ptr<TFile> sampleFile= std::make_shared<TFile>("../../ntuples_ewkino/"+ (const TString&) std::get<1>(samples[sam]),"read");
        sampleFile->cd("blackJackAndHookers");
        //determine hcounter for cross section scaling
        double hCounter;
        if(sam != 0){
            TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
            _hCounter->Read("hCounter");
            hCounter = _hCounter->GetBinContent(1);
            delete _hCounter;
        }
        std::shared_ptr<TTree> sampleTree = std::shared_ptr<TTree> ( (TTree*) (sampleFile->Get("blackJackAndHookers/blackJackAndHookersTree")) );
        Init(sampleTree.get(), sam == 0); //don't store generator info for data ( first entry )

        double scale;
        if(sam != 0) scale = std::get<2>(samples[sam])*DataLuminosity*1000/(hCounter);
        //find number of entries in sample
        long unsigned nEntries = sampleTree->GetEntries();
        std::cout<<"Entries in "<< std::get<1>(samples[sam]) << " " << nEntries << std::endl;
        double progress = 0; 	//for printing progress bar
        for(long unsigned it = 0; it < nEntries/100; ++it){
            //print progress bar	
            if(it%100 == 0 && it != 0){
                progress += (double) (100./nEntries);
                tools::printProgress(progress);
            } else if(it == nEntries -1){
                progress = 1.;
                tools::printProgress(progress);
            }
            sampleTree->GetEntry(it);
            double weight;
            if(sam == 0) weight = 1;
            else weight = scale*_weight;	
            //set flavors (temporary) 
            setFlavors();
            //vector containing good lepton indices
            std::vector<unsigned> ind(2);
            //select leptons
            const unsigned lCount = selectLep(ind);
            if(lCount < 2) continue;
            //require pt cuts (25, 20) to be passed
            if(!passPtCuts(ind)) continue;
            //determine flavor compositions
            unsigned flav = dilFlavorComb(ind) + 1; //reserve 0 for inclusive
            //determine run perios
            unsigned run = ewk::runPeriod(_runNb) + 1; //reserve 0 for inclusive
            //loop over leading leptons
            for(unsigned l = 0; l < 2; ++l){
                double fill[9] = {_3dIPSig[ind[l]], _dxy[ind[l]], _dz[ind[l]], _miniIso[ind[l]], _leptonMva[ind[l]], _ptRel[ind[l]], _ptRatio[ind[l]], _closestJetCsv[ind[l]], (double) _selectedTrackMult[ind[l]]};
                //fill histograms
                for(unsigned dist = 0; dist < 9; ++dist){
                    hists[run][flav][dist][sam]->Fill(fill[dist], weight);
                    hists[run][0][dist][sam]->Fill(fill[dist], _weight);
                    hists[0][0][dist][sam]->Fill(fill[dist], _weight);
                    hists[0][flav][dist][sam]->Fill(fill[dist], _weight);
                }
            }
            //make lorentzvectors for leptons
            TLorentzVector lepV[lCount];
            for(unsigned l = 0; l < lCount; ++l) lepV[l].SetPtEtaPhiE(_lPt[ind[l]], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
            double fill[nDist - 9] = {0, (lepV[0] + lepV[1]).M(), _lPt[ind[0]], _lPt[ind[1]], (double) _nVertex, (double) nJets(), (double) nBJets(0, false), (double) nBJets()}; //replace 0 by _met for correct trees
            for(unsigned dist = 9; dist < nDist; ++dist){
                hists[run][flav][dist][sam]->Fill(fill[dist - 9], weight);
                hists[run][0][dist][sam]->Fill(fill[dist - 9], weight);
                hists[0][0][dist][sam]->Fill(fill[dist - 9], weight);
                hists[0][flav][dist][sam]->Fill(fill[dist - 9], weight);
            }
        }
        //set histograms to 0 if negative
        for(unsigned run = 0; run < nRuns; ++run){
            for(unsigned flav = 0; flav < nFlav; ++flav){
                for(unsigned dist = 0; dist < nDist; ++dist){
                    tools::setNegativeZero(hists[run][flav][dist][sam]);
                }	
            }
        }
    }
    //merge histograms with the same physical background
    std::vector<std::string> proc = {"obs.", "DY", "TT + Jets", "WJets", "VV", "TT + X", "T + X"};
    std::vector< std::vector <std::vector < std::vector<TH1D*> > > > mergedHists(nRuns);
    for(unsigned run = 0; run < nRuns; ++run){
        for(unsigned flav = 0; flav < nFlav; ++flav){
            mergedHists[run].push_back(std::vector< std::vector <TH1D*> >() );
            for(unsigned dist = 0; dist < nDist; ++dist){
                mergedHists[run][flav].push_back(std::vector<TH1D*>(proc.size()) );
                for(size_t m = 0, sam = 0; m < proc.size(); ++m){
                    mergedHists[run][flav][dist][m] = (TH1D*) hists[run][flav][dist][sam]->Clone();
                    while(sam < samples.size() - 1 && std::get<0>(samples[sam]) == std::get<0>(samples[sam + 1]) ){
                        mergedHists[run][flav][dist][m]->Add(hists[run][flav][dist][sam + 1]);
                        ++sam;
                    }
                    ++sam;
                }
            }
        }
    }
    //plot all distributions
    for(unsigned run = 0; run < nRuns; ++run){
        for(unsigned flav = 0; flav < nFlav; ++flav){
            for(unsigned dist = 0; dist < nDist; ++dist){
                plotDataVSMC(mergedHists[run][flav][dist][0], &mergedHists[run][flav][dist][1], &proc[0], mergedHists[run][flav][dist].size() - 1, runNames[run] + "/" + flavNames[flav] + "/" + std::get<0>(histInfo[dist]) + "_" + flavNames[flav] + "_" + runNames[run], "ewkinoDilep", false, true);             //linear plots
                plotDataVSMC(mergedHists[run][flav][dist][0], &mergedHists[run][flav][dist][1], &proc[0], mergedHists[run][flav][dist].size() - 1, runNames[run] + "/" + flavNames[flav] + "/" + std::get<0>(histInfo[dist]) + "_" + flavNames[flav] + "_" + runNames[run] + "_log", "ewkinoDilep", true, true);     //log plots
            }
        }
    }
}

int main(){
    treeReader reader;
    reader.Analyze();
    return 0;
}
