//
//
//include c++ library classes

//include ROOT classes 
#include "TTree.h"

//include general parts of framework
#include "../Tools/interface/analysisTools.h"
#include "../Tools/interface/systemTools.h"
#include "../Tools/interface/stringTools.h"
#include "../Tools/interface/HistInfo.h"
#include "../weights/interface/ConcreteReweighterFactory.h"
#include "../Tools/interface/SusyScan.h"
#include "../Tools/interface/histogramTools.h"
#include "../plotting/plotCode.h"
#include "../plotting/tdrStyle.h"

#include "../Tools/interface/SampleCrossSections.h"

//include ewkino specific code
#include "interface/ttZSelection.h"
#include "interface/ttZVariables.h"
//#include "interface/ttZSearchRegions.h"

#include <functional> 

//compare floating points
bool floatEqual( const double lhs, const double rhs ){
    return ( fabs( ( lhs - rhs ) / lhs ) < 1e-6 );
}


//build histograms 
std::vector< HistInfo > makeDistributionInfo(){

    //make general plots
    std::vector< HistInfo > histInfoVec;
    histInfoVec = {

        HistInfo( "leptonPtLeading", "p_{T}^{leading lepton} (GeV)", 13, 40, 300 ),
        HistInfo( "leptonPtSubLeading", "p_{T}^{subleading lepton} (GeV)", 14, 15, 155 ),
        HistInfo( "leptonPtTrailing", "P_{T}^{trailing lepton} (GeV)", 14, 15, 155 ),

        HistInfo( "leptonEtaLeading", "|#eta|^{leading lepton}", 10, 0, 2.5 ),
        HistInfo( "leptonEtaSubLeading", "|#eta|^{subleading lepton}", 10, 0, 2.5 ),
        HistInfo( "leptonEtaTrailing", "|#eta|^{trailing lepton}", 10, 0, 2.5 ),

        HistInfo( "met", "E_{T}^{miss} (GeV)", 10, 0, 200 ),
        HistInfo( "nTI", "nTI (GeV)", 100, 0, 100 ),
        HistInfo( "mll", "M_{ll} (GeV)", 10, 75, 105 ),
        HistInfo( "ht", "H_{T} (GeV)", 10, 0, 800 ),

        HistInfo( "nJets", "number of jets", 7, 1, 8 ),
        HistInfo( "nBJets", "number of b-jets (medium deep flavor)", 5, -0.5, 4.5 ),
        HistInfo( "nVertex", "number of vertices", 30, 0, 70 ),

        HistInfo( "ttZSR", "Signal Regions", 28, 0, 28 ),
//        HistInfo( "ttZSR", "Signal Regions", 14, 0, 14 ),
        HistInfo( "ttZFlav", "Lepton flavors", 4, 0, 4 ),

        HistInfo( "ptZ", "p_{T}^{Z} [GeV]", 16, 0, 400 ),
        HistInfo( "cosThetaStar", "cos(#theta^{*})", 5, -1, 1 ),

    };
    return histInfoVec;
}


std::vector< double > buildFillingVector( Event& event, const std::string& uncertainty){
    
    auto varMap = ttZ::computeVariables( event, uncertainty );
    std::vector< double > fillValues;
    unsigned searchttZ = ttZ::SR_main( event, event.numberOfFOLeptons(), ttZ::numberOfVariedJets( event, uncertainty ), ttZ::numberOfVariedBJets( event, uncertainty ) );
    unsigned ttZFlav = ttZ::ttZFlavPlot( event );
    fillValues = {
        event.lepton( 0 ).pt(),
        event.lepton( 1 ).pt(),
        event.lepton( 2 ).pt(),

        event.lepton( 0 ).absEta(),
        event.lepton( 1 ).absEta(),
        event.lepton( 2 ).absEta(),

        varMap.at("met"),
        varMap.at("nTI"),
        varMap.at("mll"),
        varMap.at("ht"),

        varMap.at("numberOfJets"),
        varMap.at("numberOfBJets"),
        static_cast< double >( event.numberOfVertices() ),

        static_cast< double >( searchttZ ),
        static_cast< double >( ttZFlav ),

        varMap.at("ptZ"),
        varMap.at("cosThetaStar"),
    };
    
    return fillValues;
}


//auto varMap = ewkino::computeVariables( event, uncertainty );
//unsigned searchROld = ewkino::SR_EWK_3lOSSF_old( varMap.at("mtW"), varMap.at("met"), varMap.at("mll") );
//unsigned searchRNew = ewkino::SR_EWK_3lOSSF_new( varMap.at("mll"), varMap.at("mtW"), varMap.at("mt3l"), varMap.at("met"), varMap.at("ht") );



void analyze( const std::string& year, const std::string& controlRegion, const std::string& sampleDirectoryPath , const std::string& procName){


    analysisTools::checkYearString( year );

    //selection that defines the control region
    const std::map< std::string, std::function< bool (Event&, const std::string&) > > crSelectionFunctionMap{
        { "ttZ", ttZ::passSelectionTTZ},
        { "ttZclean", ttZ::passSelectionTTZclean},
        { "WZ", ttZ::passSelectionWZCR},
        { "DY", ttZ::passSelectionDYCR},
        { "ttbar", ttZ::passSelectionttbarCR},
        { "XGamma", ttZ::passSelectionXgammaCR},
        { "ZZ", ttZ::passSelectionZZCR},
        { "TTZ", ttZ::passSelectionTTZ }
    };

    auto passSelection = crSelectionFunctionMap.at( controlRegion );


    //build TreeReader and loop over samples
    std::cout << "building treeReader" << std::endl;
    TreeReader treeReader( "sampleLists/samples_ttZ_" + year + ".txt", sampleDirectoryPath );

    //build ttZ reweighter
    std::cout << "building reweighter" << std::endl;
    std::shared_ptr< ReweighterFactory >reweighterFactory( new ttZReweighterFactory() );
    CombinedReweighter reweighter = reweighterFactory->buildReweighter( "../weights/", year, treeReader.sampleVector() );

    //read FR maps
    std::cout << "building FR maps" << std::endl;
//    TFile* frFileMuons = TFile::Open( ( "frMaps/muFR_QCD_MC_Marek_" + year + ".root" ).c_str() );
//    std::shared_ptr< TH2 > frMapMuons = std::shared_ptr< TH2 >( dynamic_cast< TH2* >( frFileMuons->Get( "passed" ) ) );
    TFile* frFileMuons = TFile::Open( ( "frMaps/fakeRateMap_data_muon_" + year + "_mT.root" ).c_str() );
    std::shared_ptr< TH2 > frMapMuons = std::shared_ptr< TH2 >( dynamic_cast< TH2* >( frFileMuons->Get( ("fakeRate_muon_" + year ).c_str() ) ) );
    frMapMuons->SetDirectory( gROOT );
    frFileMuons->Close();

//    TFile* frFileElectrons = TFile::Open( ( "frMaps/elFR_QCD_MC_Marek_" + year + ".root" ).c_str() );
//    std::shared_ptr< TH2 > frMapElectrons = std::shared_ptr< TH2 >( dynamic_cast< TH2* >( frFileElectrons->Get( "passed" ) ) );
    TFile* frFileElectrons = TFile::Open( ( "frMaps/fakeRateMap_data_electron_" + year + "_mT.root" ).c_str() );
    std::shared_ptr< TH2 > frMapElectrons = std::shared_ptr< TH2 >( dynamic_cast< TH2* >( frFileElectrons->Get( ( "fakeRate_electron_" + year ).c_str() ) ) );
    frMapElectrons->SetDirectory( gROOT );
    frFileElectrons->Close();

    //histogram collection, histInfoVector will contain histinfo on all histograms defined at the top of the file withink the makeDist... function
    std::cout << "building histograms" << std::endl;
    std::vector< HistInfo > histInfoVector = makeDistributionInfo();

    //make histograms for each process, and integral signal to check shapes
    //add an additional histogram for the nonprompt prediction
    std::vector< Sample > sampleVec = treeReader.sampleVector();
    std::vector< std::vector< std::shared_ptr< TH1D > > > histograms( histInfoVector.size(), std::vector< std::shared_ptr< TH1D > >( sampleVec.size() + 1 ) );
    // loop over all histograms to be made
    for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
        // loop over all sample types considered.
        for( size_t p = 0; p < sampleVec.size() + 1; ++p ){
            if( p < sampleVec.size() ){
                histograms[ dist ][ p ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_" + sampleVec[p].uniqueName() );
            } else {
                histograms[ dist ][ p ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_nonprompt" );
            }
        }
    }

    //const std::vector< std::string > shapeUncNames = {  "JEC_" + year, "JER_" + year, "scale", "pileup", "prefire", "lepton_reco", "lepton_id" }; //, "pdf" }; //"scaleXsec", "pdfXsec" }
    const std::vector< std::string > shapeUncNames = {  "JEC_" + year, "JER_" + year, "scale", "pileup", "bTag_heavy_" + year, "bTag_light_" + year, "prefire", "lepton_reco", "lepton_id", "pdf" }; //, "scaleXsec", "pdfXsec" }
    std::map< std::string, int > shapesIntMap{{"JEC_" + year, 0}, { "JER_" + year, 1}, {"scale", 2}, {"pileup", 3}, {"bTag_heavy_" + year, 4}, {"bTag_light_" + year, 5}, {"prefire", 6}, {"lepton_reco", 7}, {"lepton_id", 8}, {"pdf" , 9}};

    unsigned numberOfPdfVariations = 100;
    // histogram for each uncertianty, dist, process
    std::map< std::string, std::vector< std::vector< std::shared_ptr< TH1D > > > > histogramsUncDown;
    std::map< std::string, std::vector< std::vector< std::shared_ptr< TH1D > > > > histogramsUncUp;
    // histogram for each uncertianty, dist
//    std::map< std::string, std::vector< std::shared_ptr< TH1D > > > histogramsUncDecompDown;
//    std::map< std::string, std::vector< std::shared_ptr< TH1D > > > histogramsUncDecompUp;
    std::vector<  std::vector< std::shared_ptr< TH1D > > > histogramsUncDecompDown( histInfoVector.size(), std::vector< std::shared_ptr< TH1D > >( shapeUncNames.size() ) );
    std::vector<  std::vector< std::shared_ptr< TH1D > > > histogramsUncDecompUp( histInfoVector.size(), std::vector< std::shared_ptr< TH1D > >( shapeUncNames.size() ) );
    // for each uncertainty
    for( const auto& unc : shapeUncNames ){
        histogramsUncDown[ unc ] = std::vector< std::vector< std::shared_ptr< TH1D > > >( histInfoVector.size(), std::vector< std::shared_ptr< TH1D > >( sampleVec.size() + 1 )  );
        histogramsUncUp[ unc ] = std::vector< std::vector< std::shared_ptr< TH1D > > >( histInfoVector.size(), std::vector< std::shared_ptr< TH1D > >( sampleVec.size() + 1 )  );

//        histogramsUncDecompDown[ unc ] = std::vector< std::shared_ptr< TH1D > >( histInfoVector.size() );
//        histogramsUncDecompUp[ unc ]   = std::vector< std::shared_ptr< TH1D > >( histInfoVector.size() );
        
        //create a set of all histograms
        for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
            // for each type of sample
            for( size_t p = 0; p < sampleVec.size() + 1; ++p ){
                if( p < sampleVec.size() ){
                    histogramsUncDown[ unc ][ dist ][ p ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_" + sampleVec[p].uniqueName() + unc + "Down" );
                    histogramsUncUp[ unc ][ dist ][ p ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_" + sampleVec[p].uniqueName() + unc + "Up" );
                } else {
                    histogramsUncDown[ unc ][ dist ][ p ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_nonprompt"  + unc + "Down" );
                    histogramsUncUp[ unc ][ dist ][ p ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_nonprompt" + unc + "Up" );
                }
            }
            histogramsUncDecompDown[ dist ][ shapesIntMap.at(unc) ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_" + unc + "Combined_Down" );
            histogramsUncDecompUp[ dist ][ shapesIntMap.at(unc) ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_" + unc + "Combined_Up" );
        }
    }

    // create histograms for each pdf variation
    std::map< unsigned, std::vector< std::vector< std::shared_ptr< TH1D > > > > histogramsPDFVars;
    // hard coded number of pdf variations!!
    for( unsigned pdf_i = 0; pdf_i < numberOfPdfVariations; ++pdf_i){
        histogramsPDFVars[ pdf_i ] = std::vector< std::vector< std::shared_ptr< TH1D > > >( histInfoVector.size(), std::vector< std::shared_ptr< TH1D > >( sampleVec.size() + 1 ) );
        //create a set of all histograms
        for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
            // for each pdf variation
            for( size_t p = 0; p < sampleVec.size() + 1; ++p ){
                if( p < sampleVec.size() ){
                    histogramsPDFVars[ pdf_i ][ dist ][ p ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_" + sampleVec[p].uniqueName() + "_pdf_" + std::to_string( pdf_i ) );
                } else {
                    histogramsPDFVars[ pdf_i ][ dist ][ p ] = histInfoVector[ dist ].makeHist( histInfoVector[ dist ].name() + "_nonprompt_pdf_" + std::to_string( pdf_i ) );
                }
            }
        }
    }

    std::cout << "event loop" << std::endl;

    bool FRfromMC = false;

    std::cout << "Fake rate from data? " << (FRfromMC ? "no":"yes") << std::endl;

    for( unsigned sampleIndex = 0; sampleIndex < treeReader.numberOfSamples(); ++sampleIndex ){
        treeReader.initSample();

        std::cout << treeReader.currentSample().fileName() << std::endl;

        bool sampleHasPdfAndScale = true;

        for( long unsigned entry = 0; entry < treeReader.numberOfEntries(); ++entry ){
            if ( sampleVec[ sampleIndex ].processName() != procName && procName != "all" ) break;
            Event event = treeReader.buildEvent( entry );
            
            //check if sample has pdf and scale information stored, try-catch blocks for every event are too slow
            if( entry == 0 && treeReader.isMC() ){
                try{
                    event.generatorInfo().relativeWeightPdfVar( 99 );
                }catch( std::out_of_range& ){
                    sampleHasPdfAndScale = false;
                }
            }

//            if(entry > treeReader.numberOfEntries()/20) break;            
//            if(entry > 1000) break;
            //apply baseline selection

            if( !ttZ::passBaselineSelection( event, true, true ) ) continue;

            //apply lepton pT cuts
            if( !ttZ::passPtCuts( event ) ) continue;

            //require triggers
            if( !ttZ::passTriggerSelection( event ) ) continue;
            if( !( event.passMetFilters() ) ) continue;

            //remove photon overlap
//            if( !ttZ::passPhotonOverlapRemoval( event ) ) continue;

            //require the right number of tight and FO(loose) leptons in 3(4) lepton events
            if( !ttZ::passSelectionLNumber( event ) ) continue;

            //require MC events to only contain prompt leptons
            size_t fillIndex = sampleIndex;
            if( event.isMC() && !ttZ::leptonsArePrompt( event ) ){
                if ( FRfromMC ) fillIndex = treeReader.numberOfSamples();
                else continue;
            }

            //apply scale-factors and reweighting
            double weight = event.weight();
            if( event.isMC() ){
                weight *= reweighter.totalWeight( event );
            }

            //apply fake-rate weight
            if( !ttZ::leptonsAreTight( event ) ){
                if( FRfromMC ) continue;
                fillIndex = treeReader.numberOfSamples();
                weight *= ttZ::fakeRateWeight( event, frMapMuons, frMapElectrons );
                if( event.isMC() ) weight *= -1.;
            }

//            // in 4L nonptrompt is calcualted from MC, should add WZ to nonprompt too.
//            if( event.isMC() && event.numberOfTightLeptons() == 4 && !ttZ::leptonsArePrompt( event ) ){
//                fillIndex = treeReader.numberOfSamples();
//            }
            
//            // for nonprompt from MC, reject events with 
//            if( event.numberOfTightLeptons() < 3 && FRfromMC ) continue;

            //fill nominal histograms
//            if( passSelection( event, "nominal" ) || ttZ::passSelectionWZCR( event, "nominal" ) ){
            if( passSelection( event, "nominal" ) ){
                auto fillValues = buildFillingVector( event, "nominal" );
                for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
                    histogram::fillValue( histograms[ dist ][ fillIndex ].get(), fillValues[ dist ], weight );
                }

                //in case of data fakes fill all uncertainties for nonprompt with nominal values
                if( event.isData() && ( fillIndex == treeReader.numberOfSamples() ) ){
                    for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
                        for( const auto& key : shapeUncNames ){
                            histogram::fillValue( histogramsUncDown[ key ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight );
                            histogram::fillValue( histogramsUncUp[ key ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight );       
                        }
                    }
                }

            }

            //no uncertainties for data
            if( event.isData() ) continue;
            
            //fill JEC down histograms
            if( passSelection( event, "JECDown" ) ){
                auto fillValues = buildFillingVector( event, "JECDown" );
                for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
                    histogram::fillValue( histogramsUncDown[ "JEC_" + year ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight );
                }
            }

            //fill JEC up histograms
            if( passSelection( event, "JECUp" ) ){
                auto fillValues = buildFillingVector( event, "JECUp" );
                for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
                    histogram::fillValue( histogramsUncUp[ "JEC_" + year ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight );
                }
            }

            //fill JER down histograms
            if( passSelection( event, "JERDown" ) ){
                auto fillValues = buildFillingVector( event, "JERDown" );
                for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
                    histogram::fillValue( histogramsUncDown[ "JER_" + year ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight );
                }
            }

            //fill JER up histograms
            if( passSelection( event, "JERUp" ) ){
                auto fillValues = buildFillingVector( event, "JERUp" );
                for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
                    histogram::fillValue( histogramsUncUp[ "JER_" + year ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight );
                }
            }

//    //            //fill unclustered down histograms
//    //            if( passSelection( event, "UnclDown" ) ){
//    //                auto fillValues = buildFillingVector( event, "UnclDown", massSplitting, nnReader );
//    //                for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
//    //                    histogram::fillValue( histogramsUncDown[ "uncl" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight );
//    //                }
//    //            }
//    //
//    //            //fill unclustered up histograms 
//    //            if( passSelection( event, "UnclUp" ) ){
//    //                auto fillValues = buildFillingVector( event, "UnclUp", massSplitting, nnReader );
//    //                for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
//    //                    histogram::fillValue( histogramsUncUp[ "uncl" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight );
//    //                }
//    //            }
            
            //apply nominal selection and compute nominal variables
            if( !passSelection( event, "nominal" ) ) continue;
            auto fillValues = buildFillingVector( event, "nominal" );

          //fill scale down histograms
          double weightScaleDown;
          if( sampleHasPdfAndScale ){
              weightScaleDown =  event.generatorInfo().relativeWeight_MuR_0p5_MuF_0p5();
          } else {
              weightScaleDown = 1.;
          }
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncDown[ "scale" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightScaleDown );
          }
      
          //fill scale up histograms
          double weightScaleUp;
          if( sampleHasPdfAndScale ){
              weightScaleUp = event.generatorInfo().relativeWeight_MuR_2_MuF_2();
          } else {
              weightScaleUp = 1.;
          }
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncUp[ "scale" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightScaleUp );
          }

          //fill pdf histograms
          for( unsigned pdf_i = 0; pdf_i < numberOfPdfVariations; ++pdf_i){
              double weightPdf = sampleHasPdfAndScale ? event.generatorInfo().relativeWeightPdfVar(pdf_i) : 1.;
              for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
                  //std::cout << " pdf_i, dist, fillIndex " << pdf_i << ", " << dist << ", " << fillIndex << std::endl;
                  histogram::fillValue( histogramsPDFVars[ pdf_i ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightPdf );
              }
          }

          //fill pileup down histograms
          double weightPileupDown = reweighter[ "pileup" ]->weightDown( event ) / reweighter[ "pileup" ]->weight( event );
          if ( std::isnan(weightPileupDown) ) weightPileupDown = 0.;
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncDown[ "pileup" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightPileupDown );
          }

          //fill pileup up histograms
          double weightPileupUp = reweighter[ "pileup" ]->weightUp( event ) / reweighter[ "pileup" ]->weight( event );
          if ( std::isnan(weightPileupUp) ) weightPileupUp = 0.;
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncUp[ "pileup" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightPileupUp );
          }

          //fill b-tag down histograms
          //WARNING : THESE SHOULD ACTUALLY BE SPLIT BETWEEN HEAVY AND LIGHT FLAVORS
          double weightBTagHeavyDown = reweighter[ "bTag_heavy" ]->weightDown( event ) / reweighter[ "bTag_heavy" ]->weight( event );
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncDown[ "bTag_heavy_" + year ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightBTagHeavyDown );
          }

          //fill b-tag up histograms
          //WARNING : THESE SHOULD ACTUALLY BE SPLIT BETWEEN HEAVY AND LIGHT FLAVORS
          double weightBTagHeavyUp = reweighter[ "bTag_heavy" ]->weightUp( event ) / reweighter[ "bTag_heavy" ]->weight( event );
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncUp[ "bTag_heavy_" + year ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightBTagHeavyUp );
          }

          //WARNING : THESE SHOULD ACTUALLY BE SPLIT BETWEEN HEAVY AND LIGHT FLAVORS
          double weightBTagLightDown = reweighter[ "bTag_light" ]->weightDown( event ) / reweighter[ "bTag_light" ]->weight( event );
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncDown[ "bTag_light_" + year ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightBTagLightDown );
          }

          //fill b-tag up histograms
          //WARNING : THESE SHOULD ACTUALLY BE SPLIT BETWEEN HEAVY AND LIGHT FLAVORS
          double weightBTagLightUp = reweighter[ "bTag_light" ]->weightUp( event ) / reweighter[ "bTag_light" ]->weight( event );
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncUp[ "bTag_light_" + year ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightBTagLightUp );
          }

        //fill prefiring down histograms
          double weightPrefireDown = reweighter[ "prefire" ]->weightDown( event ) / reweighter[ "prefire" ]->weight( event );
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncDown[ "prefire" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightPrefireDown );
          }
      
          //fill prefiring up histograms
          double weightPrefireUp = reweighter[ "prefire" ]->weightUp( event ) / reweighter[ "prefire" ]->weight( event );
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncUp[ "prefire" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * weightPrefireUp );
          }

          double recoWeightDown;
          double recoWeightUp;
//          if( !event.is2018() ){
              recoWeightDown = reweighter[ "electronReco_pTBelow20" ]->weightDown( event ) * reweighter[ "electronReco_pTAbove20" ]->weightDown( event ) / ( reweighter[ "electronReco_pTBelow20" ]->weight( event ) * reweighter[ "electronReco_pTAbove20" ]->weight( event ) );
              recoWeightUp = reweighter[ "electronReco_pTBelow20" ]->weightUp( event ) * reweighter[ "electronReco_pTAbove20" ]->weightUp( event ) / ( reweighter[ "electronReco_pTBelow20" ]->weight( event ) * reweighter[ "electronReco_pTAbove20" ]->weight( event ) );
//          } else {
//              recoWeightDown = reweighter[ "electronReco" ]->weightDown( event ) / ( reweighter[ "electronReco" ]->weight( event ) );
//              recoWeightUp = reweighter[ "electronReco" ]->weightUp( event ) / ( reweighter[ "electronReco" ]->weight( event ) );
//          }

          //fill lepton reco down histograms 
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncDown[ "lepton_reco" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * recoWeightDown );
          }

          //fill lepton reco up histograms 
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncUp[ "lepton_reco" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * recoWeightUp );
          }

          double leptonIDWeightDown = reweighter[ "muonID" ]->weightDown( event ) * reweighter[ "electronID" ]->weightDown( event ) / ( reweighter[ "muonID" ]->weight( event ) * reweighter[ "electronID" ]->weight( event ) );
          double leptonIDWeightUp = reweighter[ "muonID" ]->weightUp( event ) * reweighter[ "electronID" ]->weightUp( event ) / ( reweighter[ "muonID" ]->weight( event ) * reweighter[ "electronID" ]->weight( event ) );

          //fill lepton id down histograms
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncDown[ "lepton_id" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * leptonIDWeightDown );
          }

          //fill lepton id up histograms
          for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
              histogram::fillValue( histogramsUncUp[ "lepton_id" ][ dist ][ fillIndex ].get(), fillValues[ dist ], weight * leptonIDWeightUp );
          }

        }
    }
      
    //make SampleCrossSectionRatio objects to remove cross section effects from theory uncertainties
    std::map< std::string, SampleCrossSections > sampleCrossSectionsMap;
    for( size_t p = 1; p < sampleVec.size(); ++p ){
        sampleCrossSectionsMap[ sampleVec[ p ].uniqueName() ] = SampleCrossSections( sampleVec[p] );
    }

    //compute pdf uncertainty histograms as RMS of all pdf variations
    for( size_t p = 1; p < sampleVec.size(); ++p ){
        //track pdf cross section variations to divide them out of the shape uncertainty
        SampleCrossSections* crossSectionPtr = nullptr;
        if( p < sampleVec.size() ){
            crossSectionPtr = &sampleCrossSectionsMap[ sampleVec[ p ].uniqueName() ];
        }
        for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
            for( int bin = 1; bin < histogramsUncDown[ "pdf" ][ dist ][ p ]->GetNbinsX() + 1; ++bin ){
                double varRMS = 0.;
                double originalBin = histograms[ dist ][ p ]->GetBinContent( bin );
                for( size_t pdf = 0; pdf < numberOfPdfVariations; ++pdf ){
                    double variedBin = histogramsPDFVars[ pdf ][ dist ][ p ]->GetBinContent( bin );
                    
                    //divide out cross section effects
                    if( p < sampleVec.size() ){
                        if( crossSectionPtr->numberOfLheVariations() >= 110 ){
                            variedBin /= crossSectionPtr->crossSectionRatio_pdfVar( p );
                        }
                    }
                    double diff = ( variedBin - originalBin );
                    varRMS += ( diff * diff );
                }
                varRMS = std::sqrt( ( 1./static_cast< double >( numberOfPdfVariations ) ) * varRMS );
                histogramsUncDown[ "pdf" ][ dist ][ p ]->SetBinContent( bin, std::max( originalBin - varRMS, 0. ) );
                histogramsUncUp[ "pdf" ][ dist ][ p ]->SetBinContent( bin, std::max( originalBin + varRMS, 0. ) );
            }
        }
    }

    //set negative contributions to zero
    for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
        
        //backgrounds, data and merged signal
        for( size_t p = 0; p < sampleVec.size() + 1; ++p ){
            analysisTools::setNegativeBinsToZero( histograms[ dist ][ p ] );

            for( const auto& unc : shapeUncNames ){
                analysisTools::setNegativeBinsToZero( histogramsUncDown[ unc ][ dist ][ p ] );
                analysisTools::setNegativeBinsToZero( histogramsUncUp[ unc ][ dist ][ p ] );
            }
        }
    }
    
    //divide out cross section ratios from scale uncertainty ( can not be done for nonprompt )
    for( size_t p = 1; p < sampleVec.size(); ++p ){
        for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
            double xSecRatioScaleDown;
            double xSecRatioScaleUp;
            try {
                xSecRatioScaleDown = sampleCrossSectionsMap[ sampleVec[p].uniqueName() ].crossSectionRatio_MuR_0p5_MuF_0p5();
                xSecRatioScaleUp = sampleCrossSectionsMap[ sampleVec[p].uniqueName() ].crossSectionRatio_MuR_2_MuF_2();
            } catch( std::out_of_range& ){
                xSecRatioScaleDown = 1.;
                xSecRatioScaleUp = 1.;
            }
            histogramsUncDown[ "scale" ][ dist ][ p ]->Scale( 1./xSecRatioScaleDown );
            histogramsUncUp[ "scale" ][ dist ][ p ]->Scale( 1./xSecRatioScaleUp );
        }
      }



    //merge process histograms
//    std::vector< std::string > proc = {"Data", "ttZ", "WZ", "Xgamma", "ZZ", "Nonprompt",  };
    std::vector< std::string > proc = {"Data", "ttZ", "ttX", "WZ", "Xgamma", "ZZ", "rare", "Nonprompt" };
    std::vector< std::vector< TH1D* > > mergedHistograms( histInfoVector.size(), std::vector< TH1D* >( proc.size() ) );
    //size_t numberOfBackgrounds = 0;
    //for( const auto& s : sampleVec ){
    //    if( !s.isNewPhysicsSignal() ) ++numberOfBackgrounds; 
    //}

    // loop over all distributions
    for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
        // loop over all samples
        for( size_t m = 0, sample = 0; m < proc.size() - 1; ++m ){
            mergedHistograms[ dist ][ m ] = dynamic_cast< TH1D* >( histograms[ dist ][ sample ]->Clone() );
            //while( sample < numberOfBackgrounds - 1 && sampleVec[ sample ].processName() == sampleVec[ sample + 1 ].processName() ){
            while( sample < sampleVec.size() - 1 && sampleVec[ sample ].processName() == sampleVec[ sample + 1 ].processName() ){
                mergedHistograms[ dist ][ m ]->Add( histograms[ dist ][ sample + 1].get() );
                ++sample;
            }
            ++sample;
        }
        
        //add nonprompt histogram
        mergedHistograms[ dist ][ proc.size() -1 ] = dynamic_cast< TH1D* >( histograms[ dist ].back()->Clone() );
    }

    //merge process histograms for uncertainties 
    std::map< std::string, std::vector< std::vector< TH1D* > > > mergedHistogramsUncDown;
    std::map< std::string, std::vector< std::vector< TH1D* > > > mergedHistogramsUncUp;

    for( const auto& unc : shapeUncNames ){
        mergedHistogramsUncDown[ unc ] = std::vector< std::vector< TH1D* > >( histInfoVector.size(), std::vector< TH1D* >( proc.size() ) );
        mergedHistogramsUncUp[ unc ] = std::vector< std::vector< TH1D* > >( histInfoVector.size(), std::vector< TH1D* >( proc.size() ) );
        for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
            for( size_t m = 0, sample = 0; m < proc.size() - 1; ++m ){
                mergedHistogramsUncDown[ unc ][ dist ][ m ] = dynamic_cast< TH1D* >( histogramsUncDown[ unc ][ dist ][ sample ]->Clone() );
                mergedHistogramsUncUp[ unc ][ dist ][ m ] = dynamic_cast< TH1D* >( histogramsUncUp[ unc ][ dist ][ sample ]->Clone() );
                //while( sample < numberOfBackgrounds - 1 && sampleVec[ sample ].processName() == sampleVec[ sample + 1 ].processName() ){
                while( sample < sampleVec.size() - 1 && sampleVec[ sample ].processName() == sampleVec[ sample + 1 ].processName() ){
                    mergedHistogramsUncDown[ unc ][ dist ][ m ]->Add( histogramsUncDown[ unc ][ dist ][ sample + 1].get() );
                    mergedHistogramsUncUp[ unc ][ dist ][ m ]->Add( histogramsUncUp[ unc ][ dist ][ sample + 1].get() );
                    ++sample;
                }
                ++sample;
            }
            mergedHistogramsUncDown[ unc ][ dist ][ proc.size() - 1 ] = dynamic_cast< TH1D * >( histogramsUncDown[ unc ][ dist ].back()->Clone() );
            mergedHistogramsUncUp[ unc ][ dist ][ proc.size() - 1 ] = dynamic_cast< TH1D * >( histogramsUncUp[ unc ][ dist ].back()->Clone() );
        }
    }

    //make total uncertainty histograms for plotting
    const std::vector< std::string > uncorrelatedBetweenProcesses = {"scale", "pdf", "scaleXsec", "pdfXsec"};
    double lumiUncertainty = 1.025;
    std::vector<double> flatUnc = { lumiUncertainty, 1.02 }; //lumi, trigger
//    std::map< std::string, double > backgroundSpecificUnc;
//    if( controlRegion == "TTZ" || controlRegion == "NP" ){
//        backgroundSpecificUnc = {
//            {"Nonprompt", 1.3},
//            {"WZ", 1.1},
//            {"X + #gamma", 1.1},
//            {"ZZ/H", 1.1},
//            {"t#bar{t}/t + X", 1.15 },
//            {"Multiboson", 1.5}
//        };
//    } else {
//        backgroundSpecificUnc = {
//            {"Nonprompt", 1.3},
//            {"WZ", 1.1},
//            {"X + #gamma", 1.1},
//            {"ZZ/H", 1.1},
//            {"t#bar{t}/t + X", 1.5 },
//            {"Multiboson", 1.5}
//        };
//    }

    //const std::set< std::string > acceptedShapes = { "JEC_" + year, "JER_" + year, "uncl", "scale", "pileup", "bTag_" + year, "prefire", "lepton_reco", "lepton_id"}
    //const std::set< std::string > acceptedShapes = { "JEC_" + year, "scale", "bTag_" + year, "prefire", "lepton_reco", "lepton_id"};
    //const std::set< std::string > acceptedShapes = { "JEC_" + year };

    std::vector< TH1D* > totalSystUncertainties( histInfoVector.size() );
    for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
        totalSystUncertainties[ dist ] = dynamic_cast< TH1D* >( mergedHistograms[ dist ][ 0 ]->Clone() );
        for( int bin = 1; bin < totalSystUncertainties[ dist ]->GetNbinsX() + 1; ++bin ){
            double binUnc = 0;

            //add shape uncertainties
            for( auto& shape : shapeUncNames ){
                //if( acceptedShapes.find( shape ) == acceptedShapes.cend() ) continue;
                bool nuisanceIsUncorrelated = ( std::find( uncorrelatedBetweenProcesses.cbegin(), uncorrelatedBetweenProcesses.cend(), shape ) != uncorrelatedBetweenProcesses.cend() );

                //correlated case : linearly add up and down variations
                double varDown = 0.;
                double varUp = 0.;

                //uncorrelated case : quadratically add the maximum of the up and down variations
                double var = 0.;

                for( size_t p = 1; p < proc.size(); ++p ){
                    double nominalContent = mergedHistograms[ dist ][ p ]->GetBinContent( bin );
                    double downVariedContent = mergedHistogramsUncDown[ shape ][ dist ][ p ]->GetBinContent( bin );
                    double upVariedContent = mergedHistogramsUncUp[ shape ][ dist ][ p ]->GetBinContent( bin );
                    double down = fabs( downVariedContent - nominalContent );
                    double up = fabs( upVariedContent - nominalContent );
                    
                    //uncorrelated case : 
                    if( nuisanceIsUncorrelated ){
                        double variation = std::max( down, up );
                        var += variation*variation;
                    
                    //correlated case :     
                    } else {
                        varDown += down;
                        varUp += up;
                    }

                }
                //correlated case : 
                if( !nuisanceIsUncorrelated ){
                    var = std::max( varDown, varUp );
                    var = var*var;
                }
                
                //add (already quadratic) uncertainties 
                binUnc += var;

                //uncorrelated case : 
                if( nuisanceIsUncorrelated ){
                    histogramsUncDecompDown[ dist ][ shapesIntMap.at(shape) ]->SetBinContent( bin, sqrt( var ) * -1. );
                    histogramsUncDecompUp[   dist ][ shapesIntMap.at(shape) ]->SetBinContent(   bin, sqrt( var ) );
                    //correlated case :     
                } else {
                    histogramsUncDecompDown[ dist ][ shapesIntMap.at(shape) ]->SetBinContent( bin, varDown * -1. );
                    histogramsUncDecompUp[   dist ][ shapesIntMap.at(shape) ]->SetBinContent(   bin, varUp );
                }
            }

            //add general flat uncertainties (considered correlated among all processes)
            for( double unc : flatUnc ){
                double var = 0;
                for( size_t p = 1; p < proc.size(); ++p ){
                    if( proc[p] == "Nonprompt" ){
                        continue;
                    }
                    double binContent = mergedHistograms[ dist ][ p ]->GetBinContent( bin );
                    double variation = binContent*(unc - 1.);
                    var += variation;
                }
                binUnc += var*var;
            }

//            //add background specific uncertainties (uncorrelated between processes)
//            for(auto& uncPair : backgroundSpecificUnc){
//                for( size_t p = 1; p < proc.size(); ++p ){
//                    if(proc[p] == uncPair.first){
//                        double var = mergedHistograms[ dist ][ p ]->GetBinContent( bin )*( uncPair.second - 1. );
//                        binUnc += var*var;
//                    }
//                }
//            }

            //square root of quadratic sum is total uncertainty
            totalSystUncertainties[ dist ]->SetBinContent( bin, sqrt( binUnc ) );
        }
    }
    
    //make plots
    for( size_t dist = 0; dist < histInfoVector.size(); ++dist ){
        std::string header;
        if( year == "2016" ){
            header = "35.9 fb^{-1} (13 TeV)";
        } else if( year == "2017" ){
            header = "41.5 fb^{-1} (13 TeV)";
        } else {
            header = "59.7 fb^{-1} (13 TeV)";
        }
        std::string directoryName;
        std::string plotNameAddition;

            directoryName = stringTools::formatDirectoryName( "plots/" + year + "/" + procName + "/" + controlRegion );
            plotNameAddition = "_" + procName + "_" + controlRegion + "_" + year;
//            directoryName = stringTools::formatDirectoryName( "plots/ttZ/" + year + "/" );
//            plotNameAddition = "_" + year;

        systemTools::makeDirectory( directoryName );
        plotDataVSMC( mergedHistograms[dist][0], &mergedHistograms[dist][1], &proc[0], proc.size() - 1, directoryName + histInfoVector[ dist ].name() + plotNameAddition + ".pdf" , "ttZ", false, false, header, totalSystUncertainties[ dist ], nullptr );
        plotDataVSMC( mergedHistograms[dist][0], &mergedHistograms[dist][1], &proc[0], proc.size() - 1, directoryName + histInfoVector[ dist ].name() + plotNameAddition + "_log.pdf" , "ttZ", true, false, header, totalSystUncertainties[ dist ], nullptr );

        plotUncAll( &mergedHistograms[dist][1], proc.size() - 1, histogramsUncDecompUp[dist], histogramsUncDecompDown[dist], &shapeUncNames[0], shapeUncNames.size(), directoryName + histInfoVector[ dist ].name() + plotNameAddition + "_unc.pdf", 1.5 );
//        plotUncAll( &mergedHistograms[dist][1], proc.size() - 1, histogramsUncDecompUp[dist], histogramsUncDecompDown[dist], &shapeUncNames[0], shapeUncNames.size(), directoryName + histInfoVector[ dist ].name() + plotNameAddition + "_uncZoom.pdf", 0.5 );
    }

    // prepare root file for datacards / combine
    std::string shape_name = "shapeFile_";
    shape_name += controlRegion + "_" + year + ".root";
    //shape_name += ".root";
    TFile *file = TFile::Open( shape_name.c_str(), "RECREATE");
    std::vector< std::string > proc_names = {"data_obs", "ttZ", "ttX", "WZ", "Xgamma", "ZZ", "rare", "Nonprompt",  };
    TH1D *hist, *histUp, *histDown;
    // loop over all sample types considered.
    for( size_t p = 0; p < proc.size() ; ++p ){
        const char * name = proc_names[p].c_str();
        hist = (TH1D*)(*mergedHistograms[ 13 ][ p ]).Clone( name );
        hist->SetName( name );
        hist->Write();
        for( const auto& unc : shapeUncNames ){
            std::string nameUp = proc_names[p]+"_"+unc+"Up";

            histUp = (TH1D*)(*mergedHistogramsUncUp[ unc ][ 13 ][ p ]).Clone( nameUp.c_str() );
            histUp->SetName( nameUp.c_str() );
            histUp->Write();

            std::string nameDown = proc_names[p]+"_"+unc+"Down";

            histDown = (TH1D*)(*mergedHistogramsUncDown[ unc ][ 13 ][ p ]).Clone( nameDown.c_str() );
            histDown->SetName( nameDown.c_str() );
            histDown->Write();
        }
    }
    file->Close();
}



int main( int argc, char* argv[] ){
    setTDRStyle();
//    const std::string sampleDirectoryPath = "/pnfs/iihe/cms/store/user/wverbeke/ntuples_ewkino/";
//    const std::string sampleDirectoryPath = "/pnfs/iihe/cms/store/user/mniedzie/old_ntuples/ntuples_ttV_2017/";
    const std::string sampleDirectoryPath = "/user/mniedzie/Work/ntuples_ttz_new/";
    std::vector< std::string > argvStr( &argv[0], &argv[0] + argc );
    
    //run specific model and mass splitting and year
    // std::string model = argvStr[1];
    // std::string deltaM = argvStr[2];

    // create lists of allowed argument values. std::set would allow for faster search, 
    // but it's irrelevant here and vectors will be more convenient to handle.
    // std::set< std::string > years; years.insert("2016"); years.insert("2017"); years.insert("2018");
    std::vector< std::string > years{"2016", "2017", "2018"};
    std::vector< std::string > cRegions{"ttZ", "ttZclean", "WZ", "DY", "ttbar", "XGamma", "ZZ", "TTZ"};
    std::vector< std::string > processes{"data", "ttZ", "ttX", "WZ", "Xgam", "ZZ", "rare", "all"};

    // if not enough arguments, complain. Last argument is optional
    if(argc < 3){
        std::cerr << "please specify input arguments, usage:" << std::endl;
        std::cerr << "./ttZAnalysis year region (process)" << std::endl;
        return 1;
    }    
    std::string year = argvStr[1];
    std::string controlRegion = argvStr[2];
    std::string procName;
    if ( argc == 4 ) procName = argvStr[3];
    else procName = "all";

    if ( std::find(std::begin(years),     std::end(years),     year)          == std::end(years) || 
         std::find(std::begin(cRegions),  std::end(cRegions),  controlRegion) == std::end(cRegions) || 
         std::find(std::begin(processes), std::end(processes), procName)      == std::end(processes) 
    ){
        std::cerr << "At least one of arguments not recognized" << std::endl;
        std::cerr << "allowed arguments:" << std::endl;
        std::cerr << "year:                          2016, 2017, 2018"                      << std::endl;
        std::cerr << "control region:                ttZ, ttZclean, WZ, DY, ttbar, XGamma, ZZ, TTZ" << std::endl;
        std::cerr << "process to analyze (optional): data, ttZ, ttX, WZ, Xgam, ZZ, rare"    << std::endl;
        return 1;
    }
    
    analyze( year, controlRegion, sampleDirectoryPath, procName );

    return 0;
}
