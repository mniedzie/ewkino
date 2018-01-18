/*
Class containing data and background histograms for plotting
*/

//include c++ library classes
#include <map>
#include <string>

//include root classes


#ifndef Plot_H
#define Plot_H
class Plot{
    public:
        Plot(const std::string& fileName, std::shared_ptr<TH1D> obs, std::map< std::string, std::shared_ptr < TH1D > > back, std::map< std::string, 
            std::shared_ptr< TH1D > > unc = std::map< std::string, std::shared_ptr< TH1D > >(), std::map< std::string, std::shared_ptr< TH1D > > sig = std::map< std::string, std::shared_ptr< TH1D > >() ):
            data(obs), bkg(back), syst(unc), signal(sig) {}
        void draw(const std::string& analysis = "", bool log = false, bool normToData = false, const std::string& header = "", TH1D** bkgSyst = nullptr, const bool* isSMSignal = nullptr, const bool sigNorm = true);
    private:
        std::shared_ptr<TH1D> data;
        std::map< std::string, std::shared_ptr< TH1D > > bkg;
        std::map< std::string, std::shared_ptr< TH1D > > syst; 
        std::map< std::string, std::shared_ptr< TH1D > > signal;
        std::string fileName;
};
#endif
