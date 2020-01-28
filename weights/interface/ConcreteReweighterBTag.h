#ifndef ConcreteReweighterBTag_H
#define ConcreteReweighterBTag_H

#include "ReweighterBTagHeavyFlavor.h"
#include "ReweighterBTagLightFlavor.h"


class ReweighterBTagHeavyFlavorDeepCSV : public ReweighterBTagHeavyFlavor {

    using ReweighterBTagHeavyFlavor::ReweighterBTagHeavyFlavor;

    private:
        virtual double CSVValue( const Jet& jet ) const override{ return jet.deepCSV(); }
};



class ReweighterBTagLightFlavorDeepCSV : public ReweighterBTagLightFlavor {
    
    using ReweighterBTagLightFlavor::ReweighterBTagLightFlavor;
    
    private:
        virtual double CSVValue( const Jet& jet ) const override{ return jet.deepCSV(); }
};



class ReweighterBTagHeavyFlavorDeepFlavor : public ReweighterBTagHeavyFlavor {

    using ReweighterBTagHeavyFlavor::ReweighterBTagHeavyFlavor;
    
    private:
        virtual double CSVValue( const Jet& jet ) const override{ return jet.deepFlavor(); }
};



class ReweighterBTagLightFlavorDeepFlavor : public ReweighterBTagLightFlavor {
    
    using ReweighterBTagLightFlavor::ReweighterBTagLightFlavor;

    private:
        virtual double CSVValue( const Jet& jet ) const override{ return jet.deepFlavor(); }
};

#endif
