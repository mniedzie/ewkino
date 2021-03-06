#include "../interface/ElectronCollection.h"


ElectronCollection::ElectronCollection( const TreeReader& treeReader ){
    for( unsigned e = treeReader._nMu; e < treeReader._nLight; ++ e){
        push_back( Electron( treeReader, e ) );
    }
}
