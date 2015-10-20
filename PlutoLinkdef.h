// Pragma links for Pluto 19.8.03
#ifdef __CINT__
// General
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
// Utilities
#pragma link C++ class PUtils;
#pragma link C++ class PUtilsREngine;
#pragma link C++ class PBatch;
#pragma link C++ function makeGlobalBatch();
#pragma link C++ class PStaticData;
#pragma link C++ class PDynamicData;
#pragma link C++ function makeDataBase();
#pragma link C++ function makeStaticData();
#pragma link C++ function makeDynamicData();
#pragma link C++ class PMesh;
#pragma link C++ class PAdaptiveMesh;
#pragma link C++ class PAdaptiveMeshN;
#pragma link C++ class PValues;
#pragma link C++ class PArray;
#pragma link C++ class PSplash;
#pragma link C++ class PCommandList;
// Particle Data
#pragma link C++ class PData;
#pragma link C++ class PStdData;
#pragma link C++ function listParticle(int);
#pragma link C++ function listParticle(char *);
#pragma link C++ function listModes(int);
#pragma link C++ function listModes(char *);
#pragma link C++ function makeStdData();
#pragma link C++ class PDataBase;
// Cross sections
#pragma link C++ class PCross;
// Particle
#pragma link C++ class PParticle;
// Fireball
#pragma link C++ class PFireball;
// Dilepton
#pragma link C++ class PDiLepton;
// File IO interface
#pragma link C++ class PFileInput;
#pragma link C++ class PFileOutput;
#pragma link C++ class PHGeantOutput;
#pragma link C++ class PDebugFileOutput;

// Channel
#pragma link C++ class PChannel;
// SAID
#pragma link C++ class PSaid;
// Reaction
#pragma link C++ class PReaction;
// Filter
#pragma link C++ class PFilter;
// Thermal source
#pragma link C++ class PThermal;
// Particle decay channel
#pragma link C++ class PDecayManager;
// Particle decay list
#pragma link C++ class PDecayChannel;
// Bulk decays
#pragma link C++ class PBulkInterface;
#pragma link C++ class PPlutoBulkDecay;
#pragma link C++ class PPythiaBulkDecay;
#pragma link C++ class PEmbeddedParticles;
#pragma link C++ class PVertexFile;
#pragma link C++ class PProjector;

// Classes for distribution interface
#pragma link C++ class PDistributionManager;
#pragma link C++ class PDistributionManagerUtil;
#pragma link C++ class PDistributionCollection;

//simple distributions:
#pragma link C++ class PDistribution;
#pragma link C++ class PAngularDistribution;
#pragma link C++ class PDeltaAngularDistribution;
#pragma link C++ class PPiOmegaAngularDistribution;
#pragma link C++ class PDalitzDistribution;
#pragma link C++ class PScatterDistribution;
#pragma link C++ class PFermiMomentum;
#pragma link C++ class PFermiMomentumDD;
#pragma link C++ class PBeamSmearing;


// Channel models
#pragma link C++ class PStdModels;
#pragma link C++ class PChannelModel;
#pragma link C++ class PHadronModel;
#pragma link C++ class PBreitWigner;
#pragma link C++ class PComplexBreitWigner;
#pragma link C++ class PMassSampling;
#pragma link C++ class PDalitzDecay;
#pragma link C++ class PHadronDecay;
#pragma link C++ class PHadronDecayM1;
#pragma link C++ class PHadronDecayM1N;
#pragma link C++ class PHadronDecayM2;
#pragma link C++ class PHadronDecayM3;
#pragma link C++ class PInclusiveModel;
#pragma link C++ class PEEDirectDecay;
#pragma link C++ class PFixedDecay;
#pragma link C++ class PFixedProduction;
#pragma link C++ class PGenBod;
#pragma link C++ class PNNFSI;

#pragma link C++ class PTCrossWeight;

#pragma link C++ class PSimpleVMDFF;
#pragma link C++ class PPropagator;

#pragma link C++ class PFunction;
#pragma link C++ class PFormula;
#pragma link C++ class PF2;



// HADES-specific functions
//#pragma link C++ class PHADES;
//#pragma link C++ function PHADES::smear(PParticle*,Double_t);
#pragma link C++ function PUtils::FindIndex(Int_t, Double_t*, Double_t);
#pragma link C++ function PUtils::Tokenize(char*, char*, char **, int *);
#pragma link C++ function PUtils::remove_brackets(char **, char , char );
#pragma link C++ function makeDistributionManager();
#pragma link C++ function makePUtilsREngine();

// Track
//#pragma link C++ class PTrack;
// TrackH
//#pragma link C++ class PTrackH;
//gSaid pointer
#pragma link C++ global gSaid;
#pragma link C++ global gDM;

//#pragma link C++ function select;

#include "PluginLinkdef.h"

#endif








