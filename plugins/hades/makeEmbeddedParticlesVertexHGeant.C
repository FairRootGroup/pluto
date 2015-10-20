//TITLE <UserClass> Embedded particles with vertex info and seq-nr for HADES

//Author: J. Markert

void makeEmbeddedParticlesVertexHGeant(
                           TString vertex_ntuple = "be05278234631_dst_vertex_real.root",
                           TString outFile = "be05278234631_pluto_embedded_vertex",
			   TString outDir  = "./"
			  )
{   // This macro creates n pluto particle per sector (60 deg) with white disribution
    // in theta , phi and momentum. The particles will be created
    // at the vertex given my the input ntuple extracted from REAL data.


    //---------------CONIGURATION----------------------
    //
    Int_t asciiOut     = 0;    // write pluto ascci output for HGeant (==0 if we use HGeantOutput)
    Int_t rootOut      = 0;    // write pluto root file
    Int_t calcVertex   = 1;    // write the vertex to the ascii output for HGeant
    PChannel* channel  = 0;    // no channel needed here
    Int_t    nChannel  = 0;    // number of channels will be 0
    TString myParticle = "e+"; // particle type used for embedding
    Bool_t  writeSeqNumber = kTRUE;  // write eventSeqNumber in addition
    Bool_t  writeIndex     = kTRUE;  // write parentIndex in addition

    Float_t pmin               = 0.050;   // minimum mom  [GeV/c]
    Float_t pmax               = 1.000;   // maximum mom  [GeV/c]
    Float_t thetamin           =    10;   // minimum theta angle  [deg]
    Float_t thetamax           =    90;   // maximum theta angle  [deg]
    Float_t phimin             =     0;   // minimum phi angle    [deg]
    Float_t phimax             =    60;   // maximum phi angle    [deg]
    Float_t phistart           =     0;   // start phi angle (+ sector * 60. deg will ber rotated later)
    Int_t numParticlePerSector =     1;   // number of Particles embedded per sector

    //-------------------------------------------------


    //-------------------------------------------------
    // the output file name will be constructed
    // outFile and OutDir (.hld or .root will be stripped)
    TString filename = outFile;
    filename.ReplaceAll(".root","");
    filename.ReplaceAll(".hld" ,"");
    filename = outDir + "/" + filename;

    cout<<"Output file name : "<<filename.Data()<<endl;
    //-------------------------------------------------



    PReaction my_reaction(channel,filename.Data(),nChannel,rootOut,0,calcVertex,asciiOut);
    my_reaction.SetWriteIndex(writeIndex);

    //Construct the vertex container:
    PVertexFile *vertex = new PVertexFile();
    vertex->OpenFile(vertex_ntuple);
    //add to prologue action
    my_reaction.AddPrologueBulk(vertex);


    PHGeantOutput* output = new PHGeantOutput();
    output->SetWriteSeqNumber(writeSeqNumber);
    output->OpenFile(Form("%s.evt",filename.Data()));
    my_reaction.AddFileOutput(output);

    //Construct the embedded container:
    PEmbeddedParticles* embedded = new PEmbeddedParticles();

    //Add an e+ which we emit at a single point:
    PParticle* particle = new PParticle(myParticle.Data(),1.,2.,3.);
    //Just add the particle to the container:
    embedded->AddParticle(particle);

    // added particle will be cloned with
    // n times per sector
    if (1) {
        // nparticle / sector
	embedded->SetSamplingSector(pmin    , pmax,                   // pmin,pmax
				    thetamin, thetamax,               // theta_min, theta_max
				    phimin  , phimax,                 // phi_min,phi_max
				    phistart, numParticlePerSector ); // phiStart, numParticle/sector
    } else {
        // nparticle / full acceptance
	embedded->SetSamplingSector(pmin    , pmax,                 // pmin,pmax
				    thetamin, thetamax,             // theta_min, theta_max
				    -180.    , 180.,                // phi_min,phi_max
				    phistart, numParticlePerSector, // phiStart, numParticle/sector
				    360.    , 1);                   // Delta-phi, nSec
	
    }




    //Add our container to the reaction:
    my_reaction.AddBulk(embedded);

    // number of events
    my_reaction.Loop(10000);
    output->CloseFile();

}
