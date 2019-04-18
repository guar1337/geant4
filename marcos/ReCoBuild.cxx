#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "stdio.h"
#include "Riostream.h"
#include "ReCoBuild.hh"
#include <iostream>
#include "TLorentzVector.h"

// voron1392

void ReCoBuild() 
//TFile *inF  = new TFile(TString::Format("he6_000%d.root", iii), "READ");
{
gSystem->Load("libPhysics.so");
TFile *inF  = new TFile("gurney.root", "READ");
TFile *outF = new TFile("ReCoGu.root","recreate");
TTree *inTree = (TTree*)inF->Get("simevents");
TTree *outTree= new TTree("Reconstructed","reconstructed");
outTree->SetMakeClass(1);

inTree->SetBranchAddress("vBEAM.",  &vBEAM);
inTree->SetBranchAddress("vDeutCM.",  &vDeutCM);
inTree->SetBranchAddress("v2H.",    &v2H);
inTree->SetBranchAddress("v6He.",   &v6He);



inTree->Print();

//Creating addresses of BEAM holding branches
inTree->SetBranchAddress("Tbeam",	    &Tbeam);
inTree->SetBranchAddress("thetaCM",	    &thetaCM);
inTree->SetBranchAddress("Xpos",	    &inXpos);
inTree->SetBranchAddress("Ypos",	    &inYpos);
inTree->SetBranchAddress("Zpos",	    &inZpos);

inTree->SetBranchAddress("X6He",	    &X6He);
inTree->SetBranchAddress("Y6He",	    &Y6He);
inTree->SetBranchAddress("Z6He",	    &Z6He);

//Creating addresses of deuterium holding branches
inTree->SetBranchAddress("Tdeut",	   &Tdeut);
inTree->SetBranchAddress("CsIdeut",	 CsIdeut);
inTree->SetBranchAddress("SideutX",	 SideutX);
inTree->SetBranchAddress("SideutY",	 SideutY);
inTree->SetBranchAddress("labAng2H",	&labAng2H);

//Creating addresses of helium holding branches
inTree->SetBranchAddress("The",	     &The);
inTree->SetBranchAddress("CsIhe",	   CsIhe);
inTree->SetBranchAddress("SiheY",	   SiheY);
inTree->SetBranchAddress("SiheX",	   SiheX);
inTree->SetBranchAddress("labAngHe",	&labAngHe);


vecBEAM= new TLorentzVector();

//ReCo - detectors
outTree->Bronch("vecBEAM.",	    "TLorentzVector",   &vecBEAM);
outTree->Bronch("vec2H.",	    "TLorentzVector",   &vec2H);
outTree->Bronch("vec6He.",	    "TLorentzVector",   &vec6He);

outTree->Branch("Xpos",         &outXpos,	     "Xpos/D");
outTree->Branch("Ypos",	        &outYpos,	     "Ypos/D");
outTree->Branch("Zpos",	        &outZpos,	     "Zpos/D");

outTree->Branch("sq11mult",	    &sq11mult,	 "sq11mult/S");
outTree->Branch("sq21mult",	    &sq21mult,	 "sq21mult/S");
outTree->Branch("sq12mult",	    &sq12mult,	 "sq12mult/S");
outTree->Branch("sq22mult",	    &sq22mult,	 "sq22mult/S");
outTree->Branch("sq13mult",	    &sq13mult,	 "sq13mult/S");
outTree->Branch("sq23mult",	    &sq23mult,	 "sq23mult/S");

outTree->Branch("sq11Edep",	    sq11Edep,	  "sq11Edep[16]/D");
outTree->Branch("sq21Edep",	    sq21Edep,	  "sq21Edep[16]/D");
outTree->Branch("sq12Edep",	    sq12Edep,	  "sq12Edep[16]/D");
outTree->Branch("sq22Edep",	    sq22Edep,	  "sq22Edep[16]/D");
outTree->Branch("sq13Edep",	    sq13Edep,	  "sq13Edep[16]/D");
outTree->Branch("sq23Edep",     sq23Edep,   "sq23Edep[16]/D");

outTree->Branch("sq11strip",    sq11strip,  "sq11strip[16]/S");
outTree->Branch("sq21strip",    sq21strip,  "sq21strip[16]/S");
outTree->Branch("sq12strip",    sq12strip,  "sq12strip[16]/S");
outTree->Branch("sq22strip",    sq22strip,  "sq22strip[16]/S");
outTree->Branch("sq13strip",    sq13strip,  "sq13strip[16]/S");
outTree->Branch("sq23strip",    sq23strip,  "sq23strip[16]/S");

outTree->Branch("sq1ang",       &sq1ang,     "sq1ang/D");
outTree->Branch("sq1de",        &sq1de,      "sq1de/D");
outTree->Branch("sq1etot",      &sq1etot,    "sq1etot/D");
outTree->Branch("sq2ang",       &sq2ang,     "sq2ang/D");
outTree->Branch("sq2de",        &sq2de,      "sq2de/D");
outTree->Branch("sq2etot",      &sq2etot,    "sq2etot/D");


outTree->Branch("proba",      &proba,    "proba/D");


const double m2H=1875.6128;
const double m6He=5605.53497;
const double m4He=3727.38;
TVector3 beamDirection(0,0,1);

TLorentzVector *vec2H= new TLorentzVector(0,0,0,m2H);


//ReCo - straight from simulation
//outTree->Branch("labAngHe",	&relabAngHe,     "labAngHe/D");
//outTree->Branch("labAng2H",	&relabAng2H,     "labAng2H/D");
outTree->Branch("thetaCM",	 &rethetaCM,	"thetaCM/D");
outTree->Branch("Tbeam",	   &reTbeam,	  "Tbeam/D");




Long64_t nEntries = inTree->GetEntries();
for (Long64_t entry=0; entry<nEntries; entry++)
{
//TLorentzVector filling part
inTree->GetEntry(entry);
double Px6He=v6He->Px();
double Py6He=v6He->Py();
double Pz6He=v6He->Pz();
double Energy6He=v6He->Energy();
vec6He->SetPxPyPzE(Px6He,Py6He,Pz6He,Energy6He);

double Px2H=v2H->Px();
double Py2H=v2H->Py();
double Pz2H=v2H->Pz();
double Energy2H=v2H->Energy();
vec2H->SetPxPyPzE(Px2H,Py2H,Pz2H,Energy2H);

double PxBeam=vBEAM->Px();
double PyBeam=vBEAM->Py();
double PzBeam=vBEAM->Pz();
double EnergyBeam=vBEAM->Energy();
vecBEAM->SetPxPyPzE(PxBeam,PyBeam,PzBeam,EnergyBeam);
//End of TLorentzVector filling part

outXpos=inXpos;
outYpos=inYpos;
outZpos=inZpos;

sq11mult=0;
sq12mult=0;
sq13mult=0;
sq21mult=0;
sq22mult=0;
sq23mult=0;
//relabAngHe=180.0/3.14159*labAngHe;
//relabAng2H=180.0/3.14159*labAng2H;
rethetaCM=thetaCM;
reTbeam=Tbeam;
    //loop on every crystal
    for (int iii=0; iii<16; iii++)
    {    
    if (SideutX[iii]>0.1)
	  {
	  sq11Edep[sq11mult]=SideutX[iii];
	  sq11strip[sq11mult]=iii;
	  sq11mult++;
	  }
    if (SideutY[iii]>0.1)
	  {
	  sq12Edep[sq12mult]=SideutY[iii];
	  sq12strip[sq12mult]=iii;
	  sq12mult++;
	  }
    if (CsIdeut[iii]>0.5)
	  {
	  sq13Edep[sq13mult]=CsIdeut[iii];
	  sq13strip[sq13mult]=iii;     
	  sq13mult++;
	  }
    if (SiheX[iii]>0.1)
	  {
	  sq21Edep[sq21mult]=SiheX[iii];
	  sq21strip[sq21mult]=iii;
	  sq21mult++;
	  }
    if (SiheY[iii]>0.1)
	  {
	  sq22Edep[sq22mult]=SiheY[iii];
	  sq22strip[sq22mult]=iii;
	  sq22mult++;
	  }
    if (CsIhe[iii]>0.5)
	  {
	  sq23Edep[sq23mult]=CsIhe[iii];
	  sq23strip[sq23mult]=iii;
	  sq23mult++;
	  }
    }
    double Tbeam=v2H->Energy();
    double invariant = (m2H+m6He)*(m2H+m6He)+2.*m2H*Tbeam;
    double shorty=(invariant-m2H*m2H-m6He*m6He)*(invariant-m2H*m2H-m6He*m6He);
    double CMmom = sqrt((shorty-4.*m2H*m2H*m6He*m6He)/(4.*invariant));
    double chi = log((CMmom+sqrt(m2H*m2H+CMmom*CMmom))/m2H);
    
    
    if (sq11mult*sq12mult*sq13mult==1)
    {   
    sq1de=sq11Edep[0];
    sq1etot=sq11Edep[0]+sq13Edep[0];
    X2H=180.*sin(45./180.*3.14159)+(30.-4.*sq11strip[0])*cos(45./180.*3.14159);
    Y2H=+30.-4.* sq12strip[0];
    Z2H=180.*cos(45./180.*3.14159)-(30.-4.*sq11strip[0])*sin(45./180.*3.14159);        
    TVector3 vect2H(X2H-inXpos, Y2H-inYpos, Z2H+inXpos);
    sq1ang=vect2H.Angle(beamDirection)*180.0/3.1415927;
    

    //TLorentzVector vreco2H(vect2H,m2H);
    TLorentzVector vrecoCM(*vec2H+*vBEAM);
    //vrecoCM=vBEAM+vreco2H;
    TVector3 boostVECT =vrecoCM.BoostVector();
    //TLorentzVector labVect2H(vect2H,sq1etot+m2H);
    //vec2H->Boost(-boostVECT);
    
    double Ecm2H=vec2H->Energy();
    double newEcm=Ecm2H+sqrt(Ecm2H*Ecm2H+m6He*m6He-m2H*m2H);
    
    cout<<vec2H->Energy()-m2H<<" "<<vDeutCM->Energy()<<endl;
    proba=m6He+m2H-newEcm+(Ecm2H-m2H);
    
    
    
    
    
    
    
    
    
    
    }
    if(sq21mult*sq22mult*sq23mult==1)
    {
    sq2de=sq21Edep[0];
    sq2etot=sq21Edep[0]+sq23Edep[0];
    X6He=-180.*sin(15./180.*3.14159)+(30.-4.*sq21strip[0])*cos(15./180.*3.14159);
    Y6He=+30.-4.*sq22strip[0];
    Z6He=180.*cos(15./180.*3.14159)+(30.-4.*sq21strip[0])*sin(15./180.*3.14159);
    TVector3 vect6He(X6He-inXpos, Y6He-inYpos, Z6He+inXpos);
    sq2ang=vect6He.Angle(beamDirection)*180.0/3.1415927;   
    }
    if(sq21mult*sq22mult*sq23mult==1 || sq11mult*sq12mult*sq13mult==1)
    {
    outTree->Fill();
    }
    fill(sq11strip,sq11strip+16,0);
    fill(sq12strip,sq12strip+16,0);
    fill(sq13strip,sq13strip+16,0);
    fill(sq21strip,sq21strip+16,0);
    fill(sq22strip,sq22strip+16,0);
    fill(sq23strip,sq23strip+16,0);
    
    fill(sq11Edep,sq11Edep+16,0.);
    fill(sq12Edep,sq12Edep+16,0.);
    fill(sq13Edep,sq13Edep+16,0.);
    fill(sq21Edep,sq21Edep+16,0.);
    fill(sq22Edep,sq22Edep+16,0.);
    fill(sq23Edep,sq23Edep+16,0.);
}


outTree->Write();
outF->cd();
//outTree->Print();

cout<<"or maybe here?"<<endl;
outF->Close();
return 0;
}

