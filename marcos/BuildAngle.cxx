#include "TFile.h"
#include "TTree.h"
#include "stdio.h"
#include "Riostream.h"
#include "BuildAngle.hh"
#include <iostream>
#include "TLorentzVector.h"
#include "TVector3.h"
// voron1392

void BuildAngle()
//TFile *inF  = new TFile(TString::Format("he6_000%d.root", iii), "READ");
{
/*
ifstream plik("/home/guar/Desktop/simulation/include/params.dat");
string dummy;
float CConst;
if (!plik)
    {
    cout << "Cannot open file: params.dat\n";
    return;
    }
getline(plik, dummy);
plik>>heAng>>heDist;
plik>>deAng>>deDist;
//cout<<heDist<<"  "<<deDist<<endl; 
plik.close();
   */
TVector3 beamDirection(0,0,1);
TFile *inF  = new TFile("ReCoGu.root", "READ");
TFile *outF = new TFile("ReAngle.root","recreate");
TTree *inTree = (TTree*)inF->Get("Reconstructed");
TTree *outTree= new TTree("GettingAngles","reconstructed");
outTree->SetMakeClass(1);

inTree->Print();
v2H_N=  new TLorentzVector();
v6He_N= new TLorentzVector();
vBEAM_N= new TLorentzVector();
inTree->SetBranchAddress("vBEAM.",&vBEAM_N);
inTree->SetBranchAddress("v2H.",&v2H_N);
inTree->SetBranchAddress("v6He.",&v6He_N);

inTree->SetBranchAddress("Xpos",          &Xpos);
inTree->SetBranchAddress("Ypos",          &Ypos);
inTree->SetBranchAddress("Zpos",          &Zpos);

//ReCo - detectors
inTree->SetBranchAddress("sq11mult",      &Oldsq11mult);
inTree->SetBranchAddress("sq21mult",      &Oldsq21mult);
inTree->SetBranchAddress("sq12mult",      &Oldsq12mult);
inTree->SetBranchAddress("sq22mult",      &Oldsq22mult);
inTree->SetBranchAddress("sq13mult",      &Oldsq13mult);
inTree->SetBranchAddress("sq23mult",      &Oldsq23mult);

inTree->SetBranchAddress("sq11Edep",      Oldsq11Edep);
inTree->SetBranchAddress("sq21Edep",      Oldsq21Edep);
inTree->SetBranchAddress("sq12Edep",      Oldsq12Edep);
inTree->SetBranchAddress("sq22Edep",      Oldsq22Edep);
inTree->SetBranchAddress("sq13Edep",      Oldsq13Edep);
inTree->SetBranchAddress("sq23Edep",      Oldsq23Edep);

inTree->SetBranchAddress("sq11strip",     Oldsq11strip);
inTree->SetBranchAddress("sq21strip",     Oldsq21strip);
inTree->SetBranchAddress("sq12strip",     Oldsq12strip);
inTree->SetBranchAddress("sq22strip",     Oldsq22strip);
inTree->SetBranchAddress("sq13strip",     Oldsq13strip);
inTree->SetBranchAddress("sq23strip",     Oldsq23strip);

outTree->Bronch("vBEAM.",       "TLorentzVector",   &vBEAM_N);
outTree->Bronch("v2H.",         "TLorentzVector",   &v2H_N);
outTree->Bronch("v6He.",        "TLorentzVector",   &v6He_N);

outTree->Branch("sq1dE",         sq1dE,          "sq1dE[16]/D");
outTree->Branch("sq1Etot",       sq1Etot,        "sq1Etot[16]/D");
outTree->Branch("sq1angle",      sq1angle,       "sq1angle[16]/D");
outTree->Branch("sq2dE",         sq2dE,          "sq2dE[16]/D");
outTree->Branch("sq2Etot",       sq2Etot,        "sq2Etot[16]/D");
outTree->Branch("sq2angle",      sq2angle,       "sq2angle[16]/D");

//ReCo - straight from simulation
//inTree->SetBranchAddress("labAngHe",      &relabAngHe,     "labAngHe/D");
//inTree->SetBranchAddress("labAng2H",      &relabAng2H,     "labAng2H/D");
inTree->SetBranchAddress("thetaCM",       &OldthetaCM);
inTree->SetBranchAddress("Tbeam",         &OldTbeam);
int checker=0;
Long64_t nEntries = inTree->GetEntries();
cout<<nEntries<<endl;
for (Long64_t entry=0; entry<nEntries; entry++)
{
inTree->GetEntry(entry);
int multX=0;
int multY=0;
int newMult1=0;
int newMult2=0;
vBEAM_S=vBEAM_N;
v2H_S=v2H_N;
v6He_S=v6He_N;
int X_mult_No=0;
int Y_mult_No=0;

for (int iii=0; iii<4; iii++)//looping on multiplicity - but I want to look at CsI
    {
    int X_No=(Oldsq13strip[iii]%4)*4; //This number of the first strip ina' row     
    for(int jjj=0; jjj<4; jjj++)//for a given multi in CsI I am scanning Si
        {
        //if ((Oldsq11Edep[jjj])>0.5)
        //{cout<<X_No<<" "<<Oldsq11strip[jjj]<<" "<<X_No+3<<endl; }          
        if((Oldsq11Edep[jjj])>0.5 && (X_No+3)>=Oldsq11strip[jjj] && (X_No)<=Oldsq11strip[jjj])
            {
            //cout<<"Do you even?"<<endl;
            multX++;
            if(multX==1) 
                {
                X_mult_No=jjj;
                }
            }
        }
    for(int jjj=0; jjj<16; jjj++)
        {
        int Y_No=(Oldsq13strip[iii]/4)*4; //This number of the first strip ina' row        
        if((Oldsq12Edep[jjj])>0.5 && (Y_No)<=Oldsq12strip[jjj] && (Y_No+3)>=Oldsq12strip[jjj])
            {                
            multY++;
            if(multY==1) 
                {
                Y_mult_No=jjj;
                }
            }
        }   
    //cout<<"multX= "<<multX<<", "<<multY<<""<<endl;    
    if(multX==1 && multY==1 && Oldsq13Edep[iii]>5.0)
        {       
        newMult1++;
        sq1dE[newMult1]=Oldsq11Edep[X_mult_No];
        sq1Etot[newMult1]=Oldsq11Edep[X_mult_No]+Oldsq13Edep[iii];
        X2H=180.*sin(45./180.*3.14159)+(30.-4.*X_mult_No)*cos(45./180.*3.14159);
        Y2H=+30.-4.*Y_mult_No;
        Z2H=180.*cos(45./180.*3.14159)-(30.-4.*X_mult_No)*sin(45./180.*3.14159);
        
        TVector3 vect2H(X2H-Xpos, Y2H-Ypos, Z2H+Xpos);
        sq1angle[newMult1]=vect2H.Angle(beamDirection)*180.0/3.1415927;
//cout<<"un "<<X<<" zen "<<iii<<" du "<<sq1Etot[newMult1]<<" tre "<<Oldsq13Edep[iii]<<endl;
        }
    multX=0;
    multY=0;
    }


for (int iii=0; iii<16; iii++)//looping on multiplicity - but I want to look at CsI
    {   
    for(int jjj=0; jjj<16; jjj++)//for a given multi in CsI I am scanning Si
        {
        int X_No=(Oldsq23strip[iii]%4)*4; //This number of the first strip ina' row        
        if((Oldsq21Edep[jjj])>0.5 && (X_No)<=Oldsq21strip[jjj] && (X_No+3)>=Oldsq21strip[jjj])
            {
            multX++;
            if(multX==1) 
                {
                X_mult_No=jjj;
                }
            }
        }
    for(int jjj=0; jjj<16; jjj++)
        {
        int Y_No=(Oldsq23strip[iii]/4)*4; //This number of the first strip ina' row        
        if((Oldsq22Edep[jjj])>0.5 && (Y_No)<=Oldsq22strip[jjj] && (Y_No+3)>=Oldsq22strip[jjj])
            {                
            multY++;
            if(multY==1) 
                {
                Y_mult_No=jjj;
                }
            }
        }      
    if(multX==1 && multY==1 && Oldsq23Edep[iii]>5.0)
        {
        newMult2++;
        sq2dE[newMult2]=Oldsq21Edep[X_mult_No];
        sq2Etot[newMult2]=Oldsq21Edep[X_mult_No]+Oldsq23Edep[iii];
        X6He=-180.*sin(15./180.*3.14159)+(30.-4.*X_mult_No)*cos(15./180.*3.14159);
        Y6He=+30.-4.*Y_mult_No;
        Z6He=180.*cos(15./180.*3.14159)+(30.-4.*X_mult_No)*sin(15./180.*3.14159);
        TVector3 vect6He(X6He-Xpos, Y6He-Ypos, Z6He+Xpos);
        sq2angle[newMult2]=vect6He.Angle(beamDirection)*180.0/3.1415927;
//cout<<"un "<<X<<" zen "<<iii<<" du "<<sq1Etot[newMult1]<<" tre "<<Oldsq13Edep[iii]<<endl;
        }
    multX=0;
    multY=0;
    }



    if((newMult2+newMult1)>0)
    {
    outTree->Fill();
    }
    fill(sq1dE,sq1dE+16,0.);
    fill(sq1Etot,sq1Etot+16,0.);
    fill(sq1angle,sq1angle+16,0.);
    fill(sq2dE,sq2dE+16,0.);
    fill(sq2Etot,sq2Etot+16,0.);
    fill(sq2angle,sq2angle+16,0.);    
}
cout<<checker;

 
outF->cd();
outTree->Write();
cout<<" or maybe here?"<<endl;
outF->Close();
return 0;
}
/*

X6He=-180.*sin(15./180.*3.14159)+(30.-4.*mult2X)*cos(15./180.*3.14159);
Y6He=+30.-4.*mult2Y;
Z6He=180.*cos(15./180.*3.14159)+(30.-4.*mult2X)*sin(15./180.*3.14159);
        
        TVector3 vect6He(X6He-Xpos, Y6He-Ypos, Z6He+Xpos);
        sq2Angle[mult2]=vect6He.Angle(beamDirection);

        X2H=180.*sin(45./180.*3.14159)+(30.-4.*mult1X)*cos(45./180.*3.14159);
        Y2H=+30.-4.*mult1Y;
        Z2H=180.*cos(45./180.*3.14159)-(30.-4.*mult1X)*sin(45./180.*3.14159);
        
        TVector3 vect2H(X2H-Xpos, Y2H-Ypos, Z2H+Xpos);
        sq1Angle[mult1]=vect2H.Angle(beamDirection);
        
        
        */

