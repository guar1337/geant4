#ifndef BuildAngle_h
#define BuildAngle_h 1


TLorentzVector* v2H_N;
TLorentzVector* v6He_N;
TLorentzVector* vBEAM_N;
  
TLorentzVector* v2H_S;
TLorentzVector* v6He_S;
TLorentzVector* vBEAM_S;

double phiCM;
double thetaCM;


double BEAMenergy;
double phiBEAM;
double theBEAM;
double Tbeam;

double labAng2H;
double labAngHe;
//double exp2;
double exp1;

double mass2H;
double mass6He;
double mass4He;

double radThetaCM;
double reTheta2H;
double reTheta6He;

double Xpos;
double Ypos;
double Zpos;

double X6He;
double Y6He;
double Z6He;

double X2H;
double Y2H;
double Z2H;

//out tree branches
Short_t Oldsq11mult;
Short_t Oldsq12mult;
Short_t Oldsq13mult;
Short_t Oldsq21mult;
Short_t Oldsq22mult;
Short_t Oldsq23mult;

Double_t Oldsq11Edep[16];
Double_t Oldsq12Edep[16];
Double_t Oldsq13Edep[16];
Double_t Oldsq21Edep[16];
Double_t Oldsq22Edep[16];
Double_t Oldsq23Edep[16];

Short_t Oldsq11strip[16];
Short_t Oldsq12strip[16];
Short_t Oldsq13strip[16];
Short_t Oldsq21strip[16];
Short_t Oldsq22strip[16];
Short_t Oldsq23strip[16];

Double_t relabAngHe;
Double_t relabAng2H;
Double_t OldthetaCM;
Double_t OldTbeam;


Double_t sq1dE[16];
Double_t sq1Etot[16];
Double_t sq1angle[16];
Double_t sq2dE[16];
Double_t sq2Etot[16];
Double_t sq2angle[16];






















#endif
