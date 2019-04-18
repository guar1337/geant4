#ifndef ReCoBuild_h
#define ReCoBuild_h 1



  
TLorentzVector* v2H;
TLorentzVector* v6He;
TLorentzVector* vBEAM;
TLorentzVector* vDeutCM;

TLorentzVector* vec2H;
TLorentzVector* vec6He;
TLorentzVector* vecBEAM;

double CsIdeut[16];
double SideutX[16];
double SideutY[16];
double phiDeut;
double Tdeut;
double deutEDEP;

double proba;
double CsIhe[16];
double SiheX[16];
double SiheY[16];
double thetaDeut;
double thetaHe;
double phiHe;
double The;
double heEDEP;

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

double inXpos;
double inYpos;
double inZpos;

double outXpos;
double outYpos;
double outZpos;

double X6He;
double Y6He;
double Z6He;

double X2H;
double Y2H;
double Z2H;

//out tree branches
Short_t sq11mult;
Short_t sq12mult;
Short_t sq13mult;
Short_t sq21mult;
Short_t sq22mult;
Short_t sq23mult;

Double_t sq11Edep[16];
Double_t sq12Edep[16];
Double_t sq13Edep[16];
Double_t sq21Edep[16];
Double_t sq22Edep[16];
Double_t sq23Edep[16];

Short_t sq11strip[16];
Short_t sq12strip[16];
Short_t sq13strip[16];
Short_t sq21strip[16];
Short_t sq22strip[16];
Short_t sq23strip[16];

Double_t relabAngHe;
Double_t relabAng2H;
Double_t rethetaCM;
Double_t reTbeam;



Double_t sq1de;
Double_t sq1etot;
Double_t sq1ang;

Double_t sq2de;
Double_t sq2etot;
Double_t sq2ang;












#endif
