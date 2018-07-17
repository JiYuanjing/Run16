#include "AnaCuts.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "math.h"
#include "string.h"

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "myTree.h"
#include "THnSparse.h"
using namespace std;
#endif


  typedef struct kpipair{
    int kch;
    int pich;
//    TLorentzVector lorentz;
    float px;
    float py;
    float pz;
    float e;
  }KpiPair;

  typedef struct softpi{
    int ch;
//    TLorentzVector lorentz;
   float px;
    float py;
    float pz;
    float e;
  }SoftPi;

  typedef struct events{
    vector<KpiPair> D0vec;
    vector<softpi> softpivec;
  }Events;

  const int MixEventBufferSize=6;
  TH1F* mh1CentWg = NULL;
  TH1F* mh1Cent = NULL;

  TH3F* mh3InvariantMassVsPtVsCent = NULL;
  TH3F* mh3InvariantMassVsPtVsCentLike = NULL;  

  THnSparseF* mhnDstarD0PiMassCentCharge = NULL;
  THnSparseF* mhnDstarD0PiMassCentChargeLK = NULL;
  THnSparseF* mhnDstarD0PiMassCentChargeMix = NULL;
  THnSparseF* mhnDstarD0PiMassCentChargeMixLK = NULL;

  int getD0PtIndex(float pt)
  {
     for (int i = 0; i < anaCuts::nPtBins; i++)
     {
        if ((pt >= anaCuts::PtBinsEdge[i]) && (pt < anaCuts::PtBinsEdge[i + 1]))
           return i;
     }
     return anaCuts::nPtBins - 1;
  }

  int getD0CentIndex(int cent)
  {
     for (int i = 0; i < anaCuts::nCent; i++)
     {
        if ((cent >= anaCuts::CentEdge[i]) && (cent < anaCuts::CentEdge[i + 1]))
           return i;
     }
     return anaCuts::nCent - 1;
  }

  bool isGoodPair(float cosTheta,float piDca,float kDca,float dcaDaughters,float decayLength,float dcaD0ToPV,float pt,int cent){
    int tmpIndex = getD0PtIndex(pt);
    int centIndex = getD0CentIndex(cent);
    //leave cent bin
    return cos(cosTheta) > anaCuts::cosTheta &&
    piDca > anaCuts::pDca[centIndex][tmpIndex] && kDca > anaCuts::kDca[centIndex][tmpIndex] &&
    dcaDaughters < anaCuts::dcaDaughters[centIndex][tmpIndex] &&
    decayLength > anaCuts::decayLength[centIndex][tmpIndex] &&
    dcaD0ToPV < anaCuts::dcaV0ToPv[centIndex][tmpIndex];;
  }

  void bookHists();
  void WriteHists(TFile* out);
  void fillMixEventsHists( KpiPair*  D0,  SoftPi*  pion, Float_t weight,int cent);