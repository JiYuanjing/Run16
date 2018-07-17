#ifndef AnaCuts_H
#define AnaCuts_H

#include "Rtypes.h"
#include <string>
#include <array>

namespace anaCuts
{

   //event
   std::array<unsigned int, 6> const triggers = {
    520001,    //   VPDMB-5-p-sst_
    520011,   
    520021,  
    520031,    
    520041,
    520051 
    };    

   float const vz = 6.0;// cm.

   float const minPt = 0.3;//1.2
   // float const minPt = 0.3;//1.2
   int const nHitsFit = 20;

   //track eta cut
   float const Eta = 1.0;
   //pions
   float const nSigmaPion = 3.0;
   float const piTofBetaDiff = 0.03;
   const float spiTofBetaDiff_low[11] = {-0.03,-0.03,-0.022,-0.02,-0.02,-0.02,-0.0218,-0.0226,-0.023,-0.0235,-0.025};
   const float spiTofBetaDiff_high[11] = {0.18,0.178,0.095,0.073,0.059,0.05,0.05,0.047,0.043,0.0386,0.035};    
   //kaons
   float const nSigmaKaon = 2.0;
   float const kTofBetaDiff = 0.03;

   //kaons
   float const nSigmaProton = 2.0;
   float const pTofBetaDiff = 0.03;

   float const RapidityCut = 1.0;

  //tracking
   const int nCent = 5;
   const float CentEdge[nCent+1] = { -0.5, 1.5, 3.5, 5.5, 6.5, 8.5 };
   const int nPtBins = 6;
   const float PtBinsEdge[nPtBins+1] = {0, 0.5, 1., 2., 3., 5., 15.};
    // default
    float const cosTheta = 0.95;
    float const kDca[nCent][nPtBins] = {
        { 0.0113, 0.0103, 0.0081, 0.0066, 0.0046, 0.0038},  //60%-80%
        { 0.0110, 0.0112, 0.0081, 0.0063, 0.0064, 0.0044},  //40%-60%
        { 0.0098, 0.0089, 0.0074, 0.0085, 0.0063, 0.0049},  //20%-40%
        { 0.0111, 0.0099, 0.0091, 0.0097, 0.0062, 0.0061},  //10%-20%
        { 0.0104, 0.0099, 0.0073, 0.0086, 0.0067, 0.0056}   //0-10%
    };
    float const pDca[nCent][nPtBins] = {
        { 0.0100, 0.0096, 0.0093, 0.0094, 0.0059, 0.0050},  //60%-80%
        { 0.0107, 0.0106, 0.0097, 0.0078, 0.0063, 0.0056},  //40%-60%
        { 0.0117, 0.0106, 0.0097, 0.0066, 0.0064, 0.0056},  //20%-40%
        { 0.0098, 0.0110, 0.0101, 0.0086, 0.0070, 0.0061},  //10%-20%
        { 0.0109, 0.0106, 0.0081, 0.0092, 0.0080, 0.0057}   //0-10%
    };
    float const dcaV0ToPv[nCent][nPtBins] = {
        { 0.0075, 0.0066, 0.0064, 0.0050, 0.0058, 0.0038},  //60%-80%
        { 0.0065, 0.0064, 0.0046, 0.0049, 0.0054, 0.0057},  //40%-60%
        { 0.0045, 0.0048, 0.0042, 0.0043, 0.0052, 0.0055},  //20%-40%
        { 0.0052, 0.0049, 0.0043, 0.0042, 0.0043, 0.0050},  //10%-20%
        { 0.0061, 0.0046, 0.0042, 0.0041, 0.0037, 0.0048}   //0-10%
    };
    float const dcaDaughters[nCent][nPtBins] = {
        { 0.0073, 0.0088, 0.0092, 0.0082, 0.0083, 0.0104},  //60%-80%
        { 0.0076, 0.0087, 0.0090, 0.0082, 0.0101, 0.0093},  //40%-60%
        { 0.0078, 0.0067, 0.0069, 0.0066, 0.0073, 0.0099},  //20%-40%
        { 0.0070, 0.0070, 0.0062, 0.0067, 0.0076, 0.0085},  //10%-20%
        { 0.0066, 0.0079, 0.0061, 0.0063, 0.0076, 0.0068}   //0-10%
    };
    float const decayLength[nCent][nPtBins] = {
        { 0.0150, 0.0107, 0.0175, 0.0187, 0.0164, 0.0175},  //60%-80%
        { 0.0140, 0.0133, 0.0190, 0.0201, 0.0215, 0.0219},  //40%-60%
        { 0.0149, 0.0170, 0.0205, 0.0236, 0.0234, 0.0237},  //20%-40%
        { 0.0151, 0.0173, 0.0204, 0.0240, 0.0237, 0.0231},  //10%-20%
        { 0.0128, 0.0163, 0.0222, 0.0214, 0.0241, 0.0253}   //0-10%
    };

//------------------------------------------------------------
	float const minD0Mass = 1.6;
	float const maxD0Mass = 2.2;
//D_star cut 
float const mD0_max[9]={1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.89,1.89};
float const mD0_min[9]={1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.83,1.83};
float const cosThetaStar=0.8;
float const ptD0_min=7;
float const ptSoftPion_max=1.5;
float const ptSoftPion_min=0.15;
float const DcaSoftPion=3;
std::pair<float, float> const D0toPion_pt
{ 7, 20};
}
#endif
