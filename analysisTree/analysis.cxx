#include "AnaCuts.h"
#include "HistsAndFun.h"   //could also include other headers
#include "myTree.h"  //D0 tree

#include "TLorentzVector.h"
#include "TString.h"

bool Debug = true;
int main(int argc, char** argv){
  if (argc!=3) {
    cout << "please input two value!" <<endl;
    cout << "first one is the file list and the second one is the output root file name  "<<endl;
    return -1;
  }

  TString inputlist;
  TString outputname;
  if (argc==3){
    inputlist = argv[1];
    outputname = argv[2];
    outputname +=".hists.root";
  }
  cout <<outputname.Data() << endl;
  //add files to the chain
  TChain *chain = new TChain("T");
  int nfile = 0;
  char tmp[2000];
  ifstream readlists;
  readlists.open(inputlist.Data());
  while (readlists.good()){
    readlists.getline(tmp,2000);
    TFile *ftmp = new TFile(tmp);
    if (!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
      cout<<"Could not open this file: "<< tmp  <<endl;
    }
    else {
       if(Debug && nfile%30==0) cout<<"read in "<<nfile<<"th file: "<< tmp <<endl;
      chain->Add(tmp);
      nfile++;      
    }
  }
  
  bookHists();
  
  TFile* outfile = new TFile(outputname.Data(),"recreate");
  myTree* t = new myTree(chain);
  Long64_t nEvents = chain->GetEntries();
  cout << "total "<<nEvents <<" events"<<endl;

  //events array
  std::vector <Events> mixevents[10][9];   //10 is nvzbins 9 is n centrality bins

  for (Long64_t ievt =0;ievt<nEvents;ievt++){
    if (Debug && ievt%10000==0) cout <<"process "<<ievt <<" events" <<endl;
    t->GetEntry(ievt);
    Int_t cent = t->centrality;
    Float_t pVtx_z = t->pVtx_z;
    int const vzbin = (int)((6 + pVtx_z) / 1.2);
    Int_t nkpi = t->nkpi;
    Int_t nspi = t->nspi;
    if (nkpi>100) cout<< "kpi array size limited!" <<endl;
    Float_t weight = t->gweight;
    mh1Cent->Fill(cent);
    mh1CentWg->Fill(cent,weight);
    Events mEvent;
    bool morethanoneD0 = false;
    if (nkpi<1) {
      // cout<<"None D0 int this event"<<endl;
      continue;
    }
    for (int id0=0;id0<nkpi;id0++){
      if (t->kpi_kidx[id0]==t->kpi_piidx[id0]) continue;
      if (!(t->kpi_knFit[id0] > anaCuts::nHitsFit && t->kpi_pinFit[id0] > anaCuts::nHitsFit) ) continue;
      float kpipt = t->kpi_pt[id0];
      float d0mass = t->kpi_mass[id0];
      bool gooddaughter = (t->kpi_kPt[id0]>anaCuts::minPt) && (t->kpi_piPt[id0]>anaCuts::minPt);
      if (!gooddaughter) continue;
      bool goodkpipair = isGoodPair(t->kpi_cosTheta[id0], t->kpi_piDca[id0], t->kpi_kDca[id0], t->kpi_dcaDaughters[id0], t->kpi_decayLength[id0], t->kpi_dcaD0ToPV[id0], kpipt, cent);
      if (!goodkpipair) continue;
      bool D0unlikesign = t->kpi_kCh[id0] * t->kpi_piCh[id0] <0 ? true : false;
      TLorentzVector d0Lorentz;
      d0Lorentz.SetPxPyPzE(t->kpi_px[id0],t->kpi_py[id0],t->kpi_pz[id0],t->kpi_E[id0]);
      if (fabs(d0Lorentz.Rapidity())>=1) continue;
      if (D0unlikesign) 
        mh3InvariantMassVsPtVsCent->Fill(kpipt,cent,d0mass,weight);
      else 
        mh3InvariantMassVsPtVsCentLike->Fill(kpipt,cent,d0mass,weight);

      //reconstruct D*
      bool goodD0Mass = (d0mass<anaCuts::mD0_max[cent]) && (d0mass>anaCuts::mD0_min[cent]);
      if (!goodD0Mass) continue;
      if (!D0unlikesign) continue;
    
      //save kpi to do mix event
      //only save 4Mom and k pi charg
      morethanoneD0 = true;
      KpiPair D0tmp;
      D0tmp.kch = t->kpi_kCh[id0];
      D0tmp.pich = t->kpi_piCh[id0];
      D0tmp.px = t->kpi_px[id0];
      D0tmp.py = t->kpi_py[id0];
      D0tmp.pz = t->kpi_pz[id0];
      D0tmp.e = t->kpi_E[id0];
      mEvent.D0vec.push_back(D0tmp);
      //////////////
      //cout<<"push D0 is ok"<<endl;
      if (nspi<1) {
        // cout<<"None softpion in this event"<<endl;
        continue;
      }
      for (int ipi=0;ipi<nspi;ipi++){
        if (t->kpi_kidx[id0]==t->spi_idx[ipi] || t->kpi_piidx[id0]==t->spi_idx[ipi] ) continue;

        //tof pid cut
        float tof = t->spi_Tof[ipi];
        bool tofnotavailable = fabs(tof-(-999))<1e-6;
        if (!tofnotavailable) {
          if (fabs(tof) > anaCuts::piTofBetaDiff) continue; 
        }

        TLorentzVector spiLorentz;
        spiLorentz.SetPxPyPzE(t->spi_px[ipi],t->spi_py[ipi],t->spi_pz[ipi],t->spi_E[ipi]);
        TLorentzVector DstarLorentz = d0Lorentz + spiLorentz;
        if (fabs(DstarLorentz.Rapidity())>=1) continue;
        float dstarmass = DstarLorentz.M();
        float deltamass = dstarmass - d0mass;
        if (deltamass<0.13 || deltamass>0.2) continue;
        float dstarPt = DstarLorentz.Pt();
        int dstarcharge = t->spi_Ch[ipi];
        bool rightsign = t->kpi_piCh[id0] * t->spi_Ch[ipi]>0;
        Double_t ww = weight;
        Double_t cc = cent;
        Double_t ch = dstarcharge;
        Double_t tmpForfill[4] = {dstarPt,spiLorentz.Pt(),cc,deltamass};
        if (rightsign)
          mhnDstarD0PiMassCentCharge->Fill(tmpForfill,ww);
        else 
          mhnDstarD0PiMassCentChargeLK->Fill(tmpForfill,ww);
      }  //softpion loop
    }  //D0 loop

    //////////////////////////////
    ///add a softpion loop for mixevent
    ///the cut condition keep exactly same as Dstar loop
    //save 4Mom and charge
    if (!morethanoneD0) continue;
    bool morethanonePi=false;
    if (nspi<1) continue;
    for (int ipi=0;ipi<nspi;ipi++){
        //tof pid cut
      float tof = t->spi_Tof[ipi];
      bool tofnotavailable = fabs(tof-(-999))<1e-6;
      if (!tofnotavailable) {
        if (fabs(tof) > anaCuts::piTofBetaDiff) continue; 
      }
      morethanonePi=true;
      SoftPi pitmp;
      pitmp.ch = t->spi_Ch[ipi];
      pitmp.px = t->spi_px[ipi];
      pitmp.py = t->spi_py[ipi];
      pitmp.pz = t->spi_pz[ipi];
      pitmp.e = t->spi_E[ipi];
      mEvent.softpivec.push_back(pitmp);
      morethanonePi = true;
    }  //softpion loop
    if (!morethanonePi) continue;
    //cout<<"push pion is ok"<<endl;
    //have been pushed the D0 and softpion into the event vector 
    ///////not done ! present not save nD0 and npion
    //push the event into the vz and cent array
    //cout<<vzbin<<endl;
    mixevents[vzbin][cent].push_back(mEvent);
    //cout<<"push mix events is ok"<<endl;
    int mEventsize = mixevents[vzbin][cent].size();

    if (mEventsize==MixEventBufferSize){
   // if (mEventsize==3){ //only for debug
      //mixevent with each event in the buffer
      //the last event in the buffer is the current event!
      int ND0_1 = mixevents[vzbin][cent][mEventsize-1].D0vec.size();
      int NPion_1 = mixevents[vzbin][cent][mEventsize-1].softpivec.size();
      int ND0_i;
      int NPion_i;
      
      for (int iEvt=0; iEvt<mEventsize-1;iEvt++){
        NPion_i=mixevents[vzbin][cent][iEvt].softpivec.size();
        ND0_i=mixevents[vzbin][cent][iEvt].D0vec.size();

        for (int jD0_1=0;jD0_1<ND0_1;jD0_1++){
          for (int jPion_i=0;jPion_i<NPion_i;jPion_i++){
            //fill hists
            fillMixEventsHists(&mixevents[vzbin][cent][mEventsize-1].D0vec[jD0_1], &mixevents[vzbin][cent][iEvt].softpivec[jPion_i],weight,cent);
            }
          }
        
        for (int jPion_1=0;jPion_1<NPion_1;jPion_1++){
          for (int jD0_i=0;jD0_i<ND0_i;jD0_i++){
            fillMixEventsHists(&mixevents[vzbin][cent][iEvt].D0vec[jD0_i], &mixevents[vzbin][cent][mEventsize-1].softpivec[jPion_1],weight,cent);
          }
        }
      }//eventsbuffer loop
      //remove the oldest event
       mixevents[vzbin][cent].erase( mixevents[vzbin][cent].begin());
    } //mEventsize==MixEventBufferSize
    //cout<<"end of mix events"<<endl;
  } //event loop

  WriteHists(outfile);
  outfile->Close();
  delete chain;
  cout<<"end of program"<<endl;
  return 1;
}

void bookHists(){
  //event
  mh1Cent         = new TH1F("mh1Cent", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
  mh1CentWg         = new TH1F("mh1CentWg", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);

  //D0
  mh3InvariantMassVsPtVsCent        = new TH3F("mh3InvariantMassVsPtVsCent", "invariantMassVsPtVsCent;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
  mh3InvariantMassVsPtVsCentLike    = new TH3F("mh3InvariantMassVsPtVsCentLike", "invariantMassVsPtVsCentLike;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);

  //D*
  int const dim = 4;
  int nbins[dim] = {120,100,10,90};
  double xmin[dim] = {0.,0.,-1.5, 0.135};
  double xmax[dim] = {12.,2., 8.5, 0.18};
  mhnDstarD0PiMassCentCharge = new THnSparseF("mhnDstarD0PiMassCentCharge","mhnDstarD0PiMassCentCharge;Dstar p_{T};#pi p_{T};Centrality;m_{K#pi#pi}-m_{K#pi}(GeV/c^{2});",dim,nbins,xmin,xmax);
  mhnDstarD0PiMassCentChargeLK = new THnSparseF("mhnDstarD0PiMassCentChargeLK","mhnDstarD0PiMassCentChargeLK;Dstar p_{T};#pi p_{T};Centrality;m_{K#pi#pi}-m_{K#pi}(GeV/c^{2});",dim,nbins,xmin,xmax);
  mhnDstarD0PiMassCentChargeMix = new THnSparseF("mhnDstarD0PiMassCentChargeMix","mhnDstarD0PiMassCentChargeMix;Dstar p_{T};#pi p_{T};Centrality;m_{K#pi#pi}-m_{K#pi}(GeV/c^{2});",dim,nbins,xmin,xmax);
  mhnDstarD0PiMassCentChargeMixLK = new THnSparseF("mhnDstarD0PiMassCentChargeMixLK","mhnDstarD0PiMassCentChargeMixLK;Dstar p_{T};#pi p_{T};Centrality;m_{K#pi#pi}-m_{K#pi}(GeV/c^{2});",dim,nbins,xmin,xmax);

  mhnDstarD0PiMassCentCharge->Sumw2();
  mhnDstarD0PiMassCentChargeLK->Sumw2();
  mhnDstarD0PiMassCentChargeMix->Sumw2();
  mhnDstarD0PiMassCentChargeMixLK->Sumw2();
}

void WriteHists(TFile* out){
  out->cd();
 //event
  mh1Cent->Write();
  mh1CentWg->Write();
 
  //D0
  mh3InvariantMassVsPtVsCent->Write();
  mh3InvariantMassVsPtVsCentLike->Write();
  //D*
  mhnDstarD0PiMassCentCharge->Write();
  mhnDstarD0PiMassCentChargeLK->Write();
  mhnDstarD0PiMassCentChargeMix->Write();
  mhnDstarD0PiMassCentChargeMixLK->Write();
}

void fillMixEventsHists(KpiPair* D0, SoftPi* pion, Float_t weight,int cent){
  TLorentzVector Dstar;
  TLorentzVector pionlorentz;
  pionlorentz.SetPxPyPzE(pion->px,pion->py,pion->pz,pion->e);
  TLorentzVector D0lorentz;
  D0lorentz.SetPxPyPzE(D0->px,D0->py,D0->pz,D0->e);
  Dstar = D0lorentz + pionlorentz;
  if ( Dstar.Rapidity() < anaCuts::RapidityCut){
    float dstarmass = Dstar.M();
    float deltamass = dstarmass - D0lorentz.M();
    if (deltamass > 0.2 && deltamass<0.12) return;
    float dstarPt = Dstar.Pt();
    //int dstarcharge = pion->ch;
    bool rightsign = D0->pich * pion->ch > 0;
    Double_t ww = weight;
    Double_t cc = cent;
    //Double_t ch = dstarcharge;
    Double_t tmpForfill[4] = {dstarPt,pionlorentz.Pt(),cc,deltamass};
    if (rightsign)
      mhnDstarD0PiMassCentChargeMix->Fill(tmpForfill,ww);
    else 
      mhnDstarD0PiMassCentChargeMixLK->Fill(tmpForfill,ww);
  }
      ///////debug//
  //  cout<<D0lorentz.M()<<" "<<Dstar.M()<<" "
  //  <<pion->px<<" "<<pion->e<<" "<<pion->py<<" "
  //  <<pion->pz<<" "<<pion->e<<endl;
   ////////////////////////////
}
