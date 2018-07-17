//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 12 02:57:25 2018 by ROOT version 5.34/30
// from TTree T/T
// found on file: 3FDF867F403E69FFB759AAEEECF87085_877.picoD0.root
//////////////////////////////////////////////////////////

#ifndef myTree_h
#define myTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
#define MAXD0 100
#define MAXPI 1200

class myTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runId;
   Int_t           eventId;
   Float_t         bfield;
   Float_t         pVtx_z;
   Int_t           centrality;
   Float_t         gweight;
   Float_t         grefMultCorr;
   Int_t           nkpi;
   Float_t         kpi_mass[MAXD0];   //[nkpi]
   Float_t         kpi_pt[MAXD0];   //[nkpi]
   Float_t         kpi_dcaDaughters[MAXD0];   //[nkpi]
   Float_t         kpi_decayLength[MAXD0];   //[nkpi]
   Float_t         kpi_dcaD0ToPV[MAXD0];   //[nkpi]
   Float_t         kpi_cosTheta[MAXD0];   //[nkpi]
   Float_t         kpi_kPt[MAXD0];   //[nkpi]
   Float_t         kpi_kDca[MAXD0];   //[nkpi]
   Int_t           kpi_kCh[MAXD0];   //[nkpi]
   Int_t           kpi_kidx[MAXD0];   //[nkpi]
   Int_t           kpi_knFit[MAXD0];   //[nkpi]
   Float_t         kpi_piPt[MAXD0];   //[nkpi]
   Float_t         kpi_piDca[MAXD0];   //[nkpi]
   Int_t           kpi_piCh[MAXD0];   //[nkpi]
   Int_t           kpi_piidx[MAXD0];   //[nkpi]
   Int_t           kpi_pinFit[MAXD0];   //[nkpi]
   Float_t         kpi_px[MAXD0];   //[nkpi]
   Float_t         kpi_py[MAXD0];   //[nkpi]
   Float_t         kpi_pz[MAXD0];   //[nkpi]
   Float_t         kpi_E[MAXD0];   //[nkpi]
   Int_t           nspi;
   Float_t         spi_px[MAXPI];   //[nspi]
   Float_t         spi_py[MAXPI];   //[nspi]
   Float_t         spi_pz[MAXPI];   //[nspi]
   Float_t         spi_E[MAXPI];   //[nspi]
   Float_t         spi_Dca[MAXPI];   //[nspi]
   Int_t           spi_Ch[MAXPI];   //[nspi]
   Int_t           spi_idx[MAXPI];   //[nspi]
   Float_t         spi_Tof[MAXPI];   //[nspi]
   Int_t           spi_nFit[MAXPI];   //[nspi]

   // List of branches
   TBranch        *b_runId;   //!
   TBranch        *b_eventId;   //!
   TBranch        *b_bfield;   //!
   TBranch        *b_pVtx_z;   //!
   TBranch        *b_centrality;   //!
   TBranch        *b_gweight;   //!
   TBranch        *b_grefMultCorr;   //!
   TBranch        *b_nkpi;   //!
   TBranch        *b_kpi_mass;   //!
   TBranch        *b_kpi_pt;   //!
   TBranch        *b_kpi_dcaDaughters;   //!
   TBranch        *b_kpi_decayLength;   //!
   TBranch        *b_kpi_dcaD0ToPV;   //!
   TBranch        *b_kpi_cosTheta;   //!
   TBranch        *b_kpi_kPt;   //!
   TBranch        *b_kpi_kDca;   //!
   TBranch        *b_kpi_kCh;   //!
   TBranch        *b_kpi_kidx;   //!
   TBranch        *b_kpi_knFit;   //!
   TBranch        *b_kpi_piPt;   //!
   TBranch        *b_kpi_piDca;   //!
   TBranch        *b_kpi_piCh;   //!
   TBranch        *b_kpi_piidx;   //!
   TBranch        *b_kpi_pinFit;   //!
   TBranch        *b_kpi_px;   //!
   TBranch        *b_kpi_py;   //!
   TBranch        *b_kpi_pz;   //!
   TBranch        *b_kpi_E;   //!
   TBranch        *b_nspi;   //!
   TBranch        *b_spi_px;   //!
   TBranch        *b_spi_py;   //!
   TBranch        *b_spi_pz;   //!
   TBranch        *b_spi_E;   //!
   TBranch        *b_spi_Dca;   //!
   TBranch        *b_spi_Ch;   //!
   TBranch        *b_spi_idx;   //!
   TBranch        *b_spi_Tof;   //!
   TBranch        *b_spi_nFit;   //!

   myTree(TTree *tree=0);
   virtual ~myTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myTree_cxx
myTree::myTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("3FDF867F403E69FFB759AAEEECF87085_877.picoD0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("3FDF867F403E69FFB759AAEEECF87085_877.picoD0.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

myTree::~myTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runId", &runId, &b_runId);
   fChain->SetBranchAddress("eventId", &eventId, &b_eventId);
   fChain->SetBranchAddress("bfield", &bfield, &b_bfield);
   fChain->SetBranchAddress("pVtx_z", &pVtx_z, &b_pVtx_z);
   fChain->SetBranchAddress("centrality", &centrality, &b_centrality);
   fChain->SetBranchAddress("gweight", &gweight, &b_gweight);
   fChain->SetBranchAddress("grefMultCorr", &grefMultCorr, &b_grefMultCorr);
   fChain->SetBranchAddress("nkpi", &nkpi, &b_nkpi);
   fChain->SetBranchAddress("kpi_mass", &kpi_mass, &b_kpi_mass);
   fChain->SetBranchAddress("kpi_pt", &kpi_pt, &b_kpi_pt);
   fChain->SetBranchAddress("kpi_dcaDaughters", &kpi_dcaDaughters, &b_kpi_dcaDaughters);
   fChain->SetBranchAddress("kpi_decayLength", &kpi_decayLength, &b_kpi_decayLength);
   fChain->SetBranchAddress("kpi_dcaD0ToPV", &kpi_dcaD0ToPV, &b_kpi_dcaD0ToPV);
   fChain->SetBranchAddress("kpi_cosTheta", &kpi_cosTheta, &b_kpi_cosTheta);
   fChain->SetBranchAddress("kpi_kPt", &kpi_kPt, &b_kpi_kPt);
   fChain->SetBranchAddress("kpi_kDca", &kpi_kDca, &b_kpi_kDca);
   fChain->SetBranchAddress("kpi_kCh", &kpi_kCh, &b_kpi_kCh);
   fChain->SetBranchAddress("kpi_kidx", &kpi_kidx, &b_kpi_kidx);
   fChain->SetBranchAddress("kpi_knFit", &kpi_knFit, &b_kpi_knFit);
   fChain->SetBranchAddress("kpi_piPt", &kpi_piPt, &b_kpi_piPt);
   fChain->SetBranchAddress("kpi_piDca", &kpi_piDca, &b_kpi_piDca);
   fChain->SetBranchAddress("kpi_piCh", &kpi_piCh, &b_kpi_piCh);
   fChain->SetBranchAddress("kpi_piidx", &kpi_piidx, &b_kpi_piidx);
   fChain->SetBranchAddress("kpi_pinFit", &kpi_pinFit, &b_kpi_pinFit);
   fChain->SetBranchAddress("kpi_px", &kpi_px, &b_kpi_px);
   fChain->SetBranchAddress("kpi_py", &kpi_py, &b_kpi_py);
   fChain->SetBranchAddress("kpi_pz", &kpi_pz, &b_kpi_pz);
   fChain->SetBranchAddress("kpi_E", &kpi_E, &b_kpi_E);
   fChain->SetBranchAddress("nspi", &nspi, &b_nspi);
   fChain->SetBranchAddress("spi_px", &spi_px, &b_spi_px);
   fChain->SetBranchAddress("spi_py", &spi_py, &b_spi_py);
   fChain->SetBranchAddress("spi_pz", &spi_pz, &b_spi_pz);
   fChain->SetBranchAddress("spi_E", &spi_E, &b_spi_E);
   fChain->SetBranchAddress("spi_Dca", &spi_Dca, &b_spi_Dca);
   fChain->SetBranchAddress("spi_Ch", &spi_Ch, &b_spi_Ch);
   fChain->SetBranchAddress("spi_idx", &spi_idx, &b_spi_idx);
   fChain->SetBranchAddress("spi_Tof", &spi_Tof, &b_spi_Tof);
   fChain->SetBranchAddress("spi_nFit", &spi_nFit, &b_spi_nFit);
   Notify();
}

Bool_t myTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myTree_cxx
