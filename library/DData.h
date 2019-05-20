//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 15 10:31:53 2018 by ROOT version 6.15/01
// from TTree h100/Event
// found on file: event.root
//////////////////////////////////////////////////////////

#ifndef DData_h
#define DData_h

#include <TChain.h>
// #include <TClonesArray.h>
#include <TFile.h>

#include <map>

//#include "DParticle.h"

// Header file for the classes stored in the TTree Delphi events
// root file = event.root
// Tree      = h100 
class TH1F; 
class TH1D; 
class TH2D;
class TAxis;
class DPoolManager;
class DEvent;
class DParticle;
class DData : public TNamed {

public :
   enum Epid {kNoid = 0, kElectron = 11, kMuon = 13, kPhoton = 22, 
              kPi0 = 111, kK0L =130, kPion = 211, kK0S = 310, kK0 = 311, kK = 321,
              kN = 2121, kP = 2212, kL = 3122 }; 
   enum  Epopt  {kpNull, kAll, kCharged, kHadrons, kChargedHadrons}; 
   enum  Eopt   {kControlHisto, kControlDiff, kCorrelation, kSingleHisto};
   enum  Eparam {kpaNULL, kThrust, kMultiplicity};

   DData(const TString &name, const TString &title, const TString &dirname);
   virtual ~DData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Loop(const Eopt opt, const Epopt copt = kAll, const Eparam par = kpaNULL);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void             Control(); 
   void             Correlation(); 
   Float_t          GetEcm() const                         { fChain->GetEntry(0); return Ecm; }
   Bool_t           IsCharged(Int_t index) const; 
   Bool_t           IsChargedHadron(Int_t index) const; 
   Bool_t           IsHadron(Int_t index) const; 
   Bool_t           IsVerbose() const                      { return fVerbose; } 
   Long64_t         LoadEvent(Long64_t jentry); 
   Int_t            Leading(); 
   Double_t         Mass(Int_t pdg)                        { return fPMasses[pdg]; }
   void             Plot(const Eopt opt) const; 
   void             PrintMomentum() const; 
   void             PrintPID() const; 
   void             Run(const Eopt opt, const Epopt copt = kAll, const Eparam par = kpaNULL);
   void             SetMix(Int_t m)                        { fMixing = m; }
   void             SetMulBin(Int_t low, Int_t high)       { fMulLow = low; fMulHigh = high; }
   void             SetVerbosity(Int_t verb)               { fVerbose = verb; }
   void             SingleHisto(Eparam par); 
   void             WriteOutput(); 

   void SetPtBinning(Int_t nBins, Double_t *limits);
   void SetMultBinning(Int_t nBins, Double_t *limits);
   void SetZvtxBinning(Int_t nBins, Double_t *limits);
   Int_t GetMultBin() const;
   Bool_t RotateTracks(const DParticle *trig, const DParticle *assoc, Double_t &dphi, Double_t &dtheta, Double_t &deta) const;
   
private :
   void             CreateHistograms(Eopt opt, Eparam par = kpaNULL);
   Int_t            LoadEventMultiplicity(Long64_t jentry); 
   void             MakeMulPool(Int_t method);
   void             MakeParticlesList(Epopt opt = kAll);

   TTree                     *fChain;                 //!pointer to the analyzed TTree or TChain
   Int_t                     fCurrentTree;            //!current Tree number in a TChain
   Long64_t                  fCurrentEvent;           //!current event number
   Epopt                     fCurrentPartSel;         //!current particle tye selection
   TString                   fDataDirName;            //!path to data directory
   DEvent                    *fEvent;                 //!a Delphi event
   TList                     *fHistoList;             //!list of histograms to be filled 
   TFile                     *fHistosOutFile;         //!output file where histos are saved 
   Bool_t                    fHCreated;               //!false if no histo created for Control task
   Bool_t                    fHCorrelation;           //!false if no histo created for Correlation task
   TFile                     *fInputFile;             //!input file where ntuple and orginal histos are saved
   Long64_t                  fEvents;                 //!number of entries in the TTree
   Int_t                     fMixing;                 //!number of events used for mixing
   Int_t                     fMulHigh;                //!high multiplicity bin [fMulHigh, infinite]
   TObjArray                 fMulHighEventPool;       //!pool of high multiplicity events
   std::vector<Long64_t>     fMulHighIndexPool;       //!pool of high multiplicity events by index
   Int_t                     fMulLow;                 //!low multiplicity bin [0,fMulLow]
   TObjArray                 fMulLowEventPool;        //!pool of low multiplicity events
   std::vector<Long64_t>     fMulLowIndexPool;        //!pool of high multiplicity events by index
   Int_t                     fMultiplicity;           //!number of selected particles in the event
   std::map<Int_t, Double_t> fPMasses;                //!map key = pdg code, value = mass
   std::map<Epopt, TString>  fPopt;                   //!map key = Epopt, value = string
   TString                   fOutputFilename;         //!output file name where created histos are saved
   Long64_t fTriggersH;                               //!number of high multiplicity triggers
   Long64_t fTriggersL;                               //!number of low multiplicity triggers
   Int_t fVerbose;                                    //!level of verbosity: 0 = silent; 1 = warning; 2 = info

   // 2-particle correlation histos
   static const Int_t fNMaxBinsMult = 5;
   static const Int_t fNMaxBinsPt = 10;
   static const Int_t fNMaxBinsZvtx = 5;
   Int_t fNbinsMult;
   Int_t fNbinsTrackPt;
   Int_t fNbinsZvtx;
   TAxis *fMultAxis;
   TAxis *fTrackPtAxis;
   TAxis *fZvtxAxis;
   TH1D *fHistTrig[fNMaxBinsMult][fNMaxBinsZvtx][fNMaxBinsPt]; //!
   TH2D *fHistDPhiEta[fNMaxBinsMult][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsPt]; //!
   TH2D *fHistDPhiEtaMix[fNMaxBinsMult][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsPt]; //!   
   TH2D *fHistDPhiTheta[fNMaxBinsMult][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsPt]; //!
   TH2D *fHistDPhiThetaMix[fNMaxBinsMult][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsPt]; //!   
   DPoolManager *fPoolMgr; //!
   
   
   // Fixed size dimensions of array or collections stored in the TTree if any.

   // // Declaration of leaf types
   Float_t         Ecm;
   Float_t         Bt; 
   Float_t         Bz; 
   Float_t         Sphvec1x;
   Float_t         Sphvec1y;
   Float_t         Sphvec1z;
   Float_t         Sphvec2x;
   Float_t         Sphvec2y;
   Float_t         Sphvec2z;
   Float_t         Sphvec3x;
   Float_t         Sphvec3y;
   Float_t         Sphvec3z;
   Float_t         Sphval1;
   Float_t         Sphval2;
   Float_t         Sphval3;
   Float_t         Thrvec1x;
   Float_t         Thrvec1y;
   Float_t         Thrvec1z;
   Float_t         Thrvec2x;
   Float_t         Thrvec2y;
   Float_t         Thrvec2z;
   Float_t         Thrvec3x;
   Float_t         Thrvec3y;
   Float_t         Thrvec3z;
   Float_t         Thrval1;
   Float_t         Thrval2;
   Float_t         Thrval3;
   Int_t           Npa;
   UShort_t        Npac;
   Int_t           Paid[300];   //[Npa]
   Int_t           Pafl[300];   //[Npa]
   Int_t           Pav0d[300];   //[Npa]
   Int_t           Pav0m[300];   //[Npa]
   Float_t         Rvtx[300];   //[Npa]
   Float_t         Zvtx[300];   //[Npa]
   Float_t         Papx[300];   //[Npa]
   Float_t         Papy[300];   //[Npa]
   Float_t         Papz[300];   //[Npa]
   Int_t           Paje[300];   //[Npa]
   UChar_t         Pani[300];   //[Npa]
   Int_t           Njer;
   Float_t         Tgenr;
   Float_t         Dminr;
   Float_t         Jepx[10];   //[Njer]
   Float_t         Jepy[10];   //[Njer]
   Float_t         Jepz[10];   //[Njer]
   Float_t         Jee[10];   //[Njer]
   Float_t         Jem[10];   //[Njer]
   Int_t           Jep[10];   //[Njer]

   // // List of branches
   TBranch        *b_Ecm;   //!
   TBranch        *b_Bt;   //!
   TBranch        *b_Bz;   //!
   TBranch        *b_Sphvec1x;   //!
   TBranch        *b_Sphvec1y;   //!
   TBranch        *b_Sphvec1z;   //!
   TBranch        *b_Sphvec2x;   //!
   TBranch        *b_Sphvec2y;   //!
   TBranch        *b_Sphvec2z;   //!
   TBranch        *b_Sphvec3x;   //!
   TBranch        *b_Sphvec3y;   //!
   TBranch        *b_Sphvec3z;   //!
   TBranch        *b_Sphval1;   //!
   TBranch        *b_Sphval2;   //!
   TBranch        *b_Sphval3;   //!
   TBranch        *b_Thrvec1x;   //!
   TBranch        *b_Thrvec1y;   //!
   TBranch        *b_Thrvec1z;   //!
   TBranch        *b_Thrvec2x;   //!
   TBranch        *b_Thrvec2y;   //!
   TBranch        *b_Thrvec2z;   //!
   TBranch        *b_Thrvec3x;   //!
   TBranch        *b_Thrvec3y;   //!
   TBranch        *b_Thrvec3z;   //!
   TBranch        *b_Thrval1;   //!
   TBranch        *b_Thrval2;   //!
   TBranch        *b_Thrval3;   //!
   TBranch        *b_Npa;   //!
   TBranch        *b_Npac;   //!
   TBranch        *b_Paid;   //!
   TBranch        *b_Pafl;   //!
   TBranch        *b_Pav0d;   //!
   TBranch        *b_Pav0m;   //!
   TBranch        *b_Rvtx;   //!
   TBranch        *b_Zvtx;   //!
   TBranch        *b_Papx;   //!
   TBranch        *b_Papy;   //!
   TBranch        *b_Papz;   //!
   TBranch        *b_Paje;   //!
   TBranch        *b_Pani;   //!
   TBranch        *b_Njer;   //!
   TBranch        *b_Tgenr;   //!
   TBranch        *b_Dminr;   //!
   TBranch        *b_Jepx;   //!
   TBranch        *b_Jepy;   //!
   TBranch        *b_Jepz;   //!
   TBranch        *b_Jee;   //!
   TBranch        *b_Jem;   //!
   TBranch        *b_Jep;   //!
   ClassDef(DData,1)
};
#endif
