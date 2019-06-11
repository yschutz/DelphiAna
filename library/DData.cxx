#define DData_cxx
#include "DData.h"
#include "DEvent.h"
#include "DParticle.h"
#include "DPool.h"

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>
#include <iostream>
#include <chrono>

ClassImp(DData)

    //==========================================================================
    DData::DData(const TString &name, const TString &title, const TString &dirname) : TNamed(name, title), fChain(nullptr),
                                                                                      fDataDirName(dirname),
                                                                                      fEvent(nullptr),
                                                                                      fHCreated(kFALSE), fHCorrelation(kFALSE),
                                                                                      fHistoList(nullptr), fHistosOutFile(nullptr),
                                                                                      fMixing(10), fMultiplicity(0),
                                                                                      fTriggersH(0), fTriggersL(0),
                                                                                      fVerbose(0),
                                                                                      fNbinsMult(0),
                                                                                      fNbinsTrackPt(0),
                                                                                      fNbinsZvtx(1),
                                                                                      fMultAxis(0x0),
                                                                                      fTrackPtAxis(0x0),
                                                                                      fZvtxAxis(0x0),
                                                                                      fPoolMgr(0x0)

{
    // ctor
    fPopt[kAll] = TString("All");
    fPopt[kCharged] = TString("Charged");
    fPopt[kHadrons] = TString("Charged Hadrons");

    // where to find the root data file
    TString inputFileName("event.root");
    inputFileName.Prepend(fDataDirName);

    //where to store the output histograms
    fOutputFilename = "histos";
    fOutputFilename.Prepend(fDataDirName);
    // name o the ntuple in the root data file
    const TString kNtupleName("h100");
    fInputFile = new TFile(inputFileName, "READ");
    if (!fInputFile || !fInputFile->IsOpen())
    {
        std::cout << "ERROR: DData --> " << inputFileName << " not found!" << std::endl;
        exit(1);
    }
    fInputFile->GetObject(kNtupleName, fChain);
    Init();

    for (Int_t iMult = 0; iMult < fNMaxBinsMult; iMult++)
    {
        for (Int_t iZvtx = 0; iZvtx < fNMaxBinsZvtx; ++iZvtx)
        {
            for (Int_t iPtBin = 0; iPtBin < fNMaxBinsPt; iPtBin++)
            {
                fHistTrig[iMult][iZvtx][iPtBin] = NULL;
                for (Int_t jPtBin = 0; jPtBin < fNMaxBinsPt; jPtBin++)
                {
                    fHistDPhiEta[iMult][iZvtx][iPtBin][jPtBin] = NULL;
                    fHistDPhiEtaMix[iMult][iZvtx][iPtBin][jPtBin] = NULL;
                    fHistDPhiTheta[iMult][iZvtx][iPtBin][jPtBin] = NULL;
                    fHistDPhiThetaMix[iMult][iZvtx][iPtBin][jPtBin] = NULL;
                }
            }
        }
    }
}

//==========================================================================
    DData::DData(const TString &name, const TString &title, TChain *chain) : TNamed(name, title), fChain(chain),
                                                                                      fEvent(nullptr),
                                                                                      fHCreated(kFALSE), fHCorrelation(kFALSE),
                                                                                      fHistoList(nullptr), fHistosOutFile(nullptr),
                                                                                      fMixing(10), fMultiplicity(0),
                                                                                      fTriggersH(0), fTriggersL(0),
                                                                                      fVerbose(0),
                                                                                      fNbinsMult(0),
                                                                                      fNbinsTrackPt(0),
                                                                                      fNbinsZvtx(1),
                                                                                      fMultAxis(0x0),
                                                                                      fTrackPtAxis(0x0),
                                                                                      fZvtxAxis(0x0),
                                                                                      fPoolMgr(0x0)

{
    // ctor
    fPopt[kAll] = TString("All");
    fPopt[kCharged] = TString("Charged");
    fPopt[kHadrons] = TString("Charged Hadrons");

   
    //where to store the output histograms
    fOutputFilename = "histos";
    // name o the ntuple in the root data file
    Init();

    for (Int_t iMult = 0; iMult < fNMaxBinsMult; iMult++)
    {
        for (Int_t iZvtx = 0; iZvtx < fNMaxBinsZvtx; ++iZvtx)
        {
            for (Int_t iPtBin = 0; iPtBin < fNMaxBinsPt; iPtBin++)
            {
                fHistTrig[iMult][iZvtx][iPtBin] = NULL;
                for (Int_t jPtBin = 0; jPtBin < fNMaxBinsPt; jPtBin++)
                {
                    fHistDPhiEta[iMult][iZvtx][iPtBin][jPtBin] = NULL;
                    fHistDPhiEtaMix[iMult][iZvtx][iPtBin][jPtBin] = NULL;
                    fHistDPhiTheta[iMult][iZvtx][iPtBin][jPtBin] = NULL;
                    fHistDPhiThetaMix[iMult][iZvtx][iPtBin][jPtBin] = NULL;
                }
            }
        }
    }
}
//==========================================================================
DData::~DData()
{
    // dtor
    if (fChain)
        delete fChain->GetCurrentFile();
    if (fHistosOutFile)
        fHistosOutFile->Close();
    fHistoList->Delete();
    delete fHistoList;
    fMulLowEventPool.Delete();
    fMulHighEventPool.Delete();
    // fParticles->Delete();
    // delete fParticles;
    delete fPoolMgr;
}

//==========================================================================
void DData::Control()
{
    // Fill the histograms
    if (!fHCreated)
        CreateHistograms(kControlHisto);

    (dynamic_cast<TH1F *>(fHistoList->FindObject("ChargedMultiplicity")))->Fill(fEvent->ChargedMultiplicity());
    (dynamic_cast<TH1F *>(fHistoList->FindObject("ChargedEnergy")))->Fill(fEvent->ChargedEnergy());
    (dynamic_cast<TH1F *>(fHistoList->FindObject("Thrust")))->Fill(fEvent->Thrust());
    (dynamic_cast<TH1F *>(fHistoList->FindObject("Oblateness")))->Fill(fEvent->Oblateness());
    (dynamic_cast<TH1F *>(fHistoList->FindObject("Major")))->Fill(fEvent->ThrustMajor());
    (dynamic_cast<TH1F *>(fHistoList->FindObject("Sphericity")))->Fill(fEvent->Sphericity());
    (dynamic_cast<TH1F *>(fHistoList->FindObject("Aplanarity")))->Fill(fEvent->Aplanarity());
    TIter next(fEvent->Particles());
    Int_t index = -1;
    Double_t sump = 0, sumpz = 0.;
    while (DParticle *part = (DParticle *)next())
    {
        index++;
        (dynamic_cast<TH1F *>(fHistoList->FindObject("RPhi")))->Fill(part->VtxR());
        (dynamic_cast<TH1F *>(fHistoList->FindObject("Z")))->Fill(part->VtxZ());
	sumpz += (part->Momentum()*TMath::Cos(part->ThetaRad()));
	sump += part->Momentum();
        if (index < Npac) // charged only
        {
            (dynamic_cast<TH1F *>(fHistoList->FindObject("Momentum")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("Theta")))->Fill(part->Theta());
            if (part->Momentum() < 0.8)
                (dynamic_cast<TH1F *>(fHistoList->FindObject("Theta08")))->Fill(part->Theta());
            else if (part->Momentum() >= 0.8 && part->Momentum() < 2.0)
                (dynamic_cast<TH1F *>(fHistoList->FindObject("Theta082")))->Fill(part->Theta());
            else if (part->Momentum() >= 2.0 && part->Momentum() < 5.0)
                (dynamic_cast<TH1F *>(fHistoList->FindObject("Theta25")))->Fill(part->Theta());
            else
                (dynamic_cast<TH1F *>(fHistoList->FindObject("Theta5")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("Phi")))->Fill(part->Phi());
            TVector3 p = part->Momentum3();
            TVector3 vThrust(Thrvec1x, Thrvec1y, Thrvec1z);
            (dynamic_cast<TH1F *>(fHistoList->FindObject("pLong")))->Fill(p.Dot(vThrust));
            TVector3 vMajor(Thrvec2x, Thrvec2y, Thrvec2z);
            (dynamic_cast<TH1F *>(fHistoList->FindObject("pIn")))->Fill(p.Dot(vMajor));
            TVector3 vMinor(Thrvec3x, Thrvec3y, Thrvec3z);
            (dynamic_cast<TH1F *>(fHistoList->FindObject("pOut")))->Fill(p.Dot(vMinor));
        }
        switch (TMath::Abs(part->PdgCode()))
        {
        case kElectron:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumE")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaE")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(50.);
            break;
        case kMuon:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumMu")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaMu")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(60.);
            break;
        case kPhoton:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumPh")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaPh")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(1.);

            break;
        case kPi0:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumPi0")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaPi0")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(15.);

            break;
        case kK0L:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumK0")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaK0")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(25.);

            break;
        case kPion:
            // if (Pafl[part] != -1)  // only good identified pions
            {
                (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumPi")))->Fill(part->Momentum());
                (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaPi")))->Fill(part->Theta());
                (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(10.);
            }
            break;
        case kK0S:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumK0")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaK0")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(25.);
            break;
        case kK0:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumK0")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaK0")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(25.);
            break;
        case kK:
        {
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumK")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaK")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(20.);
        }
        break;
        case kN:
        {
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumN")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaN")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(30.);
        }
        break;
        case kP:
        {
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumN")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaN")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(30.);
        }
        break;
        case kL:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumL")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaL")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(40.);
            break;
        default:
            (dynamic_cast<TH1F *>(fHistoList->FindObject("MomentumO")))->Fill(part->Momentum());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("ThetaO")))->Fill(part->Theta());
            (dynamic_cast<TH1F *>(fHistoList->FindObject("partCom")))->Fill(99.);
            break;
        }
    }
    (dynamic_cast<TH1F *>(fHistoList->FindObject("sumpZ")))->Fill(sumpz);
    (dynamic_cast<TH2F *>(fHistoList->FindObject("sumpZp")))->Fill(sumpz,sump);
}

//==========================================================================
void DData::Correlation()
{
    // produce 2 particle correla~tion

    const Int_t kPoolMixingMethod = 2;
    if (!fHCreated)
    {
        CreateHistograms(kCorrelation);
    }
    TH2F *hSL = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiSL"));
    TH2F *hML = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiML"));
    TH2F *hCL = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiCL"));
    TH1F *hELL = dynamic_cast<TH1F *>(fHistoList->FindObject("ELeadingL"));

    TH2F *hSH = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiSH"));
    TH2F *hMH = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiMH"));
    TH2F *hCH = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiCH"));
    TH1F *hELH = dynamic_cast<TH1F *>(fHistoList->FindObject("ELeadingH"));

    Int_t max = Leading();
    (dynamic_cast<TH1F *>(fHistoList->FindObject("Mul")))->Fill(fMultiplicity);

    DParticle *leader = static_cast<DParticle *>(fEvent->Particles()->At(max));
    if (fMultiplicity < fMulLow)
    {
        hELL->Fill(leader->Energy());
        fTriggersL++;
    }
    else if (fMultiplicity > fMulHigh)
    {
        hELH->Fill(leader->Energy());
        fTriggersH++;
    }
    Double_t etaL = leader->Eta();
    Double_t phiL = leader->Phi();

    // Make signal
    Int_t index = 0;
    TIter next(fEvent->Particles());
    while (DParticle *part = (DParticle *)next())
    {
        if (index != max)
        {
            if (fMultiplicity < fMulLow)
            {
                hSL->Fill(part->Eta() - etaL, (part->Phi() - phiL) * TMath::DegToRad());
                hCL->Fill(part->Eta() - etaL, (part->Phi() - phiL) * TMath::DegToRad());
            }
            else if (fMultiplicity > fMulHigh)
            {
                hSH->Fill(part->Eta() - etaL, (part->Phi() - phiL) * TMath::DegToRad());
                hCH->Fill(part->Eta() - etaL, (part->Phi() - phiL) * TMath::DegToRad());
            }
        }
        index++;
    }

    // make 2-particle correlations including event mixing
    Int_t multBin = GetMultBin();
    if (multBin < 0)
        return;

    // event z vertex
    Double_t fzvtx = fEvent->BeamSpotL(); 
    if (fzvtx >= fZvtxAxis->GetXmax())
        return;
    if (fzvtx <= fZvtxAxis->GetXmin())
        return;
    Int_t zvtxBin = fZvtxAxis->FindBin(fzvtx) - 1;

    DPool *poolTrk = fPoolMgr->GetEventPool(-fMultiplicity, 0); // for the moment no z vtx binning
    Bool_t poolReady = poolTrk->IsReady();
    // same event
    for (Int_t iTrack = 0; iTrack < fEvent->Particles()->GetEntriesFast(); iTrack++)
    {
        DParticle *trigTrack = (DParticle *)fEvent->Particles()->UncheckedAt(iTrack);
	if (trigTrack->Eta()<=-1.7 ||
	    trigTrack->Eta()>= 1.7) continue;
        if (trigTrack->VtxR() > 0.1 ||
            trigTrack->VtxZ() > 1.0)
            continue; // only primary tracks

	//        Int_t ptBinTrk = fTrackPtAxis->FindBin(trigTrack->Pt());
        Int_t ptBinTrk = fTrackPtAxis->FindBin(trigTrack->Momentum());
        if (ptBinTrk < 1 || ptBinTrk > fNbinsTrackPt)
        {
            continue;
        }
        fHistTrig[multBin][zvtxBin][ptBinTrk - 1]->Fill(0.5);

        for (Int_t jTrack = 0; jTrack < fEvent->Particles()->GetEntriesFast(); jTrack++)
        {
            DParticle *assocTrack = (DParticle *)fEvent->Particles()->UncheckedAt(jTrack);
            if (jTrack == iTrack)
                continue;
	    if (assocTrack->Eta()<=-1.7 ||
		assocTrack->Eta()>= 1.7) continue;
            if (assocTrack->VtxR() > 0.1 ||
                assocTrack->VtxZ() > 1.0)
                continue; // only primary tracks

            Int_t ptBinAssocTrk = fTrackPtAxis->FindBin(assocTrack->Momentum());
            if (ptBinAssocTrk < 1 || ptBinAssocTrk > fNbinsTrackPt)
            {
                continue;
            }

            Double_t deltaPhi = assocTrack->PhiRad() - trigTrack->PhiRad();
            if (deltaPhi > 1.5 * TMath::Pi())
                deltaPhi -= TMath::TwoPi();
            if (deltaPhi < -0.5 * TMath::Pi())
                deltaPhi += TMath::TwoPi();

	    Double_t deltaTheta = trigTrack->ThetaRad() - assocTrack->ThetaRad();
            if (deltaTheta > 1.5 * TMath::Pi())
                deltaTheta -= TMath::TwoPi();
            if (deltaTheta < -0.5 * TMath::Pi())
                deltaTheta += TMath::TwoPi();
	    Double_t deltaEta = assocTrack->Eta() - trigTrack->Eta();
	      
	    Double_t dphi,dtheta,deta;
	    if (!RotateTracks(trigTrack,assocTrack,dphi,dtheta,deta))
	      printf("Rotation of tracks failed...\n");
	    
	    // printf("DEBUG: dphi: %f %f    dtheta: %f vs %f\n",
	    // 	     deltaPhi,dphi,deltaTheta,dtheta);
	    if (0) {
	      fHistDPhiEta[multBin][zvtxBin][ptBinTrk - 1][ptBinAssocTrk - 1]->Fill(deltaPhi,deltaEta);
	      fHistDPhiTheta[multBin][zvtxBin][ptBinTrk - 1][ptBinAssocTrk - 1]->Fill(deltaPhi,deltaTheta);
	    }
	    else {
	      fHistDPhiEta[multBin][zvtxBin][ptBinTrk - 1][ptBinAssocTrk - 1]->Fill(dphi,deta);
	      fHistDPhiTheta[multBin][zvtxBin][ptBinTrk - 1][ptBinAssocTrk - 1]->Fill(dphi,dtheta);
	    }
        }
        // mixed event
        if (poolReady)
        {
            for (Int_t jMix = 0; jMix < poolTrk->GetCurrentNEvents(); jMix++)
            {
                TClonesArray *mixedTrk = poolTrk->GetEvent(jMix);
                for (Int_t jTrack = 0; jTrack < mixedTrk->GetEntriesFast(); jTrack++)
                {
                    DParticle *assocTrack = (DParticle *)mixedTrk->UncheckedAt(jTrack);
                    if (assocTrack->Eta()<=-1.7 ||
                     	assocTrack->Eta()>= 1.7) continue;
                    if (assocTrack->VtxR() > 0.1 ||
                        assocTrack->VtxZ() > 1.0)
                        continue; // only primary tracks

                    Int_t ptBinAssocTrk = fTrackPtAxis->FindBin(assocTrack->Momentum());
                    if (ptBinAssocTrk < 1 || ptBinAssocTrk > fNbinsTrackPt)
                    {
                        continue;
                    }

                    Double_t deltaPhi = assocTrack->PhiRad() - trigTrack->PhiRad();
                    if (deltaPhi > 1.5 * TMath::Pi())
                        deltaPhi -= TMath::TwoPi();
                    if (deltaPhi < -0.5 * TMath::Pi())
                        deltaPhi += TMath::TwoPi();

		    Double_t deltaTheta = trigTrack->ThetaRad() - assocTrack->ThetaRad();
		    if (deltaTheta > 1.5 * TMath::Pi())
		      deltaTheta -= TMath::TwoPi();
		    if (deltaTheta < -0.5 * TMath::Pi())
		      deltaTheta += TMath::TwoPi();
		    Double_t deltaEta = assocTrack->Eta() - trigTrack->Eta();
		    
		    Double_t dphi,dtheta,deta;
		    if (!RotateTracks(trigTrack,assocTrack,dphi,dtheta,deta))
		      printf("Rotation of mixed tracks failed...\n");

		    if (0) {
		      fHistDPhiEtaMix[multBin][zvtxBin][ptBinTrk - 1][ptBinAssocTrk - 1]->Fill(deltaPhi,deltaEta);
		      fHistDPhiThetaMix[multBin][zvtxBin][ptBinTrk - 1][ptBinAssocTrk - 1]->Fill(deltaPhi,deltaTheta);
		    }
		    else {
		      fHistDPhiEtaMix[multBin][zvtxBin][ptBinTrk - 1][ptBinAssocTrk - 1]->Fill(dphi,deta);
		      fHistDPhiThetaMix[multBin][zvtxBin][ptBinTrk - 1][ptBinAssocTrk - 1]->Fill(dphi,dtheta);
		    }
                }
            }
        }
    }

    poolTrk->UpdatePool(fEvent->Particles());
}

//==========================================================================
void DData::CreateHistograms(Eopt opt, Eparam par)
{
    // open output file and create histograms

    if (!fHistoList)
        fHistoList = new TList();
    fHistoList->SetOwner(kTRUE);

    if (opt == kControlHisto || opt == kControlDiff)
    {
        fHCreated = kTRUE;
        fHistoList->Add(new TH1F("ChargedMultiplicity", "Charged multiplicity", 100, 0.0, 100.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h4")));

        fHistoList->Add(new TH1F("ChargedEnergy", "Charged energy / ECM", 100, 0.0, 1.2));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h5")));

        fHistoList->Add(new TH1F("Thrust", "Event Thrust", 100, 0.5, 1.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h7")));

        fHistoList->Add(new TH1F("Oblateness", "Event Oblateness", 100, 0.0, 1.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h8")));

        fHistoList->Add(new TH1F("Major", "Event Major", 100, 0.0, 1.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h9")));

        fHistoList->Add(new TH1F("Sphericity", "Event Sphericity", 100, 0.0, 1.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h10")));

        fHistoList->Add(new TH1F("Aplanarity", "Event Aplanarity", 100, 0.0, 1.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h11")));

        fHistoList->Add(new TH1F("RPhi", "Impact R/Phi(pT)", 100, -5.0, 5.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h20")));

        fHistoList->Add(new TH1F("Z", "Impact Z(pT)", 110, -11.0, 11.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h21")));

        fHistoList->Add(new TH1F("Momentum", "p[GeV]", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h22")));

        fHistoList->Add(new TH1F("MomentumPh", "p[GeV] photon", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h220")));

        fHistoList->Add(new TH1F("MomentumPi", "p[GeV] pion", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h221")));

        fHistoList->Add(new TH1F("MomentumPi0", "p[GeV] pi0", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h222")));

        fHistoList->Add(new TH1F("MomentumK", "p[GeV] K", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h223")));

        fHistoList->Add(new TH1F("MomentumK0", "p[GeV] K0", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h224")));

        fHistoList->Add(new TH1F("MomentumN", "p[GeV] N", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h225")));

        fHistoList->Add(new TH1F("MomentumL", "p[GeV] L", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h226")));

        fHistoList->Add(new TH1F("MomentumE", "p[GeV] e", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h227")));

        fHistoList->Add(new TH1F("MomentumMu", "p[GeV] mu", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h228")));

        fHistoList->Add(new TH1F("MomentumO", "p[GeV] other", 500, 0.0, 10.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h229")));

        fHistoList->Add(new TH1F("Theta", "Theta ", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h25")));

        fHistoList->Add(new TH1F("ThetaPh", "Theta photon", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h250")));

        fHistoList->Add(new TH1F("ThetaPi", "Theta pion", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h251")));

        fHistoList->Add(new TH1F("ThetaPi0", "Theta pi0", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h252")));

        fHistoList->Add(new TH1F("ThetaK", "Theta K", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h253")));

        fHistoList->Add(new TH1F("ThetaK0", "Theta K0", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h254")));

        fHistoList->Add(new TH1F("ThetaN", "Theta N", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h255")));

        fHistoList->Add(new TH1F("ThetaL", "Theta L", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h256")));

        fHistoList->Add(new TH1F("ThetaE", "Theta e", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h257")));

        fHistoList->Add(new TH1F("ThetaMu", "Theta mu", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h258")));

        fHistoList->Add(new TH1F("ThetaO", "Theta other", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h259")));

        fHistoList->Add(new TH1F("Theta08", "Theta p < 0.8 GeV", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h260")));

        fHistoList->Add(new TH1F("Theta082", "Theta 0.8 < p < 2 GeV", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h261")));

        fHistoList->Add(new TH1F("Theta25", "Theta 2 < p < 5 GeV", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h262")));

        fHistoList->Add(new TH1F("Theta5", "Theta p > 5 GeV", 180, 0.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h263")));

        fHistoList->Add(new TH1F("Phi", "Phi", 360, -180.0, 180.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h27")));

        fHistoList->Add(new TH1F("pLong", "p Long [GeV]", 200, 0.0, 50.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h28")));

        fHistoList->Add(new TH1F("pIn", "pt In [GeV]", 200, 0.0, 20.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h29")));

        fHistoList->Add(new TH1F("pOut", "pt Out [GeV]", 200, 0.0, 5.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h30")));

        fHistoList->Add(new TH1F("partCom", "particles composition", 100, 0.0, 100.0));
        fHistoList->Add(dynamic_cast<TH1F *>(fInputFile->Get("h32")));

        fHistoList->Add(new TH1F("sumpZ", "Sum(pZ) [GeV]", 1000, -200.0, 200.0));
        fHistoList->Add(new TH2F("sumpZp", "Sum(p) vs Sum(pZ) [GeV]", 100, -200.0, 200.0, 150, 0., 300.));
    }
    else if (opt == kCorrelation)
    {
        fHCreated = kTRUE;
        fHistoList->Add(new TH1F("Mul", "Particle multiplicity", 80, 0., 80.));

        // Low multiplicity
        TH2F *hsl = new TH2F("EtaPhiSL", "Signal Low M", 100, -3., 3., 100, -1.5, 4.5);
        hsl->GetXaxis()->SetTitle("#Delta#eta");
        hsl->GetYaxis()->SetTitle("#Delta#Phi");
        hsl->GetZaxis()->SetTitle("#frac{1}{N_{trigger}}#times #frac{d^2N}{d#eta d#Phi}");
        hsl->SetStats(kFALSE);
        fHistoList->Add(hsl);
        TH2F *hml = new TH2F("EtaPhiML", "Mixed Low M", 100, -3., 3., 100, -1.5, 4.5);
        hml->GetXaxis()->SetTitle("#Delta#eta");
        hml->GetYaxis()->SetTitle("#Delta#Phi");
        hml->GetZaxis()->SetTitle("#frac{1}{N_{trigger}}#times #frac{d^2N}{d#eta d#Phi}");
        hml->SetStats(kFALSE);
        fHistoList->Add(hml);
        TH2F *hcl = new TH2F("EtaPhiCL", "Correlation Low M", 100, -3., 3., 100, -1.5, 4.5);
        hcl->GetXaxis()->SetTitle("#Delta#eta");
        hcl->GetYaxis()->SetTitle("#Delta#Phi");
        hcl->GetZaxis()->SetTitle("#frac{1}{N_{trigger}}#times #frac{d^{2}N}{d#eta d#Phi}");
        hcl->SetStats(kFALSE);
        fHistoList->Add(hcl);
        fHistoList->Add(new TH1F("ELeadingL", "Leading energy Low M", 100, 0., 100.));

        // High multiplicity
        TH2F *hsh = new TH2F("EtaPhiSH", "Signal High M", 100, -3., 3., 100, -1.5, 4.5);
        hsh->GetXaxis()->SetTitle("#Delta#eta");
        hsh->GetYaxis()->SetTitle("#Delta#Phi");
        hsh->GetZaxis()->SetTitle("#frac{1}{N_{trigger}}#times #frac{d^2N}{d#eta d#Phi}");
        hsh->SetStats(kFALSE);
        fHistoList->Add(hsh);
        TH2F *hmh = new TH2F("EtaPhiMH", "Mixed High M", 100, -3., 3., 100, -1.5, 4.5);
        hmh->GetXaxis()->SetTitle("#Delta#eta");
        hmh->GetYaxis()->SetTitle("#Delta#Phi");
        hmh->GetZaxis()->SetTitle("#frac{1}{N_{trigger}}#times #frac{d^2N}{d#eta d#Phi}");
        hmh->SetStats(kFALSE);
        fHistoList->Add(hmh);
        TH2F *hch = new TH2F("EtaPhiCH", "Correlation High M", 100, -3., 3., 100, -1.5, 4.5);
        hch->GetXaxis()->SetTitle("#Delta#eta");
        hch->GetYaxis()->SetTitle("#Delta#Phi");
        hch->GetZaxis()->SetTitle("#frac{1}{N_{trigger}}#times #frac{d^{2}N}{d#eta d#Phi}");
        hch->SetStats(kFALSE);
        fHistoList->Add(hch);
        fHistoList->Add(new TH1F("ELeadingH", "Leading energy High M", 100, 0., 100.));

        fHistoList->Add(new TH1I("mulL", "mul low", 1000, 0, 100000));
        fHistoList->Add(new TH1I("mulH", "mul high", 1000, 0, 100000));

        // create 2-particle correlation histos
        for (Int_t iMult = 0; iMult < fNbinsMult; iMult++)
        {
            for (Int_t iZvtx = 0; iZvtx < fNbinsZvtx; ++iZvtx)
            {
                for (Int_t iPtBin = 0; iPtBin < fNbinsTrackPt; iPtBin++)
                {
                    fHistTrig[iMult][iZvtx][iPtBin] = new TH1D(Form("fHistTrig_Mult%02d_Z%02d_PtBin%02d", iMult, iZvtx, iPtBin),
                                                               "", 1, 0, 1);
                    fHistoList->Add(fHistTrig[iMult][iZvtx][iPtBin]);
                    for (Int_t jPtBin = 0; jPtBin < fNbinsTrackPt; jPtBin++)
                    {
                        fHistDPhiEta[iMult][iZvtx][iPtBin][jPtBin] = new TH2D(Form("fHistDPhiEta_Mult%02d_Z%02d_PtBin%02d_%02d", iMult, iZvtx, iPtBin, jPtBin),
                                                                              "",
                                                                              48, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(),
                                                                              64, -4., 4.);
                        fHistoList->Add(fHistDPhiEta[iMult][iZvtx][iPtBin][jPtBin]);
                        fHistDPhiEtaMix[iMult][iZvtx][iPtBin][jPtBin] = new TH2D(Form("fHistDPhiEtaMix_Mult%02d_Z%02d_PtBin%02d_%02d", iMult, iZvtx, iPtBin, jPtBin),
                                                                                 "",
                                                                                 48, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(),
                                                                                 64, -4., 4.);
                        fHistoList->Add(fHistDPhiEtaMix[iMult][iZvtx][iPtBin][jPtBin]);

                        fHistDPhiTheta[iMult][iZvtx][iPtBin][jPtBin] = new TH2D(Form("fHistDPhiTheta_Mult%02d_Z%02d_PtBin%02d_%02d", iMult, iZvtx, iPtBin, jPtBin),
										"",
										48, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(),
										48, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
                        fHistoList->Add(fHistDPhiTheta[iMult][iZvtx][iPtBin][jPtBin]);
                        fHistDPhiThetaMix[iMult][iZvtx][iPtBin][jPtBin] = new TH2D(Form("fHistDPhiThetaMix_Mult%02d_Z%02d_PtBin%02d_%02d", iMult, iZvtx, iPtBin, jPtBin),
                                                                                 "",
                                                                                 48, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(),
                                                                                 48, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
                        fHistoList->Add(fHistDPhiThetaMix[iMult][iZvtx][iPtBin][jPtBin]);
                    }
                }
            }
        }
        fPoolMgr = new DPoolManager("pool", 1, fMixing, fNbinsMult, (Double_t *)fMultAxis->GetXbins()->GetArray(), fNbinsZvtx, (Double_t *)fZvtxAxis->GetXbins()->GetArray());
    }
    else if (opt == kSingleHisto)
    {
        fHCreated = kTRUE;
        switch (par)
        {
        case kThrust:
        {
            TH1F *histo = new TH1F("SingleHisto", "Event Thrust", 100, 0.6, 1.0);
            histo->GetXaxis()->SetTitle("T");
            histo->GetYaxis()->SetTitle("#frac{1}{N_{trigger}}#times #frac{dN}{dT}");
            fHistoList->Add(histo);
            break;
        }
        case kMultiplicity:
        {
            TH1F *histo = new TH1F("SingleHisto", "Event Multiplicity", 100, 0., 100.);
            histo->GetXaxis()->SetTitle("M");
            histo->GetYaxis()->SetTitle("#frac{1}{N_{trigger}}#times #frac{dN}{dM}");
            fHistoList->Add(histo);
            break;
        }
        default:
            break;
        }
    }
}

//==========================================================================
Int_t DData::Cut(Long64_t /*entry*/)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}

//==========================================================================
Int_t DData::GetEntry(Long64_t entry)
{
    // Read contents of event #entry [0,fChain->GetEntriesFast()].
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}

//==========================================================================
void DData::Init()
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!fChain)
        return;
    fCurrentTree = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("Ecm", &Ecm, &b_Ecm);
    fChain->SetBranchAddress("Bt", &Bt, &b_Bt);
    fChain->SetBranchAddress("Bz", &Bz, &b_Bz);
    fChain->SetBranchAddress("Sphvec1x", &Sphvec1x, &b_Sphvec1x);
    fChain->SetBranchAddress("Sphvec1y", &Sphvec1y, &b_Sphvec1y);
    fChain->SetBranchAddress("Sphvec1z", &Sphvec1z, &b_Sphvec1z);
    fChain->SetBranchAddress("Sphvec2x", &Sphvec2x, &b_Sphvec2x);
    fChain->SetBranchAddress("Sphvec2y", &Sphvec2y, &b_Sphvec2y);
    fChain->SetBranchAddress("Sphvec2z", &Sphvec2z, &b_Sphvec2z);
    fChain->SetBranchAddress("Sphvec3x", &Sphvec3x, &b_Sphvec3x);
    fChain->SetBranchAddress("Sphvec3y", &Sphvec3y, &b_Sphvec3y);
    fChain->SetBranchAddress("Sphvec3z", &Sphvec3z, &b_Sphvec3z);
    fChain->SetBranchAddress("Sphval1", &Sphval1, &b_Sphval1);
    fChain->SetBranchAddress("Sphval2", &Sphval2, &b_Sphval2);
    fChain->SetBranchAddress("Sphval3", &Sphval3, &b_Sphval3);
    fChain->SetBranchAddress("Thrvec1x", &Thrvec1x, &b_Thrvec1x);
    fChain->SetBranchAddress("Thrvec1y", &Thrvec1y, &b_Thrvec1y);
    fChain->SetBranchAddress("Thrvec1z", &Thrvec1z, &b_Thrvec1z);
    fChain->SetBranchAddress("Thrvec2x", &Thrvec2x, &b_Thrvec2x);
    fChain->SetBranchAddress("Thrvec2y", &Thrvec2y, &b_Thrvec2y);
    fChain->SetBranchAddress("Thrvec2z", &Thrvec2z, &b_Thrvec2z);
    fChain->SetBranchAddress("Thrvec3x", &Thrvec3x, &b_Thrvec3x);
    fChain->SetBranchAddress("Thrvec3y", &Thrvec3y, &b_Thrvec3y);
    fChain->SetBranchAddress("Thrvec3z", &Thrvec3z, &b_Thrvec3z);
    fChain->SetBranchAddress("Thrval1", &Thrval1, &b_Thrval1);
    fChain->SetBranchAddress("Thrval2", &Thrval2, &b_Thrval2);
    fChain->SetBranchAddress("Thrval3", &Thrval3, &b_Thrval3);
    fChain->SetBranchAddress("Npa", &Npa, &b_Npa);
    fChain->SetBranchAddress("Npac", &Npac, &b_Npac);
    fChain->SetBranchAddress("Paid", Paid, &b_Paid);
    fChain->SetBranchAddress("Pafl", Pafl, &b_Pafl);
    fChain->SetBranchAddress("Pav0d", Pav0d, &b_Pav0d);
    fChain->SetBranchAddress("Pav0m", Pav0m, &b_Pav0m);
    fChain->SetBranchAddress("Rvtx", Rvtx, &b_Rvtx);
    fChain->SetBranchAddress("Zvtx", Zvtx, &b_Zvtx);
    fChain->SetBranchAddress("Papx", Papx, &b_Papx);
    fChain->SetBranchAddress("Papy", Papy, &b_Papy);
    fChain->SetBranchAddress("Papz", Papz, &b_Papz);
    fChain->SetBranchAddress("Paje", Paje, &b_Paje);
    fChain->SetBranchAddress("Pani", Pani, &b_Pani);
    fChain->SetBranchAddress("Njer", &Njer, &b_Njer);
    fChain->SetBranchAddress("Tgenr", &Tgenr, &b_Tgenr);
    fChain->SetBranchAddress("Dminr", &Dminr, &b_Dminr);
    fChain->SetBranchAddress("Jepx", Jepx, &b_Jepx);
    fChain->SetBranchAddress("Jepy", Jepy, &b_Jepy);
    fChain->SetBranchAddress("Jepz", Jepz, &b_Jepz);
    fChain->SetBranchAddress("Jee", Jee, &b_Jee);
    fChain->SetBranchAddress("Jem", Jem, &b_Jem);
    fChain->SetBranchAddress("Jep", Jep, &b_Jep);
    Notify();

    fPMasses[0] = -1.0;              // unidentified particle
    fPMasses[11] = 0.0005109989461;  // electron
    fPMasses[-11] = 0.0005109989461; // positron
    fPMasses[13] = 0.1056583745;     // mu+
    fPMasses[-13] = 0.1056583745;    // mu-
    fPMasses[22] = 0.0;              // photon
    fPMasses[111] = 0.1349770;       // pi0
    fPMasses[130] = 0.497611;        // K0L
    fPMasses[211] = 0.13957061;      // pi+
    fPMasses[-211] = 0.13957061;     // pi-
    fPMasses[310] = 0.497611;        // K0S
    fPMasses[311] = 0.497611;        // K0
    fPMasses[321] = 0.493677;        // K+
    fPMasses[-321] = 0.493677;       // K-
    fPMasses[2121] = 0.9395654133;   // neutron
    fPMasses[2212] = 0.9382720813;   // proton
    fPMasses[-2212] = 0.9382720813;  // anti proton
    fPMasses[3122] = 1.115683;       // lambda

    fChain->GetEntry(0);
    fEvent = new DEvent();
}

//==========================================================================
Bool_t DData::IsCharged(Int_t index) const
{
    // tells if particle is charged
    Bool_t rv = kFALSE;
    switch (TMath::Abs(Paid[index]))
    {
    case kMuon:
        rv = kTRUE;
        break;
    case kElectron:
        rv = kTRUE;
        break;
    case kPion:
        rv = kTRUE;
        break;
    case kK:
        rv = kTRUE;
        break;
    case kP:
        rv = kTRUE;
        break;
    default:
        rv = kFALSE;
        break;
    }
    return rv;
}

//==========================================================================
Bool_t DData::IsChargedHadron(Int_t index) const
{
    // tells if particle is a charged hadron
    Bool_t rv = kFALSE;
    switch (TMath::Abs(Paid[index]))
    {
    case kPion:
        rv = kTRUE;
        break;
    case kK:
        rv = kTRUE;
        break;
    case kP:
        rv = kTRUE;
        break;
    default:
        rv = kFALSE;
        break;
    }
    return rv;
}

//==========================================================================
Bool_t DData::IsHadron(Int_t index) const
{
    // tells if particle is a charged hadron
    Bool_t rv = kFALSE;

    switch (TMath::Abs(Paid[index]))
    {
    case kPi0:
        rv = kTRUE;
        break;
    case kK0L:
        rv = kTRUE;
        break;
    case kPion:
        rv = kTRUE;
        break;
    case kK0S:
        rv = kTRUE;
        break;
    case kK:
        rv = kTRUE;
        break;
    case kK0:
        rv = kTRUE;
        break;
    case kN:
        rv = kTRUE;
        break;
    case kP:
        rv = kTRUE;
        break;
    case kL:
        rv = kTRUE;
        break;
    default:
        rv = kFALSE;
        break;
    }
    return rv;
}

//==========================================================================
Int_t DData::Leading()
{
    // search index of leading charged particle
    Int_t rv = -1;
    Int_t index = 0;
    Double_t pmax = 0.;
    TIter next(fEvent->Particles());
    while (DParticle *part = (DParticle *)next())
    {
        if (part->Momentum() > pmax)
        {
            pmax = part->Momentum();
            rv = index;
        }
        index++;
    }
    return rv;
}

//==========================================================================
Long64_t DData::LoadEvent(Long64_t jentry)
{
    // fill event data
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
        return -1;
    Long64_t nb = 0;
    nb = fChain->GetEntry(jentry);
    if (Cut(ientry) < 0)
        return 0;
    fEvent->Clear();
    fEvent->SetECM(Ecm);
    fEvent->SetEventNumber(jentry);
    fEvent->SetBeamSpot(Bt, Bz);
    fEvent->SetMul(Npa, Npac);
    char name[50];
    sprintf(name, "entry # %12lld ", jentry);
    fEvent->SetName(name);
    sprintf(name, "particle selection: %s", fPopt[fCurrentPartSel].Data());
    fEvent->SetTitle(name);
    fEvent->SetSpericity(Sphval1, Sphval2, Sphval3);
    fEvent->SetThrust(Thrval1, Thrval2, Thrval3);
    fEvent->SetJet(Njer, Dminr); 
    MakeParticlesList(fCurrentPartSel);
    // if (fMultiplicity >= fMulLow && fMultiplicity <= fMulHigh)
    //    return 0;
    return fMultiplicity;
}

//==========================================================================
Int_t DData::LoadEventMultiplicity(Long64_t jentry)
{
    // fill event multiplicity
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
        return -1;
    Long64_t nb = 0;
    nb = fChain->GetEntry(jentry);
    if (Cut(ientry) < 0)
        return 0;

    Int_t mul = 0;
    Int_t study = 0;
    if (fCurrentPartSel == kAll || fCurrentPartSel == kHadrons)
        study = Npa;
    else
        study = Npac;
    Int_t pindex = 0;
    for (Int_t index = 0; index < study; index++)
    {
        if ((fCurrentPartSel == kHadrons || fCurrentPartSel == kChargedHadrons) && !IsHadron(index))
            continue;
        else
            mul++;
    }
    return mul - 1;
}

//==========================================================================
Long64_t DData::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
        return centry;
    if (fChain->GetTreeNumber() != fCurrentTree)
    {
        fCurrentTree = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

//==========================================================================
void DData::Loop(const Eopt opt, const Epopt popt, const Eparam par)
{
    //   In a ROOT session, you can do:
    //      root> .L DData.C
    //      root> DData data
    //      root> data.GetEntry(12); // Fill t data members with entry number 12
    //      root> data.Show();       // Show values of entry 12
    //      root> data.Show(16);     // Read and show values of entry 16
    //      root> data.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch

    if (fChain == 0)
        return;

    printf("INFO: Loop(%i, %i)\n", opt, popt);
    fCurrentPartSel = popt;
    fEvents = fChain->GetEntriesFast();
    Long64_t nbytes = 0;
    for (Long64_t jentry = 0; jentry < fEvents; jentry++)
    //    for (Long64_t jentry = 0; jentry < 10; jentry++)
    {
      if (jentry%10000 == 0) printf("Event %lld\n",jentry);
        fCurrentEvent = jentry;
        Long64_t rv = LoadEvent(jentry);

	// Event selection cuts
	if (((Bt-0.33)*(Bt-0.33)/0.025/0.025+(Bz+0.71)*(Bz+0.71)/0.76/0.76)>9.) continue;
	
        if (rv == -1)
            break;
        else if (rv == 0)
            continue;
        else
        {
            nbytes += rv;

            switch (opt)
            {
            case kControlHisto:
                Control();
                break;
            case kControlDiff:
                Control();
                break;
            case kCorrelation:
                Correlation();
                break;
            case kSingleHisto:
                SingleHisto(par);
                break;
            default:
                break;
            }
        }
    }
    Plot(opt);
}

//==========================================================================
void DData::MakeMulPool(Int_t method)
{
    // Make event pool for 2 multiplicity bin: < lowM and > highM
    std::cout << "INFO: MakeMulPool --> start " << std::endl;
    auto start = std::chrono::system_clock::now();
    switch (method)
    {
    case 1:
    {
        const Int_t kMix = 2;
        Int_t lowCounter = 0, highCounter = 0;
        while (lowCounter < fMixing * kMix || highCounter < fMixing * kMix)
        {
            Int_t eventIndex = gRandom->Integer(fEvents);
            Int_t mul = LoadEvent(eventIndex);
            if (mul < fMulLow && lowCounter < fMixing * kMix)
            {
                fMulLowEventPool.Add(fEvent);
                lowCounter++;
            }
            else if (mul > fMulHigh && highCounter < fMixing * kMix)
            {
                fMulHighEventPool.Add(fEvent);
                highCounter++;
            }
        }
        break;
    }
    case 2:
    {
        fChain->SetBranchStatus("*", 0); //disable all branches
        fChain->SetBranchStatus("Npa", 1);
        fChain->SetBranchStatus("Npac", 1);
        fChain->SetBranchStatus("Paid", 1);

        // store event numbers with low and high multiplicity
        Long64_t eventIndexMax = fEvents;
        for (Int_t eventIndex = 0; eventIndex < eventIndexMax; eventIndex++)
        {
            Int_t mul = LoadEventMultiplicity(eventIndex);
            if (mul < fMulLow)
                fMulLowIndexPool.push_back(eventIndex);
            else if (mul > fMulHigh)
                fMulHighIndexPool.push_back(eventIndex);
        }
        fChain->SetBranchStatus("*", 1); //enable all branches
        break;
    }
    default:
        break;
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    printf("INFO: MakeMulPool --> done in %f seconds\n", elapsed_seconds.count());
}

//==========================================================================
void DData::MakeParticlesList(Epopt opt)
{
    // makes a list of selectected particles opt=kAll, kCharged, kHadrons, kChargedHadrons
    if (IsVerbose())
        std::cout << "INFO: MakeParticlesList --> " << fEvent->GetName() << " " << fEvent->GetTitle() << std::endl;
    if (fEvent->Particles()->GetEntriesFast() != 0)
        fEvent->Particles()->Clear("C");
    Int_t study = 0;
    if (opt == kAll || opt == kHadrons)
        study = Npa;
    else
        study = Npac;
    Int_t pindex = 0;
    Int_t pindexJetVeto = 0;
    fMultiplicity = 0;
    for (Int_t index = 0; index < study; index++)
    {
        if ( (TMath::Abs(Rvtx[index]) > 0.1 || TMath::Abs(Zvtx[index]) > 1.0) || Pav0d[index] != 0 || Pani[index] == kTRUE)
            continue; 
        if ((opt == kHadrons || opt == kChargedHadrons) && !IsHadron(index))
            continue;
        DParticle *part = (DParticle *)fEvent->Particles()->ConstructedAt(pindex++);
        part->Set(Paid[index], Papx[index], Papy[index], Papz[index], Mass(Paid[index]), Pav0d[index], Rvtx[index], Zvtx[index]);
        // if (Pani[index])
        //     std::cout << index + 1 << " " << Paid[index] << " " << Rvtx[index] << " " << 
        //                      Pav0d[index] << " " << Pav0m[index] << " true" << std::endl;  
        // else 
        //     std::cout << index + 1 << " " << Paid[index] << " " << Rvtx[index] << " " << 
        //                      Pav0d[index] << " " << Pav0m[index] << " false" << std::endl;  
        if (Paje[index] == 0) pindexJetVeto++;
    }
    fMultiplicity = pindex;
    //    fMultiplicity = pindexJetVeto;
}

// //==========================================================================
Bool_t DData::Notify()
{
    //    // The Notify() function is called when a new file is opened. This~
    //    // can be either for a new TTree in a TChain or when when a new TTree
    //    // is started when using PROOF. It is normally not necessary to make changes
    //    // to the generated code, but the routine can be extended by the
    //    // user if needed. The return value is currently not used.

    return kTRUE;
}

//==========================================================================
void DData::Plot(const Eopt opt) const
{
    // plot histograms
    Int_t nHistos = TMath::Sqrt(fHistoList->GetSize() / 2.0) + 1;

    TCanvas *can = new TCanvas("can", GetTitle());
    can->SetCanvasSize(1000, 1000);
    can->SetWindowSize(500, 500);
    can->Divide(nHistos, nHistos);

    switch (opt)
    {
    case kControlHisto:
    {
        TIter next(fHistoList);
        Int_t npad = 1;
        while (TObject *h = next())
        {
            TH1F *ho = dynamic_cast<TH1F *>(h);
            TH1F *hi = dynamic_cast<TH1F *>(next());
            can->cd(npad++)->SetLogy();
            if (hi) hi->Draw("hist");
            ho->SetLineColor(2);
            ho->Draw(hi?"hist same":"hist");
        }
        break;
    }
    case kControlDiff:
    {
        TIter next(fHistoList);
        Int_t npad = 1;
        while (TObject *h = next())
        {
            TH1F *ho = dynamic_cast<TH1F *>(h);
            TH1F *hi = dynamic_cast<TH1F *>(next());
            can->cd(npad++)->SetLogy(kFALSE);
            TString hname(ho->GetName());
            hname.Prepend("diff_");
            TH1F *hio = static_cast<TH1F *>(ho->Clone(hname));
            hio->Add(ho, hi, 1.0, -1.0);
            hio->Divide(hio, hi, 1., 1.);
            hio->Draw("hist");
        }
        break;
    }
    case kCorrelation:
    {
        // Low M
        TH2F *hSL = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiSL"));
        TH2F *hML = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiML"));
        hML->Scale(1.0 / fMixing);
        TH2F *hCL = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiCL"));
        hCL->Divide(hML);
        hCL->Scale(1.0 / fTriggersL);
        can->cd(1);
        hSL->Draw();
        can->cd(2);
        hML->Draw();
        can->cd(3);
        hCL->Draw();
        can->cd(4);
        (dynamic_cast<TH1F *>(fHistoList->FindObject("ELeadingL")))->Draw();

        // High M
        TH2F *hSH = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiSH"));
        TH2F *hMH = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiMH"));
        hMH->Scale(1.0 / fMixing);
        TH2F *hCH = dynamic_cast<TH2F *>(fHistoList->FindObject("EtaPhiCH"));
        hCH->Divide(hMH);
        hCH->Scale(1.0 / fTriggersH);
        can->cd(5);
        hSH->Draw();
        can->cd(6);
        hMH->Draw();
        can->cd(7);
        hCH->Draw();
        can->cd(8);
        (dynamic_cast<TH1F *>(fHistoList->FindObject("ELeadingH")))->Draw();
        can->cd(9);
        (dynamic_cast<TH1F *>(fHistoList->FindObject("Mul")))->Draw();
        break;
    }
    case kSingleHisto:
    {
        can->cd(1);
        TH1F *histo = dynamic_cast<TH1F *>(fHistoList->FindObject("SingleHisto"));
        histo->Scale(1 / histo->Integral());
        histo->SetStats(kFALSE);
        histo->SetMarkerStyle(20);
        histo->SetMarkerSize(1.0);
        histo->SetLineColor(1);
        gPad->SetLogy();
        histo->Draw();
        break;
    }
    default:
        break;
    }
}

//==========================================================================
void DData::PrintMomentum() const
{
    // print ntuple 3-momentum components of selected particles
    TIter next(fEvent->Particles());
    Int_t index = 0;
    while (DParticle *part = (DParticle *)next())
        printf("INFO:  %4i -> px = %8.5f py = %8.5f pz = %8.5f \n", index++, part->Momentum3().X(),
               part->Momentum3().Y(),
               part->Momentum3().Z());
}

//==========================================================================
void DData::PrintPID() const
{
    // print PID of selected particles
    TIter next(fEvent->Particles());
    Int_t index = 0;
    while (DParticle *part = (DParticle *)next())
        printf("INFO:  %4i -> %8i : %8s : %s \n", index++, part->PdgCode(),
               part->GetName(),
               part->GetTitle());
}

//==========================================================================
void DData::Run(const Eopt opt, const Epopt copt, const Eparam par)
{
    // the steering method to run the process
    // static Bool_t gMulPoolExist = kFALSE;
    switch (opt)
    {
    case kControlHisto:
        Loop(opt, copt);
        break;
    case kControlDiff:
        Loop(opt, copt);
        break;
    case kCorrelation:
    {
        // if (!gMulPoolExist)
        // std::cout << "Exist " << gMulPoolExist << std::endl;
        // {
        //     Loop(kMakeMulPool, copt);
        //     gMulPoolExist = kTRUE;
        // }
        Loop(opt, copt);
        break;
    }
    case kSingleHisto:
        Loop(opt, copt, par);
        break;
    default:
        break;
    }
}

//==========================================================================
void DData::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain)
        return;
    fChain->Show(entry);
}

//==========================================================================
void DData::SingleHisto(Eparam par)
{
    //makes a single histogram of par

    if (!fHCreated)
        CreateHistograms(kSingleHisto, par);
    switch (par)
    {
    case kThrust:
        (dynamic_cast<TH1F *>(fHistoList->FindObject("SingleHisto")))->Fill(fEvent->Thrust());
        break;
    case kMultiplicity:
        (dynamic_cast<TH1F *>(fHistoList->FindObject("SingleHisto")))->Fill(fEvent->Multiplicity());
        break;
    default:
        break;
    }
}

//==========================================================================
void DData::WriteOutput(const char *opt1, const char *opt2, const char *opt3)
{
    // write to the outputfile all memory objects attached to it
    if (!fHistosOutFile)
      fHistosOutFile = new TFile(Form("%s_%s_%s%s.root",fOutputFilename.Data(),opt1,opt2,opt3), "RECREATE", "Delphi basic histograms");
    if (fHistosOutFile)
    {
        //     fHistosOutFile->Write();
        fHistoList->Write("histos", TObject::kSingleKey);
        std::cout << "INFO: WriteOutPut --> generated data saved in " << fHistosOutFile->GetName() << std::endl;
    }
    else
        std::cout << "ERROR: WriteOutPut --> file "
                  << "no outpufile  open" << std::endl;
}

//====================================================================================================================================================
void DData::SetMultBinning(Int_t nBins, Double_t *limits)
{

    if (nBins > fNMaxBinsMult)
    {
        std::cout << "WARNING : only " << fNMaxBinsMult << " centrality bins (out of the " << nBins << " proposed) will be considered" << std::endl;
        nBins = fNMaxBinsMult;
    }
    if (nBins <= 0)
    {
        std::cout << "WARNING : at least one centrality bin must be considered" << std::endl;
        nBins = 1;
    }

    fNbinsMult = nBins;
    fMultAxis = new TAxis(fNbinsMult, limits);
}

//====================================================================================================================================================
void DData::SetZvtxBinning(Int_t nBins, Double_t *limits)
{

    if (nBins > fNMaxBinsZvtx)
    {
        std::cout << "WARNING : only " << fNMaxBinsZvtx << " Zvtx bins (out of the " << nBins << " proposed) will be considered" << std::endl;
        nBins = fNMaxBinsZvtx;
    }
    if (nBins <= 0)
    {
        std::cout << "WARNING : at least one Zvtx bin must be considered" << std::endl;
        nBins = 1;
    }

    fNbinsZvtx = nBins;
    fZvtxAxis = new TAxis(fNbinsZvtx, limits);
}

//====================================================================================================================================================
void DData::SetPtBinning(Int_t nBins, Double_t *limits)
{

    if (nBins > fNMaxBinsPt)
    {
        std::cout << "WARNING : only " << fNMaxBinsPt << " pt bins (out of the " << nBins << " proposed) will be considered" << std::endl;
        nBins = fNMaxBinsPt;
    }
    if (nBins <= 0)
    {
        std::cout << "WARNING : at least one pt bin must be considered" << std::endl;
        nBins = 1;
    }

    fNbinsTrackPt = nBins;
    fTrackPtAxis = new TAxis(fNbinsTrackPt, limits);
}

//====================================================================================================================================================
Int_t DData::GetMultBin() const
{

    Int_t bin = fMultAxis->FindBin(-fMultiplicity) - 1;
    if (bin >= fNbinsMult)
        bin = -1;
    return bin;
}

//====================================================================================================================================================
Bool_t DData::RotateTracks(const DParticle *trig, const DParticle *assoc, Double_t &dphi, Double_t &dtheta, Double_t &deta) const
{
  TVector3 vtrig = trig->Momentum3();
  if (vtrig.Perp2() == 0.0) {
    dphi = dtheta = 1e6;
    return kFALSE;
  }
  TVector3 vassoc = assoc->Momentum3();
  if (vassoc.Mag2() == 0.0) {
    dphi = dtheta = 1e6;
    return kFALSE;
  }
  Double_t phi = vtrig.Phi();
  vtrig.RotateZ(-phi);
  vassoc.RotateZ(-phi);

  dphi = vassoc.Phi();
  if (dphi > 1.5 * TMath::Pi())
    dphi -= TMath::TwoPi();
  if (dphi < -0.5 * TMath::Pi())
    dphi += TMath::TwoPi();

  Double_t thetatrig  = TMath::Pi()+TMath::ATan2( -vtrig.Px(),  -vtrig.Pz());
  Double_t thetaassoc = TMath::Pi()+TMath::ATan2(-vassoc.Px(), -vassoc.Pz());

  dtheta = thetaassoc - thetatrig;
  if (dtheta > 1.5 * TMath::Pi())
    dtheta -= TMath::TwoPi();
  if (dtheta < -0.5 * TMath::Pi())
    dtheta += TMath::TwoPi();

  // vtrig.RotateY(TMath::Pi()/2.-thetatrig);
  // printf("DEBUG: %f\n",vtrig.Eta());
  
  vassoc.RotateY(TMath::Pi()/2.-thetatrig);
  deta = vassoc.Eta();
  
  return kTRUE;
}
