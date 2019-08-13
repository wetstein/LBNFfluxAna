#define dk2nuTree_cxx
#include "dk2nuTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

void dk2nuTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L dk2nuTree.C
//      Root > dk2nuTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1D* arrival_time = new TH1D("at","at",500,0.,4.);
   TH2D* at_vs_E = new TH2D("atvsE","atvsE",100,0.,10.,100,0.,4.);

   TH1D* Energy = new TH1D("Energy","Energy",100,0.,10.);
   TH1D* EnergyTb1 = new TH1D("EnergyTb1","EnergyTb1",100,0.,10.);
   TH1D* EnergyTb2 = new TH1D("EnergyTb2","EnergyTb2",100,0.,10.);
   TH1D* EnergyTb3 = new TH1D("EnergyTb3","EnergyTb3",100,0.,10.);

   // Tree variables
   int nparticles,nutype,nuparent,mpotnum;
   double nustartx,nustarty,nustartz,nustartT;
   double nuE, bwt0, bwt1, bwt2, decay_wt;
   double nupx,nupy,nupz,nupperp,nuang;
   double parentpx,parentpy,parentpz,parentpperp;
   double parentE,parentang;
   double dT,dTsm0,dTsm1,dTsm2,dTsm3,dTsm4;

   bwt0=0.; bwt1=0.; bwt2=0.; decay_wt=0.;

   double nuendx = 0.;
   double nuendy = 0.;
   double nuendz = 0.;

   TFile* tf = new TFile(ofname,"recreate");
   TTree *outtree = new TTree("ANNIEOutTree","ANNIEOutTree");


   outtree->Branch("potnum",&mpotnum);
   outtree->Branch("nutype",&nutype);
   outtree->Branch("nuparent",&nuparent);
   outtree->Branch("nparticles",&nparticles);
   outtree->Branch("nuE",&nuE);
   outtree->Branch("nupx",&nupx);
   outtree->Branch("nupy",&nupy);
   outtree->Branch("nupz",&nupz);
   outtree->Branch("nupperp",&nupperp);
   outtree->Branch("nuang",&nuang);
   outtree->Branch("nustartx",&nustartx);
   outtree->Branch("nustarty",&nustarty);
   outtree->Branch("nustartz",&nustartz);
   outtree->Branch("nustartT",&nustartT);
   outtree->Branch("nuendx",&nuendx);
   outtree->Branch("nuendy",&nuendy);
   outtree->Branch("nuendz",&nuendz);
   outtree->Branch("dT",&dT);
   outtree->Branch("dTs0",&dTsm0);
   outtree->Branch("dTs1",&dTsm1);
   outtree->Branch("dTs2",&dTsm2);
   outtree->Branch("dTs3",&dTsm3);
   outtree->Branch("dTs4",&dTsm4);
   outtree->Branch("parentE",&parentE);
   outtree->Branch("parentpx",&parentpx);
   outtree->Branch("parentpy",&parentpy);
   outtree->Branch("parentpz",&parentpz);
   outtree->Branch("parentpperp",&parentpperp);
   outtree->Branch("parentang",&parentang);

   outtree->Branch("bwt0",&bwt0);
   outtree->Branch("bwt1",&bwt1);
   outtree->Branch("bwt2",&bwt2);
   outtree->Branch("decaywt",&decay_wt);


   //int rseed=40411199;
   int rseed = theseed;
   cout<<ofname<<" "<<rseed<<endl;

   TRandom3* mrand1 = new TRandom3(rseed);
   TRandom3* mrand2 = new TRandom3(rseed+1);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry%10000==0) cout<<jentry<<" out of "<<nentries<<" "<<potnum<<" POT "<<endl;

      mpotnum = potnum;

      nparticles=ancestor_;

      nutype = ancestor_pdg[ancestor_-1];
      nuparent = ancestor_pdg[ancestor_-2];

      nustartx = ancestor_startx[ancestor_-1];
      nustarty = ancestor_starty[ancestor_-1];
      nustartz = ancestor_startz[ancestor_-1];
      nustartT = ancestor_startt[ancestor_-1];

      nupx = nuray_px[1];
      nupy = nuray_py[1];
      nupz = nuray_pz[1];
      nupperp = sqrt(nupx*nupx+nupy*nupy);
      nuang = TMath::ATan(nupperp/nupz);

      parentpx = ancestor_startpx[ancestor_-2];
      parentpy = ancestor_startpy[ancestor_-2];
      parentpz = ancestor_startpz[ancestor_-2];
      parentpperp = sqrt(parentpx*parentpx+parentpy*parentpy);
      parentE = decay_ppenergy;
      parentang = TMath::ATan(parentpperp/parentpz);

      nuE = nuray_E[1];

      decay_wt = (double)decay_nimpwt;

      bwt0 = (double)decay_nimpwt*nuray_wgt[0];
      bwt1 = (double)decay_nimpwt*nuray_wgt[1];
      bwt2 = (double)decay_nimpwt*nuray_wgt[2];

      dT = nustartT - ((nustartz+360.)/29.98);

      arrival_time->Fill(dT,bwt1);
      Energy->Fill(nuE,bwt1);

      if(dT>0.15){
          EnergyTb1->Fill(nuE,bwt1);
      }

      if(dT>0.25){
          EnergyTb2->Fill(nuE,bwt1);
      }

      if(dT>0.5){
          EnergyTb3->Fill(nuE,bwt1);
      }

      double ts0 = mrand1->Gaus(0,0.1);
      double ts1 = mrand1->Gaus(0,0.25);
      double ts2 = mrand1->Gaus(0,0.5);
      double ts3 = mrand2->Gaus(0,0.9);
      double ts4 = mrand2->Gaus(0,1.2);

      dTsm0 = dT + ts0;
      dTsm1 = dT + ts1;
      dTsm2 = dT + ts2;
      dTsm3 = dT + ts3;
      dTsm4 = dT + ts4;

      at_vs_E->Fill(nuE,dT,bwt1);
      outtree->Fill();
    }

    cout<<"done and drawing"<<endl;
    arrival_time->Draw();
    at_vs_E->Draw("colz");

/*
    Energy->Draw();

    EnergyTb1->SetLineColor(2);
    EnergyTb1->Draw("SAME");
    EnergyTb2->SetLineColor(4);
    EnergyTb2->Draw("SAME");
    EnergyTb3->SetLineColor(6);
    EnergyTb3->Draw("SAME");
*/
    cout<<"writing tree"<<endl;

    std::ofstream myoutfile;
    myoutfile.open("POTlog.txt", std::ios_base::app);
    myoutfile << _fnameid <<" "<<mpotnum << endl;;



    outtree->Write();

      // if (Cut(ientry) < 0) continue;

}
