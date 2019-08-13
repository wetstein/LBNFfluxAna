#define nuTree2GlobesFlux_cxx
#include "nuTree2GlobesFlux.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void nuTree2GlobesFlux::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L nuTree2GlobesFlux.C
//      Root > nuTree2GlobesFlux t
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

   bool dopileup;  // apply pile up from previous bunch
   bool pileuprange;
   bool smearing; // apply bunch and detector smearing
   double detectionsmearing=0.1; // size of detector smearing (in nsec)
   int bs = 0;  // which mini-bunch type: 0 = center of the original bunch; 1 = edge of th original bunch
   bool isFHC; // FHC or RHC
   bool isNear;

   // nPOT
   double nPOT = 12199789.0;
   if(!isFHC) nPOT = 9899817.0;

   // specific configuration of the script
   smearing = true;
   isFHC = false;
   dopileup=true;
   pileuprange = true;
   isNear=false;

   // how many different time-bins to divide the flux into
   const int nEbins=7;

   // binning and range of the energy spectrum in each flux
   int nhistbins = 501;
   double lowE = 0.;
   double hiE = 125.250;

   // output file name
   TString fnext;
   if(smearing){
      fnext+="smeared_";
      fnext+=detectionsmearing;
      fnext+="_";
      fnext+=bs;
    }
    else{ fnext+="unsmeared"; }

    if(dopileup) fnext+="_pileup";
    if(pileuprange) fnext+="_prange";

   TString outfilename;
   outfilename = "FluxTree_";
   if(isFHC) outfilename+="FHC_";
   else outfilename+="RHC_";
   if(isNear) outfilename+="NearDet_";
   else outfilename+="FarDet_";
   outfilename+=fnext;
   outfilename+=".root";

   // make the histograms for the different flux bins
   TH1D** BeamEspect_bin;
   TH1D** BeamEspect_vmu_bin;
   TH1D** BeamEspect_vmubar_bin;
   TH1D** BeamEspect_ve_bin;
   TH1D** BeamEspect_vebar_bin;
   BeamEspect_bin = new TH1D*[nEbins];
   BeamEspect_vmu_bin = new TH1D*[nEbins];
   BeamEspect_vmubar_bin = new TH1D*[nEbins];
   BeamEspect_ve_bin = new TH1D*[nEbins];
   BeamEspect_vebar_bin = new TH1D*[nEbins];


   TH1D* BeamEspect = new TH1D("BeamEspect","BeamEspect",nhistbins,lowE,hiE);

   for(int i = 0; i<nEbins; i++){
     TString tname,bname;
     tname+="BeamEspect";
     bname+="_bin";
     bname+=i;
     BeamEspect_bin[i] = new TH1D(tname+bname,tname+bname,nhistbins,lowE,hiE);
     BeamEspect_vmu_bin[i] = new TH1D(tname+"_vmu"+bname,tname+"_vmu"+bname,nhistbins,lowE,hiE);
     BeamEspect_vmubar_bin[i] = new TH1D(tname+"_vmubar"+bname,tname+"_vmubar"+bname,nhistbins,lowE,hiE);
     BeamEspect_ve_bin[i] = new TH1D(tname+"_ve"+bname,tname+"_ve"+bname,nhistbins,lowE,hiE);
     BeamEspect_vebar_bin[i] = new TH1D(tname+"_vebar"+bname,tname+"_vebar"+bname,nhistbins,lowE,hiE);
   }

   // Normalize the beam structure pdf
   double NormIt,NormIt2;
   NormIt = bunchsmearing->Integral();
   NormIt2 = bunchsmearing2->Integral();
   bunchsmearing->Scale(1/NormIt);
   bunchsmearing2->Scale(1/NormIt2);

   // set the time ranges for the different flux bins
   double binsize;
   double lR = -0.3;
   double hR = 1.1;
   double bc[nEbins];
   double binsize = (hR-lR)/((double)nEbins);
   for(int uu=0; uu<nEbins; uu++){bc[uu] = lR + binsize*(((double)uu)+1.0); cout<<bc[uu]<<endl;}

   TRandom3* tr = new TRandom3(234342);
   double dTr;


   Long64_t nentries = fChain->GetEntriesFast();
   //nentries = 10000;
   Long64_t nbytes = 0, nb = 0;

   // Loop over all entries
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if(jentry%100000==0) cout<<jentry<<" out of "<<nentries<<endl;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double theweight;
      if(isNear) theweight = bwt1;
      else theweight = bwt2;

      BeamEspect->Fill(nuE,theweight);

      //convolute the true time with the beam structure and gaussian detector smearing
      if(smearing==false) dTr = dT;
      else dTr = dT + this->SmearIt(detectionsmearing,bs);

      // If we don't worry about pileup, we chop the flux into bins
      // ranging from -1.0 ns to 3.0 ns
      if(!pileuprange){lR=-1.0; hR=3.0;}
      // If we *do* worry about pileup, we only chop the flux into
      // bins over the 1.88 ns spacing between the bunches

      if( dTr>lR && dTr<=hR ){
        // loop over the flux bins (no pileup applied)
        for(int aa=0; aa<nEbins; aa++){
          if(aa==0 && dTr<bc[aa]){
            BeamEspect_bin[aa]->Fill(nuE,theweight);
            if(nutype==12) BeamEspect_ve_bin[aa]->Fill(nuE,theweight);
            if(nutype==14) BeamEspect_vmu_bin[aa]->Fill(nuE,theweight);
            if(nutype==-12) BeamEspect_vebar_bin[aa]->Fill(nuE,theweight);
            if(nutype==-14) BeamEspect_vmubar_bin[aa]->Fill(nuE,theweight);
          }
          if(aa>0 && dTr>bc[aa-1] && dTr<bc[aa]){
            BeamEspect_bin[aa]->Fill(nuE,theweight);
            if(nutype==12) BeamEspect_ve_bin[aa]->Fill(nuE,theweight);
            if(nutype==14) BeamEspect_vmu_bin[aa]->Fill(nuE,theweight);
            if(nutype==-12) BeamEspect_vebar_bin[aa]->Fill(nuE,theweight);
            if(nutype==-14) BeamEspect_vmubar_bin[aa]->Fill(nuE,theweight);
          }
        }
      }

      if(dopileup){

        // apply pileup from the previous N bunches and subsequent N bunches
        for(int oo=1; oo<2; oo++){

          double dTrlow = dTr - ((double)oo)*1.88;
          double dTrhi  = dTr + ((double)oo)*1.88;

          // loop over the flux bins
          for(int aa=0; aa<nEbins; aa++){

            if( dTrlow>lR && dTrlow<=hR ){
              if(aa==0 && dTrlow<bc[aa]){
                BeamEspect_bin[aa]->Fill(nuE,theweight);
                if(nutype==12) BeamEspect_ve_bin[aa]->Fill(nuE,theweight);
                if(nutype==14) BeamEspect_vmu_bin[aa]->Fill(nuE,theweight);
                if(nutype==-12) BeamEspect_vebar_bin[aa]->Fill(nuE,theweight);
                if(nutype==-14) BeamEspect_vmubar_bin[aa]->Fill(nuE,theweight);
              }
              if(aa>0 && dTrlow>bc[aa-1] && dTrlow<bc[aa]){
                BeamEspect_bin[aa]->Fill(nuE,theweight);
                if(nutype==12) BeamEspect_ve_bin[aa]->Fill(nuE,theweight);
                if(nutype==14) BeamEspect_vmu_bin[aa]->Fill(nuE,theweight);
                if(nutype==-12) BeamEspect_vebar_bin[aa]->Fill(nuE,theweight);
                if(nutype==-14) BeamEspect_vmubar_bin[aa]->Fill(nuE,theweight);
              }
            }

            if( dTrhi>lR && dTrhi<=hR ){
              if(aa==0 && dTrhi<bc[aa]){
                BeamEspect_bin[aa]->Fill(nuE,theweight);
                if(nutype==12) BeamEspect_ve_bin[aa]->Fill(nuE,theweight);
                if(nutype==14) BeamEspect_vmu_bin[aa]->Fill(nuE,theweight);
                if(nutype==-12) BeamEspect_vebar_bin[aa]->Fill(nuE,theweight);
                if(nutype==-14) BeamEspect_vmubar_bin[aa]->Fill(nuE,theweight);
              }
              if(aa>0 && dTrhi>bc[aa-1] && dTrhi<bc[aa]){
                BeamEspect_bin[aa]->Fill(nuE,theweight);
                if(nutype==12) BeamEspect_ve_bin[aa]->Fill(nuE,theweight);
                if(nutype==14) BeamEspect_vmu_bin[aa]->Fill(nuE,theweight);
                if(nutype==-12) BeamEspect_vebar_bin[aa]->Fill(nuE,theweight);
                if(nutype==-14) BeamEspect_vmubar_bin[aa]->Fill(nuE,theweight);
              }
            }
          } // end loop over flux bins
        } // end loop over "pile-up bunches"
      }  // end if(dopileup)
   }

   // Write Out the Flux Hists
   TFile* theoutputhists = new TFile(outfilename,"RECREATE");
   //BeamTiming->Write();

   double POTnorm = (1.0/nPOT);
   BeamEspect->Scale(POTnorm);
   BeamEspect->Write();

   // POT normalize and save the flux histos to a root file
   for(int i=0; i<nEbins; i++){

     BeamEspect_bin[i]->Scale(POTnorm);
     BeamEspect_vmu_bin[i]->Scale(POTnorm);
     BeamEspect_vmubar_bin[i]->Scale(POTnorm);
     BeamEspect_ve_bin[i]->Scale(POTnorm);
     BeamEspect_vebar_bin[i]->Scale(POTnorm);

     BeamEspect_bin[i]->Write();
     BeamEspect_vmu_bin[i]->Write();
     BeamEspect_vmubar_bin[i]->Write();
     BeamEspect_ve_bin[i]->Write();
     BeamEspect_vebar_bin[i]->Write();
   }

   bunchsmearing->Write("bs");
   bunchsmearing2->Write("bs2");

   // Store the Fluxes to Text Files

   for(int i=0; i<nEbins; i++){

     TString FFileName;
     if(isFHC) FFileName+="FHC_";
     else FFileName+="RHC_";
     if(isNear) FFileName+="ND_";
     else FFileName+="FD_";
     FFileName+="Flux_bin";
     FFileName+=i;
     FFileName+=".txt";

     ofstream nuflux;
     nuflux.open(FFileName);

    for(int j=0; j<nhistbins; j++){
      double theEnergy = BeamEspect_ve_bin[i]->GetBinCenter(j+1);
      double theFlux = BeamEspect_ve_bin[i]->GetBinContent(j+1);
      nuflux<<theEnergy<<" "<<theFlux<<" ";
      theFlux = BeamEspect_vmu_bin[i]->GetBinContent(j+1);
      nuflux<<theFlux<<" 0.0 ";
      theFlux = BeamEspect_vebar_bin[i]->GetBinContent(j+1);
      nuflux<<theFlux<<" ";
      theFlux = BeamEspect_vmubar_bin[i]->GetBinContent(j+1);
      nuflux<<theFlux<<" 0.0"<<endl;
    }
    nuflux.close();
   }


   TCanvas *tc5 = new TCanvas("EnergyBins","EnergyBins");
   BeamEspect->Draw();
   for(int i=0; i<nEbins; i++){
     BeamEspect_bin[i]->SetLineColor(kBlue+i);
     BeamEspect_bin[i]->Draw("SAME");
   }

   cout<<"means: "<<bunchsmearing->GetMean()<<"    "<<bunchsmearing2->GetMean()<<"  rmss: "<<bunchsmearing->GetRMS()<<" "<<bunchsmearing2->GetRMS()<<endl;
   theoutputhists->Close();
}


// This is the routine that smears
double nuTree2GlobesFlux::SmearIt(double timeresolution, int wbs)
{
  //int wbs=0;
  double smR=0.;
  double RessmR=0.;
  double beamsmR=0.;

  // Gaussian detector smearing
  RessmR = trs->Gaus(0.,timeresolution);

  // Bunch structure smearing
  if(wbs==0){ // mini-bunch at center of original bunch
    beamsmR = bunchsmearing->GetRandom();
  }
  if(wbs==1){ // mini-bunch at edge of original bunch
    beamsmR = bunchsmearing2->GetRandom();
  }

  smR=RessmR+beamsmR;
  return smR;
}
