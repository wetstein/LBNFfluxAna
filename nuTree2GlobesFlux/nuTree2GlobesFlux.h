//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 24 13:43:36 2018 by ROOT version 5.34/36
// from TTree ANNIEOutTree/ANNIEOutTree
// found on file: outfileFHC_all.root
//////////////////////////////////////////////////////////

#ifndef nuTree2GlobesFlux_h
#define nuTree2GlobesFlux_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class nuTree2GlobesFlux {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nutype;
   Int_t           nuparent;
   Int_t           nparticles;
   Double_t        nuE;
   Double_t        nupx;
   Double_t        nupy;
   Double_t        nupz;
   Double_t        nupperp;
   Double_t        nuang;
   Double_t        nustartx;
   Double_t        nustarty;
   Double_t        nustartz;
   Double_t        nustartT;
   Double_t        nuendx;
   Double_t        nuendy;
   Double_t        nuendz;
   Double_t        dT;
   Double_t        dTs0;
   Double_t        dTs1;
   Double_t        dTs2;
   Double_t        dTs3;
   Double_t        dTs4;
   Double_t        parentE;
   Double_t        parentpx;
   Double_t        parentpy;
   Double_t        parentpz;
   Double_t        parentpperp;
   Double_t        parentang;
   Double_t        bwt0;
   Double_t        bwt1;
   Double_t        bwt2;


   TH1D* tmpbunchsmearing;
   TH1D* bunchsmearing;
   TH1D* bunchsmearing2;

   // List of branches
   TBranch        *b_nutype;   //!
   TBranch        *b_nuparent;   //!
   TBranch        *b_nparticles;   //!
   TBranch        *b_nuE;   //!
   TBranch        *b_nupx;   //!
   TBranch        *b_nupy;   //!
   TBranch        *b_nupz;   //!
   TBranch        *b_nupperp;   //!
   TBranch        *b_nuang;   //!
   TBranch        *b_nustartx;   //!
   TBranch        *b_nustarty;   //!
   TBranch        *b_nustartz;   //!
   TBranch        *b_nustartT;   //!
   TBranch        *b_nuendx;   //!
   TBranch        *b_nuendy;   //!
   TBranch        *b_nuendz;   //!
   TBranch        *b_dT;   //!
   TBranch        *b_dTs0;   //!
   TBranch        *b_dTs1;   //!
   TBranch        *b_dTs2;   //!
   TBranch        *b_dTs3;   //!
   TBranch        *b_dTs4;   //!
   TBranch        *b_parentE;   //!
   TBranch        *b_parentpx;   //!
   TBranch        *b_parentpy;   //!
   TBranch        *b_parentpz;   //!
   TBranch        *b_parentpperp;   //!
   TBranch        *b_parentang;   //!
   TBranch        *b_bwt0;   //!
   TBranch        *b_bwt1;   //!
   TBranch        *b_bwt2;   //!


   TRandom3* trs;

   TString _fname;

   nuTree2GlobesFlux(TString fname, TTree *tree=0);
   nuTree2GlobesFlux(TTree *tree=0);
   virtual ~nuTree2GlobesFlux();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
//   double SmearIt(double res, TH1D* beamshape);
   double SmearIt(double res, int wbs);
//   double SmearIt(double res);
};

#endif

#ifdef nuTree2GlobesFlux_cxx

nuTree2GlobesFlux::nuTree2GlobesFlux(TString fname,TTree *tree) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
     if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
        if (!f || !f->IsOpen()) {
          cout<<"THE NAME: "<<fname<<endl;
           f = new TFile(fname);
        }
        f->GetObject("ANNIEOutTree",tree);

     }
     Init(tree);
}


nuTree2GlobesFlux::nuTree2GlobesFlux(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("FHC_all_2019-08-01.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("input_files/FHC_all_2019-08-01.root");
      }
      f->GetObject("ANNIEOutTree",tree);
   }
   Init(tree);
}

nuTree2GlobesFlux::~nuTree2GlobesFlux()
{
   cout<<"ok, at the end"<<endl;
   if (!fChain) return;
   //delete bunchsmearing;
   //delete bunchsmearing2;
   delete fChain->GetCurrentFile();
}

Int_t nuTree2GlobesFlux::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t nuTree2GlobesFlux::LoadTree(Long64_t entry)
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

void nuTree2GlobesFlux::Init(TTree *tree)
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

   fChain->SetBranchAddress("nutype", &nutype, &b_nutype);
   fChain->SetBranchAddress("nuparent", &nuparent, &b_nuparent);
   fChain->SetBranchAddress("nparticles", &nparticles, &b_nparticles);
   fChain->SetBranchAddress("nuE", &nuE, &b_nuE);
   fChain->SetBranchAddress("nupx", &nupx, &b_nupx);
   fChain->SetBranchAddress("nupy", &nupy, &b_nupy);
   fChain->SetBranchAddress("nupz", &nupz, &b_nupz);
   fChain->SetBranchAddress("nupperp", &nupperp, &b_nupperp);
   fChain->SetBranchAddress("nuang", &nuang, &b_nuang);
   fChain->SetBranchAddress("nustartx", &nustartx, &b_nustartx);
   fChain->SetBranchAddress("nustarty", &nustarty, &b_nustarty);
   fChain->SetBranchAddress("nustartz", &nustartz, &b_nustartz);
   fChain->SetBranchAddress("nustartT", &nustartT, &b_nustartT);
   fChain->SetBranchAddress("nuendx", &nuendx, &b_nuendx);
   fChain->SetBranchAddress("nuendy", &nuendy, &b_nuendy);
   fChain->SetBranchAddress("nuendz", &nuendz, &b_nuendz);
   fChain->SetBranchAddress("dT", &dT, &b_dT);
   fChain->SetBranchAddress("dTs0", &dTs0, &b_dTs0);
   fChain->SetBranchAddress("dTs1", &dTs1, &b_dTs1);
   fChain->SetBranchAddress("dTs2", &dTs2, &b_dTs2);
   fChain->SetBranchAddress("dTs3", &dTs3, &b_dTs3);
   fChain->SetBranchAddress("dTs4", &dTs4, &b_dTs4);
   fChain->SetBranchAddress("parentE", &parentE, &b_parentE);
   fChain->SetBranchAddress("parentpx", &parentpx, &b_parentpx);
   fChain->SetBranchAddress("parentpy", &parentpy, &b_parentpy);
   fChain->SetBranchAddress("parentpz", &parentpz, &b_parentpz);
   fChain->SetBranchAddress("parentpperp", &parentpperp, &b_parentpperp);
   fChain->SetBranchAddress("parentang", &parentang, &b_parentang);
   fChain->SetBranchAddress("bwt0", &bwt0, &b_bwt0);
   fChain->SetBranchAddress("bwt1", &bwt1, &b_bwt1);
   fChain->SetBranchAddress("bwt2", &bwt2, &b_bwt2);


   bunchsmearing = new TH1D("bunchsmearing","bunchsmearing",163,-0.8366,0.7934);
   bunchsmearing2 = new TH1D("bunchsmearing2","bunchsmearing2",163,-0.8366,0.7934);

   //TFile* tbf = new TFile("BeamSimulation/531mhz_down10_up20_buff30.root","READ");
   TFile* tbf = new TFile("input_files/531mhz_fnal-real_buffertime_down30_up20_buf10_20hz_1e5.root","READ");

   //TH1D* tmpbunchsmearing = (TH1D*) tbf->Get("times_onebucket");
   tmpbunchsmearing = (TH1D*) tbf->Get("times");

   cout<<"NBins: "<<tmpbunchsmearing->GetNbinsX()<<" "<<tmpbunchsmearing->GetXaxis()->GetXmin()<<" "<<tmpbunchsmearing->GetXaxis()->GetXmax()<<endl;

   for(int aaa=0; aaa<163; aaa++){

     //int bn = tmpbunchsmearing->FindBin(-10.0-0.327+0.005+aaa*0.01);
     int bn = tmpbunchsmearing->FindBin(-0.8366+0.005+aaa*0.01);
     int bn2 = tmpbunchsmearing->FindBin(-9.42-0.8366+0.005+aaa*0.01);

     double bc = tmpbunchsmearing->GetBinContent(bn);
     double bc2 = tmpbunchsmearing->GetBinContent(bn2);

     bunchsmearing->SetBinContent(aaa+1,bc);
     bunchsmearing2->SetBinContent(aaa+1,bc2);
   }


   trs = new TRandom3(12441);

   Notify();
}

Bool_t nuTree2GlobesFlux::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void nuTree2GlobesFlux::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t nuTree2GlobesFlux::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef nuTree2GlobesFlux_cxx
