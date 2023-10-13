#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
#include "TMath.h"

using namespace std;

using std::cout;
using std::endl;


void KFPEvtMatch(){


	TFile * finZZ = new TFile("UsefulTrkInfo.root");
	finZZ->cd();

	TTree * TrkInfoforKF = (TTree*) finZZ->Get("TrkInfoforKF");
	

	TFile * finKFP = new TFile("mypipiReco/Test3.root");
	finKFP->cd();

	TTree * DecayTree = (TTree *) finKFP->Get("DecayTree");
	
	float KSPX;
	float track_1_TruthParentPx;


	TrkInfoforKF->SetBranchAddress("KSPX",&KSPX);
	DecayTree->SetBranchAddress("track_1_TruthParentPx",&track_1_TruthParentPx);


	int NEventsZZ = TrkInfoforKF->GetEntries();
	int NEventsKFP = DecayTree->GetEntries();
	

	TFile * fout = new TFile("FinalFile.root","RECREATE");
	fout->cd();

	TTree * TrkInfoforKF_new = TrkInfoforKF->CloneTree(0);
	TTree * DecayTree_new = DecayTree->CloneTree(0);


	for(int i = 0; i < NEventsKFP; i++){

		DecayTree->GetEntry(i);

		for(int j = 0; j < NEventsZZ; j++){
			
			TrkInfoforKF->GetEntry(j);

			if(track_1_TruthParentPx == KSPX){

				DecayTree_new->Fill();
				TrkInfoforKF_new->Fill();
			}

		}

	}

	DecayTree_new->Write();
	TrkInfoforKF_new->Write();
	fout->Close();
	

}
