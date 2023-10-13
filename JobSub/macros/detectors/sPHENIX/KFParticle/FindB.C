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


void FindB(int Opt){


	TString Name;

	if(Opt == 0) Name = "FindB";
	if(Opt == 1) Name = "FindBFlip";

	if(Opt == 2) Name = "FindBFlip2";

	TCanvas * c = new TCanvas("c","",600,600);
	
	c->SetLeftMargin(0.17);

	gStyle->SetOptStat(0);

	TString infile = Form("%s/KFFile.root",Name.Data());

	TFile * fin = new TFile(infile.Data());
	
	fin->cd();


	int KFIndex;


	float KFX;
	float KFY;
	float KFZ;

	float KFPx;
	float KFPy;
	float KFPz;
	float KFPt;


		
	float KFTruthX;
	float KFTruthY;
	float KFTruthZ;

	float KFTruthPx;
	float KFTruthPy;
	float KFTruthPz;


	TTree * KFTree = (TTree *) fin->Get("KFTree");
	KFTree->SetBranchAddress("KFIndex",&KFIndex);
	KFTree->SetBranchAddress("KFX",&KFX);
	KFTree->SetBranchAddress("KFY",&KFY);
	KFTree->SetBranchAddress("KFZ",&KFZ);
	KFTree->SetBranchAddress("KFPx",&KFPx);
	KFTree->SetBranchAddress("KFPy",&KFPy);
	KFTree->SetBranchAddress("KFPz",&KFPz);
	KFTree->SetBranchAddress("KFPt",&KFPt);

	KFTree->SetBranchAddress("KFTruthX",&KFTruthX);
	KFTree->SetBranchAddress("KFTruthY",&KFTruthY);
	KFTree->SetBranchAddress("KFTruthZ",&KFTruthZ);
	KFTree->SetBranchAddress("KFTruthPx",&KFTruthPx);
	KFTree->SetBranchAddress("KFTruthPy",&KFTruthPy);
	KFTree->SetBranchAddress("KFTruthPz",&KFTruthPz);

	
	int NEvents = KFTree->GetEntries();

	
	float DevPX;
	float DevPY;
	float DevPZ;
	float DevPT;

	float DevX;

	int NBins = 20000;


	TH1D * KFPxKFX = new TH1D("KFPxKFX","",NBins,-10000,10000);
	KFPxKFX->GetXaxis()->SetTitle("KFX (cm)");
	KFPxKFX->GetYaxis()->SetTitle("KFPx (GeV/c)");
	KFPxKFX->GetXaxis()->CenterTitle();
	KFPxKFX->GetYaxis()->CenterTitle();
	KFPxKFX->GetYaxis()->SetTitleOffset(1.8);



	KFPxKFX->SetMinimum(-0.2);
	KFPxKFX->SetMaximum(0.2);


	TH1D * KFPyKFX = new TH1D("KFPyKFX","",NBins,-10,10);
	KFPyKFX->GetXaxis()->SetTitle("KFX (cm)");
	KFPyKFX->GetYaxis()->SetTitle("KFPy (GeV/c)");
	KFPyKFX->GetXaxis()->CenterTitle();
	KFPyKFX->GetYaxis()->CenterTitle();
	KFPyKFX->GetYaxis()->SetTitleOffset(1.8);


	TH1D * KFPzKFX = new TH1D("KFPzKFX","",NBins,-10,10);
	KFPzKFX->GetXaxis()->SetTitle("KFX (cm)");
	KFPzKFX->GetYaxis()->SetTitle("KFPz (GeV/c)");
	KFPzKFX->GetXaxis()->CenterTitle();
	KFPzKFX->GetYaxis()->CenterTitle();
	KFPzKFX->GetYaxis()->SetTitleOffset(1.8);


	for(int i = 0; i < NEvents; i++){


		KFTree->GetEntry(i);


		KFPxKFX->SetBinContent(i+1,KFPx);
		KFPyKFX->SetBinContent(i+1,KFPy);
		KFPzKFX->SetBinContent(i+1,KFPz);

		
	}


	c->cd();


	KFPxKFX->Draw("hist");
	c->SaveAs(Form("%s/Plot/KFPxKFX.png",Name.Data()));

	KFPyKFX->Draw("hist");
	c->SaveAs(Form("%s/Plot/KFPyKFX.png",Name.Data()));


	KFPzKFX->Draw("hist");
	c->SaveAs(Form("%s/Plot/KFPzKFX.png",Name.Data()));




}
