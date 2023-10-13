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


void DrawDev(){

	TCanvas * c = new TCanvas("c","",600,600);
	
	c->SetLeftMargin(0.17);

	gStyle->SetOptStat(0);

	TString infile = "DebugMe/KFFile.root";

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

	int NBins = 200000;

	TH1D * DevPxHis = new TH1D("DevPxHis","",NBins,-1,1);
	DevPxHis->GetXaxis()->SetTitle("(KF Input X for Propagation - Truth X)/Truth X");
	DevPxHis->GetYaxis()->SetTitle("(KF Output Px After Propagation - Truth Px)/Truth Px");
	DevPxHis->GetXaxis()->CenterTitle();
	DevPxHis->GetYaxis()->CenterTitle();
	DevPxHis->GetYaxis()->SetTitleOffset(1.8);


	TH1D * DevPyHis = new TH1D("DevPyHis","",NBins,-1,1);
	DevPyHis->GetXaxis()->SetTitle("(KF Input X for Propagation - Truth X)/Truth X");
	DevPyHis->GetYaxis()->SetTitle("(KF Output Py After Propagation - Truth Py)/Truth Py");
	DevPyHis->GetXaxis()->CenterTitle();
	DevPyHis->GetYaxis()->CenterTitle();
	DevPyHis->GetYaxis()->SetTitleOffset(1.8);

	
	TH1D * DevPzHis = new TH1D("DevPzHis","",NBins,-1,1);
	DevPzHis->GetXaxis()->SetTitle("(KF Input X for Propagation - Truth X)/Truth X");
	DevPzHis->GetYaxis()->SetTitle("(KF Output Pz After Propagation - Truth Pz)/Truth Pz");
	DevPzHis->GetXaxis()->CenterTitle();
	DevPzHis->GetYaxis()->CenterTitle();
	DevPzHis->GetYaxis()->SetTitleOffset(1.8);

	TH1D * DevPtHis = new TH1D("DevPtHis","",NBins,-1,1);
	DevPtHis->GetXaxis()->SetTitle("(KF Input X for Propagation - Truth X)/Truth X");
	DevPtHis->GetYaxis()->SetTitle("(KF Output Pt After Propagation - Truth Pt)/Truth Pt");
	DevPtHis->GetXaxis()->CenterTitle();
	DevPtHis->GetYaxis()->CenterTitle();
	DevPtHis->GetYaxis()->SetTitleOffset(1.8);

	
	float CalPt;
	float TruthPt;

	for(int i = 0; i < NEvents; i++){


		KFTree->GetEntry(i);
			

		DevX = KFX - KFTruthX;

		DevPX = (KFPx - KFTruthPx)/KFTruthPx;
		DevPY = (KFPy - KFTruthPy)/KFTruthPy;
		DevPZ = (KFPz - KFTruthPz)/KFTruthPz;

		CalPt = sqrt(KFPx * KFPx + KFPy * KFPy); 
		TruthPt = sqrt(KFTruthPx * KFTruthPx + KFTruthPy * KFTruthPy); 


		DevPT = (CalPt - TruthPt)/TruthPt;
		 
		int BinX = DevPxHis->GetXaxis()->FindBin(DevX);

		DevPxHis->SetBinContent(BinX,DevPX);
		DevPyHis->SetBinContent(BinX,DevPY);
		DevPzHis->SetBinContent(BinX,DevPZ);
		DevPtHis->SetBinContent(BinX,DevPT);

	}

	c->cd();


	DevPxHis->Draw("hist");
	c->SaveAs("DebugMe/Plots/DevPxHis.png");

	DevPyHis->Draw("hist");
	c->SaveAs("DebugMe/Plots/DevPyHis.png");


	DevPzHis->Draw("hist");
	c->SaveAs("DebugMe/Plots/DevPzHis.png");

	
	DevPtHis->Draw("hist");
	c->SaveAs("DebugMe/Plots/DevPtHis.png");

}

