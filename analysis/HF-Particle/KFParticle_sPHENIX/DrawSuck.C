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


void DrawSuck(){

	TFile * fin = new TFile("KFFile.root");

	fin->cd();


	TTree * KFTree = (TTree *) fin->Get("KFTree");


	TCanvas * c = new TCanvas("c","",600,600);
	c->cd();
	//c->SetLeftMargin(0.17);

	gStyle->SetOptStat(0);


	TH1D * KFXHis = new TH1D("KFXHis","",100,18.4,18.5);
	KFXHis->GetXaxis()->SetTitle("SV Position X (cm)");
	KFXHis->GetYaxis()->SetTitle("Counts");
	KFXHis->GetXaxis()->CenterTitle();
	KFXHis->GetYaxis()->CenterTitle();
	KFXHis->GetYaxis()->SetTitleOffset(1.5);
	KFXHis->SetLineColor(1);


	TH1D * KFYHis = new TH1D("KFYHis","",100,36.9,37.2);
	KFYHis->GetXaxis()->SetTitle("SV Position Y (cm)");
	KFYHis->GetYaxis()->SetTitle("Counts");
	KFYHis->GetXaxis()->CenterTitle();
	KFYHis->GetYaxis()->CenterTitle();
	KFYHis->GetYaxis()->SetTitleOffset(1.5);
	KFYHis->SetLineColor(1);



	TH1D * KFZHis = new TH1D("KFZHis","",100,-22.3,-22.1);
	KFZHis->GetXaxis()->SetTitle("SV Position Z (cm)");
	KFZHis->GetYaxis()->SetTitle("Counts");
	KFZHis->GetXaxis()->CenterTitle();
	KFZHis->GetYaxis()->CenterTitle();
	KFZHis->GetYaxis()->SetTitleOffset(1.5);
	KFZHis->SetLineColor(1);

	TH1D * KFPxHis = new TH1D("KFPxHis","",100,2.40,2.45);
	KFPxHis->GetXaxis()->SetTitle("K_{S}^{0} Px (GeV/c)");
	KFPxHis->GetYaxis()->SetTitle("Counts");
	KFPxHis->GetXaxis()->CenterTitle();
	KFPxHis->GetYaxis()->CenterTitle();
	KFPxHis->GetYaxis()->SetTitleOffset(1.5);
	KFPxHis->SetLineColor(1);


	TH1D * KFPyHis = new TH1D("KFPyHis","",100,4.90,4.93);
	KFPyHis->GetXaxis()->SetTitle("K_{S}^{0} Py (GeV/c)");
	KFPyHis->GetYaxis()->SetTitle("Counts");
	KFPyHis->GetXaxis()->CenterTitle();
	KFPyHis->GetYaxis()->CenterTitle();
	KFPyHis->GetYaxis()->SetTitleOffset(1.5);
	KFPyHis->SetLineColor(1);


	TH1D * KFPzHis = new TH1D("KFPzHis","",100,-2.97,-2.94);
	KFPzHis->GetXaxis()->SetTitle("K_{S}^{0} Pz (GeV/c)");
	KFPzHis->GetYaxis()->SetTitle("Counts");
	KFPzHis->GetXaxis()->CenterTitle();
	KFPzHis->GetYaxis()->CenterTitle();
	KFPzHis->GetYaxis()->SetTitleOffset(1.5);
	KFPzHis->SetLineColor(1);



	TH1D * CovXXHis = new TH1D("CovXXHis","",100,0,0.03);
	CovXXHis->GetXaxis()->SetTitle("SV Covariance Matrix (0,0) (cm^{2})");
	CovXXHis->GetYaxis()->SetTitle("Counts");
	CovXXHis->GetXaxis()->CenterTitle();
	CovXXHis->GetYaxis()->CenterTitle();
	CovXXHis->GetYaxis()->SetTitleOffset(1.5);
	CovXXHis->SetLineColor(1);



	TH1D * CovYYHis = new TH1D("CovYYHis","",100,0,0.12);
	CovYYHis->GetXaxis()->SetTitle("SV Covariance Matrix (1,1) (cm^{2})");
	CovYYHis->GetYaxis()->SetTitle("Counts");
	CovYYHis->GetXaxis()->CenterTitle();
	CovYYHis->GetYaxis()->CenterTitle();
	CovYYHis->GetYaxis()->SetTitleOffset(1.5);
	CovYYHis->SetLineColor(1);



	TH1D * CovZZHis = new TH1D("CovZZHis","",100,0.0,0.05);
	CovZZHis->GetXaxis()->SetTitle("SV Covariance Matrix (2,2) (cm^{2})");
	CovZZHis->GetYaxis()->SetTitle("Counts");
	CovZZHis->GetXaxis()->CenterTitle();
	CovZZHis->GetYaxis()->CenterTitle();
	CovZZHis->GetYaxis()->SetTitleOffset(1.5);
	CovZZHis->SetLineColor(1);


	TH1D * CovPxPxHis = new TH1D("CovPxPxHis","",100,0.0002,0.0003);
	CovPxPxHis->GetXaxis()->SetTitle("SV Covariance Matrix (3,3) [(GeV/c)^{2}]");
	CovPxPxHis->GetYaxis()->SetTitle("Counts");
	CovPxPxHis->GetXaxis()->CenterTitle();
	CovPxPxHis->GetYaxis()->CenterTitle();
	CovPxPxHis->GetYaxis()->SetTitleOffset(1.5);
	CovPxPxHis->SetLineColor(1);
	

	TH1D * CovPyPyHis = new TH1D("CovPyPyHis","",100,0.001,0.002);
	CovPyPyHis->GetXaxis()->SetTitle("SV Covariance Matrix (4,4) [(GeV/c)^{2}]");
	CovPyPyHis->GetYaxis()->SetTitle("Counts");
	CovPyPyHis->GetXaxis()->CenterTitle();
	CovPyPyHis->GetYaxis()->CenterTitle();
	CovPyPyHis->GetYaxis()->SetTitleOffset(1.5);
	CovPyPyHis->SetLineColor(1);
	

	TH1D * CovPzPzHis = new TH1D("CovPzPzHis","",100,0.0003,0.0006);
	CovPzPzHis->GetXaxis()->SetTitle("SV Covariance Matrix (5,5) [(GeV/c)^{2}]");
	CovPzPzHis->GetYaxis()->SetTitle("Counts");
	CovPzPzHis->GetXaxis()->CenterTitle();
	CovPzPzHis->GetYaxis()->CenterTitle();
	CovPzPzHis->GetYaxis()->SetTitleOffset(1.5);
	CovPzPzHis->SetLineColor(1);


	TH1D * KsMassHis = new TH1D("KsMassHis","",200,-0.5,0.8);
	KsMassHis->GetXaxis()->SetTitle("Reconstructed K^{0}_{S} Invariant Mass (GeV/c^{2})");
	KsMassHis->GetYaxis()->SetTitle("Counts");
	KsMassHis->GetXaxis()->CenterTitle();
	KsMassHis->GetYaxis()->CenterTitle();
	KsMassHis->GetYaxis()->SetTitleOffset(1.5);
	KsMassHis->SetLineColor(1);


	TH1D * KsMassCalHis = new TH1D("KsMassCalHis","",200,0.2,0.8);
	KsMassCalHis->GetXaxis()->SetTitle("Reconstructed K^{0}_{S} Invariant Mass (GeV/c^{2})");
	KsMassCalHis->GetYaxis()->SetTitle("Counts");
	KsMassCalHis->GetXaxis()->CenterTitle();
	KsMassCalHis->GetYaxis()->CenterTitle();
	KsMassCalHis->GetYaxis()->SetTitleOffset(1.5);
	KsMassCalHis->SetLineColor(1);


	TH2D * KsMassCov = new TH2D("KsMassCov","",100,0,0.001,100,0.2,0.8);
	KsMassCov->GetXaxis()->SetTitle("Covariance Matrix (0,0) (cm^{2})");
	KsMassCov->GetYaxis()->SetTitle("Reconstructed K^{0}_{S} Invariant Mass (GeV/c^{2})");
	KsMassCov->GetXaxis()->CenterTitle();
	KsMassCov->GetYaxis()->CenterTitle();
	KsMassCov->GetYaxis()->SetTitleOffset(1.5);




	KFTree->Project("KFXHis","MotherX");
	KFXHis->Draw("hist");
	c->SaveAs("AnotherPlots/KFXHis.png");

	KFTree->Project("KFYHis","MotherY");
	KFYHis->Draw("hist");
	c->SaveAs("AnotherPlots/KFYHis.png");
	
	KFTree->Project("KFZHis","MotherZ");
	KFZHis->Draw("hist");
	c->SaveAs("AnotherPlots/KFZHis.png");

	KFTree->Project("KFPxHis","MotherPx");
	KFPxHis->Draw("hist");
	c->SaveAs("AnotherPlots/KFPxHis.png");

	KFTree->Project("KFPyHis","MotherPy");
	KFPyHis->Draw("hist");
	c->SaveAs("AnotherPlots/KFPyHis.png");

	KFTree->Project("KFPzHis","MotherPz");
	KFPzHis->Draw("hist");
	c->SaveAs("AnotherPlots/KFPzHis.png");


	KFTree->Project("CovXXHis","MotherCovXX");
	CovXXHis->Draw("hist");
	c->SaveAs("AnotherPlots/CovXXHis.png");

	KFTree->Project("CovYYHis","MotherCovYY");
	CovYYHis->Draw("hist");
	c->SaveAs("AnotherPlots/CovYYHis.png");

	KFTree->Project("CovZZHis","MotherCovZZ");
	CovZZHis->Draw("hist");
	c->SaveAs("AnotherPlots/CovZZHis.png");



	KFTree->Project("CovPxPxHis","MotherCovPxPx");
	CovPxPxHis->Draw("hist");
	c->SaveAs("AnotherPlots/CovPxPxHis.png");


	KFTree->Project("CovPyPyHis","MotherCovPyPy");
	CovPyPyHis->Draw("hist");
	c->SaveAs("AnotherPlots/CovPyPyHis.png");

	KFTree->Project("CovPzPzHis","MotherCovPzPz");
	CovPzPzHis->Draw("hist");
	c->SaveAs("AnotherPlots/CovPzPzHis.png");

	KFTree->Project("KsMassHis","MotherMass");
	KsMassHis->Draw("hist");
	c->SaveAs("AnotherPlots/KsMassHis.png");

	KFTree->Project("KsMassCalHis","KsMassCal");
	KsMassCalHis->Draw("hist");
	c->SaveAs("AnotherPlots/KsMassCalHis.png");

	KFTree->Project("KsMassCov","KsMassCal:MotherCovXX");
	KsMassCov->Draw("COLZ");
	c->SaveAs("AnotherPlots/KsMassCov.png");


}
