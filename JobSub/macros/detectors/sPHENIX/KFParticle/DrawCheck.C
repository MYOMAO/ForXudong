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


void DrawCheck(){

	TCanvas * c = new TCanvas("c","",600,600);

	c->SetLeftMargin(0.17);

	gStyle->SetOptStat(0);

	int KFIndex;
	float KFPropX;
	float KFPropY;
	float KFPropZ;

	float KFX;
	float KFY;
	float KFZ;

	float KFPx;
	float KFPy;
	float KFPz;
	float KFPt;


	int ACTSIndex;
	float ACTSPropX;
	float ACTSPropY;
	float ACTSPropZ;

	float ACTSX;
	float ACTSY;
	float ACTSZ;

	float ACTSPx;
	float ACTSPy;
	float ACTSPz;
	float ACTSPt;



	TFile * finKF = new TFile("Solution/MotherMody/KFFile.root");
	finKF->cd();



	TTree * KFTree = (TTree *) finKF->Get("KFTree");
	KFTree->SetBranchAddress("KFIndex",&KFIndex);
	KFTree->SetBranchAddress("KFPropX",&KFPropX);
	KFTree->SetBranchAddress("KFPropY",&KFPropY);
	KFTree->SetBranchAddress("KFPropZ",&KFPropZ);
	KFTree->SetBranchAddress("KFX",&KFX);
	KFTree->SetBranchAddress("KFY",&KFY);
	KFTree->SetBranchAddress("KFZ",&KFZ);
	KFTree->SetBranchAddress("KFPx",&KFPx);
	KFTree->SetBranchAddress("KFPy",&KFPy);
	KFTree->SetBranchAddress("KFPz",&KFPz);
	KFTree->SetBranchAddress("KFPt",&KFPt);



	TFile * finACTS = new TFile("Solution/MotherMody/ACTSFile.root");
	finACTS->cd();


	TTree * ACTSTree = (TTree *) finACTS->Get("ACTSTree");
	ACTSTree->SetBranchAddress("ACTSIndex",&ACTSIndex);
	ACTSTree->SetBranchAddress("ACTSPropX",&ACTSPropX);
	ACTSTree->SetBranchAddress("ACTSPropY",&ACTSPropY);
	ACTSTree->SetBranchAddress("ACTSPropZ",&ACTSPropZ);
	ACTSTree->SetBranchAddress("ACTSX",&ACTSX);
	ACTSTree->SetBranchAddress("ACTSY",&ACTSY);
	ACTSTree->SetBranchAddress("ACTSZ",&ACTSZ);
	ACTSTree->SetBranchAddress("ACTSPx",&ACTSPx);
	ACTSTree->SetBranchAddress("ACTSPy",&ACTSPy);
	ACTSTree->SetBranchAddress("ACTSPz",&ACTSPz);
	ACTSTree->SetBranchAddress("ACTSPt",&ACTSPt);



	float InitialX = 0;
	float InitialY = 0;
	float InitialZ = 0;


	int NBins = 10000;

	float Step = 0.01;
	float BinEnd = float(NBins) * Step;

	int NEvents = KFTree->GetEntries();

	TH1D * KFPtHis = new TH1D("KFPtHis","",NBins,0,BinEnd);
	KFPtHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	KFPtHis->GetYaxis()->SetTitle("Track p_{T} (GeV/c)");
	KFPtHis->GetXaxis()->CenterTitle();
	KFPtHis->GetYaxis()->CenterTitle();
	KFPtHis->GetYaxis()->SetTitleOffset(2.3);
	
	KFPtHis->SetLineColor(2);



	TH1D * ACTSPtHis = new TH1D("ACTSPtHis","",NBins,0,BinEnd);
	ACTSPtHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	ACTSPtHis->GetYaxis()->SetTitle("Track p_{T} (GeV/c)");
	ACTSPtHis->GetXaxis()->CenterTitle();
	ACTSPtHis->GetYaxis()->CenterTitle();
	ACTSPtHis->GetYaxis()->SetTitleOffset(2.3);

	ACTSPtHis->SetLineColor(1);


	TH1D * KFPzHis = new TH1D("KFPzHis","",NBins,0,BinEnd);
	KFPzHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	KFPzHis->GetYaxis()->SetTitle("Track p_{z} (GeV/c)");
	KFPzHis->GetXaxis()->CenterTitle();
	KFPzHis->GetYaxis()->CenterTitle();
	KFPzHis->GetYaxis()->SetTitleOffset(2.4);
	KFPzHis->SetLineColor(2);


	
	TH1D * ACTSPzHis = new TH1D("ACTSPzHis","",NBins,0,BinEnd);
	ACTSPzHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	ACTSPzHis->GetYaxis()->SetTitle("Track p_{z} (GeV/c)");
	ACTSPzHis->GetXaxis()->CenterTitle();
	ACTSPzHis->GetYaxis()->CenterTitle();
	ACTSPzHis->GetYaxis()->SetTitleOffset(2.4);
	ACTSPzHis->SetLineColor(1);

	TH1D * KFPHis = new TH1D("KFPHis","",NBins,0,BinEnd);
	KFPHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	KFPHis->GetYaxis()->SetTitle("Track Total |p| (GeV/c)");
	KFPHis->GetXaxis()->CenterTitle();
	KFPHis->GetYaxis()->CenterTitle();
	KFPHis->GetYaxis()->SetTitleOffset(1.9);
	KFPHis->SetLineColor(2);



	TH1D * ACTSPHis = new TH1D("ACTSPHis","",NBins,0,BinEnd);
	ACTSPHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	ACTSPHis->GetYaxis()->SetTitle("Track Total |p| (GeV/c)");
	ACTSPHis->GetXaxis()->CenterTitle();
	ACTSPHis->GetYaxis()->CenterTitle();
	ACTSPHis->GetYaxis()->SetTitleOffset(1.9);
	ACTSPHis->SetLineColor(1);



	TH1D * KFPxHis = new TH1D("KFPxHis","",NBins,0,BinEnd);
	KFPxHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	KFPxHis->GetYaxis()->SetTitle("Track p_{x} (GeV/c)");
	KFPxHis->GetXaxis()->CenterTitle();
	KFPxHis->GetYaxis()->CenterTitle();
	KFPxHis->GetYaxis()->SetTitleOffset(2.3);
	
	KFPxHis->SetLineColor(2);



	TH1D * ACTSPxHis = new TH1D("ACTSPxHis","",NBins,0,BinEnd);
	ACTSPxHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	ACTSPxHis->GetYaxis()->SetTitle("Track p_{x} (GeV/c)");
	ACTSPxHis->GetXaxis()->CenterTitle();
	ACTSPxHis->GetYaxis()->CenterTitle();
	ACTSPxHis->GetYaxis()->SetTitleOffset(2.3);

	ACTSPxHis->SetLineColor(1);



	TH1D * KFPyHis = new TH1D("KFPxHis","",NBins,0,BinEnd);
	KFPyHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	KFPyHis->GetYaxis()->SetTitle("Track p_{y} (GeV/c)");
	KFPyHis->GetXaxis()->CenterTitle();
	KFPyHis->GetYaxis()->CenterTitle();
	KFPyHis->GetYaxis()->SetTitleOffset(2.3);
	
	KFPyHis->SetLineColor(2);



	TH1D * ACTSPyHis = new TH1D("ACTSPyHis","",NBins,0,BinEnd);
	ACTSPyHis->GetXaxis()->SetTitle("Distance Propagated from the Starting Vertex X (cm)");
	ACTSPyHis->GetYaxis()->SetTitle("Track p_{y} (GeV/c)");
	ACTSPyHis->GetXaxis()->CenterTitle();
	ACTSPyHis->GetYaxis()->CenterTitle();
	ACTSPyHis->GetYaxis()->SetTitleOffset(2.3);

	ACTSPyHis->SetLineColor(1);





	TH2D * ACTSTrackXY = new TH2D("ACTSTrackXY","",200,0,100,100,0,200);
	ACTSTrackXY->GetXaxis()->SetTitle("Track X (cm)");
	ACTSTrackXY->GetYaxis()->SetTitle("Track Y (cm)");
	ACTSTrackXY->GetXaxis()->CenterTitle();
	ACTSTrackXY->GetYaxis()->CenterTitle();
	ACTSTrackXY->GetYaxis()->SetTitleOffset(2.3);

	ACTSTrackXY->SetMarkerColor(1);

//	TH2D * KFTrackXY = new TH2D("KFTrackXY","",100,-5,50,100,-12,5);
	TH2D * KFTrackXY = new TH2D("KFTrackXY","",100,0,100,100,0,200);
	KFTrackXY->GetXaxis()->SetTitle("Track X (cm)");
	KFTrackXY->GetYaxis()->SetTitle("Track Y (cm)");
	KFTrackXY->GetXaxis()->CenterTitle();
	KFTrackXY->GetYaxis()->CenterTitle();
	KFTrackXY->GetYaxis()->SetTitleOffset(2.3);

	KFTrackXY->SetMarkerColor(2);

	std::vector<float> ACTSXVec;
	std::vector<float> ACTSYVec;
	std::vector<float> ACTSZVec;


	std::vector<float> KFXVec;
	std::vector<float> KFYVec;
	std::vector<float> KFZVec;


	for(int i = 0; i < NEvents; i++){


		KFTree->GetEntry(i);
		ACTSTree->GetEntry(i);
	
		KFPt = sqrt(KFPx * KFPx+ KFPy * KFPy);
		ACTSPt = sqrt(ACTSPx * ACTSPx+ ACTSPy * ACTSPy);
	
		float KFP = sqrt(KFPt * KFPt + KFPz * KFPz);
		float ACTSP = sqrt(ACTSPt * ACTSPt + ACTSPz * ACTSPz);

		if(i == 0){

			InitialX = KFPropX;
			InitialY = KFPropY;
			InitialZ = KFPropZ;

		}


		KFPtHis->SetBinContent(i+1,KFPt);
		KFPzHis->SetBinContent(i+1,KFPz);
		KFPHis->SetBinContent(i+1,KFP);
		KFPxHis->SetBinContent(i+1,KFPx);
		KFPyHis->SetBinContent(i+1,KFPy);
		
		ACTSPtHis->SetBinContent(i+1,ACTSPt);
		ACTSPzHis->SetBinContent(i+1,ACTSPz);

		ACTSPHis->SetBinContent(i+1,ACTSP);
		ACTSPxHis->SetBinContent(i+1,ACTSPx);
		ACTSPyHis->SetBinContent(i+1,ACTSPy);
			

		ACTSXVec.push_back(ACTSX);
		ACTSYVec.push_back(ACTSY);
		ACTSZVec.push_back(ACTSZ);


		KFXVec.push_back(KFX);
		KFYVec.push_back(KFY);
		KFZVec.push_back(KFZ);
	

		ACTSTrackXY->Fill(ACTSX,ACTSY);
		KFTrackXY->Fill(KFX,KFY);

	
	}






	c->cd();



	KFPtHis->Draw("hist");
	ACTSPtHis->Draw("histSAME");


	TLegend * leg = new TLegend(0.32,0.65,0.70,0.80,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.034);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);

	leg->AddEntry(KFPtHis,"KFParticle","L");
	leg->AddEntry(ACTSPtHis,"ACTS","L");

	leg->Draw("SAME");

	c->SaveAs("CheckDebug/PtComp.png");




		
	KFPzHis->Draw("hist");
	ACTSPzHis->Draw("histSAME");


	TLegend * leg2 = new TLegend(0.32,0.65,0.70,0.80,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.034);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->SetLineWidth(3);

	leg2->AddEntry(KFPtHis,"KFParticle","L");
	leg2->AddEntry(ACTSPtHis,"ACTS","L");

	leg2->Draw("SAME");

	c->SaveAs("CheckDebug/PzComp.png");

	KFPHis->Draw("hist");
	ACTSPHis->Draw("histSAME");


	TLegend * leg3 = new TLegend(0.32,0.65,0.70,0.80,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.034);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->SetLineWidth(3);

	leg3->AddEntry(KFPtHis,"KFParticle","L");
	leg3->AddEntry(ACTSPtHis,"ACTS","L");

	leg3->Draw("SAME");

	c->SaveAs("CheckDebug/MomComp.png");


	TLegend * leg4 = new TLegend(0.32,0.65,0.70,0.80,NULL,"brNDC");
	leg4->SetBorderSize(0);
	leg4->SetTextSize(0.034);
	leg4->SetTextFont(42);
	leg4->SetFillStyle(0);
	leg4->SetLineWidth(3);

	leg4->AddEntry(KFPtHis,"KFParticle","L");
	leg4->AddEntry(ACTSPtHis,"ACTS","L");


	//KFPxHis->SetMaximum(2);
	//KFPxHis->SetMinimum(0);

	KFPxHis->Draw("hist");
	ACTSPxHis->Draw("histSAME");
	leg4->Draw("SAME");

	c->SaveAs("CheckDebug/PxComp.png");


	TLegend * leg5 = new TLegend(0.32,0.65,0.70,0.80,NULL,"brNDC");
	leg5->SetBorderSize(0);
	leg5->SetTextSize(0.034);
	leg5->SetTextFont(42);
	leg5->SetFillStyle(0);
	leg5->SetLineWidth(3);

	leg5->AddEntry(KFPtHis,"KFParticle","L");
	leg5->AddEntry(ACTSPtHis,"ACTS","L");


	//KFPyHis->SetMaximum(5);
	//KFPyHis->SetMinimum(-5);

	KFPyHis->Draw("hist");
	ACTSPyHis->Draw("histSAME");
	leg5->Draw("SAME");

	c->SaveAs("CheckDebug/PyComp.png");




	ACTSTrackXY->Draw("p");


	c->SaveAs("CheckDebug/ACTSTrackXY.png");



	KFTrackXY->Draw("P");

	c->SaveAs("CheckDebug/KFTrackXY.png");


	TLegend * leg6 = new TLegend(0.32,0.65,0.70,0.80,NULL,"brNDC");
	leg6->SetBorderSize(0);
	leg6->SetTextSize(0.034);
	leg6->SetTextFont(42);
	leg6->SetFillStyle(0);
	leg6->SetLineWidth(3);

	leg6->AddEntry(ACTSTrackXY,"KFParticle","P");
	leg6->AddEntry(KFTrackXY,"ACTS","P");

	KFTrackXY->Draw("p");
	ACTSTrackXY->Draw("pSAME");

	//leg6->Draw("SAME");
	
	c->SaveAs("CheckDebug/TrackXY.png");
	


	cout << "Size of vector = " << KFYVec.size() << endl;

	double KFR = -1;

	double ACTSR = -1;
	double a;
	double b;
	double cl;
	double s;

	double Area;


	float ACTSB;


	float KFB;


	TH1D * KFBHis = new TH1D("KFBHis","",100,0,2);
	KFBHis->GetYaxis()->SetTitle("Counts");
	KFBHis->GetXaxis()->SetTitle("B (T)");
	KFBHis->GetXaxis()->CenterTitle();
	KFBHis->GetYaxis()->CenterTitle();
	KFBHis->GetYaxis()->SetTitleOffset(2.3);
	
	KFBHis->SetLineColor(2);



	TH1D * ACTSBHis = new TH1D("ACTSBHis","",100,0,2);
	ACTSBHis->GetYaxis()->SetTitle("Counts");
	ACTSBHis->GetXaxis()->SetTitle("B (T)");
	ACTSBHis->GetXaxis()->CenterTitle();
	ACTSBHis->GetYaxis()->CenterTitle();
	ACTSBHis->GetYaxis()->SetTitleOffset(2.3);

	ACTSBHis->SetLineColor(1);


	int DisPic = 5000;

	int EndValue = NEvents - DisPic * 2;

	for(int i = 2; i < 103; i++){

		a = sqrt((KFXVec[i + DisPic] - KFXVec[i]) *(KFXVec[i+ DisPic] - KFXVec[i]) + (KFYVec[i + DisPic] - KFYVec[i]) *(KFYVec[i + DisPic] - KFYVec[i]) );
		b = sqrt((KFXVec[i + 2 * DisPic] - KFXVec[i+ DisPic]) *(KFXVec[i + 2 * DisPic] - KFXVec[i+DisPic]) + (KFYVec[i + 2 * DisPic] - KFYVec[i+DisPic]) *(KFYVec[i + 2 * DisPic] - KFYVec[i+ DisPic]) );
		cl = sqrt((KFXVec[i + 2 * DisPic] - KFXVec[i]) *(KFXVec[i + 2 * DisPic] - KFXVec[i]) + (KFYVec[i + 2 * DisPic] - KFYVec[i]) *(KFYVec[i + 2 * DisPic] - KFYVec[i]) );
		s = (a + b + cl) * 0.5;
		
		Area = sqrt(s * (s-a) * (s-b) * (s-cl));
		KFR = a * b * cl / (Area * 4);

		KFB = 100 * ACTSPt/(0.3 * KFR);

		//cout << "KFR = " << KFR << endl;


		KFBHis->Fill(KFB);
	}


	DisPic = 3000;

	for(int i = 5000; i < 5103; i++){

		a = sqrt((ACTSXVec[i + DisPic] - ACTSXVec[i]) *(ACTSXVec[i+ DisPic] - ACTSXVec[i]) + (ACTSYVec[i + DisPic] - ACTSYVec[i]) *(ACTSYVec[i + DisPic] - ACTSYVec[i]) );
		b = sqrt((ACTSXVec[i + 2 * DisPic] - ACTSXVec[i+ DisPic]) *(ACTSXVec[i + 2 * DisPic] - ACTSXVec[i+DisPic]) + (ACTSYVec[i + 2 * DisPic] - ACTSYVec[i+DisPic]) *(ACTSYVec[i + 2 * DisPic] - ACTSYVec[i+ DisPic]) );
		cl = sqrt((ACTSXVec[i + 2 * DisPic] - ACTSXVec[i]) *(ACTSXVec[i + 2 * DisPic] - ACTSXVec[i]) + (ACTSYVec[i + 2 * DisPic] - ACTSYVec[i]) *(ACTSYVec[i + 2 * DisPic] - ACTSYVec[i]) );
		s = (a + b + cl) * 0.5;
	
		Area = sqrt(s * (s-a) * (s-b) * (s-cl));
		ACTSR = a * b * cl / (Area * 4);
	//	cout << "ACTSR = " << ACTSR << endl;
		ACTSB = 100 * ACTSPt/(0.3 * ACTSR);
		
		ACTSBHis->Fill(ACTSB);

	}
	


	KFBHis->Draw("hist");

	c->SaveAs("CheckDebug/KFBHis.png");

	ACTSBHis->Draw("hist");
	c->SaveAs("CheckDebug/ACTSBHis.png");
	
}
