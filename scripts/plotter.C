#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include <vector>
#include "TMath.h"
#include <iostream>

#ifdef __CINT__
#pragma link C++ class std::vector<Long64_t>+;
#endif


float plotterFull(int run, bool old) {

    
    TChain *t = new TChain("tracks");
    TString fname;
    if (old)
	//fname = "../../captSummary/cmt/tree_full_golden_old.root";
	fname = "../../captSummary/cmt/tree_ASCII_old.root";

    else
	//fname = "../../captSummary/cmt/tree_full_golden.root";
	fname = "../../captSummary/cmt/tree_ASCII.root";
    
    //std::cout<<"FNAME="<<fname;
    //t->Add("../../captSummary/cmt/tracks_12135.root");
    t->Add(fname);
    
    //std::cout<<"FNAME="<<fname;
    Int_t nEntries = t->GetEntries();

    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    
    Int_t RUN = 0;
    
    std::vector<double> *first_hit_X = 0;
    std::vector<double> *first_hit_Y = 0;
    std::vector<double> *first_hit_Z = 0;

    std::vector<double> *last_hit_X = 0;
    std::vector<double> *last_hit_Y = 0;
    std::vector<double> *last_hit_Z = 0;

    Long64_t TPC_time = 0;
  
    std::vector<Long64_t> *PDS_RF_time = 0;
    std::vector<Long64_t> *PDS_delta_time = 0;
    std::vector<int> *PDS_beam_trigger = 0;
  
    t->SetBranchAddress("run", &RUN);

    t->SetBranchAddress("first_hit_X", &first_hit_X);
    t->SetBranchAddress("first_hit_Z", &first_hit_Z);

    t->SetBranchAddress("last_hit_X", &last_hit_X);
    t->SetBranchAddress("last_hit_Z", &last_hit_Z);

    t->SetBranchAddress("TPC_time", &TPC_time);
    t->SetBranchAddress("PDS_RF_time", &PDS_RF_time);
    t->SetBranchAddress("PDS_delta_time", &PDS_delta_time);
    t->SetBranchAddress("PDS_beam_trigger", &PDS_beam_trigger);

    ////////////////////////////////////////////////
    ////////////////////////////////////////////////

    
    TH1F *h_first_hit_Z = new TH1F("h_first_hit_Z","",100,-1000,0);
    TH1F *h_corrected_Zpos = new TH1F("h_corrected_Zpos","",100,-1000,0);

    int matches = 0;
  
    for (int i=0; i< nEntries; i++) {
	//for (int i=850; i< 945; i++) {
	t->GetEntry(i);
	if (RUN != run) continue;

	float Zposf = 0.0;
	float corrected_Zposf = 0.0;
	int deltaT = 0;
		
	for (size_t j = 0; j<PDS_beam_trigger->size(); j++) {
	    if (PDS_beam_trigger->at(j) == 1) {
		deltaT = PDS_delta_time->at(j);
	    }
	}

	std::vector<double> Zpositions;

	//printf("sizes x=%d z=%d\n",first_hit_X->size(),first_hit_Z->size());
	
	for (size_t j = 0; j<first_hit_X->size(); j++) {
	    if (first_hit_X->at(j) > last_hit_X->at(j))
		Zpositions.push_back(first_hit_Z->at(j));
	    else
		Zpositions.push_back(last_hit_Z->at(j));
	}

	for (size_t j = 0; j<Zpositions.size(); j++) {
	    if (Zpositions[j] < -300 && Zpositions[j] > -1000)
		{
		    Zposf = Zpositions[j];
		
		    for (size_t k = 0; k<PDS_delta_time->size(); k++) {		    
			if (Zposf + PDS_delta_time->at(k)*1.6/1000 < 100 && Zposf + PDS_delta_time->at(k)*1.6/1000 > -350) {
			    corrected_Zposf = Zposf + PDS_delta_time->at(k)*1.6/1000;
			    break;
			}		    
		    }

		    if (fabs(Zposf - corrected_Zposf) != 0.0 )
			matches++;
		
		    h_first_hit_Z->Fill(Zposf);
		    if (corrected_Zposf < 100 && corrected_Zposf > -350)
			h_corrected_Zpos->Fill(corrected_Zposf);
		
		}
	    
	}
    


    }

    float eff = float(h_corrected_Zpos->Integral()/h_first_hit_Z->Integral());
    
    h_first_hit_Z->Delete();
    h_corrected_Zpos->Delete();
    
    return eff;
    
    //return 1.0;
}

void plotterS(TString fname){


    TChain *t = new TChain("tracks");
    t->Add(fname.Data());
    
    Int_t nEntries = t->GetEntries();

    //std::cout<<nEntries<<std::endl;
  
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
  
    std::vector<double> *first_hit_X = 0;
    std::vector<double> *first_hit_Y = 0;
    std::vector<double> *first_hit_Z = 0;

    std::vector<double> *last_hit_X = 0;
    std::vector<double> *last_hit_Y = 0;
    std::vector<double> *last_hit_Z = 0;

    int evt;
    int run;
    
    Long64_t TPC_time = 0;
  
    std::vector<Long64_t> *PDS_RF_time = 0;
    std::vector<double> *PDS_delta_time = 0;
    std::vector<double> *PDS_energy = 0;
    std::vector<int> *PDS_beam_trigger = 0;
  
    t->SetBranchAddress("first_hit_X", &first_hit_X);
    t->SetBranchAddress("first_hit_Y", &first_hit_Y);
    t->SetBranchAddress("first_hit_Z", &first_hit_Z);

    t->SetBranchAddress("last_hit_X", &last_hit_X);
    t->SetBranchAddress("last_hit_Y", &last_hit_Y);
    t->SetBranchAddress("last_hit_Z", &last_hit_Z);

    t->SetBranchAddress("TPC_time", &TPC_time);
    t->SetBranchAddress("PDS_RF_time", &PDS_RF_time);
    t->SetBranchAddress("PDS_delta_time", &PDS_delta_time);
    t->SetBranchAddress("PDS_energy", &PDS_energy);
    t->SetBranchAddress("PDS_beam_trigger", &PDS_beam_trigger);

    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("event", &evt);
    

    ////////////////////////////////////////////////
    ////////////////////////////////////////////////

    TH1F *h_first_hit_X = new TH1F("h_first_hit_X","",100,-600,600);
    TH2F *h_first_hit_XY = new TH2F("h_first_hit_XY","",100,-600,600,100,-600,600);
    TH1F *h_first_hit_Z = new TH1F("h_first_hit_Z","",100,-1000,0);
    
    TH2F *h_ntracks_nPDShits = new TH2F("h_ntracks_nPDShits","",10,0,10,10,0,10);
    TH2F *h_ntracks_nPDShits2 = new TH2F("h_ntracks_nPDShits2","",10,0,10,10,0,10);

    TH2F *h_Zpos_deltaT = new TH2F("h_Zpos_deltaT","",200,-1000,0,200,-2000,4000);
    TH1F *h_corrected_Zpos = new TH1F("h_corrected_Zpos","",100,-1000,0);
    TH1F *h_deltaT = new TH1F("h_deltaT","",100,-700,-500);

    TH1F *h_first_hit_T = new TH1F("h_first_hit_T","",100,0,600);
    TH1F *h_corrected_T = new TH1F("h_corrected_T","",100,0,600);

    int matches = 0;
  

    for (int i=0; i< nEntries; i++) {
	//for (int i=850; i< 945; i++) {
	t->GetEntry(i);
	//std::cout<<"Event="<<i+1<<std::endl;
	//std::cout<<"TPC time="<<TPC_time<<std::endl;
	//std::cout<<"NTracks="<<first_hit_X->size()<<std::endl;

	float Zposf = 0.0;
	float corrected_Zposf = 0.0;
	int deltaT = 0;
	
	int nBeamPDSTrig = 0;
	int nBeamTPCTrac = 0;    
	
	std::cout<<"\nTEST"<<std::endl;
	int orig = 0;
	int corre = 0;

	std::vector<double> Zpositions;
	
	for (size_t j = 0; j<first_hit_X->size(); j++) {
	    Zpositions.push_back(first_hit_Z->at(j));
	}

	for (size_t j = 0; j<Zpositions.size(); j++) {
	    h_first_hit_T->Fill(Zpositions[j]/(-1.6));
	    h_first_hit_Z->Fill(Zpositions[j]);
	    //h_first_hit_Y->Fill(first_hit_Y->at(j));

	    if (Zpositions[j] < -300 && Zpositions[j] > -1000)
		{
		    h_first_hit_X->Fill(first_hit_X->at(j));
		    h_first_hit_XY->Fill(first_hit_X->at(j),first_hit_Y->at(j));
		    nBeamTPCTrac++;
		    orig++;
		    Zposf = Zpositions[j];
		    corrected_Zposf = Zpositions[j];
		
		    for (size_t k = 0; k<PDS_delta_time->size(); k++) {

			std::cout<<"PDS time="<<PDS_delta_time->at(k)<<" PDS energy="<< PDS_energy->at(k)<<std::endl;
			std::cout<<"Z0="<<Zposf<<" Z corr="<<Zposf + PDS_delta_time->at(k)*1.6*1.000000e-3<<std::endl;
			// std::cout<<"deltaZ0="<<deltaZ<<" new dZ="<< fabs(Zposf + PDS_delta_time->at(k)*1.6/1000 + 165)<<std::endl;

		    
			if (Zposf + PDS_delta_time->at(k)*1.6*1.000000e-3 < 100 && Zposf + PDS_delta_time->at(k)*1.6*1.000000e-3 > -350 && PDS_energy->at(k) != -9999) {
			    corrected_Zposf = Zposf + PDS_delta_time->at(k)*1.6*1.000000e-3;
			    corre++;
			    //deltaZ = fabs(Zposf + PDS_delta_time->at(k)*1.6/1000 + 165);
			    std::cout<<"Z="<<Zposf<<" Z corr="<<corrected_Zposf<<std::endl;
			    // std::cout<<"deltaZ="<<deltaZ<<std::endl;
			    break;
			}

			// std::cout<<"Z2="<<Zposf<<" Z corr="<<corrected_Zposf<<std::endl;
			// std::cout<<"deltaZ2="<<deltaZ<<std::endl;

		    
		    }
		
		    if (fabs(Zposf - corrected_Zposf) != 0.0 )
			matches++;
		
		    if (corrected_Zposf < 100 && corrected_Zposf > -350) {
			h_corrected_Zpos->Fill(corrected_Zposf);
			h_corrected_T->Fill(corrected_Zposf/(-1.6));
		    }
		
		}
	    
	}

	std::cout<<"Original2="<<orig<<" Corrected2="<<corre<<std::endl;
    
	if (orig != corre)
	    std::cout<<"FAIL"<<std::endl;
	std::cout<<"event"<<run<<" "<<evt<<std::endl;

	// if (nBeamPDSTrig == 1 && nBeamTPCTrac == 1) {
	//     //std::cout<<"Z="<<Zposf<<" deltaT="<<deltaT<<" corrected Z="<<corrected_Zposf<<std::endl;
	//     h_Zpos_deltaT->Fill(Zposf,deltaT);
	//     h_deltaT->Fill(deltaT);	    
	// }	    
	
	// h_ntracks_nPDShits->Fill(PDS_delta_time->size(),first_hit_X->size());
	// h_ntracks_nPDShits2->Fill(nBeamPDSTrig,nBeamTPCTrac);

    }

    
    gStyle->SetOptStat(0);
    std::cout<<"Matches="<<matches<<" Entries="<<nEntries<<std::endl;
    
    std::cout<<"Original="<<h_first_hit_Z->Integral()<<" Corrected="<<h_corrected_Zpos->Integral()<<std::endl;
    std::cout<<"Efficiency="<<h_corrected_Zpos->Integral()/h_first_hit_Z->Integral()<<std::endl;
    
    // TCanvas *c1 = new TCanvas("c1","",600,600);
    // h_ntracks_nPDShits->Draw("colz");
    // h_ntracks_nPDShits->GetXaxis()->SetTitle("Number PDS triggers");
    // h_ntracks_nPDShits->GetYaxis()->SetTitle("Number of tracks");

    // h_Zpos_deltaT->GetXaxis()->SetTitle("Z position");
    // h_Zpos_deltaT->GetYaxis()->SetTitle("Delta t");
    // h_Zpos_deltaT->Draw("colz");

    TCanvas *c2 = new TCanvas("c2","",600,600);
    h_corrected_Zpos->SetTitle("CAPTAIN Preliminary");
    h_corrected_Zpos->GetXaxis()->SetTitle("Z position [mm]");
    h_corrected_Zpos->Draw();
    h_corrected_Zpos->SetLineWidth(3);
    h_first_hit_X->SetLineColor(kRed);
    h_first_hit_X->SetLineWidth(3);
    
    h_first_hit_XY->Draw("colz");
    //h_deltaT->Draw();

    // TLegend *leg = new TLegend( 0.52, 0.20, 0.78, 0.40 ) ; 
    // leg->AddEntry( h_first_hit_Z, "Original" ) ;
    // leg->AddEntry( h_corrected_Zpos, "Corrected" ) ; 
    // leg->SetFillColor( kWhite ) ;
    // leg->Draw();    
    // c2->Print("corrected_Z.pdf");
    
    TCanvas *c3 = new TCanvas("c3","",600,600);
    h_corrected_T->SetTitle("CAPTAIN Preliminary");
    h_corrected_T->GetXaxis()->SetTitle("Time since first RF [#mus]");
    h_corrected_T->Draw();
    h_corrected_T->SetLineWidth(3);
    h_first_hit_T->SetLineColor(kRed);
    h_first_hit_T->SetLineWidth(3);
    
    h_first_hit_T->Draw("same");
    //h_deltaT->Draw();

    TLegend *leg = new TLegend( 0.52, 0.60, 0.78, 0.80 ) ; 
    leg->AddEntry( h_first_hit_T, "Original" ) ;
    leg->AddEntry( h_corrected_T, "Corrected" ) ; 
    leg->SetFillColor( kWhite ) ;
    leg->Draw();    
    
    c3->Print("corrected_T.pdf");
    c3->Print("corrected_T.png");

    float eff = float(h_corrected_Zpos->Integral()/h_first_hit_Z->Integral());
    
    // h_first_hit_X->Delete();
    // h_first_hit_Z->Delete();
    // h_ntracks_nPDShits->Delete();
    // h_ntracks_nPDShits2->Delete();
    // h_Zpos_deltaT->Delete();
    // h_corrected_Zpos->Delete();
    // h_deltaT->Delete();
    // c1->Close();
    // c2->Close();
    
    //}
    //plotter(argv[0]);
}

float plotterR(int run, bool old){

    //TFile *f = new TFile("../../captSummary/cmt/tracks.root");
    //TTree *t = (TTree*) f->Get("tracks");

    TChain *t = new TChain("tracks");
    TString fname;
    if (old)
	fname.Form("../../captSummary/cmt/tree_old_%d.root",run);
    else if (run < 12117)
	fname.Form("../../captSummary/cmt/tree_newPDS_%d.root",run);
    else
	fname.Form("../../captSummary/cmt/tree_%d.root",run);
    //std::cout<<"FNAME="<<fname.Data();
    //t->Add("../../captSummary/cmt/tracks_12135.root");
    t->Add(fname.Data());
    
    Int_t nEntries = t->GetEntries();

    //std::cout<<nEntries<<std::endl;
  
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
  
    std::vector<double> *first_hit_X = 0;
    std::vector<double> *first_hit_Y = 0;
    std::vector<double> *first_hit_Z = 0;

    std::vector<double> *last_hit_X = 0;
    std::vector<double> *last_hit_Y = 0;
    std::vector<double> *last_hit_Z = 0;

    Long64_t TPC_time = 0;
  
    std::vector<Long64_t> *PDS_RF_time = 0;
    std::vector<Long64_t> *PDS_delta_time = 0;
    std::vector<int> *PDS_beam_trigger = 0;
  
    t->SetBranchAddress("first_hit_X", &first_hit_X);
    t->SetBranchAddress("first_hit_Y", &first_hit_Y);
    t->SetBranchAddress("first_hit_Z", &first_hit_Z);

    t->SetBranchAddress("last_hit_X", &last_hit_X);
    t->SetBranchAddress("last_hit_Y", &last_hit_Y);
    t->SetBranchAddress("last_hit_Z", &last_hit_Z);

    t->SetBranchAddress("TPC_time", &TPC_time);
    t->SetBranchAddress("PDS_RF_time", &PDS_RF_time);
    t->SetBranchAddress("PDS_delta_time", &PDS_delta_time);
    t->SetBranchAddress("PDS_beam_trigger", &PDS_beam_trigger);

    ////////////////////////////////////////////////
    ////////////////////////////////////////////////

    TH1F *h_first_hit_X = new TH1F("h_first_hit_X","",100,-600,600);
    TH1F *h_first_hit_Z = new TH1F("h_first_hit_Z","",100,-1000,0);
    
    TH2F *h_ntracks_nPDShits = new TH2F("h_ntracks_nPDShits","",10,0,10,10,0,10);
    TH2F *h_ntracks_nPDShits2 = new TH2F("h_ntracks_nPDShits2","",10,0,10,10,0,10);

    TH2F *h_Zpos_deltaT = new TH2F("h_Zpos_deltaT","",200,-1000,0,200,-2000,4000);
    TH1F *h_corrected_Zpos = new TH1F("h_corrected_Zpos","",100,-1000,0);
    TH1F *h_deltaT = new TH1F("h_deltaT","",100,-700,-500);

    int matches = 0;
  
    for (int i=0; i< nEntries; i++) {
	//for (int i=850; i< 945; i++) {
	t->GetEntry(i);
	//std::cout<<"Event="<<i+1<<std::endl;
	//std::cout<<"TPC time="<<TPC_time<<std::endl;
	//std::cout<<"NTracks="<<first_hit_X->size()<<std::endl;

	float Zposf = 0.0;
	float corrected_Zposf = 0.0;
	int deltaT = 0;
	
	int nBeamPDSTrig = 0;
	int nBeamTPCTrac = 0;
    
	
	//std::cout<<"TEST"<<std::endl;

	for (size_t j = 0; j<PDS_beam_trigger->size(); j++) {
	    if (PDS_beam_trigger->at(j) == 1) {
		nBeamPDSTrig++;
		deltaT = PDS_delta_time->at(j);
	    }
	}

	std::vector<double> Zpositions;
	
	for (size_t j = 0; j<first_hit_X->size(); j++) {
	    h_first_hit_X->Fill(first_hit_X->at(j));
	    if (first_hit_X->at(j) > last_hit_X->at(j))
		Zpositions.push_back(first_hit_Z->at(j));
	    else
		Zpositions.push_back(last_hit_Z->at(j));
	}

	int orig = 0 ; 
	int corre = 0 ; 

	for (size_t j = 0; j<Zpositions.size(); j++) {
	    if (Zpositions[j] < -300 && Zpositions[j] > -1000)
		{
		    nBeamTPCTrac++;
		    Zposf = Zpositions[j];
		    corrected_Zposf = Zpositions[j];

		    orig++;
		
		    //if (Zposf > -200) continue;
		
		    //float deltaZ = fabs(Zposf + 165);
		
		    for (size_t k = 0; k<PDS_delta_time->size(); k++) {

			//if (PDS_beam_trigger->at(k) == 1) continue;

			// std::cout<<"PDS time="<<PDS_delta_time->at(k)<<std::endl;
			// std::cout<<"Z0="<<Zposf<<" Z corr="<<Zposf + PDS_delta_time->at(k)*1.6/1000<<std::endl;
			// std::cout<<"deltaZ0="<<deltaZ<<" new dZ="<< fabs(Zposf + PDS_delta_time->at(k)*1.6/1000 + 165)<<std::endl;

		    
			if (Zposf + PDS_delta_time->at(k)*1.6/1000 < 100 && Zposf + PDS_delta_time->at(k)*1.6/1000 > -350) {
			    corrected_Zposf = Zposf + PDS_delta_time->at(k)*1.6/1000;
			    //deltaZ = fabs(Zposf + PDS_delta_time->at(k)*1.6/1000 + 165);
			    // std::cout<<"Z="<<Zposf<<" Z corr="<<corrected_Zposf<<std::endl;
			    // std::cout<<"deltaZ="<<deltaZ<<std::endl;
			    corre++;
			    break;
			}

			// std::cout<<"Z2="<<Zposf<<" Z corr="<<corrected_Zposf<<std::endl;
			// std::cout<<"deltaZ2="<<deltaZ<<std::endl;

		    
		    }
		
		    // if (nBeamPDSTrig == 1) {
		    //     //std::cout<<"Z="<<Zposf<<" deltaT="<<deltaT<<std::endl;//" corrected Z="<<corrected_Zposf<<std::endl;		
		    //     //corrected_Zposf = Zpositions[j] + deltaT*1.6/1000;
		    // }

		    if (fabs(Zposf - corrected_Zposf) != 0.0 )
			matches++;
		
		    h_first_hit_Z->Fill(Zposf);
		    //if (corrected_Zposf < 100 && corrected_Zposf > -350)
		    h_corrected_Zpos->Fill(corrected_Zposf);
		    if (orig != corre)
			std::cout<<"FAIL"<<std::endl;
		
		}
	    
	}
    

	// if (nBeamPDSTrig == 1 && nBeamTPCTrac == 1) {
	//     //std::cout<<"Z="<<Zposf<<" deltaT="<<deltaT<<" corrected Z="<<corrected_Zposf<<std::endl;
	//     h_Zpos_deltaT->Fill(Zposf,deltaT);
	//     h_deltaT->Fill(deltaT);	    
	// }	    
	
	// h_ntracks_nPDShits->Fill(PDS_delta_time->size(),first_hit_X->size());
	// h_ntracks_nPDShits2->Fill(nBeamPDSTrig,nBeamTPCTrac);

    }

    

    //std::cout<<"Matches="<<matches<<" Entries="<<nEntries<<std::endl;
    
    //std::cout<<"Original="<<h_first_hit_Z->Integral()<<" Corrected="<<h_corrected_Zpos->Integral()<<std::endl;
    //std::cout<<"Efficiency="<<h_corrected_Zpos->Integral()/h_first_hit_Z->Integral()<<std::endl;
    
    // TCanvas *c1 = new TCanvas("c1","",600,600);
    // h_ntracks_nPDShits->Draw("colz");
    // h_ntracks_nPDShits->GetXaxis()->SetTitle("Number PDS triggers");
    // h_ntracks_nPDShits->GetYaxis()->SetTitle("Number of tracks");

    // h_Zpos_deltaT->GetXaxis()->SetTitle("Z position");
    // h_Zpos_deltaT->GetYaxis()->SetTitle("Delta t");
    // h_Zpos_deltaT->Draw("colz");

    // TCanvas *c2 = new TCanvas("c2","",600,600);
    // // h_ntracks_nPDShits2->Draw("colz");
    // // h_ntracks_nPDShits2->GetXaxis()->SetTitle("Number PDS triggers");
    // // h_ntracks_nPDShits2->GetYaxis()->SetTitle("Number of tracks");
    // h_corrected_Zpos->Draw();
    // h_first_hit_Z->SetLineColor(kRed);
    // h_first_hit_Z->Draw("same");
    // //h_deltaT->Draw();

    // TLegend *leg = new TLegend( 0.52, 0.20, 0.78, 0.40 ) ; 
    // leg->AddEntry( h_first_hit_Z, "Original" ) ;
    // leg->AddEntry( h_corrected_Zpos, "Corrected" ) ; 
    // leg->SetFillColor( kWhite ) ;
    // leg->Draw();    

    float eff = float(h_corrected_Zpos->Integral()/h_first_hit_Z->Integral());
    
    h_first_hit_X->Delete();
    h_first_hit_Z->Delete();
    h_ntracks_nPDShits->Delete();
    h_ntracks_nPDShits2->Delete();
    h_Zpos_deltaT->Delete();
    h_corrected_Zpos->Delete();
    h_deltaT->Delete();
    // c1->Close();
    // c2->Close();
    
    //}
    //plotter(argv[0]);
    
    return eff;
}


void plotter() {
    /*
      TH1F *h = new TH1F("h","",50,0,50);
      for (int i=0;i<44;i++) {
      if (i == 17) continue;
      //float eff = plotterR(12100+i,false);
      float eff = plotterFull(12100+i,false);
      std::cout<<"Eff2="<<i<<" "<<eff<<std::endl;
      h->SetBinContent(i,eff);
      }
      h->Draw();
    
      TH1F *h2 = new TH1F("h2","",50,0,50);
      for (int i=0;i<44;i++) {
      if (i == 17) continue;
      //float eff = plotterR(12100+i,true);
      float eff = plotterFull(12100+i,true);
      std::cout<<"Eff2="<<i<<" "<<eff<<std::endl;
      h2->SetBinContent(i,eff);
      }
      h2->SetLineColor(kRed);
      h2->Draw("same");


      TLegend *leg = new TLegend( 0.52, 0.20, 0.78, 0.40 ) ; 
      leg->AddEntry( h2, "No shift" ) ;
      leg->AddEntry( h, "Shifted" ) ; 
      leg->SetFillColor( kWhite ) ;
      leg->Draw();    
    */
    plotterS("outalltree.root");
    
}
