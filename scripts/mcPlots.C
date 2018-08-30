class Selector
{
public:
    enum tag_t {
	Left,
	Right,
	Center
    };
};

void mcPlots() {

    TChain *c_data = new TChain("tracks");
    TChain *c_mc = new TChain("tracks");
    c_data->Add("nb_mtpc_spl_golden_flat.root");
    c_mc->Add("tracks_proton_beam_extended_remove_large.root");
    //c_mc->Add("mc_tree.root");

    

    std::vector<double> *track_length_data = 0;
    std::vector<double> *track_length_mc = 0;

    std::vector<double> *first_hit_X_data = 0;
    std::vector<double> *first_hit_X_mc = 0;

    std::vector<double> *first_hit_Y_data = 0;
    std::vector<double> *first_hit_Y_mc = 0;

    std::vector<double> *first_hit_Z_data = 0;
    std::vector<double> *corrected_first_hit_Z_data = 0;
    std::vector<double> *first_hit_Z_mc = 0;
    std::vector<double> *truth_particle_E_mc = 0;

    c_data->SetBranchAddress("track_length",&track_length_data);
    c_data->SetBranchAddress("first_hit_X",&first_hit_X_data);
    c_data->SetBranchAddress("first_hit_Y",&first_hit_Y_data);
    c_data->SetBranchAddress("first_hit_Z",&first_hit_Z_data);
    c_data->SetBranchAddress("corrected_first_hit_Z",&corrected_first_hit_Z_data);

    c_mc->SetBranchAddress("track_length",&track_length_mc);
    c_mc->SetBranchAddress("first_hit_X",&first_hit_X_mc);
    c_mc->SetBranchAddress("first_hit_Y",&first_hit_Y_mc);
    c_mc->SetBranchAddress("first_hit_Z",&first_hit_Z_mc);
    c_mc->SetBranchAddress("truth_particle_E",&truth_particle_E_mc);
    
    TH1F* hs_mc_0 = new TH1F("hs_mc_0","",40,-500,500);
    TH1F* hs_data_0 = new TH1F("hs_data_0","",40,-500,500);

    TH1F* hs_mc_1 = new TH1F("hs_mc_1","",40,0,1000);
    TH1F* hs_data_1 = new TH1F("hs_data_1","",40,0,1000);

    TH1F* hs_mc_2 = new TH1F("hs_mc_2","",40,-500,100);
    TH1F* hs_data_2 = new TH1F("hs_data_2","",40,-500,100);

    TH1F* hs_mc_3 = new TH1F("hs_mc_3","",40,0,700);
    TH1F* hs_data_3 = new TH1F("hs_data_3","",40,0,700);

    TH1F* hs_mc_0_c = new TH1F("hs_mc_0_c","",40,-500,500);
    TH1F* hs_mc_2_c = new TH1F("hs_mc_2_c","",40,-500,100);

    // Y position
    TH1F* hs_mc_4 = new TH1F("hs_mc_4","",40,-300,100);
    TH1F* hs_data_4 = new TH1F("hs_data_4","",40,-300,100);

    TH1F* hs_mc_5 = new TH1F("hs_mc_5","",40,-300,100);
    TH1F* hs_data_5 = new TH1F("hs_data_5","",40,-300,100);

    // Z position
    TH1F* hs_mc_6 = new TH1F("hs_mc_6","",40,-400,100);
    TH1F* hs_data_6 = new TH1F("hs_data_6","",40,-400,100);

    TH1F* hs_mc_7 = new TH1F("hs_mc_7","",40,-400,100);
    TH1F* hs_data_7 = new TH1F("hs_data_7","",40,-400,100);

    // X position with E cuts
    // E > 60 MeV
    double Emin1 = 60;
    TH1F* hs_mc_8 = new TH1F("hs_mc_8","",40,-500,100);
    TH1F* hs_data_8 = new TH1F("hs_data_8","",40,-500,100);
    // E > 100 MeV
    double Emin2 = 100;
    TH1F* hs_mc_9 = new TH1F("hs_mc_9","",40,-500,100);
    TH1F* hs_data_9 = new TH1F("hs_data_9","",40,-500,100);
    // E > 200 MeV
    double Emin3 = 200;
    TH1F* hs_mc_10 = new TH1F("hs_mc_10","",40,-500,100);
    TH1F* hs_data_10 = new TH1F("hs_data_10","",40,-500,100);
    // E > 60 MeV and < 100 MeV
    TH1F* hs_mc_11 = new TH1F("hs_mc_11","",40,-500,100);
    TH1F* hs_data_11 = new TH1F("hs_data_11","",40,-500,100);
    // E > 100 MeV and < 200 MeV
    TH1F* hs_mc_12 = new TH1F("hs_mc_12","",40,-500,100);
    TH1F* hs_data_12 = new TH1F("hs_data_12","",40,-500,100);

    
    int nEntries_data = c_data->GetEntries();
    int nEntries_mc = c_mc->GetEntries();

    double zmin = -300;
    double zmax = -1;
    
    //std::cout<<"Data="<<nEntries_data<<" MC="<<nEntries_mc<<std::endl;
    printf("Data=%d MC=%d\n",nEntries_data,nEntries_mc);
    
    for (int iev = 0; iev < nEntries_data; iev++) {
	c_data->GetEntry(iev);
	for (size_t j = 0; j < track_length_data->size(); j++) {
	    if ( (corrected_first_hit_Z_data->at(j) > zmax || corrected_first_hit_Z_data->at(j) < zmin) ) continue;
	    if (first_hit_X_data->at(j) > 400 || first_hit_X_data->at(j) < -400) continue;
	    //if (!inBeamXY(first_hit_X_data->at(j),first_hit_Y_data->at(j))) continue;
	    if (first_hit_Y_data->at(j) > 0 || first_hit_Y_data->at(j) < -200) continue;

		hs_data_0->Fill(first_hit_X_data->at(j));
		hs_data_1->Fill(track_length_data->at(j));
		hs_data_4->Fill(first_hit_Y_data->at(j));
		hs_data_6->Fill(corrected_first_hit_Z_data->at(j));
		if (first_hit_X_data->at(j) > 0) continue;
		hs_data_2->Fill(first_hit_X_data->at(j));
		hs_data_3->Fill(track_length_data->at(j));
		hs_data_5->Fill(first_hit_Y_data->at(j));
		hs_data_7->Fill(corrected_first_hit_Z_data->at(j));
	    }
	
    }
    
    for (int iev = 0; iev < nEntries_mc; iev++) {
	c_mc->GetEntry(iev);
	for (size_t j = 0; j < track_length_mc->size(); j++) {
	    if (first_hit_Z_mc->at(j) > zmax || first_hit_Z_mc->at(j) < zmin) continue;
	    if (first_hit_X_mc->at(j) > 400 || first_hit_X_mc->at(j) < -400) continue;
	    //if (!inBeamXY(first_hit_X_mc->at(j),first_hit_Y_mc->at(j))) continue;
	    if (first_hit_Y_mc->at(j) > 0 || first_hit_Y_mc->at(j) < -200) continue;

	    hs_mc_0->Fill(first_hit_X_mc->at(j));
	    hs_mc_1->Fill(track_length_mc->at(j));
	    hs_mc_4->Fill(first_hit_Y_mc->at(j));
	    hs_mc_6->Fill(first_hit_Z_mc->at(j));
	    if (first_hit_X_mc->at(j) > 0) continue;
	    hs_mc_2->Fill(first_hit_X_mc->at(j));
	    hs_mc_3->Fill(track_length_mc->at(j));
	    hs_mc_5->Fill(first_hit_Y_mc->at(j));
	    hs_mc_7->Fill(first_hit_Z_mc->at(j));
	    bool Epass1 = true;
	    bool Epass2 = true;
	    bool Epass3 = true;
	    bool Epass4 = true;
	    bool Epass5 = true;
	    for (int k = 0; k < truth_particle_E_mc->size(); k++) {
		if (truth_particle_E_mc->at(k) < Emin1)
		    Epass1 = false;
		if (truth_particle_E_mc->at(k) < Emin2)
		    Epass2 = false;
		if (truth_particle_E_mc->at(k) < Emin3)
		    Epass3 = false;
		if (truth_particle_E_mc->at(k) < Emin1 || truth_particle_E_mc->at(k) > Emin2)
		    Epass4 = false;
		if (truth_particle_E_mc->at(k) < Emin2 || truth_particle_E_mc->at(k) > Emin3)
		    Epass5 = false;
	    }
	    if (Epass1)
		hs_mc_8->Fill(first_hit_X_mc->at(j));
	    if (Epass2)
		hs_mc_9->Fill(first_hit_X_mc->at(j));
	    if (Epass3)
		hs_mc_10->Fill(first_hit_X_mc->at(j));
	    if (Epass4)
		hs_mc_11->Fill(first_hit_X_mc->at(j));
	    if (Epass5)
		hs_mc_12->Fill(first_hit_X_mc->at(j));
	    
	}
    }

    TH1F* h_ratio_1 = hs_data_1->Clone();
    h_ratio_1->Divide(hs_mc_1);
    //h_ratio_1->Scale(1./h_ratio_1->Integral());
    TH1F* h_ratio_3 = hs_data_3->Clone();
    h_ratio_3->Divide(hs_mc_3);
    //h_ratio_3->Scale(1./h_ratio_3->Integral());

    drawPlots( hs_mc_0, hs_data_0, "X position (mm)",  "firsthitX", Selector::Left); 
    drawPlots( hs_mc_1, hs_data_1, "Track length (mm)",  "tracklength", Selector::Right); 
    drawPlots( hs_mc_2, hs_data_2, "X position (mm)",  "firsthitX_reduced", Selector::Left); 
    drawPlots( hs_mc_3, hs_data_3, "Track length (mm)",  "tracklength_reduced", Selector::Right); 

    drawPlots( hs_mc_4, hs_data_4, "Y position (mm)",  "firsthitY", Selector::Left); 
    drawPlots( hs_mc_6, hs_data_6, "Z position (mm)",  "firsthitZ", Selector::Left); 
    drawPlots( hs_mc_5, hs_data_5, "Y position (mm)",  "firsthitY_reduced", Selector::Left); 
    drawPlots( hs_mc_7, hs_data_7, "Z position (mm)",  "firsthitZ_reduced", Selector::Left); 


    drawPlotsNoRatio( hs_mc_0, hs_mc_0, "X position (mm)",  "firsthitX_noratio", Selector::Left); 
    drawPlotsNoRatio( hs_mc_1, hs_mc_1, "Track length (mm)",  "tracklength_noratio", Selector::Right); 
    drawPlotsNoRatio( hs_mc_4, hs_mc_4, "Y position (mm)",  "firsthitY_noratio", Selector::Left); 
    drawPlotsNoRatio( hs_mc_6, hs_mc_6, "Z position (mm)",  "firsthitZ_noratio", Selector::Left); 
    
    // drawPlotsNoRatio( hs_mc_8, hs_mc_8, "X position (mm)",  "firsthitX_noratio_60MeV", Selector::Left); 
    // drawPlotsNoRatio( hs_mc_9, hs_mc_9, "X position (mm)",  "firsthitX_noratio_100MeV", Selector::Left); 
    // drawPlotsNoRatio( hs_mc_10, hs_mc_10, "X position (mm)",  "firsthitX_noratio_200MeV", Selector::Left); 
    // drawPlotsNoRatio( hs_mc_11, hs_mc_11, "X position (mm)",  "firsthitX_noratio_60_100MeV", Selector::Left); 
    // drawPlotsNoRatio( hs_mc_12, hs_mc_12, "X position (mm)",  "firsthitX_noratio_100_200MeV", Selector::Left); 

    
    for (int iev = 0; iev < nEntries_mc; iev++) {
	c_mc->GetEntry(iev);
	for (size_t j = 0; j < track_length_mc->size(); j++) {

	    if (first_hit_Z_mc->at(j) > zmax || first_hit_Z_mc->at(j) < zmin) continue;
	    if (first_hit_X_mc->at(j) > 400 || first_hit_X_mc->at(j) < -400) continue;
	    //if (!inBeamXY(first_hit_X_mc->at(j),first_hit_Y_mc->at(j))) continue;
	    if (first_hit_Y_mc->at(j) > 0 || first_hit_Y_mc->at(j) < -200) continue;

	    float weight_1 = h_ratio_1->GetBinContent(h_ratio_1->GetXaxis()->FindBin(track_length_mc->at(j)));
	    float weight_3 = h_ratio_3->GetBinContent(h_ratio_3->GetXaxis()->FindBin(track_length_mc->at(j)));
	    
	    hs_mc_0_c->Fill(first_hit_X_mc->at(j), weight_1);
	    if (first_hit_X_mc->at(j) > 0) continue;
	    hs_mc_2_c->Fill(first_hit_X_mc->at(j), weight_3);
	}
    }

    // hs_mc_2_c->Draw();
    // hs_data_2->SetLineColor(kRed);
    // hs_data_2->Draw("same");

    
    drawPlots( hs_mc_0_c, hs_data_0, "X position (mm)",  "firsthitXcorr", Selector::Left); 
    drawPlots( hs_mc_2_c, hs_data_2, "X position (mm)",  "firsthitXcorr_reduced", Selector::Left); 

    // h_ratio_1->Draw();
    // h_ratio_3->SetLineColor(kRed);
    // h_ratio_3->Draw("same");

    
    
}

void drawPlots(TH1* hs_mc, TH1* hs_data, TString xtitle, TString fname, Selector::tag_t channel){
    
    if (channel == Selector::Left) {
	TLegend *leg = new TLegend( 0.18, 0.60, 0.52, 0.80 ) ;
	std::cout<<"LEFT"<<std::endl;
    }
    else {
	TLegend *leg = new TLegend( 0.62, 0.50, 0.98, 0.70 ) ;
	std::cout<<"NOT LEFT"<<std::endl;
    }
    
    leg->AddEntry( hs_mc, "Reweighted simulation" ) ; 
    leg->AddEntry( hs_data, "Data");
    leg->SetFillColor( kWhite ) ; 

    printf("int mc=%f int data=%f\n",hs_mc->Integral(),hs_data->Integral());
    
    hs_mc->Scale(float(hs_data->Integral())/float(hs_mc->Integral()));

    TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
    THStack* th = new THStack();
    hs_mc->SetFillColor(kYellow);
    // hs_mc->SetFillColor(kBlue);
    // hs_mc->SetFillColor(kGreen);
    hs_mc->SetLineColor(kBlack);


    th->Add(hs_mc);

    hs_data->SetMarkerStyle(20);

    Double_t eps = 0.001;
    TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
    TPad* p2 = new TPad("p2","p2",0,0,1,0.25+eps,0); p2->Draw();
    p1->SetBottomMargin(0);
    p2->SetTopMargin(0);   
    p2->SetBottomMargin(0.5);   
    p1->cd();
    gPad->SetTickx();
    gPad->SetTicky();
    hs_data->SetStats(0);
    TH1F *ratio = (TH1F*)hs_data->Clone();
    th->SetTitle("CAPTAIN Preliminary");
    hs_data->SetTitle("CAPTAIN Preliminary");
    th->Draw("histo");
    hs_data->Draw("epsame");

    // hs_data->Draw("ep");
    // th->Draw("histo same");
    // hs_data->Draw("epsame");

    TH1F *errors = (TH1F*)(th->GetStack()->Last())->Clone();
    // TFile *scale_error_fn = new TFile("data/scale_syst.root");
    // TFile *pdf_error_fn = new TFile("data/pdf_syst.root");
    // TFile *fact_renorm_error_fn = new TFile("data/fact_renorm_syst.root");

    TString hname = "";
    
    // TH1F *scale_error_hist = (TH1F*) scale_error_fn->Get(hname);
    // TH1F *pdf_error_hist = (TH1F*) pdf_error_fn->Get(hname);
    // TH1F *fact_renorm_error_hist = (TH1F*) fact_renorm_error_fn->Get(hname);

    TF1 * scale_error_funcEE = new TF1("scale_error_funcEE","2.44792e-05*x+1.45510e-01");
    TF1 * pdf_error_funcEE = new TF1("pdf_error_funcEE","1.05646e-06*x+1.33917e-01");
    TF1 * fact_renorm_error_funcEE = new TF1("fact_renorm_error_funcEE","2.97192e-06*x+4.24292e-02");
    
    TF1 * scale_error_funcMuMu = new TF1("scale_error_funcMuMu","5.01588e-05*x+5.44939e-02");
    TF1 * pdf_error_funcMuMu = new TF1("pdf_error_funcMuMu","1.21116e-06*x+1.33917e-01");
    TF1 * fact_renorm_error_funcMuMu = new TF1("fact_renorm_error_funcMuMu","5.54490e-06*x+1.13813e-01");
      
    for(int i = 0;i<errors->GetNbinsX()+1;i++){
	float errorSum = 0;
	float errorSyst = 0;
	// if(channel == Selector::EE){
	// 	errorSum = TMath::Sqrt(
	// 			       (scale_error_hist->GetBinContent(scale_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 			       (scale_error_hist->GetBinContent(scale_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
	// 			       (pdf_error_hist->GetBinContent(pdf_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 			       (pdf_error_hist->GetBinContent(pdf_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
	// 			       (fact_renorm_error_hist->GetBinContent(fact_renorm_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 			       (fact_renorm_error_hist->GetBinContent(fact_renorm_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
			       
	// 			       // (scale_error_funcEE->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))*
	// 			       // (scale_error_funcEE->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))+
	// 			       // (pdf_error_funcEE->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))*
	// 			       // (pdf_error_funcEE->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))+
	// 			       // (fact_renorm_error_funcEE->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))*
	// 			       // (fact_renorm_error_funcEE->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))+
			       
	// 			       (errors->GetBinError(i)*errors->GetBinError(i)));

	// 	errorSyst = TMath::Sqrt(
	// 				(scale_error_hist->GetBinContent(scale_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 				(scale_error_hist->GetBinContent(scale_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
	// 				(pdf_error_hist->GetBinContent(pdf_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 				(pdf_error_hist->GetBinContent(pdf_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
	// 				(fact_renorm_error_hist->GetBinContent(fact_renorm_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 				(fact_renorm_error_hist->GetBinContent(fact_renorm_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))));
	// }
	// else if(channel == Selector::MuMu){
	// 	errorSum = TMath::Sqrt(
	// 			       (scale_error_hist->GetBinContent(scale_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 			       (scale_error_hist->GetBinContent(scale_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
	// 			       (pdf_error_hist->GetBinContent(pdf_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 			       (pdf_error_hist->GetBinContent(pdf_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
	// 			       (fact_renorm_error_hist->GetBinContent(fact_renorm_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 			       (fact_renorm_error_hist->GetBinContent(fact_renorm_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +

	// 			       // (scale_error_funcMuMu->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))*
	// 			       // (scale_error_funcMuMu->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))+
	// 			       // (pdf_error_funcMuMu->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))*
	// 			       // (pdf_error_funcMuMu->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))+
	// 			       // (fact_renorm_error_funcMuMu->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))*
	// 			       // (fact_renorm_error_funcMuMu->Eval(hs_DY->GetBinCenter(i))*hs_DY->GetBinContent(i))+

	// 			       (errors->GetBinError(i)*errors->GetBinError(i)));

	// 	errorSyst = TMath::Sqrt(
	// 				(scale_error_hist->GetBinContent(scale_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 				(scale_error_hist->GetBinContent(scale_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
	// 				(pdf_error_hist->GetBinContent(pdf_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 				(pdf_error_hist->GetBinContent(pdf_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))) +
	// 				(fact_renorm_error_hist->GetBinContent(fact_renorm_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i)))*
	// 				(fact_renorm_error_hist->GetBinContent(fact_renorm_error_hist->FindBin(hs_DY->GetBinCenter(i)))*(hs_DY->GetBinContent(i))));

	// }
      
	//if(xtitle == "Mlljj [GeV]")
	{
	    //cout<<i<<" "<<hs_data->GetBinContent(i)<<" "<<errors->GetBinContent(i)<<" "<<errors->GetBinError(i)<<" "<<errorSyst<<" "<<errorSum<<endl;
	    //errors->SetBinError(i,errorSum);
	    //hs_mc->SetBinError(i,errorSum);
	    //cout<<"Bin="<<i<<" "<<hs_data->GetBinContent(i)<<endl;
	    //cout<<"Bin="<<i<<" "<<hs_mc->GetBinContent(i)<<endl;
	}
    }
    
    
    errors->SetLineColor(0);
    errors->SetFillColor(1);
    errors->SetFillStyle(3254);
    errors->Draw("E2 same");
    TString ytitle = "Tracks/(";
    ytitle += (th->GetXaxis()->GetNbins());
    ytitle += ")";
    th->GetYaxis()->SetTitle(ytitle.Data());
    th->GetXaxis()->SetTitle(xtitle.Data());
    hs_data->GetXaxis()->SetTitle(xtitle.Data());
    //hs_data->GetYaxis()->SetTitle(ytitle.Data());
    
    ratio->GetXaxis()->SetTitle(xtitle.Data());
    //ths[icanvas]->GetXaxis()->SetTickSize(1.0);
    //ths[icanvas]->GetXaxis()->SetTitleSize(0.1);
    //ratio->GetXaxis()->SetTickSize(0.40);
    //ratio->GetXaxis()->SetTitleSize(0.2);
    ratio->SetLabelSize(0.1,"x");
    leg->Draw(); 
    mycanvas->cd();
    p2->cd();
    p2->SetGridy(); 
    ratio->Sumw2();
    ratio->SetStats(0);
  
    ratio->Divide(hs_mc);
    ratio->SetMarkerStyle(21);
    ratio->SetMarkerSize(0.5);
    ratio->SetLabelSize(0.1,"y");
    ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->Draw("p");
  
    float xmax = ratio->GetXaxis()->GetXmax();
    float xmin = ratio->GetXaxis()->GetXmin();
    TF1 *f1 = new TF1("f1","1",xmin,xmax);
    ratio->Draw("p");
    f1->Draw("same");
    //std::cout<<"OFFSET="<<ratio->GetXaxis()->GetTitleOffset()<<std::endl;
    ratio->GetXaxis()->SetTitle(xtitle.Data());
    ratio->GetXaxis()->SetTitleSize(0.12);
  
    mycanvas->cd();

    TString fn = "";

    fn = "webplots/dataMC/"+fname;
    
    mycanvas->Print((fn+".pdf").Data());
    mycanvas->Print((fn+".png").Data());
    mycanvas->Print((fn+".root").Data());
    p1->SetLogy();
    mycanvas->Print((fn+"_log.pdf").Data());
    mycanvas->Print((fn+"_log.png").Data());

    mycanvas->Close();

}

void drawPlotsNoRatio(TH1* hs_mc, TH1* hs_data, TString xtitle, TString fname, Selector::tag_t channel){
    if (channel == Selector::Left) {
	TLegend *leg = new TLegend( 0.18, 0.60, 0.42, 0.70 ) ;
	std::cout<<"LEFT"<<std::endl;
    }
    else {
	TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.60 ) ;
	std::cout<<"NOT LEFT"<<std::endl;
    }
    
    leg->AddEntry( hs_mc, "Simulated protons" ) ; 
    //leg->AddEntry( hs_data, "Simulated protons");
    leg->SetFillColor( kWhite ) ; 

    printf("int mc=%f int data=%f\n",hs_mc->Integral(),hs_data->Integral());
    
    hs_mc->Scale(float(hs_data->Integral())/float(hs_mc->Integral()));

    TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
    THStack* th = new THStack();
    hs_mc->SetFillColor(kYellow);
    hs_mc->SetFillColor(kBlue);
    hs_mc->SetFillColor(kGreen);
    hs_mc->SetLineColor(kBlack);


    th->Add(hs_mc);

    hs_data->SetMarkerStyle(20);

    Double_t eps = 0.001;
    gPad->SetTickx();
    gPad->SetTicky();
    hs_data->SetStats(0);
    TH1F *ratio = (TH1F*)hs_data->Clone();
    th->SetTitle("CAPTAIN Preliminary");
    hs_data->SetTitle("CAPTAIN Preliminary");
    //th->Draw("histo sameaxis");
    th->Draw("histo");
    hs_data->Draw("sameaxis");
    
    TF1 *f1 = new TF1("f1","expo",-400,0);
    //hs_data->Fit("f1","R");
    
    // hs_data->Draw("ep");
    // th->Draw("histo same");
    // hs_data->Draw("epsame");

    TH1F *errors = (TH1F*)(th->GetStack()->Last())->Clone();
    // TFile *scale_error_fn = new TFile("data/scale_syst.root");
    // TFile *pdf_error_fn = new TFile("data/pdf_syst.root");
    // TFile *fact_renorm_error_fn = new TFile("data/fact_renorm_syst.root");

    TString hname = "";
    
    // TH1F *scale_error_hist = (TH1F*) scale_error_fn->Get(hname);
    // TH1F *pdf_error_hist = (TH1F*) pdf_error_fn->Get(hname);
    // TH1F *fact_renorm_error_hist = (TH1F*) fact_renorm_error_fn->Get(hname);

    TF1 * scale_error_funcEE = new TF1("scale_error_funcEE","2.44792e-05*x+1.45510e-01");
    TF1 * pdf_error_funcEE = new TF1("pdf_error_funcEE","1.05646e-06*x+1.33917e-01");
    TF1 * fact_renorm_error_funcEE = new TF1("fact_renorm_error_funcEE","2.97192e-06*x+4.24292e-02");
    
    TF1 * scale_error_funcMuMu = new TF1("scale_error_funcMuMu","5.01588e-05*x+5.44939e-02");
    TF1 * pdf_error_funcMuMu = new TF1("pdf_error_funcMuMu","1.21116e-06*x+1.33917e-01");
    TF1 * fact_renorm_error_funcMuMu = new TF1("fact_renorm_error_funcMuMu","5.54490e-06*x+1.13813e-01");
      
    for(int i = 0;i<errors->GetNbinsX()+1;i++){
	float errorSum = 0;
	float errorSyst = 0;
    }
    
    
    errors->SetLineColor(0);
    errors->SetFillColor(1);
    errors->SetFillStyle(3254);
    //errors->Draw("E2 same");
    TString ytitle = "Tracks/(";
    ytitle += (th->GetXaxis()->GetNbins());
    ytitle += ")";
    th->GetYaxis()->SetTitle(ytitle.Data());
    th->GetXaxis()->SetTitle(xtitle.Data());
    th->GetYaxis()->SetTitleOffset(1.5);
    hs_data->GetXaxis()->SetTitle(xtitle.Data());
    //hs_data->GetYaxis()->SetTitle(ytitle.Data());
    
    leg->Draw(); 
    mycanvas->cd();
    TString fn = "";

    fn = "webplots/dataMC/"+fname;
    
    mycanvas->Print((fn+".pdf").Data());
    mycanvas->Print((fn+".png").Data());
    mycanvas->Print((fn+".root").Data());
    mycanvas->SetLogy();
    mycanvas->Print((fn+"_log.pdf").Data());
    mycanvas->Print((fn+"_log.png").Data());

    mycanvas->Close();
}

bool inBeamXY(double x, double y){
    
    double a = -0.1188;
    double b = -54.9894;
    double width = 54/2;
    //it should be 54/2 to get good separated peaks
    double dist=fabs(a*x-y+b)/sqrt(a*a+1);
 
    if(dist<width ) return true;
    //if(dist<width)return true;
    return false;
    
}

