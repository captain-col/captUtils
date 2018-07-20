void fitXsec() {

    TChain *c_data = new TChain("tracks");
    TChain *c_mc = new TChain("tracks");
    c_data->Add("tracks_golden.root");
    c_mc->Add("tracks_proton_beam_extended_remove_more.root");

    std::vector<double> *track_length_data = 0;
    std::vector<double> *track_length_mc = 0;

    std::vector<double> *first_hit_X_data = 0;
    std::vector<double> *first_hit_X_mc = 0;

    std::vector<double> *first_hit_Y_data = 0;
    std::vector<double> *first_hit_Y_mc = 0;

    std::vector<double> *first_hit_Z_data = 0;
    std::vector<double> *first_hit_Z_mc = 0;

    c_data->SetBranchAddress("track_length",&track_length_data);
    c_data->SetBranchAddress("first_hit_X",&first_hit_X_data);
    c_data->SetBranchAddress("first_hit_Y",&first_hit_Y_data);
    c_data->SetBranchAddress("first_hit_Z",&first_hit_Z_data);

    c_mc->SetBranchAddress("track_length",&track_length_mc);
    c_mc->SetBranchAddress("first_hit_X",&first_hit_X_mc);
    c_mc->SetBranchAddress("first_hit_Y",&first_hit_Y_mc);
    c_mc->SetBranchAddress("first_hit_Z",&first_hit_Z_mc);
    
    TH1F* hs_data_0 = new TH1F("hs_data_0","",40,-500,500);
    TH1F* hs_data_1 = new TH1F("hs_data_1","",40,0,1000);
    TH1F* hs_data_2 = new TH1F("hs_data_2","",40,-500,100);
    TH1F* hs_data_3 = new TH1F("hs_data_3","",40,0,700);

    int nEntries_data = c_data->GetEntries();
    int nEntries_mc = c_mc->GetEntries();

    for (int iev = 0; iev < nEntries_data; iev++) {
	c_data->GetEntry(iev);
	for (size_t j = 0; j < track_length_data->size(); j++) {
	    if (first_hit_Z_data->at(j) > 0 || first_hit_Z_data->at(j) < -350 || first_hit_Y_data->at(j) > 0 || first_hit_Y_data->at(j) < -200) continue;
	    if (first_hit_X_data->at(j) > 400 || first_hit_X_data->at(j) < -400) continue;
	    hs_data_0->Fill(first_hit_X_data->at(j));
	    hs_data_1->Fill(track_length_data->at(j));
	    if (first_hit_X_data->at(j) > 0) continue;
	    hs_data_2->Fill(first_hit_X_data->at(j));
	    hs_data_3->Fill(track_length_data->at(j));

	}
    }

    TH1F * h_scan_0 = new TH1F("h_scan_0","",100,0,1);
    TH1F * h_scan_1 = new TH1F("h_scan_1","",100,0,1);
    TH1F * h_scan_2 = new TH1F("h_scan_2","",100,0,1);
    TH1F * h_scan_3 = new TH1F("h_scan_3","",100,0,1);
    
    for (int nScan = 0; nScan < 100; nScan++) {
	TH1F* hs_mc_0 = new TH1F("hs_mc_0","",40,-500,500);
	TH1F* hs_mc_1 = new TH1F("hs_mc_1","",40,0,1000);
        TH1F* hs_mc_2 = new TH1F("hs_mc_2","",40,-500,100);
	TH1F* hs_mc_3 = new TH1F("hs_mc_3","",40,0,700);

	float weight = TMath::Exp(-0.01*(nScan));
	//std::cout<<"weight="<<weight<<std::endl;
	//weight = 1.0;
	
	for (int iev = 0; iev < nEntries_mc; iev++) {
	    c_mc->GetEntry(iev);
	    for (size_t j = 0; j < track_length_mc->size(); j++) {
		if (first_hit_Z_mc->at(j) > 0 || first_hit_Z_mc->at(j) < -350 || first_hit_Y_mc->at(j) > 0 || first_hit_Y_mc->at(j) < -200) continue;
		if (first_hit_X_mc->at(j) > 400 || first_hit_X_mc->at(j) < -400) continue;
		    
		hs_mc_0->Fill(first_hit_X_mc->at(j), 3*weight);

		hs_mc_1->Fill(track_length_mc->at(j));	
		if (first_hit_X_mc->at(j) > 0) continue;
		hs_mc_2->Fill(first_hit_X_mc->at(j), 3*weight);
		hs_mc_3->Fill(track_length_mc->at(j));
	    }
	}

	float chi2_0 = 0.0;
	for (int b = 0; b < hs_data_0->GetXaxis()->GetNbins(); b++) {
	    //std::cout<<hs_data_0->GetBinContent(b) <<" " << hs_mc_0->GetBinContent(b)<<std::endl;
	    chi2_0 += (hs_data_0->GetBinContent(b) - hs_mc_0->GetBinContent(b))*(hs_data_0->GetBinContent(b) - hs_mc_0->GetBinContent(b));
	}
	h_scan_0->SetBinContent(nScan+1, chi2_0);
	float chi2_1 = 0.0;
	for (int b = 0; b < hs_data_1->GetXaxis()->GetNbins(); b++) {
	    chi2_1 += (hs_data_1->GetBinContent(b) - hs_mc_1->GetBinContent(b))*(hs_data_1->GetBinContent(b) - hs_mc_1->GetBinContent(b));
	}
	h_scan_1->SetBinContent(nScan+1, chi2_1);
	float chi2_2 = 0.0;
	for (int b = 0; b < hs_data_2->GetXaxis()->GetNbins(); b++) {
	    chi2_2 += (hs_data_2->GetBinContent(b) - hs_mc_2->GetBinContent(b))*(hs_data_2->GetBinContent(b) - hs_mc_2->GetBinContent(b));
	}
	h_scan_2->SetBinContent(nScan+1, chi2_2);
	float chi2_3 = 0.0;
	for (int b = 0; b < hs_data_3->GetXaxis()->GetNbins(); b++) {
	    chi2_3 += (hs_data_3->GetBinContent(b) - hs_mc_3->GetBinContent(b))*(hs_data_3->GetBinContent(b) - hs_mc_3->GetBinContent(b));
	}
	h_scan_3->SetBinContent(nScan+1, chi2_3);

	hs_mc_0->Delete();
	hs_mc_1->Delete();
	hs_mc_2->Delete();
	hs_mc_3->Delete();
    }

    gStyle->SetOptStat(0);
    h_scan_0->DrawNormalized();
    h_scan_2->SetLineColor(kRed);
    h_scan_2->DrawNormalized("same");



    TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.70 ) ;
    leg->AddEntry( h_scan_0, "Full range" ) ; 
    leg->AddEntry( h_scan_2, "Restricted range");
    leg->SetFillColor( kWhite ) ; 
    leg->Draw(); 

    
}
