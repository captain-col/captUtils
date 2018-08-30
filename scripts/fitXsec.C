bool inBeamXY(double x, double y){
    double a = -0.1188;
    double b = -54.9894;
    double width = 54/2;
    double dist=fabs(a*x-y+b)/sqrt(a*a+1); 
    if(dist<width ) return true;
    return false;    
}


void fitXsecChi2() {

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

    int maxScan = 100;
    
    TH1F * h_scan_0 = new TH1F("h_scan_0","",maxScan,0,maxScan);
    TH1F * h_scan_1 = new TH1F("h_scan_1","",maxScan,0,maxScan);
    TH1F * h_scan_2 = new TH1F("h_scan_2","",maxScan,0,maxScan);
    TH1F * h_scan_3 = new TH1F("h_scan_3","",maxScan,0,maxScan);
    TString xtitle = "TMath::Exp(-0.01*(x))";
    
    for (int nScan = 0; nScan < maxScan; nScan++) {
	TH1F* hs_mc_0 = new TH1F("hs_mc_0","",40,-500,500);
	TH1F* hs_mc_1 = new TH1F("hs_mc_1","",40,0,1000);
        TH1F* hs_mc_2 = new TH1F("hs_mc_2","",40,-500,100);
	TH1F* hs_mc_3 = new TH1F("hs_mc_3","",40,0,700);

	float weight = TMath::Exp(-float(1./maxScan)*(nScan));
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
    h_scan_0->GetXaxis()->SetTitle(xtitle);
    h_scan_0->GetYaxis()->SetTitle("#Sigma (data-mc)^{2}");
    h_scan_0->GetYaxis()->SetTitleOffset(1.4);
    h_scan_0->DrawNormalized();
    h_scan_2->SetLineColor(kRed);
    h_scan_2->DrawNormalized("same");



    TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.70 ) ;
    leg->AddEntry( h_scan_0, "Full range" ) ; 
    leg->AddEntry( h_scan_2, "Restricted range");
    leg->SetFillColor( kWhite ) ; 
    leg->Draw(); 

    
}

void fitXsecDataOnly( double Emin, double Emax) {
    TChain *c_data = new TChain("tracks");
    //c_data->Add("tracks_golden.root");
    c_data->Add("~/work/captain/software/work-area/captSummary/cmt/nb_mtpc_spl_012*_000_flat_beam*.root");

    std::vector<double> *track_length = 0;
    std::vector<double> *first_hit_X = 0;
    std::vector<double> *first_hit_Y = 0;
    std::vector<double> *first_hit_Z = 0;
    std::vector<double> *corrected_first_hit_Z = 0;
    std::vector<double> *corrected_PDS_energy = 0;
    
    c_data->SetBranchAddress("track_length",&track_length);
    c_data->SetBranchAddress("first_hit_X",&first_hit_X);
    c_data->SetBranchAddress("first_hit_Y",&first_hit_Y);
    c_data->SetBranchAddress("first_hit_Z",&first_hit_Z);
    c_data->SetBranchAddress("corrected_first_hit_Z",&corrected_first_hit_Z);
    c_data->SetBranchAddress("corrected_PDS_energy",&corrected_PDS_energy);
    
    TH1F* hs_data = new TH1F("hs_data","",60,-500,100);

    Int_t nEntries = c_data->GetEntries();

    for (int iev=0; iev<nEntries; iev++) {
	c_data->GetEntry(iev);
	for (int t=0; t<first_hit_X->size(); t++) {
	    double x = first_hit_X->at(t);
	    double y = first_hit_Y->at(t);
	    double z = corrected_first_hit_Z->at(t);
	    double e = corrected_PDS_energy->at(t);
	    if (inBeamXY(x,y) && z > -195 && z < -145 && x < 0.0 && x > -400 && e > Emin && e < Emax) {
		hs_data->Fill(x);
	    }
	}
    }

    hs_data->Draw("ep");
    TCanvas *c = new TCanvas("c");
    TF1 *func = new TF1("func","expo",-400,0);
    hs_data->Fit("func","R");
    hs_data->Draw("ep");
    TString fname;
    fname.Form("fits/fit_%.1f_%.1f.pdf",Emin,Emax); 
    c->Print(fname);
    
    hs_data->Delete();
    //c->Delete();
}

void fitXsec() {

    fitXsecDataOnly(100.0,200.0);
    fitXsecDataOnly(200.0,300.0);
    fitXsecDataOnly(300.0,400.0);
    fitXsecDataOnly(400.0,500.0);
    fitXsecDataOnly(500.0,600.0);
    fitXsecDataOnly(600.0,700.0);
    fitXsecDataOnly(700.0,800.0);
    
}
