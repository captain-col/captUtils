void truthEff() {

    TChain *c = new TChain("tracks");
    //c->Add("tracks_proton_beam_extended_remove_large.root");
    c->Add("tracks_010100.root");

    std::vector<double> *track_length = 0;
    std::vector<double> *first_hit_X = 0;
    std::vector<double> *last_hit_X = 0;
    std::vector<double> *truth_trajectory_length = 0;
    std::vector<double> *truth_trajectory_first_X = 0;
    std::vector<double> *truth_trajectory_last_X = 0;
    std::vector<double> *truth_vertex_X = 0;
    std::vector<double> *truth_particle_E = 0;

    c->SetBranchAddress("track_length",&track_length);
    c->SetBranchAddress("first_hit_X",&first_hit_X);
    c->SetBranchAddress("last_hit_X",&last_hit_X);
    c->SetBranchAddress("truth_trajectory_length",&truth_trajectory_length);
    c->SetBranchAddress("truth_trajectory_first_X",&truth_trajectory_first_X);
    c->SetBranchAddress("truth_trajectory_last_X",&truth_trajectory_last_X);
    c->SetBranchAddress("truth_vertex_X",&truth_vertex_X);
    c->SetBranchAddress("truth_particle_E",&truth_particle_E);

    Int_t nEntries = c->GetEntries();

    TH2F *h = new TH2F("h","",102,-22,1100,100,0,1100);
    TH2F *h1 = new TH2F("h1","",100,0,1100,100,0,1100);
    TH2F *h2 = new TH2F("h2","",10,0,10,100,0,1100);
    TH2F *h3 = new TH2F("h3","",100,-500,500,100,-500,500);
    TH2F *h4 = new TH2F("h4","",100,-500,500,100,-500,500);
    TH2F *h5 = new TH2F("h5","",100,-500,500,100,0,100);

    TH2F *h6 = new TH2F("h6","",100,0,1000,100,0,1000);
    TH2F *h7 = new TH2F("h7","",100,-500,500,100,0,1000);
    TH2F *h8 = new TH2F("h8","",120,-20,100,100,0,100);

    TH1F *g1 = new TH1F("g1","",100,-1,1);
    TH1F *g2 = new TH1F("g2","",100,-50,50);
    TH1F *g3 = new TH1F("g3","",100,-50,50);
    TH1F *g4 = new TH1F("g4","",100,-50,50);
    
    for (int iev=0; iev< nEntries; iev++) {
	c->GetEntry(iev);

	double x = truth_vertex_X->at(0);
	double t = truth_trajectory_length->at(0);
	double fx = truth_trajectory_first_X->at(0);
	double lx = truth_trajectory_last_X->at(0);

	double e = truth_particle_E->at(0);
	
	if (x < -500 || fx < -500 || lx < -500) continue;
	
	if (t > 1100) t = 1099.;

	h6->Fill(t,e);
	h7->Fill(x,e);
	
	//if (track_length->size() != 1) continue;
	
	if (track_length->size() == 0) {
	    h->Fill(-19,t);
	    h8->Fill(-19,t);
	    h5->Fill(x,t);
	}
	
	for (int r = 0; r < track_length->size(); r++) {
	    double tr = track_length->at(r);
	    double tx = first_hit_X->at(r);
	    double lx = last_hit_X->at(r);
	    //if (tx > 0) continue;
	    //cout<<tr<<endl;
	    if (lx < -500 ) continue;
	    h->Fill(tr,t);
	    h8->Fill(tr,t);
	    h3->Fill(tx,x);
	    g1->Fill((tx-x)/t);
	    g2->Fill(tr-t);
	    double tmin = 10.0;
	    double tmax = 2000.0;
	    if (tr < tmin || tr > tmax) continue;
	    h1->Fill(tr,t);
	    h4->Fill(tx,x);
	    g3->Fill((tx-x)/(1.0));
	    g4->Fill((tr-t));

	}

	h2->Fill(track_length->size(),t);
	
    }

    gStyle->SetOptStat(0);
    TCanvas *c0 = new TCanvas("c0");
    h->Draw("colz");
    TCanvas *c1 = new TCanvas("c1");
    h1->Draw("colz");
    TCanvas *c2 = new TCanvas("c2");
    h2->Draw("colz");
    TCanvas *c3 = new TCanvas("c3");
    h3->Draw("colz");
    TCanvas *c4 = new TCanvas("c4");
    h4->Draw("colz");
    TCanvas *c5 = new TCanvas("c5");
    g1->Draw("colz");
    TCanvas *c6 = new TCanvas("c6");
    g2->Draw("colz");
    TCanvas *c7 = new TCanvas("c7");
    g3->Draw("colz");
    TCanvas *c8 = new TCanvas("c8");
    g4->Draw("colz");
    TCanvas *c9 = new TCanvas("c9");
    h5->Draw("colz");
    TCanvas *c10 = new TCanvas("c10");
    h6->Draw("colz");
    TCanvas *c11 = new TCanvas("c11");
    h7->Draw("colz");
    TCanvas *c12 = new TCanvas("c12");
    h8->Draw("colz");
    
}
