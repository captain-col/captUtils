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
    
    TH1F* hs_mc_0 = new TH1F("hs_mc_0","",40,-500,500);
    TH1F* hs_data_0 = new TH1F("hs_data_0","",40,-500,500);

    TH1F* hs_mc_1 = new TH1F("hs_mc_1","",40,0,1000);
    TH1F* hs_data_1 = new TH1F("hs_data_1","",40,0,1000);

    TH1F* hs_mc_2 = new TH1F("hs_mc_2","",40,-500,100);
    TH1F* hs_data_2 = new TH1F("hs_data_2","",40,-500,100);

    TH1F* hs_mc_3 = new TH1F("hs_mc_3","",40,0,700);
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
    
    for (int iev = 0; iev < nEntries_mc; iev++) {
	c_mc->GetEntry(iev);
	for (size_t j = 0; j < track_length_mc->size(); j++) {
	    if (first_hit_Z_mc->at(j) > 0 || first_hit_Z_mc->at(j) < -350 || first_hit_Y_mc->at(j) > 0 || first_hit_Y_mc->at(j) < -200) continue;
	    if (first_hit_X_mc->at(j) > 400 || first_hit_X_mc->at(j) < -400) continue;
	    hs_mc_0->Fill(first_hit_X_mc->at(j));
	    hs_mc_1->Fill(track_length_mc->at(j));	
	    if (first_hit_X_mc->at(j) > 0) continue;
	    hs_mc_2->Fill(first_hit_X_mc->at(j));
	    hs_mc_3->Fill(track_length_mc->at(j));
	}
    }

    drawPlots( hs_mc_0, hs_data_0, "X position (mm)",  "firsthitX", Selector::Left); 
    drawPlots( hs_mc_1, hs_data_1, "Track length (mm)",  "tracklength", Selector::Right); 
    drawPlots( hs_mc_2, hs_data_2, "X position (mm)",  "firsthitX_reduced", Selector::Left); 
    drawPlots( hs_mc_3, hs_data_3, "Track length (mm)",  "tracklength_reduced", Selector::Right); 
    
    
}

void drawPlots(TH1* hs_mc, TH1* hs_data, TString xtitle, TString fname, Selector::tag_t channel){
    
    if (channel == Selector::Left) {
	TLegend *leg = new TLegend( 0.18, 0.60, 0.42, 0.80 ) ;
	std::cout<<"LEFT"<<std::endl;
    }
    else {
	TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.70 ) ;
	std::cout<<"NOT LEFT"<<std::endl;
    }
    
    leg->AddEntry( hs_mc, "Simulation" ) ; 
    leg->AddEntry( hs_data, "Data");
    leg->SetFillColor( kWhite ) ; 

    hs_mc->Scale(float(hs_data->Integral())/float(hs_mc->Integral()));

    TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
    THStack* th = new THStack();
    hs_mc->SetFillColor(kYellow);
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
	    cout<<i<<" "<<hs_data->GetBinContent(i)<<" "<<errors->GetBinContent(i)<<" "<<errors->GetBinError(i)<<" "<<errorSyst<<" "<<errorSum<<endl;
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


    //cout<<"Bins1="<<((TH1*)th->GetStack()->Last())->FindBin(80)<<std::endl;
    //cout<<"Bins2="<<((TH1*)th->GetStack()->Last())->FindBin(100)<<std::endl;

  
    ratio->GetXaxis()->SetTitle(xtitle.Data());
    //ths[icanvas]->GetXaxis()->SetTickSize(1.0);
    //ths[icanvas]->GetXaxis()->SetTitleSize(0.1);
    ratio->GetXaxis()->SetTickSize(0.40);
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
