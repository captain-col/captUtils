#include <algorithm>    // std::max

const Double_t m = 939.6;
const Double_t L = 23.2;
const Double_t Tgamma = 629.8e-9;
const Double_t c = 3e8;

TString modname = "PDSmod_edge";

bool inBeamXY(double x, double y){
    double a = -0.1188;
    double b = -54.9894;
    double width = 54/2;
    double dist=fabs(a*x-y+b)/sqrt(a*a+1); 
    if(dist<width ) return true;
    return false;    
}

Double_t calcKE(Double_t x) {
    Float_t xx =x;
    Double_t dtns = xx*1e-9;
    Double_t t = (dtns- Tgamma) + L/c;
    Double_t p = m*L/sqrt(t*t*c*c - L*L);
    Double_t E = sqrt(p*p + m*m) - m;
    return E;
}

Double_t sigmaKE(Double_t x) {
    Float_t xx =x;
    Double_t dtns = xx*1e-9;
    Double_t t = (dtns- Tgamma) + L/c;
    Double_t p = m*L/sqrt(t*t*c*c - L*L);

    Double_t sE = (p*p/(p*p+m*m))*(m*m*L*L*c*c*c*c*t*t/((t*t*c*c-L*L)*(t*t*c*c-L*L)*(t*t*c*c-L*L)))*1e-18*((2*8*8)+(4*4));
    return sqrt(sE);
}


Double_t reverseKE(Double_t x) {
    Float_t E =x;
    // Double_t dtns = xx*1e-9;
    // Double_t t = (dtns- Tgamma) + L/c;
    // Double_t p = m*L/sqrt(t*t*c*c - L*L);
    // Double_t E = sqrt(p*p + m*m) - m;

    // p = sqrt((E + m)*(E + m) - m*m) = m*L/(sqrt(t*t*c*c - L*L));
    // p2 = ((E + m)*(E + m) - m*m) = (m*m*L*L)/((t*t*c*c - L*L));
    // (m*m*L*L)/((E + m)*(E + m) - m*m) = (t*t*c*c - L*L);
    // t*t*c*c = (m*m*L*L)/((E + m)*(E + m) - m*m) + L*L;

    Double_t t = 1/c * sqrt((m*m*L*L)/((E + m)*(E + m) - m*m) + L*L);
    Double_t dtns = t - L/c + Tgamma;
    Double_t dt = dtns*1e9;
    return dt;
    
}

std::pair<double,double> segmentIntercept(double xmin1,double xmax1, double xmin2, double xmax2){

    double xmin = xmin1;
    double xmax = xmax1;
    if (xmin1 < xmin2)
	xmin = xmin2;
    if (xmax1 > xmax2)
	xmax = xmax2;
    return std::make_pair(xmin,xmax);
}

double calcOverlap(std::pair<double,double> overlap) {
    if (overlap.first > overlap.second)
	return 0.0;
    else
	return fabs(overlap.second - overlap.first);
	    
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

std::pair<double,double> fitXsecDataOnly( double Emin, double Emax) {
    TChain *c_data = new TChain("tracks");
    //c_data->Add("tracks_golden.root");
    c_data->Add("~/work/captain/software/work-area/captSummary/cmt/nb_mtpc_spl_012*_tree_newPDS4_*.root");

    std::vector<double> *track_length = 0;
    std::vector<double> *first_hit_X = 0;
    std::vector<double> *first_hit_Y = 0;
    std::vector<double> *first_hit_Z = 0;
    std::vector<double> *corrected_first_hit_Z = 0;
    std::vector<double> *corrected_PDS_energy = 0;

    std::vector<double> *last_hit_X = 0;
    std::vector<double> *last_hit_Y = 0;
    std::vector<double> *corrected_last_hit_Z = 0;

    c_data->SetBranchAddress("track_length",&track_length);
    c_data->SetBranchAddress("first_hit_X",&first_hit_X);
    c_data->SetBranchAddress("first_hit_Y",&first_hit_Y);
    c_data->SetBranchAddress("first_hit_Z",&first_hit_Z);
    c_data->SetBranchAddress("corrected_first_hit_Z",&corrected_first_hit_Z);
    c_data->SetBranchAddress("corrected_PDS_energy",&corrected_PDS_energy);
    c_data->SetBranchAddress("last_hit_X",&last_hit_X);
    c_data->SetBranchAddress("last_hit_Y",&last_hit_Y);
    c_data->SetBranchAddress("corrected_last_hit_Z",&corrected_last_hit_Z);
    
    TH1F* hs_data = new TH1F("hs_data","",100,-500,100);

    Int_t nEntries = c_data->GetEntries();
    std::vector<double> selected_X;
    std::vector<double> selected_Xw;

    for (int iev=0; iev<nEntries; iev++) {
	c_data->GetEntry(iev);

	int ntracks1 = 0;
	int ntracks2 = 0;
	int ntracks3 = 0;
	int ntracksC = 0;

	for(std::size_t j=0; j<first_hit_Z->size();++j){
	    if(inBeamXY((*first_hit_X)[j],(*first_hit_Y)[j])){
		if((*corrected_first_hit_Z)[j]>-195 && (*corrected_first_hit_Z)[j]<-145)//  && (*E_corr)[j]>=0)
		    {
			ntracksC++;
		    }
		if(((*first_hit_Z)[j]>-195 && (*first_hit_Z)[j]<-145))
		    {
			ntracks1++;
		    }
		if(((*first_hit_Z)[j]>-510 && (*first_hit_Z)[j]<-460))
		    {
			ntracks2++;
		    }
		if(((*first_hit_Z)[j]>-830 && (*first_hit_Z)[j]<-780))
		    {
			ntracks3++;
		    }
	    }
	}

	//std::cout<<"Tracks event="<<ntracks<<" "<<ntracksC<<std::endl;

	if (ntracks1 > 1 || ntracks2 > 1 || ntracks3 > 1)
	    {
		continue;
	    }
	
	for (int t=0; t<first_hit_X->size(); t++) {
	    double x = first_hit_X->at(t);
	    double y = first_hit_Y->at(t);
	    double z = corrected_first_hit_Z->at(t);
	    double xf = last_hit_X->at(t);
	    double yf = last_hit_Y->at(t);
	    double zf = corrected_last_hit_Z->at(t);
	    double e0 = corrected_PDS_energy->at(t);
	    double e1 = reverseKE(e0);
	    double e = calcKE(e1);
	    double tl = track_length->at(t);
	    if (inBeamXY(x,y) && z > -195 && z < -145 && x < 4.6 && x > -400.6 && tl > 15) {
		//&& e > Emin && e < Emax){// && abs(xf) < 450 && abs(yf) < 450 && zf < 0.0 && zf > -300) {				
		//if(e>0)
		if(e > Emin && e < Emax)
		    {

		    double ews[] = {0,0,0,0,0,0,0,0};
		   
		    //cout<<"E="<<e<<" Eplus="<<e+sigmaKE(e1)<<" Eminus="<<e-sigmaKE(e1)<<std::endl;
		    for (int hs = 0; hs < 700; hs += 100) {
			// std::cout<<"hs="<<hs<<","<<hs+100<<std::endl;
			// std::cout<<"E="<<e<<" Eplus="<<e+sigmaKE(e1)<<" Eminus="<<e-sigmaKE(e1)<<std::endl;
			// std::cout<<"Overlap=("<<segmentIntercept(hs,hs+100,e-sigmaKE(e1),e+sigmaKE(e1)).first<<","<<segmentIntercept(hs,hs+100,e-sigmaKE(e1),e+sigmaKE(e1)).second<<")"<<std::endl;
			ews[hs/100] = calcOverlap(segmentIntercept(hs,hs+100,e-sigmaKE(e1),e+sigmaKE(e1)))/fabs(2*sigmaKE(e1));
		    }

		    // std::cout<<"hs="<<700<<","<<900<<std::endl;
		    // std::cout<<"E="<<e<<" Eplus="<<e+sigmaKE(e1)<<" Eminus="<<e-sigmaKE(e1)<<std::endl;
		    // std::cout<<"Overlap=("<<segmentIntercept(700,900,e-sigmaKE(e1),e+sigmaKE(e1)).first<<","<<segmentIntercept(700,900,e-sigmaKE(e1),e+sigmaKE(e1)).second<<")"<<std::endl;
		    ews[7] = calcOverlap(segmentIntercept(700,900,e-sigmaKE(e1),e+sigmaKE(e1)))/fabs(2*sigmaKE(e1));

		    // for(auto ew:ews)
		    // 	cout<<"Ew="<<ew<<endl;

		    // cout<<"Emin="<<Emin<<" ew[Emin]="<<ews[int(Emin/100)]<<endl;
			
		    //hs_data->Fill(x,ews[int(Emin/100)]);
		    hs_data->Fill(x);
		    selected_X.push_back(x);
		    selected_Xw.push_back(1.0);//ews[int(Emin/100)]);

		}								

											
	    }
	}
    }

    std::cout<<"Data size="<<hs_data->GetEntries()<<" "<<selected_X.size()<<std::endl;
    
    const Int_t nq = 4;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    Double_t yq[nq];  // array to contain the quantiles
    for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
    hs_data->GetQuantiles(nq,yq,xq);
    //for (Int_t i=0;i<nq;i++) std::cout<<"i="<<xq[i]<<" Q ="<<yq[i]<<std::endl;
    // std::cout<<"IQR="<<2*abs(yq[0]-yq[2])/pow(hs_data->GetEntries(),1/3.)<<std::endl;

    double binWidth = 2*abs(yq[0]-yq[2])/pow(hs_data->GetEntries(),1/3.);
    int nBins = int(600.0/binWidth + 0.5);
    // std::cout<<"Nbins="<<nBins<<std::endl;

    TH1F* hs_xpos = new TH1F("hs_xpos","",nBins,-500,100);
    for (int ix=0; ix<selected_X.size(); ix++) hs_xpos->Fill(selected_X[ix]);//,selected_Xw[ix]);
    hs_xpos->Draw("ep");

    gStyle->SetOptFit(1);
    TCanvas *c = new TCanvas("c");
    TF1 *func = new TF1("func","expo",-391,-9);
    hs_xpos->Fit("func","RM");
    hs_xpos->Draw("ep");
    TString fname;
    fname.Form("fits/fit_%d_%d_"+modname+".pdf",int(Emin),int(Emax)); 
    c->Print(fname);
    fname.Form("fits/fit_%d_%d_"+modname+".png",int(Emin),int(Emax)); 
    c->Print(fname);
    fname.Form("fits/fit_%d_%d_"+modname+".C",int(Emin),int(Emax)); 
    c->Print(fname);
    
    hs_data->Delete();
    hs_xpos->Delete();
    
    c->Close();
    return std::make_pair(func->GetParameter(1),func->GetParError(1));
}

void fitXsec() {

    double slope = 0.0;
    double slopeE = 0.0;
    
    //slope = fitXsecDataOnly(0.0,200000.0).first;
    //slopeE = fitXsecDataOnly(0.0,2000.0).second;
    //std::cout<<"Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<1./slope<<"+/-"<< abs(1./slope)*sqrt((float(slopeE)/float(slope))*(float(slopeE)/float(slope)))<<std::endl;
    
    const Int_t nq = 6;
    Double_t xq[nq];
    Double_t yq[nq];
    Double_t xqE[nq];
    Double_t yqE[nq];
    double N_Ar = 1.3973  * 6.022E23  / 39.948; // (g/cm3 * n/mol) / g/mol = n/cm3 

    // for (int i=0;i<nq-1;i++) {
    // 	xq[i] = 150+i*100;
    // 	xqE[i] = sigmaKE(reverseKE(xq[i]));
    // 	slope = 10 * fitXsecDataOnly(100.0+100*i,200.0+100*i).first; // 1/cm
    // 	slopeE = 10 * fitXsecDataOnly(100.0+100*i,200.0+100*i).second;
    // 	yq[i] =  1e24*(slope)/(N_Ar);
    // 	yqE[i] = 1e24*slopeE/(N_Ar);
    // 	std::cout<<"Emin="<<100.0+100*i<<" Emax="<< 200.0+100*i<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq[i]<<"+/-"<< yqE[i]<<std::endl;
    // }



    // xq[0] = 150;
    // xqE[0] = 50;//std::max(50.0,sigmaKE(reverseKE(150)));
    // slope = 10 * fitXsecDataOnly(100.0,200.0).first; // 1/cm
    // slopeE = 10 * fitXsecDataOnly(100.0,200.0).second;
    // yq[0] =  1e24*(slope)/(N_Ar);
    // yqE[0] = 1e24*slopeE/(N_Ar);
    
    // xq[1] = 250;
    // xqE[1] = 50;//std::max(50.0,sigmaKE(reverseKE(250)));
    // slope = 10 * fitXsecDataOnly(200.0,300.0).first; // 1/cm
    // slopeE = 10 * fitXsecDataOnly(200.0,300.0).second;
    // yq[1] =  1e24*(slope)/(N_Ar);
    // yqE[1] = 1e24*slopeE/(N_Ar);
    
    // xq[2] = 350;
    // xqE[2] = 50;//std::max(50.0,sigmaKE(reverseKE(350)));
    // slope = 10 * fitXsecDataOnly(300.0,400.0).first; // 1/cm
    // slopeE = 10 * fitXsecDataOnly(300.0,400.0).second;
    // yq[2] =  1e24*(slope)/(N_Ar);
    // yqE[2] = 1e24*slopeE/(N_Ar);
    
    
    // xq[3] = 500;
    // xqE[3] = std::max(50.0,sigmaKE(reverseKE(500)));
    // slope = 10 * fitXsecDataOnly(500.0-xqE[3],500.0+xqE[3]).first; // 1/cm
    // slopeE = 10 * fitXsecDataOnly(500.0-xqE[3],500.0+xqE[3]).second;
    // yq[3] =  1e24*(slope)/(N_Ar);
    // yqE[3] = 1e24*slopeE/(N_Ar);
    
    // xq[4] = 800;
    // xqE[4] = std::max(50.0,sigmaKE(reverseKE(800)));
    // slope = 10 * fitXsecDataOnly(800.0-xqE[3],800.0+xqE[3]).first; // 1/cm
    // slopeE = 10 * fitXsecDataOnly(800.0-xqE[3],800.0+xqE[3]).second;
    // yq[4] =  1e24*(slope)/(N_Ar);
    // yqE[4] = 1e24*slopeE/(N_Ar);
    

    double xminE = 71.6802;
    xminE = 100.0;
    double xmaxE = 199;    
    xq[0] = xmaxE - fabs(xminE-xmaxE)/2;
    xqE[0] = fabs(xminE-xmaxE)/2;
    slope = 10 * fitXsecDataOnly(xminE,xmaxE).first; // 1/cm
    slopeE = 10 * fitXsecDataOnly(xminE,xmaxE).second;
    yq[0] =  1e24*(slope)/(N_Ar);
    yqE[0] = 1e24*slopeE/(N_Ar);

    std::cout<<"Emin="<<xminE<<" Emax="<< xmaxE<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq[0]<<"+/-"<< yqE[0]<<std::endl;
    
    xminE = xmaxE;
    xmaxE = 296;
    
    xq[1] = xmaxE - fabs(xminE-xmaxE)/2;
    xqE[1] = fabs(xminE-xmaxE)/2;
    slope = 10 * fitXsecDataOnly(xminE,xmaxE).first; // 1/cm
    slopeE = 10 * fitXsecDataOnly(xminE,xmaxE).second;
    yq[1] =  1e24*(slope)/(N_Ar);
    yqE[1] = 1e24*slopeE/(N_Ar);
    std::cout<<"Emin="<<xminE<<" Emax="<< xmaxE<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq[1]<<"+/-"<< yqE[1]<<std::endl;

    xminE = xmaxE;
    xmaxE = 369;
    
    xq[2] = xmaxE - fabs(xminE-xmaxE)/2;
    xqE[2] = fabs(xminE-xmaxE)/2;
    slope = 10 * fitXsecDataOnly(xminE,xmaxE).first; // 1/cm
    slopeE = 10 * fitXsecDataOnly(xminE,xmaxE).second;
    yq[2] =  1e24*(slope)/(N_Ar);
    yqE[2] = 1e24*slopeE/(N_Ar);
    std::cout<<"Emin="<<xminE<<" Emax="<< xmaxE<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq[2]<<"+/-"<< yqE[2]<<std::endl;

    xminE = xmaxE;
    xmaxE = 481;
    
    xq[3] = xmaxE - fabs(xminE-xmaxE)/2;
    xqE[3] = fabs(xminE-xmaxE)/2;
    slope = 10 * fitXsecDataOnly(xminE,xmaxE).first; // 1/cm
    slopeE = 10 * fitXsecDataOnly(xminE,xmaxE).second;
    yq[3] =  1e24*(slope)/(N_Ar);
    yqE[3] = 1e24*slopeE/(N_Ar);
    std::cout<<"Emin="<<xminE<<" Emax="<< xmaxE<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq[3]<<"+/-"<< yqE[3]<<std::endl;

    xminE = xmaxE;
    xmaxE = 674;
    
    xq[4] = xmaxE - fabs(xminE-xmaxE)/2;
    xqE[4] = fabs(xminE-xmaxE)/2;
    slope = 10 * fitXsecDataOnly(xminE,xmaxE).first; // 1/cm
    slopeE = 10 * fitXsecDataOnly(xminE,xmaxE).second;
    yq[4] =  1e24*(slope)/(N_Ar);
    yqE[4] = 1e24*slopeE/(N_Ar);
    std::cout<<"Emin="<<xminE<<" Emax="<< xmaxE<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq[4]<<"+/-"<< yqE[4]<<std::endl;

    xminE = xmaxE;
    xmaxE = 836.4;
    xmaxE = 900.0;
    
    xq[5] = xmaxE - fabs(xminE-xmaxE)/2;
    xqE[5] = fabs(xminE-xmaxE)/2;
    slope = 10 * fitXsecDataOnly(xminE,xmaxE).first; // 1/cm
    slopeE = 10 * fitXsecDataOnly(xminE,xmaxE).second;
    yq[5] =  1e24*(slope)/(N_Ar);
    yqE[5] = 1e24*slopeE/(N_Ar);
    std::cout<<"Emin="<<xminE<<" Emax="<< xmaxE<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq[5]<<"+/-"<< yqE[5]<<std::endl;



    xminE = 100;
    xmaxE = 836.4;
    xmaxE = 900.0;
    
    xq5 = xmaxE - fabs(xminE-xmaxE)/2;
    xqE5 = fabs(xminE-xmaxE)/2;
    slope = 10 * fitXsecDataOnly(xminE,xmaxE).first; // 1/cm
    slopeE = 10 * fitXsecDataOnly(xminE,xmaxE).second;
    yq5 =  1e24*(slope)/(N_Ar);
    yqE5 = 1e24*slopeE/(N_Ar);
    std::cout<<"Emin="<<xminE<<" Emax="<< xmaxE<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq5<<"+/-"<< yqE5<<std::endl;
    
    xminE = xmaxE;
    xmaxE = 1620.9;
    
    // xq[6] = xmaxE - fabs(xminE-xmaxE)/2;
    // xqE[6] = fabs(xminE-xmaxE)/2;
    // slope = 10 * fitXsecDataOnly(xminE,xmaxE).first; // 1/cm
    // slopeE = 10 * fitXsecDataOnly(xminE,xmaxE).second;
    // yq[6] =  1e24*(slope)/(N_Ar);
    // yqE[6] = 1e24*slopeE/(N_Ar);


    
    // pair <double,double> fitRes = fitXsecDataOnly(700,900);

    // xq[6] = 800;
    // xqE[6] = sigmaKE(reverseKE(800));
    // slope = 10 * fitRes.first; // 1/cm
    // slopeE = 10 * fitRes.second;
    // //double N_Ar = 1.3973  * 6.022E23  / 39.948; // (g/cm3 * n/mol) / g/mol = n/cm3 
    // yq[6] =  1e24*(slope)/(N_Ar);
    // yqE[6] = 1e24*slopeE/(N_Ar);
    // std::cout<<"Emin="<< xq[6]- xqE[6]<<" Emax="<< xq[6] + xqE[6]<<" Slope="<<slope<<"+/-"<<slopeE<<" Xsec="<<yq[6]<<"+/-"<< yqE[6]<<std::endl;

    gStyle->SetPadLeftMargin(0.12); gStyle->SetPadRightMargin(0.05);    
    gStyle->SetPadBottomMargin(0.15); gStyle->SetPadTopMargin(0.05);    
    //TCanvas *c = new TCanvas("c");
    
    TString captainText     = "CAPTAIN";
    float captainTextFont   = 61;  // default is helvetic-bold

    // ratio of "CAPTAIN" and extra text size
    float extraOverCaptainTextSize  = 0.76;

    float captainTextSize      = 1.35;
    float captainTextOffset    = 0.1;  // only used in outOfFrame version

    int W = 800;
    int H = 600;

    // Macro taken from CMS example 
    // Initiated by: Gautier Hamel de Monchenault (Saclay)
    // Updated by:   Dinko Ferencek (Rutgers)
    //
    int H_ref = 600; 
    int W_ref = 800; 

    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref; 
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;
    TString canvName = "Xsec_final";

    TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
    // canv->SetFillColor(0);
    // canv->SetBorderMode(0);
    // canv->SetFrameFillStyle(0);
    // canv->SetFrameBorderMode(0);
    // canv->SetLeftMargin( L/W );
    // canv->SetRightMargin( R/W );
    // canv->SetTopMargin( T/H );
    // canv->SetBottomMargin( B/H );
    // canv->SetTickx(0);
    // canv->SetTicky(0);

    float l = canv->GetLeftMargin();
    float t = canv->GetTopMargin();
    float r = canv->GetRightMargin();
    float b = canv->GetBottomMargin();

    int alignY_=3;
    int alignX_=2;
    int align_ = 10*alignX_ + alignY_;

    float relPosX    = 0.155;
    float relPosY    = 0.035;

    float posX_= l + relPosX*(1-l-r);
    float posY_ = 1-t - relPosY*(1-t-b);

    gPad->SetTickx();//, y()
    gPad->SetTicky();
    TGraphErrors *gr = new TGraphErrors(nq,xq,yq,xqE,yqE);
    gr->SetTitle("");
    //gr->SetMarkerColor(4);
    gr->SetMarkerStyle(20);
    TString xtitle = "Energy[MeV]";
    TString ytitle = "#sigma[barn]";
    gr->GetYaxis()->SetTitle(ytitle.Data());
    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->SetRangeUser(0.1,1.8);
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetNdivisions(512);
    gr->GetXaxis()->SetTitle(xtitle.Data());
    gr->GetXaxis()->SetTitleSize(0.06);
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetNdivisions(505);
    //gr->GetYaxis()->SetTitleOffset(1.5);

    gr->Draw("AP");

    Double_t xqG[nq];
    Double_t yqG[nq];
    Double_t xqGE[nq];
    Double_t yqGE[nq];

    xminE = 71.6802;
    xminE = 100.0;
    xmaxE = 199;    
    xqG[0] = xmaxE - fabs(xminE-xmaxE)/2;
    xqGE[0] = fabs(xminE-xmaxE)/2;
    yqG[0] =  0.70;
    yqGE[0] = 0.030;
    
    xminE = xmaxE;
    xmaxE = 296;
    
    xqG[1] = xmaxE - fabs(xminE-xmaxE)/2;
    xqGE[1] = fabs(xminE-xmaxE)/2;
    yqG[1] =  0.54;
    yqGE[1] = 0.022;

    xminE = xmaxE;
    xmaxE = 369;
    
    xqG[2] = xmaxE - fabs(xminE-xmaxE)/2;
    xqGE[2] = fabs(xminE-xmaxE)/2;
    yqG[2] =  0.54;
    yqGE[2] = 0.015;

    xminE = xmaxE;
    xmaxE = 481;
    
    xqG[3] = xmaxE - fabs(xminE-xmaxE)/2;
    xqGE[3] = fabs(xminE-xmaxE)/2;
    yqG[3] =  0.59;
    yqGE[3] = 0.013;

    xminE = xmaxE;
    xmaxE = 674;
    
    xqG[4] = xmaxE - fabs(xminE-xmaxE)/2;
    xqGE[4] = fabs(xminE-xmaxE)/2;
    yqG[4] =  0.57;
    yqGE[4] = 0.012;

    xminE = xmaxE;
    xmaxE = 836.4;
    xmaxE = 900.0;
    
    xqG[5] = xmaxE - fabs(xminE-xmaxE)/2;
    xqGE[5] = fabs(xminE-xmaxE)/2;
    yqG[5] =  0.56;
    yqGE[5] = 0.011;

    TGraphErrors *grG = new TGraphErrors(nq,xqG,yqG,xqGE,yqGE);
    grG->SetLineColor(kRed);
    grG->Draw("P");
    gr->SetFillColor(kWhite);
    grG->SetFillColor(kWhite);
    
    Double_t xqF[nq];
    Double_t yqF[nq];
    Double_t xqFE[nq];
    Double_t yqFE[nq];

    xminE = 71.6802;
    xminE = 100.0;
    xmaxE = 199;    
    xqF[0] = xmaxE - fabs(xminE-xmaxE)/2;
    xqFE[0] = fabs(xminE-xmaxE)/2;
    yqF[0] =  0.741;
    yqFE[0] = 0.015;
    
    xminE = xmaxE;
    xmaxE = 296;
    
    xqF[1] = xmaxE - fabs(xminE-xmaxE)/2;
    xqFE[1] = fabs(xminE-xmaxE)/2;
    yqF[1] =  0.635;
    yqFE[1] = 0.007;

    xminE = xmaxE;
    xmaxE = 369;
    
    xqF[2] = xmaxE - fabs(xminE-xmaxE)/2;
    xqFE[2] = fabs(xminE-xmaxE)/2;
    yqF[2] =  0.655;
    yqFE[2] = 0.004;

    xminE = xmaxE;
    xmaxE = 481;
    
    xqF[3] = xmaxE - fabs(xminE-xmaxE)/2;
    xqFE[3] = fabs(xminE-xmaxE)/2;
    yqF[3] =  0.658;
    yqFE[3] = 0.004;

    xminE = xmaxE;
    xmaxE = 674;
    
    xqF[4] = xmaxE - fabs(xminE-xmaxE)/2;
    xqFE[4] = fabs(xminE-xmaxE)/2;
    yqF[4] =  0.656;
    yqFE[4] = 0.003;

    xminE = xmaxE;
    xmaxE = 836.4;
    xmaxE = 900.0;
    
    xqF[5] = xmaxE - fabs(xminE-xmaxE)/2;
    xqFE[5] = fabs(xminE-xmaxE)/2;
    yqF[5] =  0.649;
    yqFE[5] = 0.003;

    TGraphErrors *grF = new TGraphErrors(nq,xqF,yqF,xqFE,yqFE);
    grF->SetLineColor(kBlue);
    grF->Draw("P");
    grF->SetFillColor(kWhite);
    
    TLegend *leg = new TLegend( 0.70, 0.75, 0.90, 0.90 ) ;
    leg->AddEntry( gr, "Data" ) ;
    leg->AddEntry( grG, "GEANT4" ) ;
    leg->AddEntry( grF, "FLUKA" ) ;
    //Float_t sizeleg = leg->GetTextSize()*1.2;
    //leg->SetTextSize( sizeleg ) ;
    leg->SetFillColor( kWhite ) ;
    leg->Draw();


  
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);    
    
    float extraTextSize = extraOverCaptainTextSize*captainTextSize;

    latex.SetTextFont(captainTextFont);
    latex.SetTextSize(captainTextSize*t);
    latex.SetTextAlign(align_);
    latex.DrawLatex(posX_, posY_, captainText);
    
    canv->Print("fits/Xsec_"+modname+".png");																		  
    canv->Print("fits/Xsec_"+modname+".pdf");												
    canv->Print("fits/Xsec_"+modname+".C");												
    
}
