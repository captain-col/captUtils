bool inBeamXY(double x, double y){
    double a = -0.1188;
    double b = -54.9894;
    double width = 54/2;
    double dist=fabs(a*x-y+b)/sqrt(a*a+1); 
    if(dist<width ) return true;
    return false;    
}


void efficiency() {

    TChain *tree = new TChain("tracks");
    //tree->Add("~/work/captain/software/work-area/captSummary/cmt/nb_mtpc_spl_012*_000_flat_*.root");
    tree->Add("~/work/captain/software/work-area/captSummary/cmt/nb_mtpc_spl_012*_000_flat_beam*.root");
    //tree->Add("~/work/captain/software/work-area/captSummary/cmt/test139.root");
    //tree->Add("~/work/captain/software/work-area/captSummary/cmt/nb_mtpc_spl_nb_mtpc_spl_012100_000_reco_80761f0a6a.root_recon_PDS_flat.root");
    //tree->Add("rec_PDS_shift_newcorrection_rightenergy.root");
    std::vector<double> * minZ = 0;
    std::vector<double> * maxZ =0;
    std::vector<double> * minX = 0;
    std::vector<double> * maxX =0;
    std::vector<double> * minY = 0;
    std::vector<double> * maxY =0;
    std::vector<double> * Energy =0;
    std::vector<double> * deltaT = 0;
    std::vector<double> * trackL = 0;
    std::vector<double> * qsum = 0;
    std::vector<double> * qmax = 0;
    std::vector<int> * trigg = 0;
    std::vector<double> * minZ_corr = 0;
    std::vector<double> * E_corr = 0;
    int event = 0;
    int run = 0;
    long long int tpc_time =0;

    tree->SetBranchAddress("first_hit_Z",&minZ);
    tree->SetBranchAddress("last_hit_Z",&maxZ);
    tree->SetBranchAddress("first_hit_X",&minX);
    tree->SetBranchAddress("last_hit_X",&maxX);
    tree->SetBranchAddress("first_hit_Y",&minY);
    tree->SetBranchAddress("last_hit_Y",&maxY);
    tree->SetBranchAddress("track_length",&trackL);
    tree->SetBranchAddress("PDS_energy",&Energy);
    tree->SetBranchAddress("PDS_delta_time",&deltaT);
    tree->SetBranchAddress("PDS_trigger_type",&trigg);
    //tree->SetBranchAddress("PDS_qSum",&qsum);
    //tree->SetBranchAddress("PDS_qMax",&qmax);
    tree->SetBranchAddress("event",&event);
    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("TPC_time",&tpc_time);
    tree->SetBranchAddress("corrected_first_hit_Z",&minZ_corr);
    tree->SetBranchAddress("corrected_PDS_energy",&E_corr);


    TH1F *h_event = new TH1F("h_event","",1e5,0,1e5);
    TH1F *h_eventC = new TH1F("h_eventC","",1e5,0,1e5);


    Float_t bins[] = { 0,20,40,60,80,
		       100,120,140,160,180,
		       200,220,240,260,
		       300,340,360,
		       400,
		       500,560,
		       640,		       
		       800,
		       1000};
    Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
    Float_t bins2[] = {0,20,40,60,80,
		       100,120,140,160,180,
		       200,220,240,260,280,
		       300,320,340,360,380,
		       400,420,440,460,480,
		       500,520,540,560,580,
		       600,700,800,900,1000};
    Int_t  binnum2 = sizeof(bins2)/sizeof(Float_t) - 1;

    Float_t bins3[] = {	
	-480,-460,-440,-420,	
	-400,-380,-360,-340,-320,	
	-300,-280,-260,-240,-220,	
	-200,-180,-160,-140,-120,	
	-100,-80,-60,-40,-20,	
	0,20,40,60,80,
	100,120,140,160,180,
	200,220,240,260,280,
	300,320,340,360,380,
	400,420,440,460,480};
    Int_t  binnum3 = sizeof(bins3)/sizeof(Float_t) - 1;
    
    TH2F *h_tL_tE = new TH2F("h_tL_tE","",binnum,bins,binnum2,bins2);
    //TH2F *h_tL_tE = new TH2F("h_tL_tE","",100,0,2000,100,0,1000);

    TH2F *h_tX_tE = new TH2F("h_tX_tE","",binnum,bins,binnum3,bins3);

    int nentries = tree->GetEntries();

    int ntracks_beamC=0;
    int ntracks_beam=0;
    
    for(int i=0;i<nentries;++i){
	tree->GetEntry(i);
	//if (run != 12139) continue;
	//std::cout<<"Run="<<run<<" Event="<<event<<std::endl;
	
	int ntracks1 = 0;
	int ntracks2 = 0;
	int ntracks3 = 0;
	int ntracksC = 0;

	//if (minZ->size() != 1) continue;

	
	
	for(std::size_t j=0; j<minZ->size();++j){
	    if(inBeamXY((*minX)[j],(*minY)[j])){
		if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145)//  && (*E_corr)[j]>=0)
		    {
			ntracksC++;
		    }
		if(((*minZ)[j]>-195 && (*minZ)[j]<-145))
		    {
			ntracks1++;
		    }
		if(((*minZ)[j]>-510 && (*minZ)[j]<-460))
		    {
			ntracks2++;
		    }
		if(((*minZ)[j]>-830 && (*minZ)[j]<-780))
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
	for(std::size_t j=0; j<minZ->size();++j){

	    int rev = 0;
	    int revC = 0;
	    
	    
	    if(inBeamXY((*minX)[j],(*minY)[j])){
		//if( ((*minZ)[j]>-510 && (*minZ)[j]<-460) || ((*minZ)[j]>-830 && (*minZ)[j]<-780))
		if(((*minZ)[j]>-195 && (*minZ)[j]<-145) || ((*minZ)[j]>-510 && (*minZ)[j]<-460) || ((*minZ)[j]>-830 && (*minZ)[j]<-780))

		    {
			ntracks_beam++;
			//std::cout<<(*minZ)[j]<<std::endl;
			h_event->Fill((run-12090)*2000+(event));
			rev = (run-12090)*2000+(event);
		    }		
		if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145 && (*E_corr)[j]>0 )
		    {
			//if (((*minZ)[j]>-510 && (*minZ)[j]<-460) || ((*minZ)[j]>-830 && (*minZ)[j]<-780))
			    if(((*minZ)[j]>-195 && (*minZ)[j]<-145) || ((*minZ)[j]>-510 && (*minZ)[j]<-460) || ((*minZ)[j]>-830 && (*minZ)[j]<-780))					    
			    {
				ntracks_beamC++;
				//std::cout<<(*minZ_corr)[j]<<std::endl;
				h_eventC->Fill((run-12090)*2000+(event));
				revC = (run-12090)*2000+(event);
				if (abs((*minX)[j]) < 400 && (*minX)[j] < 0) {
				    h_tL_tE->Fill((*E_corr)[j],(*trackL)[j]);
				    h_tX_tE->Fill((*E_corr)[j],(*minX)[j]);
				}
			    }
		    }
	    }
	    
	    if (rev != revC) std::cout<<"run="<<run<<" ev="<<event<<" rev="<<rev<<" "<<revC<<std::endl;

	}


    }

    std::cout<<"Ntracks="<<ntracks_beam<<" NtracksC="<<ntracks_beamC<<std::endl;

    // h_event->Draw();
    // h_eventC->SetLineColor(kRed);
    // h_eventC->Draw("same");
    TCanvas *c1 = new TCanvas("c1");
    h_tX_tE->Draw("colz");
    TCanvas *c2 = new TCanvas("c2");
    h_tL_tE->Draw("colz");
}
