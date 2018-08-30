
//gROOT->Reset();
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <array>
#include <utility>


#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLine.h"

///just matching

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



int Analysis8()
{
    gROOT->Reset();
    gROOT->ProcessLine("#include <vector>");
    TCanvas* c0 = new TCanvas("c0","A Simple Graph Example0",200,10,700,500);
    TCanvas* c1 = new TCanvas("c1","A Simple Graph Example1",200,10,700,500);
    TCanvas* c2 = new TCanvas("c2","A Simple Graph Example2",200,10,700,500);
    TCanvas* c3 = new TCanvas("c3","A Simple Graph Example3",200,10,700,500);
    TCanvas* c4 = new TCanvas("c4","A Simple Graph Example4",200,10,700,500);
    TCanvas* c5 = new TCanvas("c5","A Simple Graph Example5",200,10,700,500);
    TCanvas* c6 = new TCanvas("c6","A Simple Graph Example6",200,10,700,500);
    TCanvas* c7 = new TCanvas("c7","A Simple Graph Example7",200,10,700,500);
    TCanvas* c8 = new TCanvas("c8","A Simple Graph Example8",200,10,700,500);
    
    TFile *hfile;
    //TTree *tree;
    std::cout<<"Run"<<std::setw(7)<<" Event"<<std::setw(15)<<" TPC_time"<<std::setw(15)<<" nTrack"<<std::setw(7)<<" nPDSTrig"<<std::setw(7)<<std::endl;
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
    // hfile = new TFile("rec_PDS_shift.root","READ");
    hfile = new TFile("rec_PDS_shift_newcorrection_rightenergy.root","READ");
    //tree = (TTree*)hfile->Get("tracks");

    TChain *tree = new TChain("tracks");
    tree->Add("~/work/captain/software/work-area/captSummary/cmt/nb_mtpc_spl_012*_000_flat_beam*.root");
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
    
    int nentries = tree->GetEntries();
    
    int number=0;
    
    TH1F* FZ = new TH1F("fz_old","fz_old",1000,-3500,3500);
    FZ->GetXaxis()->SetTitle("Z,mm");
    FZ->GetYaxis()->SetTitle("number of events");
    
    TH1F* FZ_new = new TH1F("fz_new","fz_new",1000,-3500,3500);
    FZ->GetXaxis()->SetTitle("Z,mm");
    FZ->GetYaxis()->SetTitle("number of events");
    
    TH1F* FZ_new_1 = new TH1F("fz_new_1ene/track","fz_new_1ene/track",1000,-3500,3500);
    FZ->GetXaxis()->SetTitle("Z,mm");
    FZ->GetYaxis()->SetTitle("number of events");
    
    TH1F* FZ_new_n = new TH1F("fz_new_nene/track","fz_new_nene/track",1000,-3500,3500);
    FZ->GetXaxis()->SetTitle("Z,mm");
    FZ->GetYaxis()->SetTitle("number of events");
    
    TH1F* fEnergy = new TH1F("fenergy_old","fenergy_old",999,1,1000);
    fEnergy->GetXaxis()->SetTitle("E,MeV");
    fEnergy->GetYaxis()->SetTitle("number of events");
    TH1F* fEnergy_new = new TH1F("fenergy_new","fenergy_new",999,1,1000);
    fEnergy_new->GetXaxis()->SetTitle("E,MeV");
    fEnergy_new->GetYaxis()->SetTitle("number of events");
    
    TH2F* fEn_L = new TH2F("fenergy_new","fenergy_new",50,-100,1000,100,1,1000);
    fEn_L->GetXaxis()->SetTitle("E_neutron,MeV");
    fEn_L->GetYaxis()->SetTitle("Tracklength,mm");
    
    TH2F* fEn_L_1 = new TH2F("fenergy_new_1ene/track","fenergy_new_1ene/track",50,-100,1000,100,1,1000);
    fEn_L_1->GetXaxis()->SetTitle("E_neutron,MeV");
    fEn_L_1->GetYaxis()->SetTitle("Tracklength,mm");
    
    TH2F* fEn_L_n = new TH2F("fenergy_newf_nene/track","fenergy_new_nene/track",50,-100,1000,100,1,1000);
    fEn_L_n->GetXaxis()->SetTitle("E_neutron,MeV");
    fEn_L_n->GetYaxis()->SetTitle("Tracklength,mm");
    
    TH1F* pds = new TH1F("pds","pds",100,0,100);
    
   
    int thd=0;
    int scd=0;
    int fst=0;
    int ntracks=0;
    int ntracks_beam=0;

    for(int i=0;i<nentries;++i){
        tree->GetEntry(i);
        std::set<double> dups;
        int energies=0;
        ntracks+=(int)minZ->size();
        for(std::size_t j=0; j<minZ->size();++j){
            if((*E_corr)[j]>0){
                energies++;
            dups.insert((*E_corr)[j]);
            }
        }
        
        bool dup=0;
        if((int)dups.size()<energies)dup=1;
        
        for(std::size_t j=0; j<minZ->size();++j){
            
            if(inBeamXY((*minX)[j],(*minY)[j])){
                if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145)ntracks_beam++;
            if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145 && (*E_corr)[j]>=0)
            fEn_L->Fill((*E_corr)[j],(*trackL)[j]);
            if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145 && (*E_corr)[j]==-9999)
            fEn_L->Fill(-100,(*trackL)[j]);
            if((*E_corr)[j]>=0)
            FZ_new->Fill((*minZ_corr)[j]);
            
            if(!dup){
                if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145 && (*E_corr)[j]>=0)
                    fEn_L_1->Fill((*E_corr)[j],(*trackL)[j]);
                if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145 && (*E_corr)[j]==-9999)
                    fEn_L_1->Fill(-100,(*trackL)[j]);
                if((*E_corr)[j]>=0)
                    FZ_new_1->Fill((*minZ_corr)[j]);
                
            }
            
            if(dup){
                if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145 && (*E_corr)[j]>=0)
                    fEn_L_n->Fill((*E_corr)[j],(*trackL)[j]);
                if((*minZ_corr)[j]>-195 && (*minZ_corr)[j]<-145 && (*E_corr)[j]==-9999)
                    fEn_L_n->Fill(-100,(*trackL)[j]);
                if((*E_corr)[j]>=0)
                    FZ_new_n->Fill((*minZ_corr)[j]);
                
            }
            
            }
            
           /* if((*minZ)[j]>-1000 && (*minZ)[j]<-700)thd++;
             if((*minZ)[j]>-700 && (*minZ)[j]<-350)scd++;
            if((*minZ)[j]>-350 && (*minZ)[j]<-0)fst++;;
            FZ->Fill((*minZ)[j]);*/
        }
        
       
        int npds=deltaT->size();
        pds->Fill(npds);
            
            
        }
    
    std::cout<<"1="<<fst<<" ;2="<<scd<<" ;3="<<thd<<std::endl;
    
    c0->cd();
    FZ_new_1->Draw();
    c1->cd();
    FZ_new_n->Draw();
    c2->cd();
    FZ_new->Draw();
    c3->cd();
    fEn_L_1->Draw("COLZ");
    gPad->Print("plotCheck/En_L_1.C");
    c4->cd();
    fEn_L_n->Draw("COLZ");
    gPad->Print("plotCheck/En_L_n.C");
    c5->cd();
    fEn_L->Draw("COLZ");
    gPad->Print("plotCheck/En_L.C");
    
  //  std::cout<<"fEn_L_1.entries="<<fEn_L_1->GetEntries()<<"fEn_L_n.entries="<<fEn_L_n->GetEntries()<<"fEn_L.entries="<<fEn_L->GetEntries()<<
    
    double x[(int)fEn_L->GetXaxis()->GetNbins()];
    double y[(int)fEn_L->GetXaxis()->GetNbins()];
    double errx[(int)fEn_L->GetXaxis()->GetNbins()];
    double erry[(int)fEn_L->GetXaxis()->GetNbins()];
    for(int i=1;i<fEn_L->GetXaxis()->GetNbins()+1;++i){
        double mean=0;
        double rms=0;
        double nz=0;
        for(int j=1;j<fEn_L->GetYaxis()->GetNbins()+1;++j){
            int ibin = fEn_L->GetBin(i,j);
            if(fEn_L->GetBinContent(ibin)>0){
                nz++;
                mean+=fEn_L->GetBinContent(ibin);
            }
        }
        if(nz>0)mean;
        if(mean>0){
            for(int j=1;j<fEn_L->GetYaxis()->GetNbins()+1;++j){
                int ibin = fEn_L->GetBin(i,j);
                if(fEn_L->GetBinContent(ibin)>0){
                    rms=pow(mean-fEn_L->GetBinContent(ibin),2);
                }
            }
            rms=sqrt(rms/nz);
        }
        x[i]=fEn_L->GetXaxis()->GetBinCenter(i);
        y[i]=mean;
        errx[i]=0.01;
        erry[i]=rms;
        // std::cout<<"mean+-rms = "<<mean<<"+-"<<rms<<std::endl;
        
    }
    c6->cd();
    TGraphErrors* errplot = new TGraphErrors((int)fEn_L->GetXaxis()->GetNbins(),x,y,errx,erry);
    errplot->GetXaxis()->SetTitle("E_neutron,MeV");
    errplot->GetYaxis()->SetTitle("Mean_track_length,mm");
    errplot->Draw("AP");
    gPad->Print("plotCheck/meanLength.C");
    
    double x1[(int)fEn_L_1->GetXaxis()->GetNbins()];
    double y1[(int)fEn_L_1->GetXaxis()->GetNbins()];
    double errx1[(int)fEn_L_1->GetXaxis()->GetNbins()];
    double erry1[(int)fEn_L_1->GetXaxis()->GetNbins()];
    for(int i=1;i<fEn_L_1->GetXaxis()->GetNbins()+1;++i){
        double mean=0;
        double rms=0;
        double nz=0;
        for(int j=1;j<fEn_L_1->GetYaxis()->GetNbins()+1;++j){
            int ibin = fEn_L_1->GetBin(i,j);
            if(fEn_L_1->GetBinContent(ibin)>0){
                nz++;
                mean+=fEn_L_1->GetBinContent(ibin);
            }
        }
        if(nz>0)mean;
        if(mean>0){
            for(int j=1;j<fEn_L_1->GetYaxis()->GetNbins()+1;++j){
                int ibin = fEn_L_1->GetBin(i,j);
                if(fEn_L_1->GetBinContent(ibin)>0){
                    rms=pow(mean-fEn_L_1->GetBinContent(ibin),2);
                }
            }
            rms=sqrt(rms/nz);
        }
        x1[i]=fEn_L_1->GetXaxis()->GetBinCenter(i);
        y1[i]=mean;
        errx1[i]=0.01;
        erry1[i]=rms;
        // std::cout<<"mean+-rms = "<<mean<<"+-"<<rms<<std::endl;
        
    }
    c7->cd();
    TGraphErrors* errplot1 = new TGraphErrors((int)fEn_L_1->GetXaxis()->GetNbins(),x1,y1,errx1,erry1);
    errplot1->GetXaxis()->SetTitle("E_neutron,MeV");
    errplot1->GetYaxis()->SetTitle("Mean_track_length,mm");
    errplot1->Draw("AP");
    gPad->Print("plotCheck/meanLength_1.C");
    
    double xn[(int)fEn_L_n->GetXaxis()->GetNbins()];
    double yn[(int)fEn_L_n->GetXaxis()->GetNbins()];
    double errxn[(int)fEn_L_n->GetXaxis()->GetNbins()];
    double erryn[(int)fEn_L_n->GetXaxis()->GetNbins()];
    for(int i=1;i<fEn_L_n->GetXaxis()->GetNbins()+1;++i){
        double mean=0;
        double rms=0;
        double nz=0;
        for(int j=1;j<fEn_L_n->GetYaxis()->GetNbins()+1;++j){
            int ibin = fEn_L_n->GetBin(i,j);
            if(fEn_L_n->GetBinContent(ibin)>0){
                nz++;
                mean+=fEn_L_n->GetBinContent(ibin);
            }
        }
        if(nz>0)mean;
        if(mean>0){
            for(int j=1;j<fEn_L_n->GetYaxis()->GetNbins()+1;++j){
                int ibin = fEn_L_n->GetBin(i,j);
                if(fEn_L_n->GetBinContent(ibin)>0){
                    rms=pow(mean-fEn_L_n->GetBinContent(ibin),2);
                }
            }
            rms=sqrt(rms/nz);
        }
        xn[i]=fEn_L_n->GetXaxis()->GetBinCenter(i);
        yn[i]=mean;
        errxn[i]=0.01;
        erryn[i]=rms;
        // std::cout<<"mean+-rms = "<<mean<<"+-"<<rms<<std::endl;
        
    }
    c8->cd();
    TGraphErrors* errplotn = new TGraphErrors((int)fEn_L_n->GetXaxis()->GetNbins(),xn,yn,errxn,erryn);
    errplotn->GetXaxis()->SetTitle("E_neutron,MeV");
    errplotn->GetYaxis()->SetTitle("Mean_track_length,mm");
    errplotn->Draw("AP");
    gPad->Print("plotCheck/meanLength.C");
    
    
    std::cout<<nentries<<std::endl;
    std::cout<<"hello World!"<<std::endl;
    return 0;
    
    
}

