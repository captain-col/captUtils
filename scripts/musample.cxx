#include <memory>
#include <iostream>
#include <fstream>

// Include to get Event Loop.
#include <eventLoop.hxx>
#include <TReconTrack.hxx>
#include <TReconHit.hxx>
#include <TChannelId.hxx>
#include <TDigitContainer.hxx>
#include <TDigit.hxx>
#include <TPulseDigit.hxx>
#include <TDigitProxy.hxx>
#include <TChannelInfo.hxx>
#include <CaptGeomId.hxx>
#include <TGeometryId.hxx>
// Includes for ROOT classes
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TPad.h>
#include <HEPUnits.hxx>
#include <set>


class TMuSampleLoop: public CP::TEventLoopFunction {
public:
    TMuSampleLoop() {}
    virtual ~TMuSampleLoop() {}
    void Initialize() {
flifeTime = new TH1F("lifetime","lifeTime",10000,0,10000);
fcharges =  new TH2F("charges", "charges",50,0,36,10000,100000);
 fChargeColl =  new TH2F("chargeColl", "chargeColl",350,0,350,10000,0,100000);
fchargeslog =  new TH2F("Collected Charge X Plane", "Collected Charge X Plane",350,0,350,62,7,13); 
fchargeslog->GetXaxis()->SetTitle("Wire#");
fchargeslog->GetYaxis()->SetTitle("ln(charge,e)");  
fAngLenZ =  new TH2F("AngLenZ", "AnglLeanZ",10,0,600,36,0,180);
 fNumberTracks = new TH1F("NTracks","NTracks",10,0,10);

  }
    bool operator () (CP::TEvent& evt) {
      int ntracks=0;

bool muon=false;

  CP::THandle<CP::TReconObjectContainer> RecCont = evt.Get<CP::TReconObjectContainer>("~/fits/TCaptainRecon/final");
//CP::THandle<CP::TReconObjectContainer> RecCont = evt.Get<CP::TReconObjectContainer>("~/fits/TCaptainRecon/TTracking2D/final");
std::set<CP::THandle<CP::THit>> hitsSet;

if(RecCont){

   for(std::size_t i=0; i < RecCont->size();i++){
 std::unique_ptr<TH1F> chargeTime( new TH1F("charge","charge",32,0,320));       
        CP::THandle<CP::TReconTrack> track =(*RecCont)[i];
	if(track)ntracks++;
	TLorentzVector startL = track->GetFront()->GetPosition();
	TLorentzVector stopL = track->GetBack()->GetPosition();
	TVector3 start;
	TVector3 stop;
	if(startL.X()>stopL.X()){	
	TVector3 startch(startL.X(),startL.Y(),startL.Z());
	TVector3 stopch(stopL.X(),stopL.Y(),stopL.Z());
	start=startch;
	stop=stopch;
	}else{
	TVector3 stopch(startL.X(),startL.Y(),startL.Z());
        TVector3 startch(stopL.X(),stopL.Y(),stopL.Z());
        start=startch;
        stop=stopch;	
	}	

	if(start.Z()>-1350 && start.Z()<350)continue;
	if(stop.Z()>-1350 && stop.Z()<350)continue;
	
	
	TVector3 dir=stop-start;
	TVector3 dirZ(0,0,1);
	TVector3 dirX(1,0,0);
	TVector3 dirU(1,1.732,0);
	TVector3 dirV(1,-1.732,0);
	double	angle = dir.Angle(dirZ);
	double angleX=dir.Angle(dirX);
	double angleU=dir.Angle(dirU);
	double angleV=dir.Angle(dirV);
	angle = angle*180/3.14;
	double cosX = cos(angleX);
	if(cosX<0)cosX*=-1;
//	cosX=1;	
	double length = dir.Mag();
	double lengthZ = fabs(start.Z()-stop.Z());
//	if((angle>30 && angle<60) || (angle>120 && angle<150)){
	if(lengthZ>300){
	fAngLenZ->Fill(lengthZ,angle);
std::cout<<"event="<<evt.GetContext().GetEvent()<<"; run="<<evt.GetContext().GetRun()<<std::endl;
muon=true;
CP::THandle<CP::THitSelection> hits=track->GetHits();

for(CP::THitSelection::iterator h = hits->begin();h!=hits->end();++h){

for(int i=0;i<(*h)->GetConstituentCount();++i){

CP::THandle<CP::THit> constit = (*h)->GetConstituent(i);
int plane =  CP::GeomId::Captain::GetWirePlane(constit->GetGeomId());

if(plane==CP::GeomId::Captain::kXPlane){
hitsSet.insert(constit);
}
}
}
}
//}

if(hitsSet.size()>6){

CP::THitSelection hitsSelect;
for(std::set<CP::THandle<CP::THit>>::iterator it=hitsSet.begin();it!=hitsSet.end();++it){

hitsSelect.push_back(*it);

 fChargeColl->Fill(CP::GeomId::Captain::GetWirePlane((*it)->GetGeomId()),(*it)->GetCharge());

}

std::sort(hitsSelect.begin(),hitsSelect.end(),[](CP::THandle<CP::THit> lh, CP::THandle<CP::THit> rh){
double ltime = lh->GetTime();
double rtime = rh->GetTime();
return ltime<rtime;
});

double setSize = hitsSet.size();
std::cout<<setSize<<std::endl;


double chargeS = hitsSelect.front()->GetCharge();

chargeS=(chargeS+hitsSelect[1]->GetCharge()+hitsSelect[2]->GetCharge())/3;

double timeS = hitsSelect.front()->GetTime();
double chargeF = hitsSelect.back()->GetCharge();

chargeF=(chargeF+hitsSelect[(int)setSize-2]->GetCharge()+hitsSelect[(int)setSize-3]->GetCharge())/3;

double timeF= hitsSelect.back()->GetTime();

double step = 200/setSize;
int count = 0;
for(CP::THitSelection::iterator it=hitsSelect.begin();it!=hitsSelect.end();++it){
if(chargeS>chargeF){
fcharges->Fill(count*step,(*it)->GetCharge()/cosX);
fchargeslog->Fill(count*step,log((*it)->GetCharge()/cosX));
chargeTime->Fill(count*step,(*it)->GetCharge());
}

if(chargeS<chargeF){
fcharges->Fill(step*(setSize-count),(*it)->GetCharge()/cosX);
std::cout<<(*it)->GetCharge()/cosX<< " ; "<<log((*it)->GetCharge()/cosX)<<std::endl;
fchargeslog->Fill(step*(setSize-count),log((*it)->GetCharge()/cosX));
chargeTime->Fill(count*step,(*it)->GetCharge());
}
count++;
}

TF1* exp1 = new TF1("exp1","[2]+[1]*exp(-x/[0])",0,320);
exp1->SetParameter(0,8000);
exp1->SetParameter(1,50);
exp1->SetParameter(2,10);
//chargeTime->Fit(exp1,"R");
//std::cout<<exp1->GetParameter(0)<<std::endl;
//flifeTime->Fill(exp1->GetParameter(0));
delete exp1;
}


hitsSet.clear();

}


 }
 fNumberTracks->Fill(ntracks);
	if(muon){
      return true;}
else{
return false;
}
      }
    // Called at least once.  If multiple file are open, it will be called
    // for each one.   Notice there are two forms...
    void Finalize(CP::TRootOutput * const output) {
TF1 *exp = new TF1("exp","[1]*exp(-x/[0])",5000,70000);
exp->SetParameter(1,50);
exp->SetParameter(0,8000);
fcharges->Fit(exp,"R");
fcharges->Draw("COLZ");
gPad->Print("Charges.C");

TF1 *line = new TF1("line","[0]*x+[1]",0,200);
line->FixParameter(0,-0.012);
line->FixParameter(1,10.75);
fchargeslog->Draw("COLZ");
line->Draw("same");
gPad->Print("ChargesLog.C");
fAngLenZ->Draw("COLZ");
gPad->Print("AngLenZ.C");

flifeTime->Draw();
gPad->Print("LifeTime.C");

 fChargeColl->Draw("COLZ");
 gPad->Print("ChargeColl.C");

 fNumberTracks->Draw();
 gPad->Print("NumberTracks.C");
 

}


private:
TH1F* flifeTime;
TH2F* fcharges;
TH2F* fchargeslog;
TH2F* fAngLenZ;
  TH2F* fChargeColl;
  TH1F* fNumberTracks;
};

int main(int argc, char **argv) {
    TMuSampleLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
