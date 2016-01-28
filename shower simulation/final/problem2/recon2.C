#include "TTree.h"
const double pi = TMath::Pi();

const int MAX_N = 1000;
int evt;
int n_mu;
Float_t theta;
Float_t phi;
Float_t e;
Float_t x[MAX_N];
Float_t y[MAX_N];
Float_t z[MAX_N];
Float_t t[MAX_N];
int i_mu_ref;

float a, b, a_err, b_err, etheta;

TRandom *r;

double Dt_Point_Plane(double x, double y, double z, double theta, double phi, double d){

  double a=cos(phi)*sin(theta);
  double b=sin(phi)*sin(theta);
  double c=cos(theta);

  if (a==0 && b==0 && c==0) {
    a=1.0000e-8;
    b=1.0000e-8;
    c=1.0000e-8;
  }
  return (
  a*x+b*y+c*z+d)/0.0003; //0.0003 is C in km/ns
}
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
  
  float theta=TMath::Pi()-par[0]*TMath::DegToRad();
  float phi=par[1]*TMath::DegToRad()+TMath::Pi();
  float d=0;
  float t_ref=t[i_mu_ref];
  float x_ref=x[i_mu_ref];
  float y_ref=y[i_mu_ref];
  float z_ref=z[i_mu_ref];
  float t_exp_ref=Dt_Point_Plane(x_ref, y_ref, z_ref, theta, phi, d);
  
  double chisq = 0;
  for (int i=0; i<n_mu ; i++){
    float dt_err=20;
    float dt_obs=t[i]-t_ref;
    float dt_exp=Dt_Point_Plane(x[i], y[i], z[i], theta, phi, d)-t_exp_ref;
    float delta=dt_obs-dt_exp;
    float chi=delta/dt_err;
    chisq=chisq+chi*chi;
  }
  f = chisq;
  
  
}


void recon2(){
  TFile *f1 = new TFile("skymap.root", "recreate");
  TFile *file  = new TFile("mnu_e_theta.root");
  TFile *file1  = new TFile("Muon.root");

  TNtuple *e_ntuple = (TNtuple*) file->Get("e_ntuple");
  e_ntuple->SetBranchAddress("theta", &etheta);
  e_ntuple->SetBranchAddress("a", &a);
  e_ntuple->SetBranchAddress("a_err", &a_err);
  e_ntuple->SetBranchAddress("b", &b);
  e_ntuple->SetBranchAddress("b_err", &b_err);
  // e_ntuple->SetBranchAddress("chi2", &chi2);


  TTree *Tree_Muon = (TTree*) file1->Get("muon_tree");
  Tree_Muon->SetBranchAddress("evt",&evt);
  Tree_Muon->SetBranchAddress("n_mu",&n_mu);
  Tree_Muon->SetBranchAddress("theta",&theta);
  Tree_Muon->SetBranchAddress("phi",&phi);
  Tree_Muon->SetBranchAddress("e",&e);
  Tree_Muon->SetBranchAddress("x",x);
  Tree_Muon->SetBranchAddress("y",y);
  Tree_Muon->SetBranchAddress("z",z);
  Tree_Muon->SetBranchAddress("t",t);

  TNtuple *recon = new TNtuple("recon","recon","theta:phi:n_mu:chi2:e_recon");
  int nentries = (int)Tree_Muon->GetEntries();


  
  TMinuit *gMinuit = new TMinuit(); 
  gMinuit->SetFCN(fcn);

  r= new TRandom();

  for (int i=0 ; i< nentries ; i++){
    Tree_Muon->GetEntry(i); 
    
    if (n_mu<5) continue;

    double minchi2 = 1e12;
    double theta_fit, theta_fit_err, phi_fit, phi_fit_err;
    i_mu_ref=TMath::LocMin(n_mu, t);
    for (int j = 0; j < 100; ++j)
    {
      //Run Minuit
      //Parameter initialization, start, step, low, high
      int tmp1=0, tmp2=0, tmp3=0;
      double theta_start=acos(r ->Rndm())*180./pi;
      double phi_start=360*(r->Rndm());

      gMinuit->mnparm(0, "theta", theta_start, 0.01, 0, 90, tmp1); 
      gMinuit->mnparm(1, "phi", phi_start, 0.01, 0, 360, tmp2);      
      
      double arglist[1]={500}; //arglist[0] is number of trials
      gMinuit->mnexcm("MIGRAD", arglist,3,tmp3);  
      double chi2,edm,errdef;
      int nvpar,nparx,icstat;  

      gMinuit->mnstat(chi2,edm,errdef,nvpar,nparx,icstat);            
      //Obtain Fit parameter
      if (chi2 < minchi2){
        minchi2 = chi2;
        gMinuit->GetParameter(0, theta_fit, theta_fit_err);
        gMinuit->GetParameter(1, phi_fit, phi_fit_err);
      }
    }

    int inttheta = (int)round(theta_fit); //Elevation angle
    double e_recon = 0;
    e_ntuple ->GetEntry(90-inttheta); 
    if (90-inttheta <= 60)
    {
      // modify to zenith angle
      e_recon = (n_mu-a)/b;
    }
    cout << "theta:" << inttheta <<"e:" << e_recon << " a:" << a << " b:" << b << " n_mu:" << n_mu << endl;
     
    
    recon->Fill(theta_fit,phi_fit,n_mu,minchi2,e_recon);


  }
  TCanvas * c1 = new TCanvas("c1", "c1", 1200, 600);
  
  c1 ->Divide(2, 1);
  c1 ->cd(1);
  recon->SetMarkerStyle(2);
  recon->Draw("e_recon:n_mu");
  TH2F * skymap = new TH2F("skymap", "skymap", 360, 0, 360, 90, 0, 90); 
  c1 ->cd(2);
  skymap->SetMarkerStyle(2);

  recon->Draw("theta:phi>>skymap");
  skymap->Draw("surf1");

  int bin_x, bin_y, bin_z;
  int max_bin = skymap ->GetMaximumBin();
  skymap ->GetBinXYZ(max_bin, bin_x, bin_y, bin_z);
  cout << "x:" << bin_x << " y:" << bin_y << " z:" << bin_z << endl;
  

  f1->cd();

  recon->Write();
  skymap->Write();
  f1->Close();

}

 
  
