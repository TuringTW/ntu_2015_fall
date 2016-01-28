TProfile *profile1;
const int MAX_N = 1000;
int i_mu_ref;
float theta, phi, e;
float x[MAX_N], y[MAX_N], z[MAX_N], t[MAX_N];
int evt, n_mu;

double Dt_Point_Plain(double x, double y, double z, double theta, double phi, double d){
	double a = cos(phi)*sin(theta);
	double b = sin(phi)*sin(theta);
	double c = cos(theta);

	if (a == 0 && b==0 && c==0)
	{
		a = 1.0000e-8;
		b = 1.0000e-8;
		c = 1.0000e-8;
	}
	return (a*x+b*y+c*z+d)/0.0003; //0.0003KM/ns
}
void fcn(int &npr, double *gin, double &f, double *par, int iflag){

	double theta = TMath::Pi() - par[0];
	double phi = par[1] + TMath::Pi();
	double d = par[2];
	double t_ref = t[i_mu_ref];
	double x_ref = x[i_mu_ref];
	double y_ref = y[i_mu_ref];
	double z_ref = z[i_mu_ref];

	double t_exp_ref = Dt_Point_Plain(x_ref, y_ref, z_ref, theta, phi, d);

	double chi2 = 0;
	for (int i = 0; i < n_mu; ++i)
	{
		double dt_error = 20*sqrt(2);
		double dt_obs = t[i] - t_ref;
		double dt_exp = Dt_Point_Plain(x[i], y[i], z[i], theta, phi, d) - t_exp_ref;
		double chi = (dt_obs - dt_exp)/dt_error;
		chi2 = chi2 + chi*chi;

	}
	f = chi2;

}

void theta_phi_grid(){
	TCanvas * c1 = new TCanvas("c1", "c1", 1200, 500);	
	TFile *f1 = TFile::Open("muon.root");
	TTree *muon_tree = (TTree*)gDirectory->Get("muon_tree");
	
	muon_tree->SetBranchAddress("evt", &evt);
	muon_tree->SetBranchAddress("n_mu", &n_mu);
	muon_tree->SetBranchAddress("e", &e);
	muon_tree->SetBranchAddress("theta", &theta);
	muon_tree->SetBranchAddress("phi", &phi);
	muon_tree->SetBranchAddress("x", x);
	muon_tree->SetBranchAddress("y", y);
	muon_tree->SetBranchAddress("z", z);
	muon_tree->SetBranchAddress("t", t);
	
	int nentries = (int) muon_tree ->GetEntries();
	TH2F * h2 = new TH2F("h2", "h2", 360, 0, 360, 90, 0, 90);	

	for (int i = 0; i < 1; ++i)
	{
		muon_tree ->GetEntry(i);
		i_mu_ref = TMath::LocMin(n_mu, t);

		double dummy[100];
		int iflag = 0;
		int npar = 3;
		double par[npar];
		for (int itheta = 0; itheta < 180; ++itheta)
		{
			for (int iphi = 0; iphi < 360; ++iphi)
			{
				double chi2;
				par[0] = double(itheta)*TMath::DegToRad();
				par[1] = double(iphi)*TMath::DegToRad();
				par[2] = 0;
				fcn(npar, dummy, chi2, par, iflag);
				h2->Fill(iphi, itheta, chi2);
			}
		}
	}
	// h2 ->SetAxisRange(20, 70, "X");
	// h2 ->SetAxisRange(20, 70, "Y");
	h2 ->SetMaximum(500);

	h2 ->Draw("surf3");
}