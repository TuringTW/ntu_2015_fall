const int MAX_N = 1000;
int i_mu_ref;
float otheta, phi, e;
float x[MAX_N], y[MAX_N], z[MAX_N], t[MAX_N];
int evt, n_mu;

float theta, a, a_err, b, b_err, chi2;
typedef struct energy_t
{
	float theta;
	float a;
	float a_err;
	float b;
	float b_err;
	float chi2;
};
void energy_fit(){
	TFile *file  = new TFile("energy.root");
	TFile *file1 = new TFile("mnu_e_theta.root", "recreate");
	TTree *muon_tree = (TTree*) file->Get("muon_tree");
	
	muon_tree->SetBranchAddress("evt", &evt);
	muon_tree->SetBranchAddress("n_mu", &n_mu);
	muon_tree->SetBranchAddress("e", &e);
	muon_tree->SetBranchAddress("theta", &otheta);
	muon_tree->SetBranchAddress("phi", &phi);
	muon_tree->SetBranchAddress("x", x);
	muon_tree->SetBranchAddress("y", y);
	muon_tree->SetBranchAddress("z", z);
	muon_tree->SetBranchAddress("t", t);

	TNtuple * e_ntuple = new TNtuple("e_ntuple", "e_ntuple", "theta:a:a_err:b:b_err:chi2");
	// p0 = 3.661
	// p1 = 0.04453
	TProfile *p1 = new TProfile("p1", "p1", 20, 100, 10000, 0, 2000);
	muon_tree ->GetEntry(0);
	TCanvas * c1 = new TCanvas("c1", "c1", 800, 500);
	int nentries = muon_tree->GetEntries();
	for (int itheta = 0; itheta < 61; ++itheta)
	{
		c1 ->Divide(1, 1);
		c1 ->cd(1);


		char temp0[20];
		sprintf(temp0, "theta==%d", itheta);
		p1->SetTitle(temp0);
		muon_tree ->Draw("n_mu:e>>p1", temp0);

		p1->SetXTitle("Energy(GeV)");
		p1->SetYTitle("n_mu(#)");
		gStyle -> SetOptFit(1);
		gStyle -> SetOptStat(0);
		gStyle -> SetOptTitle(1);
		gStyle->SetOptStat(1110);
		
		p1->Fit("pol1");
		TF1 *f1 = (TF1*) p1->GetFunction("pol1");

		float chi2 = f1->GetChisquare();
		float a = f1->GetParameter(0);
		float a_err = f1->GetParError(0);
		float b = f1->GetParameter(1);
		float b_err = f1->GetParError(1);

		e_ntuple -> Fill(itheta,a,a_err,b,b_err,chi2);
		cout << itheta << " " << a << " " << b << endl;
		char temp[100];
		sprintf(temp, "./result/%d.jpg", (itheta+1000));
		c1->Print(temp);
		c1->Clear();
		p1->Reset();
	}
	c1 ->Divide(2, 1);
	c1 ->cd(1);
	e_ntuple -> SetMarkerStyle(2);
	e_ntuple -> Draw("theta:a");
	
	c1 ->cd(2);
	e_ntuple -> SetMarkerStyle(2);
	e_ntuple -> Draw("theta:b");

	c1->Print("./result/theta2ab.jpg");

	// file1->cd();
	// file1 ->Write();
	// file1->Close();
}