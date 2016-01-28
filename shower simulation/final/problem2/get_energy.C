const int MAX_N = 1000;
int evt;
float n_mu;
float theta, phi, e_recon;

int i_mu_ref;
float etheta, a, a_err, b, b_err, chi2;
void get_energy(){
	// TFile *file  = new TFile("mnu_e_theta.root");
	TFile *file1  = new TFile("skymap.root");
	TNtuple *recon = (TNtuple*) file1->Get("recon");
	recon->SetBranchAddress("n_mu",&n_mu);
	recon->SetBranchAddress("theta",&theta);
	recon->SetBranchAddress("phi",&phi);
	recon->SetBranchAddress("e_recon",&e_recon);


	TH1F * energy_bg = new TH1F("energy_bg", "energy_bg", 100, 100, 10000);
	TH1F * energy_sg = new TH1F("energy_sg", "energy_sg", 100, 100, 10000);
	TH1F * energy_total = new TH1F("energy_total", "energy_total", 100, 100, 10000);

	TCanvas * c1 = new TCanvas("c1", "c1", 1200, 400);

	
	c1 ->Divide(2, 1);
	c1 ->cd(1);
	c1->SetLogx();
	recon->Draw("e_recon>>energy_bg", "((phi>115|phi<105)&theta<60&theta>50)");
	energy_bg->SetXTitle("energy(GeV)");
	energy_bg->SetYTitle("count(#)");
	c1->Print("./energy_distribution_bg.jpg");
	c1 ->cd(2);
	c1->SetLogx();
	recon->Draw("e_recon>>energy_sg", "((phi<115&phi>105&theta<60&theta>50))");
	energy_sg->SetXTitle("energy(GeV)");
	energy_sg->SetYTitle("count(#)");
	c1->Print("./energy_distribution_sg.jpg");

	// c1 ->cd(3);
	// c1->SetLogx();	
	// recon->Draw("e_recon>>energy_total");
	// energy_total->SetXTitle("energy(GeV)");
	// energy_total->SetYTitle("count(#)");
	// c1->Print("./energy_distribution_total.jpg");
	c1->Clear();
	c1 ->Divide(2, 1);
	c1 ->cd(1);
	c1->SetLogx();
	recon->Draw("e_recon>>energy_bg", "((phi>115|phi<105)&theta<60&theta>50)*3.41349/120", "HIST");
	energy_bg->SetXTitle("energy(GeV)");
	energy_bg->SetYTitle("count(#)");
	c1->Print("./energy_distribution_bg_scaled.jpg");
	c1 ->cd(2);
	c1->SetLogx();
	recon->Draw("e_recon>>energy_sg", "(((phi<115&phi>105&theta<60&theta>50))-((phi>115|phi<105)&theta<60&theta>50)*3.41349/120)>0", "HIST");
	energy_sg->SetXTitle("energy(GeV)");
	energy_sg->SetYTitle("count(#)");
	c1->Print("./energy_distribution_sg.jpg_scaled");


}