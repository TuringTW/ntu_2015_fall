TProfile *profile1;
void fcn(int &npr, double *gin, double &f, double *par, int iflag){
	double p0 = par[0];
	double p1 = par[1];
	float chi2 = 0;
	for (int ibin = 0; ibin < profile1->GetNbinsX(); ++ibin)
	{
		float x = profile1->GetBinCenter(ibin+1);
		float nevt_obs = profile1->GetBinContent(ibin+1);

		if (nevt_obs == 0)continue;
		float error_bar = profile1 ->GetBinError(ibin+1);
		float nevt_exp = p0+p1*x;
		float chi = (nevt_obs - nevt_exp)/error_bar;
		chi2 = chi2 + chi*chi;
	}
	f = chi2;

}
void profile_muon_press(){
	TCanvas * c1 = new TCanvas("c1", "c1", 1200, 500);	
	TFile *f1 = TFile::Open("muon.root");
	TNtuple *ntuple = (TNtuple*)gDirectory->Get("ntuple");
	// TNtuple *ntuple = (TNtuple*)f1->Get("ntuple");
	profile1  = new TProfile("profile1", "profile1", 30, 0.9, 1.1, 0, 2000);

	ntuple ->Draw("nobs:press*1000>>profile1");
	profile1->Fit("pol1", "", "same");


	int temp1 = 0;
	int temp2 = 0;
	int temp3 = 0;
	TMinuit *gM = new TMinuit();
	gM->SetFCN(fcn);

	gM->mnparm(0, "p0", 0, 0.01, 0, 200, temp1);
	gM->mnparm(1, "p1", 0, 0.01, -100, 200, temp2);

	double arglist[1]= {500};
	gM->mnexcm("MIGRAD", arglist, 2, temp3);

	double p0_fit, p0_error, p1_fit, p1_error;
	gM->GetParameter(0, p0_fit, p0_error);
	gM->GetParameter(1, p1_fit, p1_error);
// draw
	profile1->Draw();
	profile1->Fit("pol1", "", "same");
	TF1 *func1 = new TF1("func1", "[0]+[1]*x", 0.9, 1.1);
	func1 -> SetParameter(0, p0_fit);
	func1 -> SetParameter(1, p1_fit);
	func1 ->Draw("same");

	func1->SetLineColor(4);
	func1->SetLineWidth(5);
	func1->SetLineStyle(2);





}