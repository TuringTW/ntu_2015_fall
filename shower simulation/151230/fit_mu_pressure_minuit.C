TProfile * profile0;

void set_style()
{
	gStyle -> SetOptFit(1);
	gStyle -> SetOptStat(0);
	gStyle -> SetOptTitle(0);
	// gStyle -> Set0ptLogx(1);
	// gStyle -> Set0ptLogy(1);
	gStyle -> SetPadTickX(1);
	gStyle -> SetPadTickY(1);
	gStyle -> SetPadGridX(1);
	gStyle -> SetPadGridY(1);
	gStyle -> SetPadColor(0);

	gStyle -> SetCanvasColor(0);
	gStyle -> SetFrameBorderMode(0);
	gStyle -> SetCanvasBorderMode(0);
	gStyle -> SetPadBorderMode(0);
	gStyle -> SetPadColor(0);
	gStyle -> SetTitleColor(1);
	gStyle -> SetTitleFillColor(0);
	gStyle -> SetStatColor(1);
	gStyle -> SetOptTitle(0);

}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	float p0 = par[0];
	float p1 = par[1];
	float chi2 = 0;
	cout << "wetdg" << endl;
	cout << profile0 -> GetNbinsX() << endl;
	cout << "wetdsdkjhfjlasdfkjlg" << endl;
	for(int ibin=0; ibin<(profile0 -> GetNbinsX()); ibin++)
	{
		cout << "wet33123131dg" << endl;
		float x = profile0 -> GetBinCenter(ibin+1);
		cout << "wetdg123123" << endl;
		float nevt_obs = profile0 -> GetBinContent(ibin+1);
		if(nevt_obs == 0)
		{
			continue;
		}
		float error_obs = profile0 -> GetBinError(ibin+1);
		float nevt_exp = p0 + p1*x;
		// if(nevt_obs == 0) 
		// {
		// 	continue;
		// }
		cout << error_obs << endl;

		float chi = (nevt_obs - nevt_exp)/error_obs;
		chi2 = chi2 + chi*chi;
	}
	cout << "sdfgsdgsdfstdg" << endl;

	f = chi2;
	// printf("p0=%g, p1=%g, chi2=%g\n", p0, p1, chi2);
}

void fit_mu_pressure_minuit()
{
	set_style();

	TCanvas * c1 = new TCanvas("c1", "c1", 100, 100, 500, 500);

	TFile * file0 = TFile::Open("muon.root");

	profile0 = new TProfile("profile0", "profile0", 20, 0.9, 1.1, 0, 2000);

	TNtuple * ntuple = (TNtuple*) file0 -> Get("ntuple");
	ntuple -> Draw("nobs:press*1000>>profile0");
	profile0 -> SetMarkerStyle(20);
	profile0 -> SetMarkerColor(2);
	profile0 -> SetXTitle("Pressure (Atm)");
	profile0 -> SetYTitle("Muon Counts");
	profile0 ->Draw();
	
	TMinuit * minuit = new TMinuit();
	minuit -> SetFCN(fcn);

	// Parameter initialization, start, step, low, high
	int temp1 = 0;
	int temp2 = 0;
	int temp3 = 0;
	minuit -> mnparm(0, "p0", 0, 0.01, -1, 100, temp1);
	minuit -> mnparm(1, "p1", 0, 0.01, -200, 200, temp2);
	// minuit -> FixParameter(1);

	// Run Minuit
	Double_t arglist[1] = {500};	// arglist[0] is number of trials 
	minuit -> mnexcm("MIGRAD", arglist, 2, temp3);

	// Obtain fit parameter
	double p0_fit, p0_fit_err, p1_fit, p1_fit_err;
	minuit -> GetParameter(0, p0_fit, p0_fit_err);
	minuit -> GetParameter(1, p1_fit, p1_fit_err);

	// Draw fit line
	profile0 -> Draw();
	TF1 * function = new TF1("function", "[0]+[1]*x", 0.9, 1.1);
	function -> SetParameter(0, p0_fit);
	function -> SetParameter(1, p1_fit);
	function -> SetLineColor(4);
	function -> SetLineColor(6);
	function -> Draw("same");

	profile0 -> Fit("pol1", "", "same");
	
}