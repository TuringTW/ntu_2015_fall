void profile_muon_theta(){
	TCanvas * c1 = new TCanvas("c1", "c1", 1200, 500);	
	TFile *f1 = TFile::Open("muon.root");
	// TNtuple *ntuple = (TNtuple*)gDirectory->Get("ntuple");
	TNtuple *ntuple = (TNtuple*)f1->Get("ntuple");

	c1 ->Divide(2, 0);
	c1 ->cd(1);
	
	ntuple ->SetMarkerSize(25);
	ntuple ->Draw("nobs:theta");

	c1 ->cd(2);
	TProfile *p1 = new TProfile("p1", "p1", 30, 0, 90, 0, 2000);
	ntuple ->Draw("nobs:theta*180/3.14>>p1");
	p1 ->SetMarkerColor(2);
	// p1->Fit("pol2");

	TF1 *fn1 = new TF1("fn1", "[0]*cos(x/180*3.14)", 0, 90)	;
	float a_high = 100;
	float a_low = 0;
	int nstep = 100;
	float a_step = (a_high-a_low)/nstep;
	TH1F * chi2_hist = new TH1F("chi2_hist", "chi2_hist", nstep, a_low, a_high);
	for (int i = 0; i < nstep; ++i)
	{
		float a = a_low + a_step*i;
		fn1 ->SetParameter(0, a);
		
		float chi2 = 0;
		for (int ibin = 0; ibin < p1->GetNbinsX(); ++ibin)
		{
			float theta = p1->GetBinCenter(ibin+1);
			float nevt_obs = p1->GetBinContent(ibin+1);

			if (nevt_obs == 0)continue;
			float error_bar = p1 ->GetBinError(ibin+1);
			float nevt_exp = fn1->Eval(theta);
			float chi = (nevt_obs - nevt_exp)/error_bar;
			chi2 = chi2 + chi*chi;
		}
		chi2_hist ->SetBinContent(i+1, chi2);

	}
	c1 ->Clear();
	c1 ->Divide(1, 2);
	c1 ->cd(1);
	chi2_hist ->Draw();

	c1 ->cd(2);
	float a_fit_result = chi2_hist->GetMinimumBin()*a_step+a_low;
	fn1 ->SetParameter(0, a_fit_result);

	p1->Draw();
	fn1 ->Draw("same");
}
