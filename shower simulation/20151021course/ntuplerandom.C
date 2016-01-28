void ntuplerandom(){
	TFile *f = new TFile("ntuple.root", "RECREATE");[google-chrome]
name=google-chrome
baseurl=http://dl.google.com/linux/chrome/rpm/stable/$basearch
enabled=1
gpgcheck=1
gpgkey=https://dl-ssl.google.com/linux/linux_signing_key.pub
	TCanvas* c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
	TNtuple* ntuple = new TNtuple("test", "ntuple", "x:y:z");
	

	TRandom *r = new TRandom();
	r->SetSeed(212312413);
	float mean = 0, sigma = 1;
	for (int i = 0; i < 100000; ++i){

		// float x, y;
		float x = r->Gaus(mean, sigma);
		float y = r->Gaus(mean, sigma);
		float z = r->Gaus(mean, sigma);
		// r->Rannor(x, y);
		// float y = pow(3*x, 1.0/3);
		
		ntuple->Fill(x, y, z);
	}
	c1->Divide(3, 1);
	
	c1->cd(1);
	ntuple->Draw("x");
	
	c1->cd(2);
	ntuple->Draw("x:y");
	ntuple->SetMarkerColor(2);
	ntuple->Draw("x:y", "x**2+y**2<1", "same");
	
	c1->cd(3);
	ntuple->Draw("x:y:z");
	ntuple->SetMarkerColor(1);
	ntuple->Draw("x:y:z", "x**2+y**2+z**2<1", "same");

	f->Write();
	// hist->Draw("box");
	// hist->Fit("");

}