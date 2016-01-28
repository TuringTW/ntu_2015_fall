void signal_box(){
	TFile *file  = new TFile("skymap.root");
	TNtuple *recon = (TNtuple*) file->Get("recon");
	
	TCanvas * c1 = new TCanvas("c1", "c1", 1200, 600);
  
	c1 ->Divide(2, 1);
	c1 ->cd(1);
	recon ->SetMarkerStyle(2);
	recon ->Draw("theta:phi");
	
	TH2F * skymap = new TH2F("skymap", "skymap", 360, 0, 360, 90, 0, 90); 
	c1 ->cd(2);
	skymap->SetMarkerStyle(2);
	recon ->Draw("theta:phi>>skymap");	
	skymap->Draw("surf1");
	skymap->SetXTitle("phi(degree)");
	skymap->SetYTitle("theta(degree)");
	skymap->SetZTitle("count");

	c1->Print("result/skymap.jpg");

	int bin_x, bin_y, bin_z;
	int max_bin = skymap ->GetMaximumBin();
	skymap ->GetBinXYZ(max_bin, bin_x, bin_y, bin_z);
	max_value = skymap->GetMaximum();
	cout << "MAX value Point: (" << bin_x << ", " << bin_y << ") value:" << max_value << endl;

// draww scatter plot

	c1->Clear();
	c1->Divide(2, 1);
	c1 ->cd(1);
	TH1F * phi_dist = new TH1F("phi_dist", "phi_dist", 100, 0, 360);	
	recon ->Draw("phi>>phi_dist");
	phi_dist->SetXTitle("phi(degree)");
	phi_dist->SetYTitle("count");
	c1 ->cd(2);
	TH1F * theta_dist = new TH1F("theta_dist", "theta_dist", 100, 0, 90);	
	recon ->Draw("theta>>theta_dist");
	theta_dist->SetXTitle("theta(degree)");
	theta_dist->SetYTitle("count");

	c1->Print("result/theta_and_phi.jpg");

	
	

// calculate background property

	double sum = 0;
	double sumsq =0;
	double counter = 0;
	double signal_counter = 0;
	
	TH1F * bin_curve = new TH1F("bin_curve", "bin_curve", 100, 0, 10);

	for (int i = 0; i <= 360; ++i)
	{
		TH1F * temp_box = new TH1F("temp_box", "temp_box", 1, i, i+10);
		recon->Draw("phi>>temp_box", "theta>=50&theta<60");
		double bin_value = temp_box->GetBinContent(1);
		if (!(i>=95&&i<115))
		{
			
			counter++;
			sum += bin_value;
			sumsq +=bin_value*bin_value;
			bin_curve->Fill(bin_value);
			temp_box->Reset();
		}
		if (i==105)
		{
			signal_counter = bin_value;
		}
	}
	
	cout << sum << " " << sumsq << " "<< counter << endl;
	double mean = sum/counter;
	double stdev = sqrt((sumsq/counter)-mean*mean);
	cout << "mean of bg: " << mean << endl;
	cout << "stdev of bg: " << stdev << endl;
	cout << "evt in signal box: "<< signal_counter << endl;

	c1->Clear();
	c1->Divide(2, 1);
	c1->cd(1);
	recon->Draw("theta:phi");
	TBox *b1 = new TBox(0, 50, 105, 60);
	b1->SetFillStyle(0);
	b1->SetLineColor(4);
	b1->SetLineWidth(2);
	b1->Draw();
	TBox *b2 = new TBox(115, 50, 360, 60);
	b2->SetFillStyle(0);
	b2->SetLineColor(4);
	b2->SetLineWidth(2);
	b2->Draw();
	
	TBox *b_s = new TBox(105, 50, 115, 60);
	b_s->SetFillStyle(0);
	b_s->SetLineColor(2);
	b_s->SetLineWidth(1);
	b_s->Draw();

	c1->cd(2);
	bin_curve->Draw();
	bin_curve->SetXTitle("Bin value of BG");
	bin_curve->SetTitle("Bin value distribution of BG BOX");
	c1->Print("result/box_and_hist.jpg");



}