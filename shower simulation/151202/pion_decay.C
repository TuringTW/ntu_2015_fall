void pion_decay(){
	double M_pion = 0.13977; //gev
	double M_muon = 0.105658;
	double M_neu = 0;

	double Pi_lifetime = 0.26033e-6; //s
	double c = 3e8;//m/s
	TCanvas *c1 = new TCanvas("test", "test");
	TRandom *r = new TRandom();
	TH1F *hist_theta = new TH1F("hist_theta", "hist_theta", 100, 175, 180);
	TH1F *hist_energy = new TH1F("hist_energy", "hist_energy", 100, 0, 1);
	
	// r->SetSeed(392);


	for (int i = 0; i < 1000; ++i)
	{
		double px = 0;
		double py = 0;
		double pz = -1; // GeV
		double energy = sqrt(px*px+py*py+pz*pz+M_pion*M_pion);

		TLorentzVector pion(px, py, pz, energy);

		double decay_particle_mass[2] = {M_muon, M_neu};
		TGenPhaseSpace event;
		event.SetDecay(pion, 2, decay_particle_mass);
		event.Generate();

		TLorentzVector Muon = *(event.GetDecay(0));
		TLorentzVector Neu = *(event.GetDecay(1));

		double theta_Muon = Muon.Theta()*TMath::RadToDeg();
		double energy_Muon = Muon.E();
		cout << theta_Muon << endl;
		hist_theta->Fill(theta_Muon);
		hist_energy->Fill(energy_Muon);


	}
	
	c1->Divide(2,1);
	c1->cd(1);
	hist_theta->Draw();
	c1->cd(2);
	hist_energy->Draw();

}