void drawing_particle(){
	double M_pion = 0.13977; //gev
	double M_muon = 0.105658;
	double M_neu = 0;

	double Pi_lifetime = 0.26033e-6; //s
	double c = 3e8;//m/s
	TCanvas *c1 = new TCanvas("test", "test");
	TRandom *r = new TRandom();
	
	TView *view = TView::CreateView(1, 0, 0);
	view->ShowAxis();
	view->SetRange(-5, -5, 0, 5, 5, 10);

	for (int i = 0; i < 100; ++i)
	{
		double theta = acos(r->Rndm());
		double phi = r->Rndm()*TMath::Pi()*2.;

		double Mu_p = 1;//Gev
		double Mu_px = -Mu_p*sin(theta)*sin(phi);
		double Mu_py = -Mu_p*sin(theta)*cos(phi);
		double Mu_pz = -Mu_p*cos(theta); // GeV
		double Mu_E = sqrt(Mu_p*Mu_p+M_muon*M_muon);

		TLorentzVector Muon(Mu_px, Mu_py, Mu_pz, Mu_E);
		TParticle *particle_Muon = new TParticle();

		particle_Muon->SetMomentum(Muon);

		double vx = 0;
		double vy = 0;
		double vz = 10;//km

		particle_Muon->SetProductionVertex(vx, vy, vz, 0);
		particle_Muon->SetLineColor(2);
		particle_Muon->Draw();

	}


}