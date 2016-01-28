TRandom *r;
int pdg_proton_p = 2112;
int pdg_proton_n = -2112;
int pdg_pion_p = 211;
int pdg_pion_0 = 111;
int pdg_pion_n = -211;
int pdg_muon_p = -13;
int pdg_muon_n = 13;
int pdg_neutrino = 14;

double tau_pion = 0.26033e-7; //s
double tau_muon = 2.2e-6;//s
const int n_branches = 4;

double M_proton = 0.938;
double M_pi = 0.13977;
double M_pi0 = 0.13498;
double M_mu = 0.105658;
double M_neu = 0;
double M_e = 511e-06; //GeV

double rho0 = 1.021/1000; //g/cm3
double A=8; //g/cm3   ==>air density

const double pi = TMath::Pi();
TNtuple *ntuple1 = new TNtuple("ntuple1", "ntuple1", "ea:eb");

struct particle
{
	int id;
	int flag;
	double e;
	double px;
	double py;
	double pz;
	double vertex_x;
	double vertex_y;
	double vertex_z;
	double fx;
	double fy;
	double fz;
};
const int bank_size = 1e5;
particle Ptcl_bank[bank_size];
int n_particle = 0;

TF1 *fn_muon_dxde;

double X2H(double x, double h0, double theta){
	double costheta = fabs(cos(theta));
	double h;
	if (costheta!=0)
	{
		h = -A*log(x*costheta/(rho0*1e5*A)+exp(-h0/A));
	}else{
		h = h0;
	}
	return h;
}



double dedx_muon(double *x, double *par){
	double E = x[0];
	int material_code = par[0];
	// cout << material_code << endl;
 
	double K = 0.307e-3; // Gev/(g/cm2)
	double I = 85.7e-9; //for air, Gev
	double ZA = 0.49919; // for air
	double delta = 0; // for air and rock
	if (material_code == 1)
	{
		I = 136.3e-9;
		ZA = 0.50;
	}

	double z = 1;
	double p = sqrt(E*E - M_mu*M_mu);
	double gamma = sqrt(1.+pow(p/M_mu, 2));
	double beta = sqrt(1.-1./(gamma*gamma));
	double g_b = p/M_mu;
	double tmax = 2*M_e*p*p/(M_mu*M_mu + M_e*M_e + 2*M_e*E);
	double dedx = K*z*z*ZA/(beta*beta)*(1/2.*log(2*M_e*g_b*g_b*tmax/(I*I)) - beta*beta -delta/2.);
	return dedx;
}
double dxde_muon(double *x, double *par){
	double dedx = dedx_muon(x, par);
	double dxde = 1./dedx;
	return dxde;
}
double get_muon_range(double E, int material_code){
	fn_muon_dxde->SetParameters(0, material_code);
	if (E<0.2){return 0;} // ignore low energy
	else{
		double R = fn_muon_dxde->Integral(0.2, E);
		return R;
	}
}
void get_muon_stop_position(double mu_E, double mu_theta, double mu_phi, double vertex_x, double vertex_y, double vertex_z, double &mu_stop_x, double &mu_stop_y, double &mu_stop_z){
	int material_code = 0; //air
	double range_muon = get_muon_range(mu_E, material_code);
	double X = r ->Exp(range_muon);

	double h0 = vertex_z;
	mu_stop_z = X2H(X, h0, mu_theta);
	double radius = mu_stop_z*tan(mu_theta);
	mu_stop_x = radius*cos(mu_phi);
	mu_stop_y = radius*sin(mu_phi);

}
void Ptcl_Register(int id, int flag, double px, double py, double pz, double x, double y, double z){
	double mass = 0;
	
	if (id == pdg_pion_0) mass = M_pi0;
	if (id == pdg_pion_p) mass = M_pi;
	if (id == pdg_proton_p ) mass = M_proton;


	Ptcl_bank[n_particle].id = id;
	Ptcl_bank[n_particle].flag = flag;
	Ptcl_bank[n_particle].e = sqrt(px*px+py*py+pz*pz+mass*mass);
	Ptcl_bank[n_particle].px = px;
	Ptcl_bank[n_particle].py = py;
	Ptcl_bank[n_particle].pz = pz;
	Ptcl_bank[n_particle].vertex_x = x;
	Ptcl_bank[n_particle].vertex_y = y;
	Ptcl_bank[n_particle].vertex_z = z;
	

	n_particle++;

}

void Get_Int_Posi(int id, double theta, double phi,double x, double y, double z, double &vertex_x, double &vertex_y, double &vertex_z){
	double lambda = 0;
	if (id == pdg_proton_p){lambda = 87;} //cm
	else if (id == pdg_pion_p){lambda = 116;} //cm
	else{
		cout << id << endl;
	}
	double X = r->Exp(lambda);
	double h0 = z;
	vertex_z = X2H(X, h0, theta);
	double radius = vertex_z*tan(theta);
	vertex_x = radius*cos(phi);
	vertex_y = radius*sin(phi);
}

void Nsplit(double all, double branch[], int branches){
	double frac[branches];
	double frac_sum = 0;
	for (int i = 0; i < branches; ++i)
	{
		frac[i] = r->Rndm();
		frac_sum += frac[i];
	}

	for (int i = 0; i < branches; ++i)
	{
		branch[i] = all*frac[i]/frac_sum;
	}
}
void twosplit(double all, double &a, double &b){
	double frac = r->Rndm();
	a = all*frac;
	b = all*(1-frac);
}

void Hillas_Split(double id, double E, double theta, double phi, double vertex_x, double vertex_y, double vertex_z){

	double E_init_proton = E; // GeV
	double E_proton = 0;
	double E_init_pion = 0;

	double mass = 0;
	if (id == pdg_pion_p) mass = M_pi;
	else mass = M_proton;
	
	while(E_proton<mass){
		twosplit(E_init_proton, E_proton, E_init_pion);
	}
	
	double p = sqrt(E_proton*E_proton-mass*mass);
	double px = -p*sin(theta)*cos(phi);
	double py = -p*sin(theta)*sin(phi);
	double pz = -p*cos(theta);


	Ptcl_Register(id, 1, px, py, pz, vertex_x, vertex_y, vertex_z);

	double pion_branch[n_branches];
	
	Nsplit(E_init_pion, pion_branch, n_branches);
	for (int i = 0; i < n_branches; ++i)
	{
		double energy_in = pion_branch[i];
		while(energy_in > 0){
			double E_a, E_b;
			twosplit(energy_in, E_a, E_b);
			
			energy_in = E_b - M_pi;
			if (energy_in>0)
			{
				// ?
				double px = -E_a*sin(theta)*cos(phi);
				double py = -E_a*sin(theta)*sin(phi);
				double pz = -E_a*cos(theta);

				Ptcl_Register(pdg_pion_p, 1, px, py, pz, vertex_x, vertex_y, vertex_z);
			}
		}
	}

	
}
double get_decay_length(int particle, double Pi_lifetime, double energy){
		double M = 0;
		if (particle == pdg_pion_p)
		{
			M = M_pi;
		}else{
			M = M_mu;
		}

		double t = r->Exp(Pi_lifetime);
		double c = 3e8;
		double gamma = energy/M;
		double length = c*t*gamma/1e3;
		return length;
}


void recursive_inter_with_air(){
	TFile *f1 = new TFile("muon.root", "recreate");
	TCanvas *c1 = new TCanvas("test", "test", 600, 700);

	TH1F * mu_c_1 = new TH1F("mu_c_1", "mu_c_1", 100, 0, 500);
	TH1F * mu_c_2 = new TH1F("mu_c_2", "mu_c_2", 100, 0, 500);
	TH1F * mu_c_3 = new TH1F("mu_c_3", "mu_c_3", 100, 0, 500);

	fn_muon_dxde = new TF1("f1", dxde_muon, M_mu, 1000, 1);
	TNtuple * ntuple = new TNtuple("ntuple", "ntuple", "nobs:detz:theta:press");
	// TNtuple *ntuple = new TNtuple("ntuple", "ntuple", "id:flag:e:px:py:pz:x:y:z");
	r = new TRandom();
	// r ->SetSeed(12232);
	for (int ievertex_nt = 0; ievertex_nt < 500; ievertex_nt++)
	{
		n_particle = 0;

		// init muon conunter
		int mu_counter1 = 0;
		int mu_counter2 = 0;
		int mu_counter3 = 0;

		double E = 1000; // GeV
		double theta = 0;
		// double theta = acos(r->Rndm());
		double phi = r->Rndm()*pi*2;
		double x = 0, y = 0, z = 1.0e10; //infinity
		double vertex_z, vertex_x, vertex_y; //km

		

		Get_Int_Posi(pdg_proton_p, theta, phi, x, y, z, vertex_x, vertex_y, vertex_z);

		Hillas_Split(pdg_proton_p, E, theta, phi, vertex_x, vertex_y, vertex_z);

		double proton_threshold = 10;
		double pion_threshold = 10;
		int i_particle = 0;
		while(i_particle < n_particle){
			int id = Ptcl_bank[i_particle].id;
			int flag = Ptcl_bank[i_particle].flag;
			double e = Ptcl_bank[i_particle].e;
			double x = Ptcl_bank[i_particle].vertex_x;
			double y = Ptcl_bank[i_particle].vertex_y;
			double z = Ptcl_bank[i_particle].vertex_z;
			if (id == pdg_proton_p && flag == 1 && e > proton_threshold)
			{
				Get_Int_Posi(id, theta, phi, x, y, z, vertex_x, vertex_y, vertex_z);
				if (vertex_z > 0)
				{
					Hillas_Split(id, e, theta, phi, vertex_x, vertex_y, vertex_z);
					Ptcl_bank[i_particle].flag = 0;
				}
			}
			if (id == pdg_pion_p &&flag == 1){

				Get_Int_Posi(id, theta, phi, x, y, z, vertex_x, vertex_y, vertex_z)	;			
				double l_to_interaction = sqrt((vertex_x-x)*(vertex_x-x)+(vertex_y-y)*(vertex_y-y)+(vertex_z-z)*(vertex_z-z));
				double l_to_decay = get_decay_length(pdg_pion_p, tau_pion, e);
				
				if (l_to_decay>l_to_interaction){
					if (vertex_z > 0){
						Hillas_Split(id, e, theta, phi, vertex_x, vertex_y, vertex_z);
						Ptcl_bank[i_particle].flag = 0;
					}
				}else{ // pion decay
					Ptcl_bank[i_particle].flag = 0;
					double px = Ptcl_bank[i_particle].px;
					double py = Ptcl_bank[i_particle].py;
					double pz = Ptcl_bank[i_particle].pz; // GeV
					TLorentzVector pion(px, py, pz, e);
					double decay_particle_mass[2] = {M_mu, M_neu};
					TGenPhaseSpace event;
					event.SetDecay(pion, 2, decay_particle_mass);
					event.Generate();
					TLorentzVector Muon = *(event.GetDecay(0));
					TLorentzVector Neu = *(event.GetDecay(1));

					double mu_theta = Muon.Theta();
					double mu_phi = Muon.Phi();
					
					if (mu_theta>pi/2.)	{
						double mu_E = Muon.E();
						double muon_decay_time = r ->Exp(tau_muon);
						double muon_decay_length = get_decay_length(pdg_muon_n, tau_muon, mu_E);

						float det_x = 0;
						float det_y = 0;
						float det_z = 0;
						mu_counter1++;
						double l_mu_vtx_to_det = sqrt(pow(vertex_x-det_x,2)+pow(vertex_y-det_y,2)+pow(vertex_z-det_z,2));
						// cout << "MuE" << mu_E << " Mass" <<M_mu<< " dl" << muon_decay_length << " lvd" << l_mu_vtx_to_det << endl;
						if (muon_decay_length > l_mu_vtx_to_det)	
						{
							mu_counter2++;

							double mu_stop_x;
							double mu_stop_y;
							double mu_stop_z;

							get_muon_stop_position(mu_E, mu_theta, mu_phi, vertex_x, vertex_y, vertex_z, mu_stop_x, mu_stop_y, mu_stop_z);
							double l_muon_stop = sqrt(pow(vertex_x-mu_stop_x,2)+pow(vertex_y-mu_stop_y,2)+pow(vertex_z-mu_stop_z,2));
							// cout << l_muon_stop << " " << l_mu_vtx_to_det <<endl;
							if (l_muon_stop > l_mu_vtx_to_det)
							{
								mu_counter3++;
							}
							
							
						}
					}
				}
			}
			i_particle++;
		}
		mu_c_1 ->Fill(mu_counter1);
		mu_c_2 ->Fill(mu_counter2);
		mu_c_3 ->Fill(mu_counter3);

		// for (int i_particle = 0; i_particle < n_particle; ++i_particle){
		// 	int id = Ptcl_bank[i_particle].id;
		// 	int flag = Ptcl_bank[i_particle].flag;
		// 	double e = Ptcl_bank[i_particle].e;
		// 	double px = Ptcl_bank[i_particle].px;
		// 	double py = Ptcl_bank[i_particle].py;
		// 	double pz = Ptcl_bank[i_particle].pz;
		// 	double x = Ptcl_bank[i_particle].vertex_x;
		// 	double y = Ptcl_bank[i_particle].vertex_y;
		// 	double z = Ptcl_bank[i_particle].vertex_z;
		// 	ntuple->Fill(id, flag, e, px, py, pz, x, y, z);
		// }
	}

	mu_c_3 ->Draw();
	mu_c_2 ->SetLineColor(2);
	mu_c_2 ->Draw("same");
	mu_c_1 ->SetLineColor(4);
	mu_c_1 ->Draw("same");
	// ntuple->Draw("z");
} 
