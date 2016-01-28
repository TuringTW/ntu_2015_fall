void gnthm(){

	// data from Seaborg, et al, Phys Rev 77 26 (1950) 
	double h = 3600;
	double d = 86400;
	double m = 2.6298E6; 
	double y = 3.15576E7;

	double z[24];
	double q[24];
	double t[24];
	// // z: atomic number
	// // q: alpha energy
	// // t: measured half-life

	z[1] = z[2] = 94; // Cm
	z[3] = z[4] = z[5] = 92; // Pu
	z[6] = z[7] = z[8] = z[9] = z[10] = 90;  // U 
	z[11] = z[12] = z[13] = z[14] = 88; z[15] = z[16] = z[17] = 86; // Th
	z[18] = z[19] = z[20] = 84; // Ra
	z[21] = z[22] = z[23] = z[24] = 82; // Po

	q[1] = 6.18; t[1] = 150*d; // Cm242 

	q[2] = 6.37; t[2] = 26.8* d;// Cm240

	q[3] = 5.60; t[3] = 92* y;// Pu238

	q[4] = 5.85; t[4] = 2.7* y;// Pu236

	q[5] = 6.26; t[5] = 8.5* h;// Pu234 

	q[6] = 4.25; t[6] = 4.51E9* y;// U238

	q[7] = 4.84; t[7] = 2.35E5* y; // U234

	q[8] = 5.40; t[8] = 70* y; // U232

	q[9] = 5.96; t[9] = 20.8* d;// U230 

	q[10] = 6.83; t[10] = 9.3*60;// U228
	    
	q[11] = 4.05; t[11] = 1.39*pow(10,10)*y;// Th232 

	q[12] = 4.76; t[12] = 8.0*pow(10,4) *y;// Th230

	q[13] = 5.52; t[13] = 1.90* y;// Th228

	q[14] = 6.41; t[14] = 30.9*60;// Th226

	q[15] = 4.88; t[15] = 1622*y;// Ra226

	q[16] = 5.78; t[16] = 3.64*d;// Ra224

	q[17] = 6.62; t[17] = 38;// Ra222

	q[18] = 5.59; t[18] = 3.83* d;// Rn222

	q[19] = 6.39; t[19] = 54.5;// Rn220

	q[20] = 7.25; t[20] = 0.019;// Rn218

	q[21] = 6.12; t[21] = 3.05*60;// Po218

	q[22] = 6.89; t[22] = 0.158;// Po216

	q[23] = 7.83; t[23] = 1.5E-4;// Po214

	q[24] = 8.95; t[24] = 3.0E-7;// Po212

	double x[24], y1[24];

	for (int i = 1; i < 25; ++i)
	{
		x[i-1] = log10(t[i]/log(2));
		y1[i-1] = z[i]/sqrt(q[i]);
		cout << i << " " << x[i-1] << "," << y1[i-1] <<endl;
	}

	n = 24;
	TGraph *g = new TGraph(n,y1,x);
	g->Draw("A*");
	g->GetYaxis()->SetTitle("log10(tau/log(2))");
	g->GetXaxis()->SetTitle("z/sqrt(q)");
	g->SetTitle("Recitation5 : 1-1 Geiger-Nuttall Law ");
	g->Fit("pol1", "B");

}