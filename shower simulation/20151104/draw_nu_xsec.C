void draw_nu_xsec(){
	TF1 *sigma_cc = new TF1("test", "5.53E-36*pow(x, 0.363)", 1.0E01, 1.0E11);
	TF1 *sigma_nc = new TF1("test1", "2.51E-36*pow(x, 0.363)", 1.0E01, 1.0E11);

	sigma_nc->Draw();
	sigma_cc->SetLineColor(3);
	sigma_cc->Draw("same");
}