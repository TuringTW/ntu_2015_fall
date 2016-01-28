double myfun(double *x, double *par){
  if (x[0]>=0)
  {
    return par[0];
  }else{
    return par[1];
  }
}
void function2(){
  TF1 *test1 = new TF1("myfun", myfun, -5, 5, 2);
  TF1 *test2 = new TF1("myfun2", myfun, -5, 5, 2);
  // TF1 *test3 = new TF1("test", "sin(x)+5", 0, 10);

  test1->SetParameters(2, 0);
  test1->SetLineStyle(1);
  test1->SetLineColor(2);
  test1->Draw();

  test2->SetParameters(1,0);
  test2->
  test2->SetLineStyle(2);
  test2->SetLineColor(4);
  test2->Draw("same");

  // test3->SetLineStyle(2);
  // test3->SetLineColor(4);
  // test3->Draw("same");
}
