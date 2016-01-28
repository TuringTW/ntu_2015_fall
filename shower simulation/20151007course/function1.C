cp void function1(){
  TF1 *test1 = new TF1("test", "x+1", 0, 10);
  TF1 *test2 = new TF1("test", "x*x", 0, 10);
  TF1 *test3 = new TF1("test", "sin(x)+5", 0, 10);

  test1->SetLineStyle(1);
  test1->SetLineColor(2);
  test1->Draw();

  test2->SetLineStyle(2);
  test2->SetLineColor(4);
  test2->Draw("same");

  test3->SetLineStyle(2);
  test3->SetLineColor(4);
  test3->Draw("same");
}
