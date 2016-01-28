void h1f_test(){
  TH1F* hist = new TH1F("test", "test", 100, 0, 100);
  int sum = 0;
    for(int i = 0; i <= 100; i++){
      sum = sum + i;
      hist->Fill(sum);
    }
    hist->Draw();
}
