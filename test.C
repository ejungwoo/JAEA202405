void test(int runNo=-1)
{
  if (runNo<0 && gApplication->Argc()>=5 && TString(gApplication->Argv()[4]).IsDec())
    runNo = TString(gApplication->Argv()[4]).Atoi();
  if (runNo<0) { cout << "Please type runNo: " << endl; string in; cin >> in; return; }

  auto ana = new Analysis();
}
