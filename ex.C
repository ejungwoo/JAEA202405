void ex()
{
  double mp = 938.791;
  double m50cr = 46524.8;
  double mt = 2809.45;
  double m48cr = 44669.3;
  double tp = 30.; // beam energy
  double tt; // triton energy (input)

  double val1 = mp*mp + m50cr*m50cr + mt+mt + 2*(mp+tp)*m50cr;
  double val2 = 2*(mt+tt)*m50cr - 2*(mp*tp)*(mt*tt) 
  double val3 = 2*sqrt((mp+tp)*(mp+tp)-mp*mp) * sqrt((mt+tt)*(mt+tt)-mt*mt)*cos(theta)
  double ex = sqrt(val1-val2+val3) - m48cr;
}
