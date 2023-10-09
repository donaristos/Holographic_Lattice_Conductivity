
elem[0][0]=
0
;
elem[0][1]=
0
;
elem[0][2]=
0
;
elem[0][3]=
0
;
elem[0][4]=
0
;
elem[0][5]=
0
;
elem[0][6]=
0
;
elem[1][0]=
0
;
elem[1][1]=
0
;
elem[1][2]=
0
;
elem[1][3]=
0
;
elem[1][4]=
0
;
elem[1][5]=
0
;
elem[1][6]=
0
;
elem[2][0]=
0
;
elem[2][1]=
0
;
elem[2][2]=
0
;
elem[2][3]=
0
;
elem[2][4]=
0
;
elem[2][5]=
0
;
elem[2][6]=
0
;
elem[3][0]=
0
;
elem[3][1]=
0
;
elem[3][2]=
0
;
elem[3][3]=
0
;
elem[3][4]=
0
;
elem[3][5]=
0
;
elem[3][6]=
0
;
elem[4][0]=
0
;
elem[4][1]=
0
;
elem[4][2]=
0
;
elem[4][3]=
0
;
elem[4][4]=
0
;
elem[4][5]=
0
;
elem[4][6]=
0
;
elem[5][0]=
0
;
elem[5][1]=
0
;
elem[5][2]=
0
;
elem[5][3]=
0
;
elem[5][4]=
0
;
elem[5][5]=
0
;
elem[5][6]=
0
;
elem[6][0]=
0
;
elem[6][1]=
0
;
elem[6][2]=
0
;
elem[6][3]=
0
;
elem[6][4]=
0
;
elem[6][5]=
0
;
elem[6][6]=
0
;
if(j==j1){
elem[0][0]+=
0
;
elem[0][1]+=
0
;
elem[0][2]+=
0
;
elem[0][3]+=
0
;
elem[0][4]+=
0
;
elem[0][5]+=
0
;
elem[0][6]+=
0
;
elem[1][0]+=
D[-1 + d10](i,j,i1,j1)/2.
;
elem[1][1]+=
D[-1 + d10](i,j,i1,j1)/2.
;
elem[1][2]+=
0
;
elem[1][3]+=
0
;
elem[1][4]+=
0
;
elem[1][5]+=
0
;
elem[1][6]+=
0
;
elem[2][0]+=
0
;
elem[2][1]+=
0
;
elem[2][2]+=
D[-1 + d10](i,j,i1,j1)
;
elem[2][3]+=
0
;
elem[2][4]+=
0
;
elem[2][5]+=
0
;
elem[2][6]+=
0
;
elem[3][0]+=
0
;
elem[3][1]+=
0
;
elem[3][2]+=
0
;
elem[3][3]+=
0
;
elem[3][4]+=
D[-1 + d10](i,j,i1,j1)
;
elem[3][5]+=
0
;
elem[3][6]+=
0
;
elem[4][0]+=
((2*Qtt.diff(d01,i,j))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) - 
    Qr1(i,j)/(2.*(1 + Qtt(i,j))))*D[-1 + d10](i,j,i1,j1)
;
elem[4][1]+=
((-2*Qtt.diff(d01,i,j))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) - 
    Qr1(i,j)/(2.*(1 + Qtt(i,j))))*D[-1 + d10](i,j,i1,j1)
;
elem[4][2]+=
((2*Qtt.diff(d01,i,j))/(L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
    (Q11.diff(d01,i,j)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)))*D[-1 + d10](i,j,i1,j1)
;
elem[4][3]+=
2*D[-1 + d10](i,j,i1,j1)
;
elem[4][4]+=
(Q22.diff(d01,i,j)*(1 + Qtt(i,j))*D[-1 + d10](i,j,i1,j1))/
  (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2))
;
elem[4][5]+=
0
;
elem[4][6]+=
(4*h.diff(d01,i,j)*pow(mu,2)*(1 + Qtt(i,j))*D[-1 + d10](i,j,i1,j1))/
  (L*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[5][0]+=
-(a0(i,j)*D[-1 + d10](i,j,i1,j1))/(2.*(1 + Qtt(i,j)))
;
elem[5][1]+=
-(a0(i,j)*D[-1 + d10](i,j,i1,j1))/(2.*(1 + Qtt(i,j)))
;
elem[5][2]+=
(a0(i,j)*D[-1 + d10](i,j,i1,j1))/(2 + 2*Q11(i,j))
;
elem[5][3]+=
0
;
elem[5][4]+=
(a0(i,j)*D[-1 + d10](i,j,i1,j1))/(2 + 2*Q22(i,j))
;
elem[5][5]+=
2*D[-1 + d10](i,j,i1,j1)
;
elem[5][6]+=
B1*mu*a0(i,j)*D[-1 + d10](i,j,i1,j1)
;
elem[6][0]+=
0
;
elem[6][1]+=
0
;
elem[6][2]+=
0
;
elem[6][3]+=
0
;
elem[6][4]+=
0
;
elem[6][5]+=
0
;
elem[6][6]+=
D[-1 + d10](i,j,i1,j1)
;
}
if(i==i1){
elem[0][0]+=
0
;
elem[0][1]+=
0
;
elem[0][2]+=
0
;
elem[0][3]+=
0
;
elem[0][4]+=
0
;
elem[0][5]+=
0
;
elem[0][6]+=
0
;
elem[1][0]+=
(-4*Qtt.diff(d01,i,j)*D[-1 + d01](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (2*(1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[1][1]+=
0
;
elem[1][2]+=
0
;
elem[1][3]+=
0
;
elem[1][4]+=
0
;
elem[1][5]+=
0
;
elem[1][6]+=
0
;
elem[2][0]+=
((4*Qtt.diff(d01,i,j))/(pow(L,2)*(-12 + pow(mu,2))*(1 + Qtt(i,j))) - 
    (2*(1 + Q11(i,j))*Qr1(i,j))/(L*(1 + Qtt(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[2][1]+=
0
;
elem[2][2]+=
(-6*Q11.diff(d01,i,j)*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
  (2*(1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[2][3]+=
0
;
elem[2][4]+=
(2*Q22.diff(d01,i,j)*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
  ((-12 + pow(mu,2))*pow(L + L*Q22(i,j),2))
;
elem[2][5]+=
0
;
elem[2][6]+=
(8*h.diff(d01,i,j)*pow(mu,2)*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*(-12 + pow(mu,2)))
;
elem[3][0]+=
0
;
elem[3][1]+=
0
;
elem[3][2]+=
0
;
elem[3][3]+=
0
;
elem[3][4]+=
(-4*Q22.diff(d01,i,j)*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))) + 
  (2*(1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[3][5]+=
0
;
elem[3][6]+=
0
;
elem[4][0]+=
((2*Q11.diff(d10,i,j))/(L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
    (4*Qr1.diff(d01,i,j))/
     (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
    (2*Q11.diff(d01,i,j)*Qr1(i,j))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
    (2*Qrr.diff(d10,i,j))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) + 
    (2*Qtt.diff(d10,i,j))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) + 
    ((-12 + pow(mu,2))*pow(Qr1(i,j),2) + 
       (-12 + pow(mu,2))*pow(Q11(i,j),2)*pow(Qr1(i,j),2) + 
       2*Q11(i,j)*(-1 + (-12 + pow(mu,2))*pow(Qr1(i,j),2) - 
          Qtt(i,j)) - 4*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Qtt(i,j))))*
  D[-1 + d01](i,j,i1,j1)
;
elem[4][1]+=
0
;
elem[4][2]+=
((-2*Qtt.diff(d01,i,j)*Qr1(i,j))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
    (Q11.diff(d10,i,j)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) + 
    (2*Qr1.diff(d01,i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
    (2*Q11.diff(d01,i,j)*Qr1(i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) + 
    ((1 + Qtt(i,j))*(4 + Q11(i,j) + 2*Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)))*D[-1 + d01](i,j,i1,j1)
;
elem[4][3]+=
((-4*Qtt.diff(d01,i,j))/
      (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - Qr1(i,j)/L + 
     (2*Q11.diff(d01,i,j)*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))*
   D[-1 + d01](i,j,i1,j1) + (2*(1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[4][4]+=
((Q22.diff(d10,i,j)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) - 
    (2*Q22.diff(d01,i,j)*Qr1(i,j)*(1 + Qtt(i,j)))/
     ((-12 + pow(mu,2))*(1 + Q11(i,j))*pow(L + L*Q22(i,j),2)) - 
    ((Q22(i,j) - 2*Qtt(i,j))*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)))*
  D[-1 + d01](i,j,i1,j1)
;
elem[4][5]+=
(-8*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j))*D[-1 + d01](i,j,i1,j1))/
  (L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j)))
;
elem[4][6]+=
((4*h.diff(d10,i,j)*pow(mu,2)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
    (4*pow(mu,2)*h(i,j)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
    (8*h.diff(d01,i,j)*pow(mu,2)*Qr1(i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[5][0]+=
(a0(i,j)*Qr1(i,j)*D[-1 + d01](i,j,i1,j1))/(L + L*Qtt(i,j))
;
elem[5][1]+=
0
;
elem[5][2]+=
(-((a0(i,j)*Qr1(i,j))/(2*L + 2*L*Q11(i,j))) - 
    (a0.diff(d01,i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))*
  D[-1 + d01](i,j,i1,j1)
;
elem[5][3]+=
-((a0(i,j)*D[-1 + d01](i,j,i1,j1))/L)
;
elem[5][4]+=
(-((a0(i,j)*Qr1(i,j))/(2*L + 2*L*Q22(i,j))) + 
    (a0.diff(d01,i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))))*
  D[-1 + d01](i,j,i1,j1)
;
elem[5][5]+=
((-2*Qr1(i,j))/L - (Q11.diff(d01,i,j)*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
     (2*h.diff(d01,i,j)*B1*mu*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
     (Q22.diff(d01,i,j)*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))))*
   D[-1 + d01](i,j,i1,j1) + (2*(1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[5][6]+=
(-((B1*mu*a0(i,j)*Qr1(i,j))/L) + 
    (2*a0.diff(d01,i,j)*B1*mu*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[6][0]+=
(2*h.diff(d01,i,j)*D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[6][1]+=
0
;
elem[6][2]+=
-((h.diff(d01,i,j)*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
    (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))
;
elem[6][3]+=
0
;
elem[6][4]+=
(h.diff(d01,i,j)*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j)))
;
elem[6][5]+=
0
;
elem[6][6]+=
((2*Qtt.diff(d01,i,j))/(pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
     Qr1(i,j)/L - (Q11.diff(d01,i,j)*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
     (Q22.diff(d01,i,j)*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))))*
   D[-1 + d01](i,j,i1,j1) + (2*(1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
}
if(j==j1 && i==i1){
elem[0][0]+=
1
;
elem[0][1]+=
-1
;
elem[0][2]+=
0
;
elem[0][3]+=
0
;
elem[0][4]+=
0
;
elem[0][5]+=
0
;
elem[0][6]+=
0
;
elem[1][0]+=
(2*Qtt.diff(d02,i,j))/(pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (72 + 108*pow(A1,2) + 2*pow(mu,2) + 3*pow(A1,2)*pow(mu,2) - 
     144*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1)) - 96*exp(A1*mu*h(i,j)) + 
     96*Qtt(i,j) + 144*pow(A1,2)*Qtt(i,j) - 8*pow(mu,2)*Qtt(i,j) - 
     12*pow(A1,2)*pow(mu,2)*Qtt(i,j) - 
     144*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1))*Qtt(i,j) - 
     96*exp(A1*mu*h(i,j))*Qtt(i,j) + 
     Q22(i,j)*(9*pow(A1,2)*(4 + pow(mu,2) - 
           16*exp((-2*mu*h(i,j))/(3.*A1))) + 
        6*(4 + pow(mu,2) - 16*exp(A1*mu*h(i,j))) + 
        (pow(A1,2)*(72 - 6*pow(mu,2) - 
              144*exp((-2*mu*h(i,j))/(3.*A1))) - 
           4*(-12 + pow(mu,2) + 24*exp(A1*mu*h(i,j))))*Qtt(i,j)) + 
     Q11(i,j)*(24 + 36*pow(A1,2) + 6*pow(mu,2) + 
        9*pow(A1,2)*pow(mu,2) - 
        144*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1)) - 
        96*exp(A1*mu*h(i,j)) + 
        (pow(A1,2)*(72 - 6*pow(mu,2) - 
              144*exp((-2*mu*h(i,j))/(3.*A1))) - 
           4*(-12 + pow(mu,2) + 24*exp(A1*mu*h(i,j))))*Qtt(i,j) + 
        Q22(i,j)*(-24 + 10*pow(mu,2) + 
           pow(A1,2)*(-36 + 15*pow(mu,2) - 
              144*exp((-2*mu*h(i,j))/(3.*A1))) - 96*exp(A1*mu*h(i,j)) + 
           (-144*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1)) - 
              96*exp(A1*mu*h(i,j)))*Qtt(i,j))))/
   ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j)))
;
elem[1][1]+=
0
;
elem[1][2]+=
(2*pow(Qtt.diff(d01,i,j),2))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
  (2*Qtt.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
  pow(1 + Qtt(i,j),2)/pow(1 + Q11(i,j),2)
;
elem[1][3]+=
0
;
elem[1][4]+=
pow(1 + Qtt(i,j),2)/pow(1 + Q22(i,j),2)
;
elem[1][5]+=
(-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)))/(-12 + pow(mu,2))
;
elem[1][6]+=
(-2*mu*exp((-2*mu*h(i,j))/(3.*A1))*
    ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
       exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) + 
      24*A1*(-1 + exp((2*mu*h(i,j))/(3.*A1) + A1*mu*h(i,j)))*
       pow(1 + Qtt(i,j),2)))/((2 + 3*pow(A1,2))*(-12 + pow(mu,2)))
;
elem[2][0]+=
(4*pow(h.diff(d01,i,j),2)*pow(mu,2))/
   (pow(L,2)*(-12 + pow(mu,2))) - 
  (3*pow(Q11.diff(d01,i,j),2))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
  (2*Q11.diff(d02,i,j))/(pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  pow(Q22.diff(d01,i,j),2)/
   ((-12 + pow(mu,2))*pow(L + L*Q22(i,j),2)) - 
  (2*pow(Qtt.diff(d01,i,j),2))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Qtt(i,j),2)) + 
  (2*Qtt.diff(d01,i,j)*(1 + Q11(i,j))*Qr1(i,j))/
   (L*pow(1 + Qtt(i,j),2)) - 
  (2*exp((-2*mu*h(i,j))/(3.*A1))*(1 + Q11(i,j))*
     ((2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) + 
       12*(3*pow(A1,2) + 2*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*
        pow(1 + Qtt(i,j),2)))/
   ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*pow(1 + Qtt(i,j),2))
;
elem[2][1]+=
0
;
elem[2][2]+=
(6*pow(Q11.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) - 
  (2*Q11.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
  (2*Qtt.diff(d01,i,j)*Qr1(i,j))/(L + L*Qtt(i,j)) + 
  (exp((-2*mu*h(i,j))/(3.*A1))*
     (2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) - 
       (1 + Qtt(i,j))*(3*pow(A1,2)*
           (24 + (-12 + pow(mu,2))*exp((2*mu*h(i,j))/(3.*A1))) + 
          2*exp((2*mu*h(i,j))/(3.*A1))*
           (-12 + pow(mu,2) + 24*exp(A1*mu*h(i,j))) + 
          24*(3*pow(A1,2) + 2*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*
           Qtt(i,j))))/
   ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[2][3]+=
-((Qtt.diff(d01,i,j)*(2 + 2*Q11(i,j)))/(L + L*Qtt(i,j)))
;
elem[2][4]+=
(-2*pow(Q22.diff(d01,i,j),2)*L*(1 + Qtt(i,j)))/
  ((-12 + pow(mu,2))*pow(L + L*Q22(i,j),3))
;
elem[2][5]+=
(4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j))*(1 + Q11(i,j)))/
  ((-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[2][6]+=
(2*mu*exp((-2*mu*h(i,j))/(3.*A1))*(1 + Q11(i,j))*
    ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
       exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) - 
      24*A1*(-1 + exp((2*mu*h(i,j))/(3.*A1) + A1*mu*h(i,j)))*
       pow(1 + Qtt(i,j),2)))/
  ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[3][0]+=
(2*Q22.diff(d02,i,j))/(pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (2*pow(Q22.diff(d01,i,j),2))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))) - 
  (2*exp((-2*mu*h(i,j))/(3.*A1))*(1 + Q22(i,j))*
     ((2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) + 
       12*(3*pow(A1,2) + 2*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*
        pow(1 + Qtt(i,j),2)))/
   ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*pow(1 + Qtt(i,j),2))
;
elem[3][1]+=
0
;
elem[3][2]+=
(-2*Q22.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
  (2*pow(Q22.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[3][3]+=
0
;
elem[3][4]+=
(2*pow(Q22.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
  (exp((-2*mu*h(i,j))/(3.*A1))*
     (2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) - 
       (1 + Qtt(i,j))*(3*pow(A1,2)*
           (24 + (-12 + pow(mu,2))*exp((2*mu*h(i,j))/(3.*A1))) + 
          2*exp((2*mu*h(i,j))/(3.*A1))*
           (-12 + pow(mu,2) + 24*exp(A1*mu*h(i,j))) + 
          24*(3*pow(A1,2) + 2*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*
           Qtt(i,j))))/
   ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[3][5]+=
(4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j))*(1 + Q22(i,j)))/
  ((-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[3][6]+=
(2*mu*exp((-2*mu*h(i,j))/(3.*A1))*(1 + Q22(i,j))*
    ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
       exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) - 
      24*A1*(-1 + exp((2*mu*h(i,j))/(3.*A1) + A1*mu*h(i,j)))*
       pow(1 + Qtt(i,j),2)))/
  ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[4][0]+=
-((Q11.diff(d01,i,j)*Q11.diff(d10,i,j))/
     (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3))) + 
  (2*Qr1.diff(d02,i,j))/(pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (4*h.diff(d01,i,j)*h.diff(d10,i,j)*pow(mu,2))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (4*h.diff(d01,i,j)*pow(mu,2)*h(i,j))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (Q22.diff(d01,i,j)*Q22.diff(d10,i,j))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
  (pow(Q11.diff(d01,i,j),2)*Qr1(i,j))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) - 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*Qr1(i,j))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (pow(Q22.diff(d01,i,j),2)*Qr1(i,j))/
   ((-12 + pow(mu,2))*(1 + Q11(i,j))*pow(L + L*Q22(i,j),2)) + 
  (Qrr.diff(d10,i,j)*Qr1(i,j))/(2.*pow(1 + Qtt(i,j),2)) + 
  (Q22.diff(d01,i,j)*(2 - Q22(i,j) + 4*Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
  Qtt.diff(d10,i,j)*((-2*Qtt.diff(d01,i,j))/
      (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2)) + 
     Qr1(i,j)/(2.*pow(1 + Qtt(i,j),2))) + 
  Qtt.diff(d01,i,j)*((2*Qrr.diff(d10,i,j))/
      (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2)) - 
     pow(Qr1(i,j),2)/(L*pow(1 + Qtt(i,j),2))) + 
  Q11.diff(d01,i,j)*((2*Qr1.diff(d01,i,j))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
     (6 + Q11(i,j) + 4*Qtt(i,j))/
      (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3))) - 
  (Qr1(i,j)*(3 + 4*Qtt(i,j) + 2*pow(Qtt(i,j),2) + 
       Q22(i,j)*(2 + 2*Qtt(i,j) + pow(Qtt(i,j),2)) + 
       Q11(i,j)*(2 + Q22(i,j) + 2*Qtt(i,j) + pow(Qtt(i,j),2))))/
   ((1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2))
;
elem[4][1]+=
0
;
elem[4][2]+=
(8*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)))/
   (L*pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)) - 
  (2*Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j))/
   (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Qtt(i,j))) + 
  (3*Q11.diff(d01,i,j)*Q11.diff(d10,i,j)*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),4)) - 
  (2*Qr1.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
  (4*h.diff(d01,i,j)*h.diff(d10,i,j)*pow(mu,2)*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
  (4*h.diff(d01,i,j)*pow(mu,2)*h(i,j)*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
  (Q22.diff(d01,i,j)*Q22.diff(d10,i,j)*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),2)) - 
  (3*pow(Q11.diff(d01,i,j),2)*Qr1(i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),4)) + 
  (Qr1(i,j)*(1 + Qtt(i,j)))/pow(1 + Q11(i,j),2) + 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*Qr1(i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
  (pow(Q22.diff(d01,i,j),2)*Qr1(i,j)*(1 + Qtt(i,j)))/
   ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*pow(L + L*Q22(i,j),2)) + 
  (Q22.diff(d01,i,j)*(Q22(i,j) - 2*Qtt(i,j))*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),2)) + 
  Qtt.diff(d01,i,j)*((-4*Q11.diff(d10,i,j))/
      (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) + 
     (4*Qr1.diff(d01,i,j))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
     (2*(3 + Q11(i,j)))/(L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) + 
     (4*Q11.diff(d01,i,j)*Qr1(i,j))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) + 
     (2*Qrr.diff(d10,i,j))/
      (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Qtt(i,j)))) + 
  Q11.diff(d01,i,j)*((-4*Qr1.diff(d01,i,j)*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) - 
     ((1 + Qtt(i,j))*(11 + 2*Q11(i,j) + 6*Qtt(i,j)))/
      (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),4)))
;
elem[4][3]+=
-(Qr1.diff(d01,i,j)/L) - Qrr.diff(d10,i,j)/(2.*(1 + Qtt(i,j))) - 
  Qtt.diff(d10,i,j)/(2.*(1 + Qtt(i,j))) + 
  (pow(Q11.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) - 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (pow(Q22.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
   ((-12 + pow(mu,2))*(1 + Q11(i,j))*pow(L + L*Q22(i,j),2)) + 
  Qtt.diff(d01,i,j)*((-2*Q11.diff(d01,i,j))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
     (2*Qr1(i,j))/(L + L*Qtt(i,j))) + 
  (-24 + 6*pow(mu,2) + 12*Qtt(i,j) + 3*pow(mu,2)*Qtt(i,j) + 
     24*pow(Qtt(i,j),2) - 2*pow(mu,2)*pow(Qtt(i,j),2) + 
     Q22(i,j)*(-36 + 7*pow(mu,2) + (-12 + 5*pow(mu,2))*Qtt(i,j) - 
        (-12 + pow(mu,2))*pow(Qtt(i,j),2)) + 
     Q11(i,j)*(-36 + 7*pow(mu,2) + (-12 + 5*pow(mu,2))*Qtt(i,j) - 
        (-12 + pow(mu,2))*pow(Qtt(i,j),2) + 
        Q22(i,j)*(8*(-6 + pow(mu,2)) + (-36 + 7*pow(mu,2))*Qtt(i,j))))/
   ((-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j)))
;
elem[4][4]+=
(-2*Q22.diff(d01,i,j)*Q22.diff(d10,i,j)*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),3)) + 
  (Qr1(i,j)*(1 + Qtt(i,j)))/pow(1 + Q22(i,j),2) + 
  (2*pow(Q22.diff(d01,i,j),2)*L*Qr1(i,j)*(1 + Qtt(i,j)))/
   ((-12 + pow(mu,2))*(1 + Q11(i,j))*pow(L + L*Q22(i,j),3)) + 
  (Q22.diff(d01,i,j)*(-1 + Q22(i,j) - 4*Qtt(i,j))*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),3))
;
elem[4][5]+=
(-8*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)))/
  (L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j)))
;
elem[4][6]+=
(-8*a0.diff(d01,i,j)*B1*pow(mu,3)*a0(i,j)*exp(B1*mu*h(i,j)))/
   (L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))) + 
  (4*h.diff(d01,i,j)*pow(mu,2)*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[5][0]+=
-((a0.diff(d01,i,j)*Q11.diff(d01,i,j))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2))) + 
  (2*a0.diff(d02,i,j))/(pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (2*a0.diff(d01,i,j)*h.diff(d01,i,j)*B1*mu)/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (a0.diff(d01,i,j)*Q22.diff(d01,i,j))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))) - 
  a0(i,j)/pow(1 + Qtt(i,j),2) + 
  (Qrr.diff(d10,i,j)*a0(i,j))/(2.*pow(1 + Qtt(i,j),2)) + 
  (Qtt.diff(d10,i,j)*a0(i,j))/(2.*pow(1 + Qtt(i,j),2)) - 
  (Qtt.diff(d01,i,j)*a0(i,j)*Qr1(i,j))/(L*pow(1 + Qtt(i,j),2))
;
elem[5][1]+=
0
;
elem[5][2]+=
a0(i,j)/(2.*pow(1 + Q11(i,j),2)) - 
  (Q11.diff(d10,i,j)*a0(i,j))/(2.*pow(1 + Q11(i,j),2)) - 
  (2*a0.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
  (2*a0.diff(d01,i,j)*h.diff(d01,i,j)*B1*mu*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
  (a0.diff(d01,i,j)*Q22.diff(d01,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))) + 
  Q11.diff(d01,i,j)*((a0(i,j)*Qr1(i,j))/(2.*L*pow(1 + Q11(i,j),2)) + 
     (2*a0.diff(d01,i,j)*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)))
;
elem[5][3]+=
(-2*a0.diff(d01,i,j))/L - (h.diff(d01,i,j)*B1*mu*a0(i,j))/L - 
  (Q11.diff(d01,i,j)*a0(i,j))/(2*L + 2*L*Q11(i,j)) - 
  (Q22.diff(d01,i,j)*a0(i,j))/(2*L + 2*L*Q22(i,j)) + 
  (Qtt.diff(d01,i,j)*a0(i,j))/(L + L*Qtt(i,j))
;
elem[5][4]+=
a0(i,j)/(2.*pow(1 + Q22(i,j),2)) - 
  (Q22.diff(d10,i,j)*a0(i,j))/(2.*pow(1 + Q22(i,j),2)) + 
  Q22.diff(d01,i,j)*((a0(i,j)*Qr1(i,j))/(2.*L*pow(1 + Q22(i,j),2)) - 
     (a0.diff(d01,i,j)*(1 + Qtt(i,j)))/
      (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)))
;
elem[5][5]+=
-(Qr1.diff(d01,i,j)/L) + h.diff(d10,i,j)*B1*mu + 
  Q11.diff(d10,i,j)/(2 + 2*Q11(i,j)) + 
  Q22.diff(d10,i,j)/(2 + 2*Q22(i,j)) - 
  (h.diff(d01,i,j)*B1*mu*Qr1(i,j))/L - 
  (Q11.diff(d01,i,j)*Qr1(i,j))/(2*L + 2*L*Q11(i,j)) - 
  (Q22.diff(d01,i,j)*Qr1(i,j))/(2*L + 2*L*Q22(i,j)) - 
  Qrr.diff(d10,i,j)/(2.*(1 + Qtt(i,j))) - 
  Qtt.diff(d10,i,j)/(2.*(1 + Qtt(i,j))) + 
  (Qtt.diff(d01,i,j)*Qr1(i,j))/(L + L*Qtt(i,j)) + 
  (Q22(i,j) + Q11(i,j)*(1 + 2*Q22(i,j) - Qtt(i,j)) - 2*Qtt(i,j) - 
     Q22(i,j)*Qtt(i,j) + 2*B1*mu*h(i,j)*(1 + Q11(i,j))*(1 + Q22(i,j))*
      (1 + Qtt(i,j)))/(2.*(1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j)))
;
elem[5][6]+=
B1*mu*a0(i,j)
;
elem[6][0]+=
-((h.diff(d01,i,j)*Q11.diff(d01,i,j))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2))) + 
  (2*h.diff(d02,i,j))/(pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (h.diff(d01,i,j)*Q22.diff(d01,i,j))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))) + 
  (exp((-2*mu*h(i,j))/(3.*A1))*
     (-((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
          exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1))) + 
       24*A1*(-1 + exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*
        pow(1 + Qtt(i,j),2)))/
   ((2 + 3*pow(A1,2))*mu*(-12 + pow(mu,2))*pow(1 + Qtt(i,j),2))
;
elem[6][1]+=
0
;
elem[6][2]+=
(-2*h.diff(d01,i,j)*Qtt.diff(d01,i,j))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
  (2*h.diff(d01,i,j)*Q11.diff(d01,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) - 
  (2*h.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
  (h.diff(d01,i,j)*Q22.diff(d01,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[6][3]+=
-(h.diff(d01,i,j)/L)
;
elem[6][4]+=
-((h.diff(d01,i,j)*Q22.diff(d01,i,j)*(1 + Qtt(i,j)))/
    (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)))
;
elem[6][5]+=
(2*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)))/((-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[6][6]+=
(exp((-2*mu*h(i,j))/(3.*A1))*((2 + 3*pow(A1,2))*pow(B1,2)*pow(mu,2)*
       pow(a0(i,j),2)*exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) + 
      (1 + Qtt(i,j))*(16 + (2 + 3*pow(A1,2))*(-12 + pow(mu,2))*
          exp((2*mu*h(i,j))/(3.*A1)) + 
         24*pow(A1,2)*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)) + 
         8*(2 + 3*pow(A1,2)*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*Qtt(i,j)\
)))/((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
}
