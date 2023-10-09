
if(rgrid.IsRightBoundary(i)){
elem[0][0]=0;
elem[0][1]=0;
elem[0][2]=0;
elem[0][3]=0;
elem[0][4]=0;
elem[0][5]=0;
elem[0][6]=0;
elem[1][0]=0;
elem[1][1]=0;
elem[1][2]=0;
elem[1][3]=0;
elem[1][4]=0;
elem[1][5]=0;
elem[1][6]=0;
elem[2][0]=0;
elem[2][1]=0;
elem[2][2]=0;
elem[2][3]=0;
elem[2][4]=0;
elem[2][5]=0;
elem[2][6]=0;
elem[3][0]=0;
elem[3][1]=0;
elem[3][2]=0;
elem[3][3]=0;
elem[3][4]=0;
elem[3][5]=0;
elem[3][6]=0;
elem[4][0]=0;
elem[4][1]=0;
elem[4][2]=0;
elem[4][3]=0;
elem[4][4]=0;
elem[4][5]=0;
elem[4][6]=0;
elem[5][0]=0;
elem[5][1]=0;
elem[5][2]=0;
elem[5][3]=0;
elem[5][4]=0;
elem[5][5]=0;
elem[5][6]=0;
elem[6][0]=0;
elem[6][1]=0;
elem[6][2]=0;
elem[6][3]=0;
elem[6][4]=0;
elem[6][5]=0;
elem[6][6]=0;
if(i==i1 && j==j1){
elem[0][0]=1;
elem[1][1]=1;
elem[2][2]=1;
elem[3][3]=1;
elem[4][4]=1;
elem[5][5]=1;
elem[6][6]=1;
}
if(i+1==i1 && j==j1){
elem[0][0]=-1;
elem[1][1]=-1;
elem[2][2]=-1;
elem[3][3]=-1;
elem[4][4]=-1;
elem[5][5]=-1;
elem[6][6]=-1;
}
}
else if(rgrid.IsLeftBoundary(i)){
elem[0][0]=0;
elem[0][1]=0;
elem[0][2]=0;
elem[0][3]=0;
elem[0][4]=0;
elem[0][5]=0;
elem[0][6]=0;
elem[1][0]=0;
elem[1][1]=0;
elem[1][2]=0;
elem[1][3]=0;
elem[1][4]=0;
elem[1][5]=0;
elem[1][6]=0;
elem[2][0]=0;
elem[2][1]=0;
elem[2][2]=0;
elem[2][3]=0;
elem[2][4]=0;
elem[2][5]=0;
elem[2][6]=0;
elem[3][0]=0;
elem[3][1]=0;
elem[3][2]=0;
elem[3][3]=0;
elem[3][4]=0;
elem[3][5]=0;
elem[3][6]=0;
elem[4][0]=0;
elem[4][1]=0;
elem[4][2]=0;
elem[4][3]=0;
elem[4][4]=0;
elem[4][5]=0;
elem[4][6]=0;
elem[5][0]=0;
elem[5][1]=0;
elem[5][2]=0;
elem[5][3]=0;
elem[5][4]=0;
elem[5][5]=0;
elem[5][6]=0;
elem[6][0]=0;
elem[6][1]=0;
elem[6][2]=0;
elem[6][3]=0;
elem[6][4]=0;
elem[6][5]=0;
elem[6][6]=0;
if(j==j1){
elem[0][0]=D[-1+d10](i-1,j,i1,j1)-D[-1+d10](i,j,i1,j1);
elem[1][1]=D[-1+d10](i-1,j,i1,j1)-D[-1+d10](i,j,i1,j1);
elem[2][2]=D[-1+d10](i-1,j,i1,j1)-D[-1+d10](i,j,i1,j1);
elem[3][3]=D[-1+d10](i-1,j,i1,j1)-D[-1+d10](i,j,i1,j1);
elem[4][4]=D[-1+d10](i-1,j,i1,j1)-D[-1+d10](i,j,i1,j1);
elem[5][5]=D[-1+d10](i-1,j,i1,j1)-D[-1+d10](i,j,i1,j1);
elem[6][6]=D[-1+d10](i-1,j,i1,j1)-D[-1+d10](i,j,i1,j1);
}
}
else{
elem[0][0]=
(-2*rgrid[i]*Qr1(i,j)*D[-1 + d11](i,j,i1,j1))/L
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
(-2*rgrid[i]*Qr1(i,j)*D[-1 + d11](i,j,i1,j1))/L
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
(-2*rgrid[i]*Qr1(i,j)*D[-1 + d11](i,j,i1,j1))/L
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
(-2*rgrid[i]*Qr1(i,j)*D[-1 + d11](i,j,i1,j1))/L
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
(-2*rgrid[i]*Qr1(i,j)*D[-1 + d11](i,j,i1,j1))/L
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
(-2*rgrid[i]*Qr1(i,j)*D[-1 + d11](i,j,i1,j1))/L
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
(-2*rgrid[i]*Qr1(i,j)*D[-1 + d11](i,j,i1,j1))/L
;
if(j==j1){
elem[0][0]+=
((-2*Qtt.diff(d10,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j)) + 
     (2*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qtt(i,j)) + 
     (-(pow(-1 + rgrid[i],2)*rgrid[i]*
           pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
             pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
           (1 + rgrid[i]*Qtt(i,j))) + 
        pow(rgrid[i],2)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
         (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
        2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (3*Qrr(i,j) + Qtt(i,j) + 2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           Q22(i,j)*(-1 + rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
           Q11(i,j)*(-1 + rgrid[i]*Q22(i,j)*
               (-2 + rgrid[i]*Qrr(i,j) - rgrid[i]*Qtt(i,j)) + 
              rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j)))))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))))*
   D[-1 + d10](i,j,i1,j1) + D[-1 + d20](i,j,i1,j1)
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
((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i])/
     (4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - pow(mu,2)*pow(rgrid[i],3)) \
- (4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*rgrid[i])/
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3)) + 
    (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))))*D[-1 + d10](i,j,i1,j1)
;
elem[0][6]+=
0
;
elem[1][0]+=
((Qtt.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     pow(1 + rgrid[i]*Qtt(i,j),2) - 
    (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
    ((1 + rgrid[i]*Qrr(i,j))*((-8 - 
            (4 + pow(mu,2))*pow(rgrid[i],3) + 
            2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
         (4 + 2*(4 + pow(mu,2))*pow(rgrid[i],3) - 
            3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j)))/
     ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)))*
  D[-1 + d10](i,j,i1,j1)
;
elem[1][1]+=
((-3*Qrr.diff(d10,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) + 
     (3*Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qrr(i,j)) + 
     (pow(-1 + rgrid[i],2)*rgrid[i]*
         pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
           pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
         (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
         (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
        pow(rgrid[i],2)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
         (1 + rgrid[i]*Qrr(i,j))*
         (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
        2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (6*Qrr(i,j) + 3*rgrid[i]*pow(Qrr(i,j),2) - Qtt(i,j) + 
           4*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           2*pow(rgrid[i],2)*pow(Qrr(i,j),2)*Qtt(i,j) + 
           Q22(i,j)*(-1 - 2*rgrid[i]*Qtt(i,j) + 
              2*rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
              pow(rgrid[i],2)*pow(Qrr(i,j),2)*(2 + rgrid[i]*Qtt(i,j))\
) + Q11(i,j)*(-1 - 2*rgrid[i]*Qtt(i,j) + 
              rgrid[i]*Q22(i,j)*
               (-2 + 2*rgrid[i]*Qrr(i,j) + 
                 pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
                 3*rgrid[i]*Qtt(i,j)) + 
              2*rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
              pow(rgrid[i],2)*pow(Qrr(i,j),2)*(2 + rgrid[i]*Qtt(i,j)))\
))/(2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j))))*D[-1 + d10](i,j,i1,j1) + 
  D[-1 + d20](i,j,i1,j1)
;
elem[1][2]+=
((Q11.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     pow(1 + rgrid[i]*Q11(i,j),2) - 
    (2*Qr1.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L + L*rgrid[i]*Q11(i,j)) - 
    (Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q11(i,j),2)) - 
    ((Q11(i,j) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
     pow(1 + rgrid[i]*Q11(i,j),2))*D[-1 + d10](i,j,i1,j1)
;
elem[1][3]+=
-2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
    pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
  D[-1 + d10](i,j,i1,j1)
;
elem[1][4]+=
((Q22.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     pow(1 + rgrid[i]*Q22(i,j),2) - 
    (Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),2)) - 
    ((Q22(i,j) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
     pow(1 + rgrid[i]*Q22(i,j),2))*D[-1 + d10](i,j,i1,j1)
;
elem[1][5]+=
((-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (1 + rgrid[i]*Qrr(i,j)))/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
    (4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
    (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))))*
  D[-1 + d10](i,j,i1,j1)
;
elem[1][6]+=
(4*pow(mu,2)*h(i,j)*(1 + rgrid[i]*Qrr(i,j)) + 
    4*h.diff(d10,i,j)*pow(mu,2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)) - 
    (4*h.diff(d01,i,j)*pow(mu,2)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/L)*D[-1 + d10](i,j,i1,j1)
;
elem[2][0]+=
0
;
elem[2][1]+=
0
;
elem[2][2]+=
((-2*Q11.diff(d10,i,j)*rgrid[i])/(1 + rgrid[i]*Q11(i,j)) + 
     (2*Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Q11(i,j)) + 
     (-(pow(-1 + rgrid[i],2)*rgrid[i]*
           pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
             pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
           (1 + rgrid[i]*Qtt(i,j))) + 
        pow(rgrid[i],2)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
         (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
        2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (3*Qrr(i,j) - Qtt(i,j) + 2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           Q22(i,j)*(-1 - 2*rgrid[i]*Qtt(i,j) + 
              rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
           Q11(i,j)*(1 - pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
              rgrid[i]*Qrr(i,j)*
               (2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)))))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))))*
   D[-1 + d10](i,j,i1,j1) + D[-1 + d20](i,j,i1,j1)
;
elem[2][3]+=
(-2*Qrr.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Q11(i,j))*
    D[-1 + d10](i,j,i1,j1))/(L + L*rgrid[i]*Qrr(i,j))
;
elem[2][4]+=
0
;
elem[2][5]+=
((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (1 + rgrid[i]*Q11(i,j)))/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
    (4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*rgrid[i]*(1 + rgrid[i]*Q11(i,j)))/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
    (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))))*
  D[-1 + d10](i,j,i1,j1)
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
((-2*Q22.diff(d10,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) + 
     (2*Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Q22(i,j)) + 
     (-(pow(-1 + rgrid[i],2)*rgrid[i]*
           pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
             pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
           (1 + rgrid[i]*Qtt(i,j))) + 
        pow(rgrid[i],2)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
         (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
        2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (3*Qrr(i,j) - Qtt(i,j) + 2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           Q22(i,j)*(1 + rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
           Q11(i,j)*(-1 - rgrid[i]*(2 + rgrid[i]*Q22(i,j))*Qtt(i,j) + 
              rgrid[i]*Qrr(i,j)*
               (2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)))))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))))*
   D[-1 + d10](i,j,i1,j1) + D[-1 + d20](i,j,i1,j1)
;
elem[3][5]+=
((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (1 + rgrid[i]*Q22(i,j)))/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
    (4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*rgrid[i]*(1 + rgrid[i]*Q22(i,j)))/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
    (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))))*
  D[-1 + d10](i,j,i1,j1)
;
elem[3][6]+=
0
;
elem[4][0]+=
((Qtt.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       pow(1 + rgrid[i]*Qtt(i,j),2)) - 
    ((-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         2*pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
     (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)))*
  D[-1 + d10](i,j,i1,j1)
;
elem[4][1]+=
(-((Qr1.diff(d10,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j))) - 
    (Qrr.diff(d01,i,j)*rgrid[i])/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qrr(i,j))) + 
    (pow(rgrid[i],2)*Qr1(i,j)*
       ((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
         pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
            pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)))/
     (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) + 
    (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (L + L*rgrid[i]*Qrr(i,j)))*D[-1 + d10](i,j,i1,j1)
;
elem[4][2]+=
((2*Qrr.diff(d01,i,j)*rgrid[i])/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) - 
    (Q11.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3)) + 
    (Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/pow(1 + rgrid[i]*Q11(i,j),2))*
  D[-1 + d10](i,j,i1,j1)
;
elem[4][3]+=
(-((Qrr.diff(d10,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j))) + 
     (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qrr(i,j)) + 
     (-3*pow(-1 + rgrid[i],2)*rgrid[i]*
         pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
           pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
         (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
         (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
        pow(rgrid[i],2)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
         (1 + rgrid[i]*Qrr(i,j))*
         (4 + rgrid[i]*Qrr(i,j) + 3*rgrid[i]*Qtt(i,j)) - 
        2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (4*Qrr(i,j) + 3*rgrid[i]*pow(Qrr(i,j),2) - Qtt(i,j) + 
           2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           2*pow(rgrid[i],2)*pow(Qrr(i,j),2)*Qtt(i,j) + 
           Q22(i,j)*(-1 + 2*rgrid[i]*Qrr(i,j) - 2*rgrid[i]*Qtt(i,j) + 
              pow(rgrid[i],2)*pow(Qrr(i,j),2)*(2 + rgrid[i]*Qtt(i,j))\
) + Q11(i,j)*(-1 + 2*rgrid[i]*Qrr(i,j) - 2*rgrid[i]*Qtt(i,j) + 
              pow(rgrid[i],2)*pow(Qrr(i,j),2)*
               (2 + rgrid[i]*Qtt(i,j)) + 
              rgrid[i]*Q22(i,j)*
               (-2 + pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
                 3*rgrid[i]*Qtt(i,j) - 
                 2*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j)))))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j))))*D[-1 + d10](i,j,i1,j1) + 
  D[-1 + d20](i,j,i1,j1)
;
elem[4][4]+=
((Q22.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       pow(1 + rgrid[i]*Q22(i,j),2)) + 
    (Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/pow(1 + rgrid[i]*Q22(i,j),2))*
  D[-1 + d10](i,j,i1,j1)
;
elem[4][5]+=
(-8*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
    (1 + rgrid[i]*Qrr(i,j))*D[-1 + d10](i,j,i1,j1))/
  (L*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
      pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
    (1 + rgrid[i]*Qtt(i,j)))
;
elem[4][6]+=
(4*h.diff(d01,i,j)*pow(mu,2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j))*
    D[-1 + d10](i,j,i1,j1))/
  (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
      pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j)))
;
elem[5][0]+=
(-(a0(i,j)*rgrid[i])/(2.*(-1 + rgrid[i])*(1 + rgrid[i]*Qtt(i,j))) - 
    (a0.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qtt(i,j)) + 
    (a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2*L + 2*L*rgrid[i]*Qtt(i,j)))*D[-1 + d10](i,j,i1,j1)
;
elem[5][1]+=
(-(a0(i,j)*rgrid[i])/(2.*(-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j))) - 
    (a0.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qrr(i,j)) + 
    (a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2*L + 2*L*rgrid[i]*Qrr(i,j)))*D[-1 + d10](i,j,i1,j1)
;
elem[5][2]+=
((a0(i,j)*rgrid[i])/(2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))) + 
    (a0.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q11(i,j)) - 
    (a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Q11(i,j))))*D[-1 + d10](i,j,i1,j1)
;
elem[5][3]+=
-((a0.diff(d01,i,j)*rgrid[i]*D[-1 + d10](i,j,i1,j1))/L)
;
elem[5][4]+=
((a0(i,j)*rgrid[i])/(2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q22(i,j))) + 
    (a0.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q22(i,j)) - 
    (a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Q22(i,j))))*D[-1 + d10](i,j,i1,j1)
;
elem[5][5]+=
(-((Qr1.diff(d01,i,j)*rgrid[i])/L) + h.diff(d10,i,j)*B1*mu*rgrid[i] + 
     (Q11.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q11(i,j)) + 
     (Q22.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q22(i,j)) - 
     (h.diff(d01,i,j)*B1*mu*pow(rgrid[i],2)*Qr1(i,j))/L - 
     (Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) - 
     (Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q22(i,j))) - 
     (Qrr.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qrr(i,j)) + 
     (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j)) - 
     (Qtt.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qtt(i,j)) + 
     (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qtt(i,j)) + 
     (4 - Q22(i,j) + 5*rgrid[i]*Q22(i,j) + Qrr(i,j) + 
        3*rgrid[i]*Qrr(i,j) + 4*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + 
        Qtt(i,j) + 3*rgrid[i]*Qtt(i,j) + 
        4*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
        2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
        2*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
        pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
        3*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
        2*B1*mu*h(i,j)*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))*
         (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
         (1 + rgrid[i]*Qtt(i,j)) + 
        Q11(i,j)*(-1 + 5*rgrid[i] + 4*pow(rgrid[i],2)*Qtt(i,j) + 
           pow(rgrid[i],2)*Qrr(i,j)*(4 + (1 + 3*rgrid[i])*Qtt(i,j)) + 
           rgrid[i]*Q22(i,j)*(-2 + 6*rgrid[i] + 
              rgrid[i]*(-1 + 5*rgrid[i])*Qtt(i,j) + 
              rgrid[i]*Qrr(i,j)*
               (-1 + 5*rgrid[i] + 4*pow(rgrid[i],2)*Qtt(i,j)))))/
      (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))*
   D[-1 + d10](i,j,i1,j1) + D[-1 + d20](i,j,i1,j1)
;
elem[5][6]+=
(a0.diff(d10,i,j)*B1*mu*rgrid[i] + 
    (B1*mu*a0(i,j)*rgrid[i])/(-1 + rgrid[i]) - 
    (a0.diff(d01,i,j)*B1*mu*pow(rgrid[i],2)*Qr1(i,j))/L)*
  D[-1 + d10](i,j,i1,j1)
;
elem[6][0]+=
(h(i,j)/(2 + 2*rgrid[i]*Qtt(i,j)) + 
    (h.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qtt(i,j)) - 
    (h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Qtt(i,j))))*D[-1 + d10](i,j,i1,j1)
;
elem[6][1]+=
(-(h(i,j)/(2 + 2*rgrid[i]*Qrr(i,j))) - 
    (h.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qrr(i,j)) + 
    (h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2*L + 2*L*rgrid[i]*Qrr(i,j)))*D[-1 + d10](i,j,i1,j1)
;
elem[6][2]+=
(h(i,j)/(2 + 2*rgrid[i]*Q11(i,j)) + 
    (h.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q11(i,j)) - 
    (h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Q11(i,j))))*D[-1 + d10](i,j,i1,j1)
;
elem[6][3]+=
-((h.diff(d01,i,j)*rgrid[i]*D[-1 + d10](i,j,i1,j1))/L)
;
elem[6][4]+=
(h(i,j)/(2 + 2*rgrid[i]*Q22(i,j)) + 
    (h.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q22(i,j)) - 
    (h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Q22(i,j))))*D[-1 + d10](i,j,i1,j1)
;
elem[6][5]+=
((2*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i])/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
    (2*a0.diff(d10,i,j)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*rgrid[i])/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
    (2*a0.diff(d01,i,j)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
       pow(rgrid[i],2)*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))))*
  D[-1 + d10](i,j,i1,j1)
;
elem[6][6]+=
(-((Qr1.diff(d01,i,j)*rgrid[i])/L) + 
     (Q11.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q11(i,j)) + 
     (Q22.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q22(i,j)) - 
     (Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) - 
     (Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q22(i,j))) - 
     (Qrr.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qrr(i,j)) + 
     (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j)) + 
     (Qtt.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qtt(i,j)) - 
     (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Qtt(i,j))) + 
     (pow(rgrid[i],2)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
         (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*
           (-Qrr(i,j) + Qtt(i,j) + 
             Q22(i,j)*(1 + rgrid[i]*(2 + rgrid[i]*Qrr(i,j))*Qtt(i,j)) + 
             Q11(i,j)*(1 + rgrid[i]*(2 + rgrid[i]*Qrr(i,j))*Qtt(i,j) + 
                rgrid[i]*Q22(i,j)*
                 (2 + 3*rgrid[i]*Qtt(i,j) + 
                   rgrid[i]*Qrr(i,j)*(1 + 2*rgrid[i]*Qtt(i,j))))))/2.)/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j))))*D[-1 + d10](i,j,i1,j1) + 
  D[-1 + d20](i,j,i1,j1)
;
}
if(i==i1){
elem[0][0]+=
((-2*Qtt.diff(d01,i,j)*rgrid[i]*
        (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qtt(i,j))) + 
     (2*Qtt.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qtt(i,j)) - 
     (2*Qr1(i,j)*(2 + rgrid[i]*Qtt(i,j)))/(L + L*rgrid[i]*Qtt(i,j)))*
   D[-1 + d01](i,j,i1,j1) + ((2 + 
       (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j)))
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
((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
       Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))) + 
    (4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))) - 
    (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))))*
  D[-1 + d01](i,j,i1,j1)
;
elem[0][6]+=
0
;
elem[1][0]+=
(-((Qtt.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
         (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
    (Qtt.diff(d01,i,j)*pow(rgrid[i],3)*pow(Qr1(i,j),2)*
       (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Qtt(i,j),2) + 
    (rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
       ((-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
            2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
         (4 + 2*(4 + pow(mu,2))*pow(rgrid[i],3) - 
            3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)))*
  D[-1 + d01](i,j,i1,j1)
;
elem[1][1]+=
(-((Qrr.diff(d01,i,j)*rgrid[i]*
          (4 + 3*(-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2) + 4*rgrid[i]*Qrr(i,j)))/
        (pow(L,2)*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Qrr(i,j)))) + 
     (3*Qrr.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qrr(i,j)) - 
     (Qr1(i,j)*(4 + rgrid[i]*Qrr(i,j) + 
          (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
           (1 + rgrid[i]*Qrr(i,j))))/(L + L*rgrid[i]*Qrr(i,j)))*
   D[-1 + d01](i,j,i1,j1) + ((2 + 
       (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j)))
;
elem[1][2]+=
(-((Q11.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
         (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q11(i,j),2))) + 
    (2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/(pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
    (Q11.diff(d01,i,j)*pow(rgrid[i],3)*pow(Qr1(i,j),2)*
       (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Q11(i,j),2) + 
    (rgrid[i]*Qr1(i,j)*(Q11(i,j) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
     (L*pow(1 + rgrid[i]*Q11(i,j),2)))*D[-1 + d01](i,j,i1,j1)
;
elem[1][3]+=
((4*Qr1.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/pow(L,2) - 
    (2*Q11.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L + L*rgrid[i]*Q11(i,j)) + 
    (2*Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/(pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
    (2*(2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2) + 
         Q11(i,j)*(rgrid[i] + 
            (-1 + rgrid[i])*pow(rgrid[i],3)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)))*
       (1 + rgrid[i]*Qrr(i,j)))/(L*rgrid[i]*(1 + rgrid[i]*Q11(i,j))))*
  D[-1 + d01](i,j,i1,j1)
;
elem[1][4]+=
(-((Q22.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
         (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),2))) + 
    (Q22.diff(d01,i,j)*pow(rgrid[i],3)*pow(Qr1(i,j),2)*
       (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Q22(i,j),2) + 
    (rgrid[i]*Qr1(i,j)*(Q22(i,j) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
     (L*pow(1 + rgrid[i]*Q22(i,j),2)))*D[-1 + d01](i,j,i1,j1)
;
elem[1][5]+=
((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
       Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
    (4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
    (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (1 + rgrid[i]*Qrr(i,j))*
       ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) - 2*(1 + rgrid[i]*Qrr(i,j))))/
     (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qtt(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[1][6]+=
((-4*pow(mu,2)*h(i,j)*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/L - 
    (4*h.diff(d10,i,j)*pow(mu,2)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/L + 
    (4*h.diff(d01,i,j)*pow(mu,2)*pow(rgrid[i],3)*pow(Qr1(i,j),2)*
       (1 + rgrid[i]*Qrr(i,j)))/pow(L,2))*D[-1 + d01](i,j,i1,j1)
;
elem[2][0]+=
(-(((-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
           2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Q11(i,j))*
         Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
       (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2))) \
+ (2*Qtt.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Qtt(i,j),2)))*
  D[-1 + d01](i,j,i1,j1)
;
elem[2][1]+=
((2*Qrr.diff(d01,i,j)*rgrid[i])/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) + 
    (2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*
       Qr1(i,j))/(pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
    (pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
       ((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
         pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
            pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) - 
    (2*Qr1.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Q11(i,j)))/
     (L + L*rgrid[i]*Qrr(i,j)))*D[-1 + d01](i,j,i1,j1)
;
elem[2][2]+=
((2*Q11.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Q11(i,j)) - 
     (2*Qr1(i,j)*(1 + rgrid[i]*Q11(i,j) - rgrid[i]*Qrr(i,j)))/
      (L + L*rgrid[i]*Q11(i,j)) - 
     (2*Q11.diff(d01,i,j)*rgrid[i]*
        (3 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 3*rgrid[i]*Qrr(i,j)))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)))*
   D[-1 + d01](i,j,i1,j1) + ((2 + 
       (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j)))
;
elem[2][3]+=
((-4*Qr1.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Q11(i,j)))/pow(L,2) - 
    (2*(1 + rgrid[i]*Q11(i,j))*
       (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)))/(L*rgrid[i]) \
+ (2*Qrr.diff(d01,i,j)*pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*
       Qr1(i,j))/(pow(L,2)*(1 + rgrid[i]*Qrr(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[2][4]+=
((2*Q22.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q22(i,j),2)) + 
    (2*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
     (L*pow(1 + rgrid[i]*Q22(i,j),2)))*D[-1 + d01](i,j,i1,j1)
;
elem[2][5]+=
((-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
       (1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
    (4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
    (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) - 2*(1 + rgrid[i]*Qrr(i,j))))/
     (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))))*
  D[-1 + d01](i,j,i1,j1)
;
elem[2][6]+=
(8*h.diff(d01,i,j)*pow(mu,2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j))*
    D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
      pow(mu,2)*pow(rgrid[i],4)))
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
((2*Q22.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Q22(i,j)) - 
     (2*(2 + rgrid[i]*Q22(i,j))*Qr1(i,j))/(L + L*rgrid[i]*Q22(i,j)) - 
     (2*Q22.diff(d01,i,j)*rgrid[i]*
        (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))))*D[-1 + d01](i,j,i1,j1) + 
  ((2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j)))
;
elem[3][5]+=
((-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
       (1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
    (4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
    (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (1 + rgrid[i]*Q22(i,j))*
       (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qtt(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[3][6]+=
0
;
elem[4][0]+=
((Qtt.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       pow(1 + rgrid[i]*Qtt(i,j),2)) - 
    (2*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/
     ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
    ((1 + rgrid[i]*Qrr(i,j))*((-1 + rgrid[i])*rgrid[i]*
          (32 + 32*rgrid[i] + 32*pow(rgrid[i],2) - 
            4*(-4 + pow(mu,2))*pow(rgrid[i],3) - 
            4*(-4 + pow(mu,2))*pow(rgrid[i],4) - 
            4*(-4 + pow(mu,2))*pow(rgrid[i],5) - 
            pow(mu,2)*(12 + pow(mu,2))*pow(rgrid[i],6) + 
            2*pow(mu,4)*pow(rgrid[i],7))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) + 2*
          ((8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
               2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
            (-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
               3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j))))/
     (2.*pow(-1 + rgrid[i],2)*
       pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(L + L*rgrid[i]*Q11(i,j))*
       pow(1 + rgrid[i]*Qtt(i,j),2)))*D[-1 + d01](i,j,i1,j1)
;
elem[4][1]+=
((2*Q11.diff(d10,i,j)*rgrid[i])/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) - 
    (2*Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) - 
    (Qrr.diff(d10,i,j)*rgrid[i])/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qrr(i,j))) + 
    (2*Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qrr(i,j))) - 
    (Qr1.diff(d01,i,j)*rgrid[i]*
       (4 + (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) + 4*rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qrr(i,j))) + 
    (Qr1.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (L + L*rgrid[i]*Qrr(i,j)) - 
    (8 + pow(rgrid[i],4)*pow(4 - 
          (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),4) + 
       pow(rgrid[i],6)*pow(Q11(i,j),2)*pow(Qr1(i,j),2)*
        ((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
          pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)) + 
       10*rgrid[i]*Qrr(i,j) - pow(rgrid[i],2)*pow(Qr1(i,j),2)*
        (8 - 5*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          6*pow(mu,2)*pow(rgrid[i],4) + 
          2*(4*rgrid[i] - (4 + pow(mu,2))*pow(rgrid[i],4) + 
             pow(mu,2)*pow(rgrid[i],5))*Qrr(i,j)) + 
       2*rgrid[i]*Q11(i,j)*(2 + 
          pow(rgrid[i],4)*pow(4 - 
             (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),4) + 
          3*rgrid[i]*Qrr(i,j) + 
          pow(rgrid[i],2)*pow(Qr1(i,j),2)*
           (-4 + 4*(4 + pow(mu,2))*pow(rgrid[i],3) - 
             5*pow(mu,2)*pow(rgrid[i],4) + 
             (-4*rgrid[i] + (4 + pow(mu,2))*pow(rgrid[i],4) - 
                pow(mu,2)*pow(rgrid[i],5))*Qrr(i,j))))/
     (2.*L*(-1 + rgrid[i])*rgrid[i]*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
       (1 + rgrid[i]*Qrr(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[4][2]+=
((-2*Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) - 
    (Q11.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3)) + 
    (2*Qr1.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
    (2*Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),3)) - 
    ((-4 + (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2) + 
         rgrid[i]*Q11(i,j)*(-1 + 
            (-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)) - 
         2*rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*rgrid[i]*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3)))*
  D[-1 + d01](i,j,i1,j1)
;
elem[4][3]+=
((2*Q11.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) \
- (Qrr.diff(d01,i,j)*rgrid[i]*
        (4 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 4*rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j))) + 
     (Qrr.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qrr(i,j)) + 
     (Qr1(i,j)*(-(pow(rgrid[i],3)*
             (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
             (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
          pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
           pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
             pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j)) - 
          (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*
           (6 + 7*rgrid[i]*Qrr(i,j) + 
             2*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
             rgrid[i]*Q11(i,j)*(4 + 3*rgrid[i]*Qrr(i,j)))))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j))))*D[-1 + d01](i,j,i1,j1) + 
  ((2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j)))
;
elem[4][4]+=
((Q22.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       pow(1 + rgrid[i]*Q22(i,j),2)) - 
    (2*Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/
     ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       pow(L + L*rgrid[i]*Q22(i,j),2)) - 
    ((Q22(i,j) + (-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       pow(1 + rgrid[i]*Q22(i,j),2)))*D[-1 + d01](i,j,i1,j1)
;
elem[4][5]+=
((-8*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (1 + rgrid[i]*Qrr(i,j)))/
     (L*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qtt(i,j))) - 
    (8*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qtt(i,j))) + 
    (16*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
       pow(rgrid[i],2)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qtt(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[4][6]+=
((4*pow(mu,2)*h(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
    (4*h.diff(d10,i,j)*pow(mu,2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
     (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
    (8*h.diff(d01,i,j)*pow(mu,2)*pow(rgrid[i],2)*Qr1(i,j)*
       (1 + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))))*
  D[-1 + d01](i,j,i1,j1)
;
elem[5][0]+=
((a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Qtt(i,j))) - 
    (a0.diff(d01,i,j)*rgrid[i]*
       (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
            pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qtt(i,j))) + 
    (a0.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2*L + 2*L*rgrid[i]*Qtt(i,j)))*D[-1 + d01](i,j,i1,j1)
;
elem[5][1]+=
((a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j))) + 
    (a0.diff(d01,i,j)*rgrid[i]*
       (1 - ((-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
            pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qrr(i,j))) + 
    (a0.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2*L + 2*L*rgrid[i]*Qrr(i,j)))*D[-1 + d01](i,j,i1,j1)
;
elem[5][2]+=
(-(a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))) - 
    (a0.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Q11(i,j))) + 
    (a0.diff(d01,i,j)*rgrid[i]*
       (-1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
            pow(Qr1(i,j),2))/2. - rgrid[i]*Qrr(i,j)))/
     ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)))*
  D[-1 + d01](i,j,i1,j1)
;
elem[5][3]+=
(-((a0.diff(d10,i,j)*rgrid[i])/L) + (a0(i,j)*rgrid[i])/(L - L*rgrid[i]) + 
    (2*a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2))*
  D[-1 + d01](i,j,i1,j1)
;
elem[5][4]+=
(-(a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q22(i,j))) - 
    (a0.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Q22(i,j))) + 
    (a0.diff(d01,i,j)*rgrid[i]*
       (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
            pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Q22(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[5][5]+=
(-((Qr1.diff(d10,i,j)*rgrid[i])/L) + 
     (2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) - 
     (h.diff(d10,i,j)*B1*mu*pow(rgrid[i],2)*Qr1(i,j))/L - 
     (Q11.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) - 
     (Q22.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q22(i,j))) + 
     (Q11.diff(d01,i,j)*rgrid[i]*
        (-1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. - rgrid[i]*Qrr(i,j)))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) \
+ (Qrr.diff(d01,i,j)*rgrid[i]*
        (1 - ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j))) + 
     (Q22.diff(d01,i,j)*rgrid[i]*
        (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))) + 
     (h.diff(d01,i,j)*B1*mu*rgrid[i]*
        (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
     (Qrr.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j)) - 
     (Qtt.diff(d01,i,j)*rgrid[i]*
        (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qtt(i,j))) + 
     (Qtt.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qtt(i,j)) + 
     (Qr1(i,j)*(2 - 6*rgrid[i] + 3*rgrid[i]*Q22(i,j) - 
          7*pow(rgrid[i],2)*Q22(i,j) + rgrid[i]*Qrr(i,j) - 
          5*pow(rgrid[i],2)*Qrr(i,j) + 
          2*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) - 
          6*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j) + rgrid[i]*Qtt(i,j) - 
          5*pow(rgrid[i],2)*Qtt(i,j) + 
          2*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
          6*pow(rgrid[i],3)*Q22(i,j)*Qtt(i,j) - 
          4*pow(rgrid[i],3)*Qrr(i,j)*Qtt(i,j) + 
          pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
          5*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
          2*B1*mu*h(i,j)*(-1 + rgrid[i])*rgrid[i]*
           (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
           (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q11(i,j)*(3 - 7*rgrid[i] + 2*rgrid[i]*Qtt(i,j) - 
             6*pow(rgrid[i],2)*Qtt(i,j) + 
             rgrid[i]*Qrr(i,j)*
              (2 - 6*rgrid[i] + 
                (rgrid[i] - 5*pow(rgrid[i],2))*Qtt(i,j)) + 
             rgrid[i]*Q22(i,j)*
              (4 - 8*rgrid[i] + (3 - 7*rgrid[i])*rgrid[i]*Qtt(i,j) + 
                rgrid[i]*Qrr(i,j)*
                 (3 - 7*rgrid[i] + 2*(1 - 3*rgrid[i])*rgrid[i]*Qtt(i,j))))\
))/(2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))*
   D[-1 + d01](i,j,i1,j1) + ((2 + 
       (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j)))
;
elem[5][6]+=
(-((a0.diff(d10,i,j)*B1*mu*pow(rgrid[i],2)*Qr1(i,j))/L) + 
    (B1*mu*a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/(L - L*rgrid[i]) + 
    (a0.diff(d01,i,j)*B1*mu*rgrid[i]*
       (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))))*
  D[-1 + d01](i,j,i1,j1)
;
elem[6][0]+=
((h.diff(d01,i,j)*rgrid[i]*(1 + 
         ((-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
            pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qtt(i,j))) - 
    (h(i,j)*rgrid[i]*Qr1(i,j))/(2.*(L + L*rgrid[i]*Qtt(i,j))) - 
    (h.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Qtt(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[6][1]+=
((h.diff(d01,i,j)*rgrid[i]*(1 - 
         ((-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
            pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qrr(i,j))) + 
    (h(i,j)*rgrid[i]*Qr1(i,j))/(2*L + 2*L*rgrid[i]*Qrr(i,j)) + 
    (h.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2*L + 2*L*rgrid[i]*Qrr(i,j)))*D[-1 + d01](i,j,i1,j1)
;
elem[6][2]+=
(-(h(i,j)*rgrid[i]*Qr1(i,j))/(2.*(L + L*rgrid[i]*Q11(i,j))) - 
    (h.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Q11(i,j))) + 
    (h.diff(d01,i,j)*rgrid[i]*
       (-1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
            pow(Qr1(i,j),2))/2. - rgrid[i]*Qrr(i,j)))/
     ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
         pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)))*
  D[-1 + d01](i,j,i1,j1)
;
elem[6][3]+=
(-(h(i,j)/L) - (h.diff(d10,i,j)*rgrid[i])/L + 
    (2*h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2))*
  D[-1 + d01](i,j,i1,j1)
;
elem[6][4]+=
(-(h(i,j)*rgrid[i]*Qr1(i,j))/(2.*(L + L*rgrid[i]*Q22(i,j))) - 
    (h.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
     (2.*(L + L*rgrid[i]*Q22(i,j))) + 
    (h.diff(d01,i,j)*rgrid[i]*
       (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
            pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*(-1 + rgrid[i])*
       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Q22(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[6][5]+=
((-2*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
    (2*a0.diff(d10,i,j)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2)*Qr1(i,j))/
     (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
    (2*a0.diff(d01,i,j)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
       (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
     (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
         pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
       (1 + rgrid[i]*Qtt(i,j))))*D[-1 + d01](i,j,i1,j1)
;
elem[6][6]+=
(-((Qr1.diff(d10,i,j)*rgrid[i])/L) + 
     (2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) - 
     (Q11.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) - 
     (Q22.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q22(i,j))) + 
     (Q11.diff(d01,i,j)*rgrid[i]*
        (-1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. - rgrid[i]*Qrr(i,j)))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) \
+ (Qrr.diff(d01,i,j)*rgrid[i]*
        (1 - ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j))) + 
     (Q22.diff(d01,i,j)*rgrid[i]*
        (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))) + 
     (Qrr.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j)) + 
     (Qtt.diff(d01,i,j)*rgrid[i]*
        (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qtt(i,j))) - 
     (Qtt.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Qtt(i,j))) - 
     (Qr1(i,j)*(pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
           (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
          ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             (2 + rgrid[i]*Qrr(i,j) + 3*rgrid[i]*Qtt(i,j) + 
               2*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
               rgrid[i]*Q22(i,j)*
                (3 + 4*rgrid[i]*Qtt(i,j) + 
                  rgrid[i]*Qrr(i,j)*(2 + 3*rgrid[i]*Qtt(i,j))) + 
               rgrid[i]*Q11(i,j)*
                (3 + 4*rgrid[i]*Qtt(i,j) + 
                  rgrid[i]*Qrr(i,j)*(2 + 3*rgrid[i]*Qtt(i,j)) + 
                  rgrid[i]*Q22(i,j)*
                   (4 + 5*rgrid[i]*Qtt(i,j) + 
                     rgrid[i]*Qrr(i,j)*(3 + 4*rgrid[i]*Qtt(i,j))))))/2.))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j))))*D[-1 + d01](i,j,i1,j1) + 
  ((2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j)))
;
}
if(j==j1 && i==i1){
elem[0][0]+=
(pow(Qtt.diff(d10,i,j),2)*pow(rgrid[i],2))/
   pow(1 + rgrid[i]*Qtt(i,j),2) + 
  (2*Qtt.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (L*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],2)*
     (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
  Qtt.diff(d10,i,j)*((-8 + 5*(4 + pow(mu,2))*pow(rgrid[i],3) - 
        6*pow(mu,2)*pow(rgrid[i],4) + 
        (8*rgrid[i] + (4 + pow(mu,2))*pow(rgrid[i],4) - 
           2*pow(mu,2)*pow(rgrid[i],5))*Qrr(i,j))/
      (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
     (2*Qtt.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((1 + rgrid[i]*Qtt(i,j))*((2 + 3*pow(A1,2))*
           exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],8)*
           pow(-12 + pow(mu,2)*(-3 + 4*rgrid[i]),2)*
           (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
           (Qrr(i,j) - Qtt(i,j)) - 
          (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],7)*pow(-12 + pow(mu,2)*(-3 + 4*rgrid[i]),2)*
           (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
           (1 + rgrid[i]*Qtt(i,j)) + 
          (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(-1 + rgrid[i],3)*pow(rgrid[i],3)*
           pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3),3)*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
           (3 + 2*rgrid[i]*Qtt(i,j)) - 
          2*(2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*
           (-8 - 12*pow(rgrid[i],5)*pow(Qr1(i,j),2) - 
             3*pow(mu,2)*pow(rgrid[i],5)*pow(Qr1(i,j),2) + 
             4*pow(mu,2)*pow(rgrid[i],6)*pow(Qr1(i,j),2) - 
             7*rgrid[i]*Qrr(i,j) - 4*rgrid[i]*Qtt(i,j) - 
             12*pow(rgrid[i],6)*pow(Qr1(i,j),2)*Qtt(i,j) - 
             3*pow(mu,2)*pow(rgrid[i],6)*pow(Qr1(i,j),2)*Qtt(i,j) + 
             4*pow(mu,2)*pow(rgrid[i],7)*pow(Qr1(i,j),2)*Qtt(i,j) - 
             4*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
             rgrid[i]*Q22(i,j)*
              (-5 - 2*rgrid[i]*Qtt(i,j) + 
                pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)) - 
                2*rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
             rgrid[i]*Q11(i,j)*
              (-5 - 4*rgrid[i]*Qrr(i,j) - 2*rgrid[i]*Qtt(i,j) - 
                2*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
                pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)) + 
                rgrid[i]*Q22(i,j)*
                 (-2 - rgrid[i]*Qrr(i,j) + 
                   pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j))))) - 
          (-1 + rgrid[i])*rgrid[i]*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*
           (288*pow(A1,2) + 192*
              exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) - 
             48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3) - 
             72*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3) - 
             12*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3) - 
             18*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3) - 
             16*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4) - 
             24*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4) + 
             288*pow(A1,2)*rgrid[i]*Q22(i,j) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Q22(i,j) + 
             48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Q22(i,j) + 72*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Q22(i,j) + 12*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Q22(i,j) + 18*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Q22(i,j) - 48*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j) - 72*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j) + 4*(2 + 3*pow(A1,2))*pow(mu,2)*
              pow(a0(i,j),2)*
              exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
              (1 + rgrid[i]*Q22(i,j)) + 
             288*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Qrr(i,j) - 
             264*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qrr(i,j) - 396*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qrr(i,j) - 66*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qrr(i,j) - 99*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qrr(i,j) + 88*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qrr(i,j) + 132*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qrr(i,j) + 288*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*
              Qrr(i,j) + 192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                  rgrid[i])/(3.*A1))*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) - 
             168*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j)*Qrr(i,j) - 
             252*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*Qrr(i,j) - 
             42*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*Qrr(i,j) - 
             63*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j)*Qrr(i,j) + 
             56*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qrr(i,j) + 
             84*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Q22(i,j)*Qrr(i,j) + 288*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Qtt(i,j) + 
             48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qtt(i,j) + 72*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qtt(i,j) + 12*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qtt(i,j) + 18*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qtt(i,j) - 48*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qtt(i,j) - 72*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qtt(i,j) + 288*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*
              Qtt(i,j) + 192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                  rgrid[i])/(3.*A1))*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
             144*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j)*Qtt(i,j) + 
             216*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*Qtt(i,j) + 
             36*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*Qtt(i,j) + 
             54*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j)*Qtt(i,j) - 
             80*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qtt(i,j) - 
             120*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Q22(i,j)*Qtt(i,j) + 
             288*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
             192*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qrr(i,j)*Qtt(i,j) - 
             288*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
             48*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
             72*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qrr(i,j)*Qtt(i,j) + 
             64*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) + 
             96*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Qrr(i,j)*Qtt(i,j) + 
             288*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*
              Qtt(i,j) + 192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                  rgrid[i])/(3.*A1))*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*
              Qtt(i,j) - 96*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             144*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             24*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             36*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
             32*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
             48*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
             rgrid[i]*Q11(i,j)*
              (rgrid[i]*Qrr(i,j)*
                 (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                      7*pow(rgrid[i],3)*
                       (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                   3*pow(A1,2)*
                    (96 + 7*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       pow(rgrid[i],3)*
                       (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                   4*rgrid[i]*
                    (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                      3*pow(A1,2)*
                       (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j)) \
+ 2*(3*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (16*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (4 + pow(mu,2) - 4*pow(mu,2)*rgrid[i])) + 
                      3*pow(A1,2)*
                       (16 - exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-4 + pow(mu,2)*(-1 + 4*rgrid[i])))) + 
                   rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                      3*pow(A1,2)*
                       (48 - exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i]))))*Qtt(i,j)\
) + rgrid[i]*Q22(i,j)*(3*rgrid[i]*Qrr(i,j)*
                    (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (32*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                      3*pow(A1,2)*
                       (32 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                      32*(3*pow(A1,2) + 
                         2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j)) + 
                   2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                      3*pow(A1,2)*
                       (48 - exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i]))) + 
                      rgrid[i]*
                       (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (60 + pow(mu,2)*(15 - 28*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                         pow(A1,2)*
                         (144 - 
                         3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-60 + pow(mu,2)*(-15 + 28*rgrid[i]))))*Qtt(i,j)\
))))) - rgrid[i]*((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],7)*pow(-12 + pow(mu,2)*(-3 + 4*rgrid[i]),2)*
           (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
           (Qrr(i,j) - Qtt(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
          (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*pow(4 - 
             (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),3)*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
           (2 + 3*rgrid[i]*Qtt(i,j) + pow(rgrid[i],2)*pow(Qtt(i,j),2)) \
- 2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*
           (72*pow(A1,2) + 48*
              exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i])) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i])) + 
             8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
             12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
             72*pow(A1,2)*rgrid[i]*Q22(i,j) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Q22(i,j) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j) + 6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j) + 9*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j) + 
             72*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Qrr(i,j) + 
             8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j) + 
             12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Qrr(i,j) + 72*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*
              Qrr(i,j) + 48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                (3.*A1))*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + 
             6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j) + 
             9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*Qrr(i,j) + 144*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Qtt(i,j) - 
             24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qtt(i,j) - 
             36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Qtt(i,j) + 10*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Qtt(i,j) + 15*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qtt(i,j) + 
             144*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
             24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qtt(i,j) - 
             36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j)*Qtt(i,j) + 
             6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qtt(i,j) + 
             9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*Qtt(i,j) + 
             144*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
             11*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j)*Qtt(i,j) + 
             (33*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                Qrr(i,j)*Qtt(i,j))/2. + 
             144*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*
              Qtt(i,j) + 96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                (3.*A1))*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
             7*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
              Qtt(i,j) + (21*pow(A1,2)*
                exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
                (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
                Qtt(i,j))/2. + 
             72*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*pow(Qtt(i,j),2) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              pow(Qtt(i,j),2) + 
             3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qtt(i,j),2) + 
             (9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                pow(Qtt(i,j),2))/2. + 
             72*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*
              pow(Qtt(i,j),2) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j)*pow(Qtt(i,j),2) + 
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
              pow(Qtt(i,j),2) + 
             (3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                Q22(i,j)*pow(Qtt(i,j),2))/2. + 
             72*pow(A1,2)*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Qrr(i,j)*pow(Qtt(i,j),2) + 
             72*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
              exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
              (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
             (rgrid[i]*Q11(i,j)*
                (144*pow(A1,2) + 
                  96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                     (3.*A1)) - 
                  48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],3) - 
                  72*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],3) - 
                  12*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],3) - 
                  18*pow(A1,2)*pow(mu,2)*
                   exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3) \
+ 288*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
                  192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                     (3.*A1))*rgrid[i]*Qtt(i,j) + 
                  48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],4)*Qtt(i,j) + 
                  72*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],4)*Qtt(i,j) + 
                  12*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],4)*Qtt(i,j) + 
                  18*pow(A1,2)*pow(mu,2)*
                   exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                   Qtt(i,j) - 
                  48*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],5)*Qtt(i,j) - 
                  72*pow(A1,2)*pow(mu,2)*
                   exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                   Qtt(i,j) + 
                  144*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                  96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                     (3.*A1))*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                  72*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
                  108*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
                  18*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
                  27*pow(A1,2)*pow(mu,2)*
                   exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                   pow(Qtt(i,j),2) - 
                  40*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],6)*pow(Qtt(i,j),2) - 
                  60*pow(A1,2)*pow(mu,2)*
                   exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
                   pow(Qtt(i,j),2) + 
                  rgrid[i]*Qrr(i,j)*
                   (6*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (8*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (8 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                     rgrid[i]*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        7*pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                        (96 + 
                        7*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ 2*pow(rgrid[i],2)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (24 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                      pow(Qtt(i,j),2)) + 
                  rgrid[i]*Q22(i,j)*
                   (8*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (6*exp(A1*mu*h(i,j)*rgrid[i]) - 
                         pow(mu,2)*pow(rgrid[i],4)) - 
                        3*pow(A1,2)*
                         (-6 + 
                         pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],4))) + 
                     2*rgrid[i]*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                        pow(rgrid[i],3)) + 
                        3*pow(A1,2)*
                         (48 - 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i]))))*
                      Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (60 + pow(mu,2)*(15 - 28*rgrid[i]))*
                        pow(rgrid[i],3)) + 
                        pow(A1,2)*
                         (144 - 
                         3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-60 + pow(mu,2)*(-15 + 28*rgrid[i]))))*
                      pow(Qtt(i,j),2) + 
                     rgrid[i]*Qrr(i,j)*
                      (4*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (12*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                         3*pow(A1,2)*
                         (12 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                        3*rgrid[i]*
                         (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (32*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                         3*pow(A1,2)*
                         (32 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ 48*(3*pow(A1,2) + 2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*pow(rgrid[i],2)*pow(Qtt(i,j),2)))))/
              2.) - (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*
           (-12 + pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2) - 12*rgrid[i]*Qrr(i,j) - 
             16*rgrid[i]*Qtt(i,j) + 
             2*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2)*Qtt(i,j) - 
             14*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
             4*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
             4*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             rgrid[i]*Q22(i,j)*
              (pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) - 
                2*(4 + 5*rgrid[i]*Qtt(i,j) + 
                   pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   rgrid[i]*Qrr(i,j)*pow(2 + rgrid[i]*Qtt(i,j),2))) + 
             rgrid[i]*Q11(i,j)*
              (pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) - 
                2*(4 + 5*rgrid[i]*Qtt(i,j) + 
                   pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   rgrid[i]*Qrr(i,j)*pow(2 + rgrid[i]*Qtt(i,j),2)) + 
                rgrid[i]*Q22(i,j)*
                 (pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) - 
                   2*(2 + 2*rgrid[i]*Qtt(i,j) + 
                      rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j)))))))))/
   (2.*(2 + 3*pow(A1,2))*pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
     pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qtt(i,j),2))
;
elem[0][1]+=
(-4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))) + 
  (2*Qtt.diff(d02,i,j)*rgrid[i])/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (2*pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],2))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (Qtt.diff(d10,i,j)*((pow(rgrid[i],3)*
          (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j)))/2. - 
       (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        (3 + 2*rgrid[i]*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q11(i,j)*(2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)))))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        pow(rgrid[i],6)*pow(-12 + pow(mu,2)*(-3 + 4*rgrid[i]),2)*
        (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (1 + rgrid[i]*Qtt(i,j)) + 
       2*(2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4),2)*(2 + rgrid[i]*Qtt(i,j))*
        (3 + 2*rgrid[i]*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q11(i,j)*(2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j))) - 
       (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        (144*pow(A1,2) + 96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
             (3.*A1)) - 192*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3) - 
          288*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3) - 
          48*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3) - 
          72*pow(A1,2)*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3) + 
          64*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4) + 
          96*pow(A1,2)*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4) + 288*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
          192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*Qtt(i,j) - 
          264*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
           Qtt(i,j) - 396*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*Qtt(i,j) - 
          66*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*Qtt(i,j) - 
          99*pow(A1,2)*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*Qtt(i,j) + 
          88*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*Qtt(i,j) + 
          132*pow(A1,2)*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*Qtt(i,j) + 
          144*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
          96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
          96*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           pow(Qtt(i,j),2) - 
          144*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
          24*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
          36*pow(A1,2)*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
          32*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
          48*pow(A1,2)*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
          rgrid[i]*Q22(i,j)*(6*
              (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 (8*exp(A1*mu*h(i,j)*rgrid[i]) + 
                   pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                3*pow(A1,2)*
                 (8 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
             rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                   7*pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                3*pow(A1,2)*
                 (96 + 7*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) + 
             2*pow(rgrid[i],2)*
              (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                   pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                3*pow(A1,2)*
                 (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
              pow(Qtt(i,j),2)) + 
          rgrid[i]*Q11(i,j)*(6*
              (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 (8*exp(A1*mu*h(i,j)*rgrid[i]) + 
                   pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))\
) + 3*pow(A1,2)*(8 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))\
)) + rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                   7*pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                3*pow(A1,2)*
                 (96 + 7*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) + 
             2*pow(rgrid[i],2)*
              (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                   pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                3*pow(A1,2)*
                 (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    pow(rgrid[i],3)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
              pow(Qtt(i,j),2) + 
             rgrid[i]*Q22(i,j)*
              (4*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    (12*exp(A1*mu*h(i,j)*rgrid[i]) + 
                      pow(rgrid[i],3)*
                       (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                   3*pow(A1,2)*
                    (12 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       pow(rgrid[i],3)*
                       (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                3*rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                    (32*exp(A1*mu*h(i,j)*rgrid[i]) + 
                      pow(rgrid[i],3)*
                       (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                   3*pow(A1,2)*
                    (32 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       pow(rgrid[i],3)*
                       (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) + 
                48*(3*pow(A1,2) + 
                   2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)))*
                 pow(rgrid[i],2)*pow(Qtt(i,j),2))))))/
   (2.*(2 + 3*pow(A1,2))*pow(rgrid[i],2)*
     pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[0][2]+=
(Qtt.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q11(i,j),2) + 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (2*Qtt.diff(d02,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (2*pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j))) + 
  ((1 + rgrid[i]*Qrr(i,j))*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4) + 
       rgrid[i]*(-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j)))/
   (pow(rgrid[i],2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),2))
;
elem[0][3]+=
(-2*Qtt.diff(d11,i,j)*rgrid[i])/L + 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))) + 
  (4*a0.diff(d01,i,j)*a0.diff(d10,i,j)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))) + 
  (2*Qtt.diff(d02,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) - 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))) - 
  (2*pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*(1 + rgrid[i]*Qtt(i,j))) - 
  (2*Qtt.diff(d01,i,j)*(2 + rgrid[i]*Qtt(i,j)))/
   (L + L*rgrid[i]*Qtt(i,j)) + 
  (Qr1(i,j)*(8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
       2*pow(mu,2)*pow(rgrid[i],4) + 
       (4*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
          3*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j)))/rgrid[i] + 
  Qtt.diff(d10,i,j)*(-((-1 + rgrid[i])*rgrid[i]*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)) + 
     (2*Qtt.diff(d01,i,j)*pow(rgrid[i],2))/(L + L*rgrid[i]*Qtt(i,j)))
;
elem[0][4]+=
(Qtt.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q22(i,j),2) + 
  ((1 + rgrid[i]*Qrr(i,j))*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4) + 
       rgrid[i]*(-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j)))/
   (pow(rgrid[i],2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q22(i,j),2))
;
elem[0][5]+=
(4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i])/
   (4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - pow(mu,2)*pow(rgrid[i],3)) - 
  (4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i])/
   (4 - (4 + pow(mu,2))*pow(rgrid[i],3) + pow(mu,2)*pow(rgrid[i],4)) \
+ (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3)))
;
elem[0][6]+=
(-2*pow(a0.diff(d10,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2))/
   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3)) + 
  (4*a0.diff(d01,i,j)*B1*pow(mu,3)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))) + 
  a0.diff(d10,i,j)*((4*B1*pow(mu,3)*a0(i,j)*
        exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2))/
      (4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
        pow(mu,2)*pow(rgrid[i],3)) + 
     (4*a0.diff(d01,i,j)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3)))) - 
  (2*pow(a0.diff(d01,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))) - 
  (2*mu*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4) + 
       24*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3)))
;
elem[1][0]+=
-((pow(Qtt.diff(d10,i,j),2)*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
     pow(1 + rgrid[i]*Qtt(i,j),3)) + 
  (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (Qrr.diff(d10,i,j)*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Qrr(i,j)))/
   (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j))*
     ((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) - 2*(1 + rgrid[i]*Qrr(i,j))))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (pow(Qtt.diff(d01,i,j),2)*L*pow(rgrid[i],4)*pow(Qr1(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Qtt(i,j),3) - 
  (Qtt.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
     (2*rgrid[i]*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) - 
       (-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          3*pow(mu,2)*pow(rgrid[i],4))*(-1 + rgrid[i]*Qtt(i,j))))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),3)) + 
  a0.diff(d10,i,j)*((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
  Qtt.diff(d10,i,j)*((2*Qtt.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j)*
        (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Qtt(i,j),3)) + 
     ((1 + rgrid[i]*Qrr(i,j))*(2*rgrid[i]*
           (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) - 
          (-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
             3*pow(mu,2)*pow(rgrid[i],4))*(-1 + rgrid[i]*Qtt(i,j))))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),3))) - 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*(1 + rgrid[i]*Qrr(i,j))*
     (-4*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*(-1 + rgrid[i])*
        pow(rgrid[i],3)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)) + 
       exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        (pow(rgrid[i],2)*(-48*(4 + pow(mu,2)) + 
             96*pow(mu,2)*rgrid[i] + 
             3*pow(4 + pow(mu,2),2)*pow(rgrid[i],3) - 
             12*pow(mu,2)*(4 + pow(mu,2))*pow(rgrid[i],4) + 
             8*pow(mu,4)*pow(rgrid[i],5)) - 
          (-96 + 48*(4 + pow(mu,2))*pow(rgrid[i],3) - 
             80*pow(mu,2)*pow(rgrid[i],4) + 
             3*pow(4 + pow(mu,2),2)*pow(rgrid[i],6) - 
             4*pow(mu,2)*(4 + pow(mu,2))*pow(rgrid[i],7) + 
             2*pow(mu,4)*pow(rgrid[i],8))*Qtt(i,j) + 
          Qrr(i,j)*(-96 - 36*(4 + pow(mu,2))*pow(rgrid[i],3) + 
             96*pow(mu,2)*pow(rgrid[i],4) + 
             6*pow(4 + pow(mu,2),2)*pow(rgrid[i],6) - 
             21*pow(mu,2)*(4 + pow(mu,2))*pow(rgrid[i],7) + 
             14*pow(mu,4)*pow(rgrid[i],8) + 
             rgrid[i]*(32 - 52*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                96*pow(mu,2)*pow(rgrid[i],4) + 
                2*pow(4 + pow(mu,2),2)*pow(rgrid[i],6) - 
                9*pow(mu,2)*(4 + pow(mu,2))*pow(rgrid[i],7) + 
                6*pow(mu,4)*pow(rgrid[i],8))*Qtt(i,j)))))/
   (2.*rgrid[i]*pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4),2)*pow(1 + rgrid[i]*Qtt(i,j),3))
;
elem[1][1]+=
(2*pow(Qr1.diff(d01,i,j),2)*pow(rgrid[i],2))/pow(L,2) + 
  2*pow(h.diff(d10,i,j),2)*pow(mu,2)*pow(rgrid[i],2) + 
  (pow(Q11.diff(d10,i,j),2)*pow(rgrid[i],2))/
   (2.*pow(1 + rgrid[i]*Q11(i,j),2)) + 
  (2*Qrr.diff(d02,i,j)*rgrid[i])/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (pow(Q22.diff(d10,i,j),2)*pow(rgrid[i],2))/
   (2.*pow(1 + rgrid[i]*Q22(i,j),2)) - 
  (4*h.diff(d01,i,j)*pow(mu,2)*h(i,j)*pow(rgrid[i],2)*Qr1(i,j))/L - 
  2*Qr1.diff(d10,i,j)*(-1 + rgrid[i])*rgrid[i]*
   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
   Qr1(i,j) + (2*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(rgrid[i],4)*
     pow(Qr1(i,j),2))/pow(L,2) + 
  (pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],4)*pow(Qr1(i,j),2))/
   (2.*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],4)*pow(Qr1(i,j),2))/
   (2.*pow(L + L*rgrid[i]*Q22(i,j),2)) + 
  h.diff(d10,i,j)*(4*pow(mu,2)*h(i,j)*rgrid[i] - 
     (4*h.diff(d01,i,j)*pow(mu,2)*pow(rgrid[i],3)*Qr1(i,j))/L) + 
  (2*Qr1.diff(d01,i,j)*(2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2) + 
       Q11(i,j)*(rgrid[i] + (-1 + rgrid[i])*pow(rgrid[i],3)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2))))/
   (L + L*rgrid[i]*Q11(i,j)) + 
  (Q22.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*
     (-2 + rgrid[i]*Q22(i,j) - 4*rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (3*Qrr.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (L*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (6*pow(Qrr.diff(d10,i,j),2)*pow(rgrid[i],2))/
   pow(2 + 2*rgrid[i]*Qrr(i,j),2) + 
  (3*pow(Qrr.diff(d01,i,j),2)*pow(rgrid[i],4)*pow(Qr1(i,j),2))/
   (2.*pow(L + L*rgrid[i]*Qrr(i,j),2)) + 
  Q11.diff(d01,i,j)*((2*Qr1.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
     (rgrid[i]*Qr1(i,j)*(-2 + rgrid[i]*Q11(i,j) - 4*rgrid[i]*Qrr(i,j)))/
      (L*pow(1 + rgrid[i]*Q11(i,j),2))) + 
  Q11.diff(d10,i,j)*((-2*Qr1.diff(d01,i,j)*pow(rgrid[i],2))/
      (L + L*rgrid[i]*Q11(i,j)) - 
     (Q11.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (L*pow(1 + rgrid[i]*Q11(i,j),2)) + 
     (2 - rgrid[i]*Q11(i,j) + 4*rgrid[i]*Qrr(i,j))/
      pow(1 + rgrid[i]*Q11(i,j),2)) + 
  Q22.diff(d10,i,j)*(-((Q22.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
        (L*pow(1 + rgrid[i]*Q22(i,j),2))) + 
     (2 - rgrid[i]*Q22(i,j) + 4*rgrid[i]*Qrr(i,j))/
      pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (pow(Qtt.diff(d10,i,j),2)*pow(rgrid[i],2))/
   (2.*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) - 4*(1 + rgrid[i]*Qrr(i,j))))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],4)*pow(Qr1(i,j),2))/
   (2.*pow(L + L*rgrid[i]*Qtt(i,j),2)) - 
  (Qtt.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*
     (8 + 4*pow(rgrid[i],3) + pow(mu,2)*pow(rgrid[i],3) - 
       2*pow(mu,2)*pow(rgrid[i],4) + 
       2*rgrid[i]*(8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
          2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
       rgrid[i]*(-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  a0.diff(d10,i,j)*((-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        pow(rgrid[i],2))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  Qtt.diff(d10,i,j)*(-((Qtt.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
        (L*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
     (8 + 4*pow(rgrid[i],3) + pow(mu,2)*pow(rgrid[i],3) - 
        2*pow(mu,2)*pow(rgrid[i],4) + 
        2*rgrid[i]*(8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
           2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
        rgrid[i]*(-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
           3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
  Qrr.diff(d10,i,j)*((-3*Qrr.diff(d01,i,j)*L*pow(rgrid[i],3)*
        Qr1(i,j))/pow(L + L*rgrid[i]*Qrr(i,j),2) + 
     (-48 + 36*pow(rgrid[i],3) + 9*pow(mu,2)*pow(rgrid[i],3) - 
        8*pow(mu,2)*pow(rgrid[i],4) - 48*rgrid[i]*Qrr(i,j) + 
        24*pow(rgrid[i],4)*Qrr(i,j) + 
        6*pow(mu,2)*pow(rgrid[i],4)*Qrr(i,j) - 
        4*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j) - 
        24*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
        12*pow(rgrid[i],5)*pow(Qrr(i,j),2) + 
        3*pow(mu,2)*pow(rgrid[i],5)*pow(Qrr(i,j),2) - 
        2*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2) - 
        40*rgrid[i]*Qtt(i,j) + 40*pow(rgrid[i],4)*Qtt(i,j) + 
        10*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
        10*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
        32*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
        32*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) + 
        8*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
        8*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) - 
        16*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) + 
        16*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) + 
        4*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) - 
        4*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*Qtt(i,j) + 
        rgrid[i]*Q22(i,j)*(-40 + 28*pow(rgrid[i],3) + 
           7*pow(mu,2)*pow(rgrid[i],3) - 
           6*pow(mu,2)*pow(rgrid[i],4) - 
           8*(4*rgrid[i] - (4 + pow(mu,2))*pow(rgrid[i],4) + 
              pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j) + 
           2*rgrid[i]*Qrr(i,j)*
            (-16 + (4 + pow(mu,2))*pow(rgrid[i],3) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j)) + 
           pow(rgrid[i],2)*pow(Qrr(i,j),2)*
            (-16 + (4 + pow(mu,2))*pow(rgrid[i],3) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j))) + 
        rgrid[i]*Q11(i,j)*(-40 + 28*pow(rgrid[i],3) + 
           7*pow(mu,2)*pow(rgrid[i],3) - 
           6*pow(mu,2)*pow(rgrid[i],4) - 32*rgrid[i]*Qtt(i,j) + 
           32*pow(rgrid[i],4)*Qtt(i,j) + 
           8*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
           8*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) + 
           2*rgrid[i]*Qrr(i,j)*
            (-16 + (4 + pow(mu,2))*pow(rgrid[i],3) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j)) + 
           pow(rgrid[i],2)*pow(Qrr(i,j),2)*
            (-16 + (4 + pow(mu,2))*pow(rgrid[i],3) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j)) + 
           rgrid[i]*Q22(i,j)*(-32 + 20*pow(rgrid[i],3) + 
              5*pow(mu,2)*pow(rgrid[i],3) - 
              4*pow(mu,2)*pow(rgrid[i],4) + 
              2*rgrid[i]*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                 2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
              pow(rgrid[i],2)*
               (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                 2*pow(mu,2)*pow(rgrid[i],4))*pow(Qrr(i,j),2) - 
              6*(4*rgrid[i] - (4 + pow(mu,2))*pow(rgrid[i],4) + 
                 pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j))))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
        (1 + rgrid[i]*Qtt(i,j)))) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        pow(rgrid[i],7)*pow(-12 + pow(mu,2)*(-3 + 4*rgrid[i]),2)*
        pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
        pow(1 + rgrid[i]*Qrr(i,j),2)*(Qrr(i,j) - Qtt(i,j))*
        (1 + rgrid[i]*Qtt(i,j)) + 
       (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        pow(rgrid[i],2)*pow(4 - 
          (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4),3)*pow(1 + rgrid[i]*Q11(i,j),2)*
        pow(1 + rgrid[i]*Q22(i,j),2)*pow(Qr1(i,j),2)*
        (4 + 7*rgrid[i]*Qrr(i,j) + 3*pow(rgrid[i],2)*pow(Qrr(i,j),2))*
        pow(1 + rgrid[i]*Qtt(i,j),2) + 
       2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (72*pow(A1,2) + 48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
             (3.*A1)) - 12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i])) - 
          18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i])) + 
          8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
          12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
          72*pow(A1,2)*rgrid[i]*Q22(i,j) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*Q22(i,j) - 
          12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j) - 
          18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j) \
+ 6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j) + 
          9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j) + 144*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
          96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*Qrr(i,j) - 
          24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qrr(i,j) - 
          36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qrr(i,j) \
+ 10*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j) + 
          15*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Qrr(i,j) + 144*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*
           Qrr(i,j) + 96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
             (3.*A1))*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) - 
          24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qrr(i,j) - 
          36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
           Q22(i,j)*Qrr(i,j) + 
          6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j) + 
          9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*Qrr(i,j) + 
          72*pow(A1,2)*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
          12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*pow(Qrr(i,j),2) - 
          18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
           pow(Qrr(i,j),2) + 
          3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qrr(i,j),2) + 
          (9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
             pow(Qrr(i,j),2))/2. + 
          72*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) - 
          12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*pow(Qrr(i,j),2) \
- 18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*pow(Qrr(i,j),2) \
+ exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
           pow(Qrr(i,j),2) + 
          (3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
             Q22(i,j)*pow(Qrr(i,j),2))/2. + 
          144*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
          96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*Qtt(i,j) - 
          12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qtt(i,j) - 
          18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qtt(i,j) \
+ 16*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qtt(i,j) + 
          24*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Qtt(i,j) + 144*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*
           Qtt(i,j) + 96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
             (3.*A1))*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
          12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qtt(i,j) - 
          18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
           Q22(i,j)*Qtt(i,j) + 
          12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qtt(i,j) + 
          18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*Qtt(i,j) + 
          288*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
          24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qrr(i,j)*Qtt(i,j) - 
          36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
           Qrr(i,j)*Qtt(i,j) + 
          23*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j)*Qtt(i,j) + 
          (69*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
             Qrr(i,j)*Qtt(i,j))/2. + 
          288*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
          192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
          24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qrr(i,j)*Qtt(i,j) \
- 36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qrr(i,j)*Qtt(i,j) \
+ 15*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
           Qtt(i,j) + (45*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
             (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
             Qtt(i,j))/2. + 144*pow(A1,2)*pow(rgrid[i],3)*
           pow(Qrr(i,j),2)*Qtt(i,j) + 
          96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) - 
          12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*pow(Qrr(i,j),2)*Qtt(i,j) \
- 18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*pow(Qrr(i,j),2)*Qtt(i,j) \
+ 9*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qrr(i,j),2)*
           Qtt(i,j) + (27*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
             (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qrr(i,j),2)*
             Qtt(i,j))/2. + 144*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*
           pow(Qrr(i,j),2)*Qtt(i,j) + 
          96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) - 
          12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
           (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*
           pow(Qrr(i,j),2)*Qtt(i,j) - 
          18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],7)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
           Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) + 
          5*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
           pow(Qrr(i,j),2)*Qtt(i,j) + 
          (15*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
             Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j))/2. + 
          72*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
          6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qtt(i,j),2) + 
          9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qtt(i,j),2) + 
          72*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) + 
          4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
           pow(Qtt(i,j),2) + 
          6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qtt(i,j),2) + 
          144*pow(A1,2)*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          9*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j)*
           pow(Qtt(i,j),2) + 
          (27*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
             Qrr(i,j)*pow(Qtt(i,j),2))/2. + 
          144*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*
           pow(Qtt(i,j),2) + 
          96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          5*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
           pow(Qtt(i,j),2) + 
          (15*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
             Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2))/2. + 
          72*pow(A1,2)*pow(rgrid[i],4)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) + 
          6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          72*pow(A1,2)*pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) \
+ 2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],8)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
           pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
           exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
           (1 + rgrid[i]*Qtt(i,j)) + 
          (rgrid[i]*Q11(i,j)*(2*
                (6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   (8*exp(A1*mu*h(i,j)*rgrid[i]) - 
                     (4 + pow(mu,2))*pow(rgrid[i],3)) - 
                  9*pow(A1,2)*
                   (-8 + (4 + pow(mu,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],3)) + 
                  12*rgrid[i]*
                   (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (4*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (-4 + pow(mu,2)*(-1 + rgrid[i]))*
                        pow(rgrid[i],3)) + 
                     3*pow(A1,2)*
                      (4 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (-4 + pow(mu,2)*(-1 + rgrid[i]))*
                        pow(rgrid[i],3)))*Qtt(i,j) + 
                  2*pow(rgrid[i],2)*
                   (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (12*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (12 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                   pow(Qtt(i,j),2)) + 
               pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                     (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                      pow(rgrid[i],3)) + 
                  3*pow(A1,2)*
                   (48 - exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],3)*
                      (-36 + pow(mu,2)*(-9 + 20*rgrid[i]))) - 
                  rgrid[i]*(-2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (96*exp(A1*mu*h(i,j)*rgrid[i]) - 
                        pow(rgrid[i],3)*
                        (12 + pow(mu,2)*(3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (-96 + 
                        exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (12 + pow(mu,2)*(3 + 4*rgrid[i]))))*Qtt(i,j) + 
                  2*pow(rgrid[i],2)*
                   (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                   pow(Qtt(i,j),2)) + 
               rgrid[i]*Qrr(i,j)*
                (6*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (16*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (4 + pow(mu,2) - 4*pow(mu,2)*rgrid[i])) + 
                     3*pow(A1,2)*
                      (16 - exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-4 + pow(mu,2)*(-1 + 4*rgrid[i])))) + 
                  3*rgrid[i]*
                   (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (64*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-28 + pow(mu,2)*(-7 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (64 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (-28 + pow(mu,2)*(-7 + 4*rgrid[i]))))*Qtt(i,j) \
+ pow(rgrid[i],2)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        5*pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (96 + 5*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                   pow(Qtt(i,j),2)) + 
               rgrid[i]*Q22(i,j)*
                (pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                   (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (60 + pow(mu,2)*(15 - 28*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                     pow(A1,2)*
                      (144 - 
                        3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-60 + pow(mu,2)*(-15 + 28*rgrid[i]))) + 
                     rgrid[i]*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                        pow(rgrid[i],3)) + 
                        3*pow(A1,2)*
                        (96 - 
                        exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (-36 + pow(mu,2)*(-9 + 20*rgrid[i]))))*
                      Qtt(i,j) + 
                     48*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
                  2*(8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (6*exp(A1*mu*h(i,j)*rgrid[i]) - 
                        pow(mu,2)*pow(rgrid[i],4)) - 
                     12*pow(A1,2)*
                      (-6 + pow(mu,2)*
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],4)) + 
                     2*rgrid[i]*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 2*rgrid[i]))) + 
                        pow(A1,2)*
                        (72 + 
                        3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 2*rgrid[i]))))*Qtt(i,j) \
+ pow(rgrid[i],2)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (24 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                      pow(Qtt(i,j),2)) + 
                  rgrid[i]*Qrr(i,j)*
                   (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                        3*pow(A1,2)*
                         (48 - 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i])))) + 
                     rgrid[i]*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (192*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (12 + pow(mu,2)*(3 - 20*rgrid[i]))*
                        pow(rgrid[i],3)) + 
                        pow(A1,2)*
                         (576 - 
                         3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 20*rgrid[i]))))*
                      Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (96 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                      pow(Qtt(i,j),2)))))/2.) + 
       (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4),2)*
        (-12 - 24*rgrid[i]*Q22(i,j) - 
          9*pow(rgrid[i],2)*pow(Q22(i,j),2) + 
          3*pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qr1(i,j),2) + 
          6*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2) + 
          3*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2) - 20*rgrid[i]*Qrr(i,j) - 
          46*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) - 
          16*pow(rgrid[i],3)*pow(Q22(i,j),2)*Qrr(i,j) + 
          6*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qr1(i,j),2)*Qrr(i,j) + 
          12*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2)*Qrr(i,j) + 
          6*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2)*Qrr(i,j) + 
          pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
          10*pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) + 
          3*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qr1(i,j),2)*pow(Qrr(i,j),2) + 
          6*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2)*pow(Qrr(i,j),2) + 
          3*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2)*pow(Qrr(i,j),2) + 
          6*pow(rgrid[i],3)*pow(Qrr(i,j),3) + 
          6*pow(rgrid[i],4)*Q22(i,j)*pow(Qrr(i,j),3) + 
          4*pow(rgrid[i],5)*pow(Q22(i,j),2)*pow(Qrr(i,j),3) - 
          24*rgrid[i]*Qtt(i,j) - 48*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
          18*pow(rgrid[i],3)*pow(Q22(i,j),2)*Qtt(i,j) + 
          6*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qr1(i,j),2)*Qtt(i,j) + 
          12*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2)*Qtt(i,j) + 
          6*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2)*Qtt(i,j) - 
          46*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
          104*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
          38*pow(rgrid[i],4)*pow(Q22(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
          12*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qr1(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
          24*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
          12*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2)*Qrr(i,j)*Qtt(i,j) - 
          10*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) - 
          44*pow(rgrid[i],4)*Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) - 
          12*pow(rgrid[i],5)*pow(Q22(i,j),2)*pow(Qrr(i,j),2)*
           Qtt(i,j) + 6*pow(rgrid[i],8)*
           (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qr1(i,j),2)*
           pow(Qrr(i,j),2)*Qtt(i,j) + 
          12*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*Qtt(i,j) + 
          6*pow(rgrid[i],10)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*Qtt(i,j) \
+ 6*pow(rgrid[i],4)*pow(Qrr(i,j),3)*Qtt(i,j) + 
          2*pow(rgrid[i],6)*pow(Q22(i,j),2)*pow(Qrr(i,j),3)*
           Qtt(i,j) - 9*pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
          18*pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) - 
          6*pow(rgrid[i],4)*pow(Q22(i,j),2)*pow(Qtt(i,j),2) + 
          3*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
          6*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
          3*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
          16*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) - 
          38*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) - 
          12*pow(rgrid[i],5)*pow(Q22(i,j),2)*Qrr(i,j)*
           pow(Qtt(i,j),2) + 
          6*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qr1(i,j),2)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          12*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          6*pow(rgrid[i],10)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2)*Qrr(i,j)*pow(Qtt(i,j),2) \
- 12*pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) - 
          pow(rgrid[i],6)*pow(Q22(i,j),2)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) + 
          3*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          6*pow(rgrid[i],10)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           Q22(i,j)*pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) \
+ 3*pow(rgrid[i],11)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           pow(Q22(i,j),2)*pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) + 
          4*pow(rgrid[i],5)*pow(Qrr(i,j),3)*pow(Qtt(i,j),2) + 
          2*pow(rgrid[i],6)*Q22(i,j)*pow(Qrr(i,j),3)*
           pow(Qtt(i,j),2) + 
          2*pow(rgrid[i],7)*pow(Q22(i,j),2)*pow(Qrr(i,j),3)*
           pow(Qtt(i,j),2) - 
          4*pow(mu,2)*pow(h(i,j),2)*pow(rgrid[i],2)*
           pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
           pow(1 + rgrid[i]*Qrr(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) + 
          pow(rgrid[i],2)*pow(Q11(i,j),2)*
           (2*pow(rgrid[i],3)*pow(Qrr(i,j),3)*
              (2 + rgrid[i]*Qtt(i,j) + 
                pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
             3*(-3 - 6*rgrid[i]*Qtt(i,j) - 
                2*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
             2*rgrid[i]*Qrr(i,j)*
              (-8 - 19*rgrid[i]*Qtt(i,j) - 
                6*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                3*pow(rgrid[i],5)*
                 (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qr1(i,j),2)*
                 pow(1 + rgrid[i]*Qtt(i,j),2)) + 
             pow(rgrid[i],3)*pow(Qrr(i,j),2)*
              (3*pow(rgrid[i],4)*
                 (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) - 
                Qtt(i,j)*(12 + rgrid[i]*Qtt(i,j))) + 
             pow(rgrid[i],2)*pow(Q22(i,j),2)*
              (-2*pow(rgrid[i],3)*pow(Qrr(i,j),3)*
                 (-1 + rgrid[i]*Qtt(i,j)) + 
                3*(-2 - 4*rgrid[i]*Qtt(i,j) - 
                   pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                2*rgrid[i]*Qrr(i,j)*
                 (-6 - 15*rgrid[i]*Qtt(i,j) - 
                   4*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                 (-1 - 14*rgrid[i]*Qtt(i,j) - 
                   2*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
             2*rgrid[i]*Q22(i,j)*
              (-(pow(rgrid[i],3)*pow(Qrr(i,j),3)*
                   (-1 + 4*rgrid[i]*Qtt(i,j) + 
                     pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
                3*(-3 - 6*rgrid[i]*Qtt(i,j) - 
                   2*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                 (-6 - 24*rgrid[i]*Qtt(i,j) - 
                   7*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                rgrid[i]*Qrr(i,j)*
                 (-19 - 44*rgrid[i]*Qtt(i,j) - 
                   15*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   6*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)))) + 
          2*rgrid[i]*Q11(i,j)*(pow(rgrid[i],3)*pow(Qrr(i,j),3)*
              (3 + pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
             3*(-4 - 8*rgrid[i]*Qtt(i,j) - 
                3*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
             pow(rgrid[i],2)*pow(Qrr(i,j),2)*
              (-5 - 22*rgrid[i]*Qtt(i,j) - 
                6*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                3*pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
             rgrid[i]*Qrr(i,j)*
              (-23 - 52*rgrid[i]*Qtt(i,j) - 
                19*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                6*pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                 pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
             2*rgrid[i]*Q22(i,j)*
              (-2*pow(rgrid[i],4)*pow(Qrr(i,j),3)*Qtt(i,j)*
                 (3 + rgrid[i]*Qtt(i,j)) + 
                3*(-4 - 8*rgrid[i]*Qtt(i,j) - 
                   3*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                 (-11 - 34*rgrid[i]*Qtt(i,j) - 
                   12*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                2*rgrid[i]*Qrr(i,j)*
                 (-13 - 29*rgrid[i]*Qtt(i,j) - 
                   11*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
             pow(rgrid[i],2)*pow(Q22(i,j),2)*
              (-(pow(rgrid[i],3)*pow(Qrr(i,j),3)*
                   (-1 + 4*rgrid[i]*Qtt(i,j) + 
                     pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
                3*(-3 - 6*rgrid[i]*Qtt(i,j) - 
                   2*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                 (-6 - 24*rgrid[i]*Qtt(i,j) - 
                   7*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                rgrid[i]*Qrr(i,j)*
                 (-19 - 44*rgrid[i]*Qtt(i,j) - 
                   15*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   6*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2))))) + 
       (1 + rgrid[i]*Qrr(i,j))*
        (-((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],6)*
             pow(-12 + pow(mu,2)*(-3 + 4*rgrid[i]),2)*
             pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
             pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j))) - 
          2*(2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],7)*pow(-12 + pow(mu,2)*(-3 + 4*rgrid[i]),2)*
           pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
           (1 + rgrid[i]*Qrr(i,j))*(Qrr(i,j) - Qtt(i,j))*
           (1 + rgrid[i]*Qtt(i,j)) - 
          (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*pow(4 - 
             (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),3)*
           pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
           pow(Qr1(i,j),2)*(7 + 6*rgrid[i]*Qrr(i,j))*
           pow(1 + rgrid[i]*Qtt(i,j),2) - 
          2*(2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*
           (-10 - 23*rgrid[i]*Q22(i,j) - 
             8*pow(rgrid[i],2)*pow(Q22(i,j),2) + 
             3*pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2) + 
             6*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*pow(Qr1(i,j),2) + 
             3*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Q22(i,j),2)*pow(Qr1(i,j),2) + rgrid[i]*Qrr(i,j) - 
             10*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + 
             3*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2)*Qrr(i,j) + 
             6*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*pow(Qr1(i,j),2)*Qrr(i,j) + 
             3*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Q22(i,j),2)*pow(Qr1(i,j),2)*Qrr(i,j) + 
             9*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
             9*pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) + 
             6*pow(rgrid[i],4)*pow(Q22(i,j),2)*pow(Qrr(i,j),2) - 
             23*rgrid[i]*Qtt(i,j) - 
             52*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
             19*pow(rgrid[i],3)*pow(Q22(i,j),2)*Qtt(i,j) + 
             6*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2)*Qtt(i,j) + 
             12*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*pow(Qr1(i,j),2)*Qtt(i,j) + 
             6*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Q22(i,j),2)*pow(Qr1(i,j),2)*Qtt(i,j) - 
             10*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
             44*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             12*pow(rgrid[i],4)*pow(Q22(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
             6*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
             12*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*pow(Qr1(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
             6*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Q22(i,j),2)*pow(Qr1(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
             9*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) + 
             3*pow(rgrid[i],5)*pow(Q22(i,j),2)*pow(Qrr(i,j),2)*
              Qtt(i,j) - 8*pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
             19*pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) - 
             6*pow(rgrid[i],4)*pow(Q22(i,j),2)*pow(Qtt(i,j),2) + 
             3*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
             6*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
             3*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Q22(i,j),2)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
             12*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) - 
             pow(rgrid[i],5)*pow(Q22(i,j),2)*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             3*pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qr1(i,j),2)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             6*pow(rgrid[i],9)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*pow(Qr1(i,j),2)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             3*pow(rgrid[i],10)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Q22(i,j),2)*pow(Qr1(i,j),2)*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             6*pow(rgrid[i],4)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
             3*pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*
              pow(Qtt(i,j),2) + 
             3*pow(rgrid[i],6)*pow(Q22(i,j),2)*pow(Qrr(i,j),2)*
              pow(Qtt(i,j),2) - 
             4*pow(mu,2)*pow(h(i,j),2)*pow(rgrid[i],2)*
              pow(1 + rgrid[i]*Q11(i,j),2)*
              pow(1 + rgrid[i]*Q22(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
              pow(1 + rgrid[i]*Qtt(i,j),2) + 
             pow(rgrid[i],2)*pow(Q11(i,j),2)*
              (-8 - 19*rgrid[i]*Qtt(i,j) - 
                6*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                3*pow(rgrid[i],5)*
                 (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qr1(i,j),2)*
                 pow(1 + rgrid[i]*Qtt(i,j),2) + 
                3*pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                 (2 + rgrid[i]*Qtt(i,j) + 
                   pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
                pow(rgrid[i],2)*Qrr(i,j)*
                 (3*pow(rgrid[i],4)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) - 
                   Qtt(i,j)*(12 + rgrid[i]*Qtt(i,j))) + 
                pow(rgrid[i],2)*pow(Q22(i,j),2)*
                 (-6 - 15*rgrid[i]*Qtt(i,j) - 
                   4*pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
                   3*pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                    (-1 + rgrid[i]*Qtt(i,j)) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
                    pow(1 + rgrid[i]*Qtt(i,j),2) - 
                   rgrid[i]*Qrr(i,j)*
                    (1 + 14*rgrid[i]*Qtt(i,j) + 
                      2*pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
                rgrid[i]*Q22(i,j)*
                 (-19 - 44*rgrid[i]*Qtt(i,j) - 
                   15*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   6*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
                    pow(1 + rgrid[i]*Qtt(i,j),2) - 
                   3*pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                    (-1 + 4*rgrid[i]*Qtt(i,j) + 
                      pow(rgrid[i],2)*pow(Qtt(i,j),2)) - 
                   2*rgrid[i]*Qrr(i,j)*
                    (6 + 24*rgrid[i]*Qtt(i,j) + 
                      7*pow(rgrid[i],2)*pow(Qtt(i,j),2)))) + 
             rgrid[i]*Q11(i,j)*
              (-23 + 6*pow(rgrid[i],5)*
                 (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qr1(i,j),2) + 
                4*rgrid[i]*(-13 + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2))*Qtt(i,j) + 
                pow(rgrid[i],2)*
                 (-19 + 6*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2))*pow(Qtt(i,j),2) + 
                3*pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                 (3 + pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
                2*rgrid[i]*Qrr(i,j)*
                 (-5 - 22*rgrid[i]*Qtt(i,j) - 
                   6*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                pow(rgrid[i],2)*pow(Q22(i,j),2)*
                 (-19 - 44*rgrid[i]*Qtt(i,j) - 
                   15*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   6*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
                    pow(1 + rgrid[i]*Qtt(i,j),2) - 
                   3*pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                    (-1 + 4*rgrid[i]*Qtt(i,j) + 
                      pow(rgrid[i],2)*pow(Qtt(i,j),2)) - 
                   2*rgrid[i]*Qrr(i,j)*
                    (6 + 24*rgrid[i]*Qtt(i,j) + 
                      7*pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
                4*rgrid[i]*Q22(i,j)*
                 (-13 - 29*rgrid[i]*Qtt(i,j) - 
                   11*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   3*pow(rgrid[i],5)*
                    (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                    pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
                    pow(1 + rgrid[i]*Qtt(i,j),2) - 
                   3*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j)*
                    (3 + rgrid[i]*Qtt(i,j)) - 
                   rgrid[i]*Qrr(i,j)*
                    (11 + 34*rgrid[i]*Qtt(i,j) + 
                      12*pow(rgrid[i],2)*pow(Qtt(i,j),2))))) - 
          (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
           (288*pow(A1,2) + 192*
              exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) - 
             48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3) - 
             72*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3) - 
             12*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3) - 
             18*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3) - 
             16*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4) - 
             24*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4) + 
             288*pow(A1,2)*rgrid[i]*Q22(i,j) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Q22(i,j) + 
             48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Q22(i,j) + 72*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Q22(i,j) + 12*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Q22(i,j) + 18*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Q22(i,j) - 48*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j) - 72*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j) + 288*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Qrr(i,j) + 
             48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qrr(i,j) + 72*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qrr(i,j) + 12*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qrr(i,j) + 18*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qrr(i,j) - 48*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qrr(i,j) - 72*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qrr(i,j) + 288*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*
              Qrr(i,j) + 192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                  rgrid[i])/(3.*A1))*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + 
             144*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j)*Qrr(i,j) + 
             216*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*Qrr(i,j) + 
             36*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*Qrr(i,j) + 
             54*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j)*Qrr(i,j) - 
             80*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qrr(i,j) - 
             120*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Q22(i,j)*Qrr(i,j) + 576*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
             384*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Qtt(i,j) - 
             360*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qtt(i,j) - 540*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qtt(i,j) - 90*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qtt(i,j) - 135*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              Qtt(i,j) + 88*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qtt(i,j) + 132*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qtt(i,j) + 576*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*
              Qtt(i,j) + 384*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                  rgrid[i])/(3.*A1))*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
             168*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j)*Qtt(i,j) - 
             252*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*Qtt(i,j) - 
             42*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*Qtt(i,j) - 
             63*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Q22(i,j)*Qtt(i,j) + 
             24*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qtt(i,j) + 
             36*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Q22(i,j)*Qtt(i,j) + 
             576*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
             384*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
             240*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qrr(i,j)*Qtt(i,j) - 
             360*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
             60*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
             90*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              Qrr(i,j)*Qtt(i,j) + 
             48*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) + 
             72*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Qrr(i,j)*Qtt(i,j) + 
             576*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*
              Qtt(i,j) + 384*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                  rgrid[i])/(3.*A1))*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*
              Qtt(i,j) - 48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             72*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             12*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             18*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             16*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             24*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
             288*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
             216*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              pow(Qtt(i,j),2) - 
             324*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
             54*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
             81*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              pow(Qtt(i,j),2) + 
             72*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
             108*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              pow(Qtt(i,j),2) + 
             288*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*
              pow(Qtt(i,j),2) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) - 
             120*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Q22(i,j)*pow(Qtt(i,j),2) - 
             180*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*pow(Qtt(i,j),2) - 
             30*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Q22(i,j)*pow(Qtt(i,j),2) - 
             45*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Q22(i,j)*pow(Qtt(i,j),2) + 
             40*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*Q22(i,j)*pow(Qtt(i,j),2) + 
             60*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              Q22(i,j)*pow(Qtt(i,j),2) + 
             288*pow(A1,2)*pow(rgrid[i],3)*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) - 
             192*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Qrr(i,j)*pow(Qtt(i,j),2) - 
             288*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) - 
             48*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) - 
             72*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              Qrr(i,j)*pow(Qtt(i,j),2) + 
             64*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             96*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              Qrr(i,j)*pow(Qtt(i,j),2) + 
             288*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) - 
             96*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) - 
             144*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) - 
             24*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) - 
             36*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             32*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],8)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             48*pow(A1,2)*pow(mu,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],8)*
              Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             4*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
              exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
              (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
             rgrid[i]*Q11(i,j)*
              (288*pow(A1,2) + 
                192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                   (3.*A1)) + 48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],3) + 
                72*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],3) + 
                12*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],3) + 
                18*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3) - 
                48*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],4) - 
                72*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4) + 
                288*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
                192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 rgrid[i]*Qrr(i,j) + 
                144*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                 Qrr(i,j) + 216*pow(A1,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                 Qrr(i,j) + 36*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                 Qrr(i,j) + 54*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                 Qrr(i,j) - 80*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                 Qrr(i,j) - 120*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                 Qrr(i,j) + 576*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
                384*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 rgrid[i]*Qtt(i,j) - 
                168*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                 Qtt(i,j) - 252*pow(A1,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                 Qtt(i,j) - 42*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                 Qtt(i,j) - 63*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                 Qtt(i,j) + 24*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                 Qtt(i,j) + 36*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                 Qtt(i,j) + 576*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*
                 Qtt(i,j) + 384*
                 exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
                48*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                 Qrr(i,j)*Qtt(i,j) - 
                72*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
                12*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
                18*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                 Qrr(i,j)*Qtt(i,j) - 
                16*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) - 
                24*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
                 Qrr(i,j)*Qtt(i,j) + 
                288*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
                120*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                 pow(Qtt(i,j),2) - 
                180*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
                30*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
                45*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
                 pow(Qtt(i,j),2) + 
                40*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
                60*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
                 pow(Qtt(i,j),2) + 
                288*pow(A1,2)*pow(rgrid[i],3)*Qrr(i,j)*
                 pow(Qtt(i,j),2) + 
                192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) - 
                96*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
                 Qrr(i,j)*pow(Qtt(i,j),2) - 
                144*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) - 
                24*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) - 
                36*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
                 Qrr(i,j)*pow(Qtt(i,j),2) + 
                32*pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 pow(rgrid[i],7)*Qrr(i,j)*pow(Qtt(i,j),2) + 
                48*pow(A1,2)*pow(mu,2)*
                 exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
                 Qrr(i,j)*pow(Qtt(i,j),2) + 
                rgrid[i]*Q22(i,j)*
                 (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                      3*pow(A1,2)*
                       (48 - exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i])))) + 
                   rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (192*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (12 + pow(mu,2)*(3 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                      pow(A1,2)*
                       (576 - 
                         3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 20*rgrid[i]))))*Qtt(i,j) \
+ pow(rgrid[i],2)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                      3*pow(A1,2)*
                       (96 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                    pow(Qtt(i,j),2) + 
                   2*rgrid[i]*Qrr(i,j)*
                    (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                       (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (60 + pow(mu,2)*(15 - 28*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                      pow(A1,2)*
                       (144 - 
                         3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-60 + pow(mu,2)*(-15 + 28*rgrid[i]))) + 
                      rgrid[i]*
                       (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                         3*pow(A1,2)*
                         (96 - 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i]))))*Qtt(i,j) \
+ 48*(3*pow(A1,2) + 2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*pow(rgrid[i],2)*pow(Qtt(i,j),2))))) \
- 2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*
           (72*pow(A1,2) + 48*
              exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i])) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i])) + 
             8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
             12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
             72*pow(A1,2)*rgrid[i]*Q22(i,j) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Q22(i,j) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j) + 6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j) + 9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j) + 144*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Qrr(i,j) - 
             24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qrr(i,j) - 
             36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Qrr(i,j) + 10*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Qrr(i,j) + 15*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j) + 
             144*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) - 
             24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qrr(i,j) - 
             36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j)*Qrr(i,j) + 
             6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j) + 
             9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*Qrr(i,j) + 
             72*pow(A1,2)*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*pow(Qrr(i,j),2) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              pow(Qrr(i,j),2) + 
             3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qrr(i,j),2) + 
             (9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                pow(Qrr(i,j),2))/2. + 
             72*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*
              pow(Qrr(i,j),2) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j)*pow(Qrr(i,j),2) + 
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
              pow(Qrr(i,j),2) + 
             (3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                Q22(i,j)*pow(Qrr(i,j),2))/2. + 
             144*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              rgrid[i]*Qtt(i,j) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qtt(i,j) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Qtt(i,j) + 16*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Qtt(i,j) + 24*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qtt(i,j) + 
             144*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qtt(i,j) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j)*Qtt(i,j) + 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qtt(i,j) + 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*Qtt(i,j) + 
             288*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
             24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qrr(i,j)*Qtt(i,j) - 
             36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Qrr(i,j)*Qtt(i,j) + 
             23*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j)*Qtt(i,j) + 
             (69*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                Qrr(i,j)*Qtt(i,j))/2. + 
             288*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
             192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
             24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
              Qtt(i,j) - 36*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
              Qtt(i,j) + 15*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
             (45*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                Q22(i,j)*Qrr(i,j)*Qtt(i,j))/2. + 
             144*pow(A1,2)*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) - 
             12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*pow(Qrr(i,j),2)*
              Qtt(i,j) - 18*pow(A1,2)*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*pow(Qrr(i,j),2)*
              Qtt(i,j) + 9*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qrr(i,j),2)*Qtt(i,j) + 
             (27*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                pow(Qrr(i,j),2)*Qtt(i,j))/2. + 
             144*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*pow(Qrr(i,j),2)*
              Qtt(i,j) + 96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                (3.*A1))*pow(rgrid[i],4)*Q22(i,j)*pow(Qrr(i,j),2)*
              Qtt(i,j) - 12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) - 
             18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*
              Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) + 
             5*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
              pow(Qrr(i,j),2)*Qtt(i,j) + 
             (15*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j))/2. + 
             72*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qtt(i,j),2) + 
             9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qtt(i,j),2) + 
             72*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) + 
             4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
              pow(Qtt(i,j),2) + 
             6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*pow(Qtt(i,j),2) + 
             144*pow(A1,2)*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             9*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             (27*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                Qrr(i,j)*pow(Qtt(i,j),2))/2. + 
             144*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             5*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*Qrr(i,j)*
              pow(Qtt(i,j),2) + 
             (15*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
                Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2))/2. + 
             72*pow(A1,2)*pow(rgrid[i],4)*pow(Qrr(i,j),2)*
              pow(Qtt(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
             4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qrr(i,j),2)*
              pow(Qtt(i,j),2) + 
             6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
             72*pow(A1,2)*pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*
              pow(Qtt(i,j),2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*
              pow(Qtt(i,j),2) + 
             2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],8)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*Q22(i,j)*
              pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
             3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],8)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
              Q22(i,j)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
             2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
              exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
              pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
              (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
              (1 + rgrid[i]*Qtt(i,j)) + 
             (rgrid[i]*Q11(i,j)*
                (2*(6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (8*exp(A1*mu*h(i,j)*rgrid[i]) - 
                        (4 + pow(mu,2))*pow(rgrid[i],3)) - 
                     9*pow(A1,2)*
                      (-8 + (4 + pow(mu,2))*
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)) + 
                     12*rgrid[i]*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (4*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (-4 + pow(mu,2)*(-1 + rgrid[i]))*
                        pow(rgrid[i],3)) + 
                        3*pow(A1,2)*
                         (4 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (-4 + pow(mu,2)*(-1 + rgrid[i]))*
                         pow(rgrid[i],3)))*Qtt(i,j) + 
                     2*pow(rgrid[i],2)*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (12*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (12 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                      pow(Qtt(i,j),2)) + 
                  pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                   (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                     3*pow(A1,2)*
                      (48 - exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i]))) - 
                     rgrid[i]*
                      (-2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (96*exp(A1*mu*h(i,j)*rgrid[i]) - 
                        pow(rgrid[i],3)*
                        (12 + pow(mu,2)*(3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (-96 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (12 + pow(mu,2)*(3 + 4*rgrid[i]))))*Qtt(i,j) + 
                     2*pow(rgrid[i],2)*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (24 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                      pow(Qtt(i,j),2)) + 
                  rgrid[i]*Qrr(i,j)*
                   (6*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (16*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (4 + pow(mu,2) - 4*pow(mu,2)*rgrid[i])) + 
                        3*pow(A1,2)*
                         (16 - 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-4 + pow(mu,2)*(-1 + 4*rgrid[i])))) + 
                     3*rgrid[i]*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (64*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-28 + pow(mu,2)*(-7 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (64 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-28 + pow(mu,2)*(-7 + 4*rgrid[i]))))*Qtt(i,j) \
+ pow(rgrid[i],2)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         5*pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                        3*pow(A1,2)*
                         (96 + 
                         5*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                      pow(Qtt(i,j),2)) + 
                  rgrid[i]*Q22(i,j)*
                   (pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                      (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (60 + pow(mu,2)*(15 - 28*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                        pow(A1,2)*
                         (144 - 
                         3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-60 + pow(mu,2)*(-15 + 28*rgrid[i]))) + 
                        rgrid[i]*
                         (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                        pow(rgrid[i],3)) + 
                         3*pow(A1,2)*
                         (96 - 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i]))))*
                         Qtt(i,j) + 
                        48*(3*pow(A1,2) + 
                         2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
                     2*(8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (6*exp(A1*mu*h(i,j)*rgrid[i]) - 
                         pow(mu,2)*pow(rgrid[i],4)) - 
                        12*pow(A1,2)*
                         (-6 + 
                         pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],4)) + 
                        2*rgrid[i]*
                         (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 2*rgrid[i]))) + 
                         pow(A1,2)*
                         (72 + 
                         3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 2*rgrid[i]))))*Qtt(i,j) \
+ pow(rgrid[i],2)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                         3*pow(A1,2)*
                         (24 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                         pow(Qtt(i,j),2)) + 
                     rgrid[i]*Qrr(i,j)*
                      (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (36 + pow(mu,2)*(9 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                         3*pow(A1,2)*
                         (48 - 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-36 + pow(mu,2)*(-9 + 20*rgrid[i])))) + 
                        rgrid[i]*
                         (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (192*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         (12 + pow(mu,2)*(3 - 20*rgrid[i]))*
                         pow(rgrid[i],3)) + 
                         pow(A1,2)*
                         (576 - 
                         3*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 20*rgrid[i]))))*Qtt(i,j) \
+ pow(rgrid[i],2)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         (96*exp(A1*mu*h(i,j)*rgrid[i]) + 
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                         3*pow(A1,2)*
                         (96 + 
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*
                         pow(Qtt(i,j),2)))))/2.))))/
   (2.*(2 + 3*pow(A1,2))*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
     pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Q22(i,j),2)*pow(1 + rgrid[i]*Qrr(i,j),2)*
     pow(1 + rgrid[i]*Qtt(i,j),2))
;
elem[1][2]+=
(2*pow(Qrr.diff(d01,i,j),2)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) - 
  (pow(Q11.diff(d10,i,j),2)*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q11(i,j),3) + 
  (Qrr.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q11(i,j),2) - 
  (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (2*Qrr.diff(d02,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) - 
  (pow(Q11.diff(d01,i,j),2)*L*pow(rgrid[i],4)*pow(Qr1(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Q11(i,j),3) + 
  Q11.diff(d10,i,j)*((2*Qr1.diff(d01,i,j)*L*pow(rgrid[i],2)*
        (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Q11(i,j),2) + 
     (2*Q11.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j)*
        (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q11(i,j),3)) + 
     ((-1 + rgrid[i]*Q11(i,j) - 4*rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qrr(i,j)))/pow(1 + rgrid[i]*Q11(i,j),3)) + 
  Q11.diff(d01,i,j)*((-2*Qr1.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j)*
        (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Q11(i,j),2) + 
     (rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
        (1 - rgrid[i]*Q11(i,j) + 4*rgrid[i]*Qrr(i,j)))/
      (L*pow(1 + rgrid[i]*Q11(i,j),3))) + 
  ((1 + rgrid[i]*Qrr(i,j))*(pow(rgrid[i],2)*
        (-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
       (12 - 6*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          7*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
       Q11(i,j)*(-12 + pow(mu,2)*pow(rgrid[i],4) + 
          rgrid[i]*(-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
             3*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j))))/
   ((-1 + rgrid[i])*rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3)) - 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*pow(1 + rgrid[i]*Qrr(i,j),2))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[1][3]+=
(-2*Qrr.diff(d11,i,j)*rgrid[i])/L + 
  (2*Qrr.diff(d02,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) - 
  (3*pow(Qrr.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) - 
  (4*h.diff(d01,i,j)*pow(mu,2)*h(i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   L - (4*h.diff(d01,i,j)*h.diff(d10,i,j)*pow(mu,2)*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/L - 
  2*Qr1.diff(d10,i,j)*(-1 + rgrid[i])*
   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
   (1 + rgrid[i]*Qrr(i,j)) - (Q11.diff(d01,i,j)*Q11.diff(d10,i,j)*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (Q22.diff(d01,i,j)*Q22.diff(d10,i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/pow(L,2) + 
  (4*Qr1.diff(d01,i,j)*rgrid[i]*
     (4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/L + 
  (pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Q11(i,j),2) + 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Q22(i,j),2) + 
  (Q22.diff(d01,i,j)*rgrid[i]*(Q22(i,j) - 2*Qrr(i,j))*
     (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (Qr1(i,j)*(-16 + 13*(4 + pow(mu,2))*pow(rgrid[i],3) - 
       16*pow(mu,2)*pow(rgrid[i],4) - 
       3*rgrid[i]*(4 - 4*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          5*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j)))/rgrid[i] - 
  (Qrr.diff(d01,i,j)*(4 + rgrid[i]*Qrr(i,j) + 
       3*(-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
        (1 + rgrid[i]*Qrr(i,j))))/(L + L*rgrid[i]*Qrr(i,j)) + 
  Q11.diff(d01,i,j)*((2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
        (1 + rgrid[i]*Qrr(i,j)))/(pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
     (rgrid[i]*(Q11(i,j) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
      (L*pow(1 + rgrid[i]*Q11(i,j),2))) + 
  Qrr.diff(d10,i,j)*((-1 + rgrid[i])*rgrid[i]*
      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j) + 
     (3*Qrr.diff(d01,i,j)*pow(rgrid[i],2))/(L + L*rgrid[i]*Qrr(i,j))) - 
  (Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*a0.diff(d01,i,j)*a0.diff(d10,i,j)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Qtt(i,j),2) + 
  (Qtt.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j))*
     ((-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
       (4 + 2*(4 + pow(mu,2))*pow(rgrid[i],3) - 
          3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2))
;
elem[1][4]+=
-((pow(Q22.diff(d10,i,j),2)*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
     pow(1 + rgrid[i]*Q22(i,j),3)) + 
  (Qrr.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q22(i,j),2) - 
  (pow(Q22.diff(d01,i,j),2)*L*pow(rgrid[i],4)*pow(Qr1(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j)))/pow(L + L*rgrid[i]*Q22(i,j),3) + 
  (Q22.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
     (1 - rgrid[i]*Q22(i,j) + 4*rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q22(i,j),3)) + 
  Q22.diff(d10,i,j)*((2*Q22.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j)*
        (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),3)) + 
     ((-1 + rgrid[i]*Q22(i,j) - 4*rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qrr(i,j)))/pow(1 + rgrid[i]*Q22(i,j),3)) + 
  ((1 + rgrid[i]*Qrr(i,j))*(pow(rgrid[i],2)*
        (-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
       (12 - 6*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          7*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
       Q22(i,j)*(-12 + pow(mu,2)*pow(rgrid[i],4) + 
          rgrid[i]*(-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
             3*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j))))/
   ((-1 + rgrid[i])*rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q22(i,j),3))
;
elem[1][5]+=
(-4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[1][6]+=
4*h.diff(d10,i,j)*pow(mu,2)*(1 + rgrid[i]*Qrr(i,j)) - 
  (4*h.diff(d01,i,j)*pow(mu,2)*rgrid[i]*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/L - 
  (2*pow(a0.diff(d10,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (4*a0.diff(d01,i,j)*B1*pow(mu,3)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (2*pow(a0.diff(d01,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j))*
     ((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) - 2*(1 + rgrid[i]*Qrr(i,j))))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  a0.diff(d10,i,j)*((-4*B1*pow(mu,3)*a0(i,j)*
        exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
        (1 + rgrid[i]*Qrr(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
     (4*a0.diff(d01,i,j)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (2*mu*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*(1 + rgrid[i]*Qrr(i,j))*
     (-((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
          exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)\
) + 2*((2 + 3*pow(A1,2))*mu*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*h(i,j)*
           (-1 + rgrid[i])*rgrid[i]*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3)) - 
          12*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
               (3.*A1)))*(1 + rgrid[i]*Qrr(i,j)))*(1 + rgrid[i]*Qtt(i,j))))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[2][0]+=
(2*Qtt.diff(d01,i,j)*rgrid[i]*
     (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),3)) - 
  (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (Q11.diff(d10,i,j)*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Qrr(i,j)))/
   (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) - 2*(1 + rgrid[i]*Qrr(i,j))))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (-4*pow(mu,2)*pow(a0(i,j),2)*exp(B1*mu*h(i,j)*rgrid[i])*
      pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j)) + 
     (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
        2*pow(mu,2)*pow(rgrid[i],4))*(2 + rgrid[i]*Q11(i,j))*
      (1 + rgrid[i]*Qrr(i,j)))/
   (2.*pow(rgrid[i],2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (2*pow(Qtt.diff(d01,i,j),2)*L*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Qtt(i,j),3)) + 
  a0.diff(d10,i,j)*((-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)))
;
elem[2][1]+=
(4*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(rgrid[i],2))/
   (pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))) + 
  (2*Q11.diff(d02,i,j)*rgrid[i])/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (3*pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q22(i,j),2)) + 
  (2*Q11.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/(L + L*rgrid[i]*Q11(i,j)) + 
  (2*Q22.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (L*pow(1 + rgrid[i]*Q22(i,j),2)) - 
  (pow(Qrr.diff(d01,i,j),2)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Qrr(i,j),2)) + 
  (2*Qr1.diff(d10,i,j)*Qrr.diff(d01,i,j)*L*pow(rgrid[i],2)*
     (1 + rgrid[i]*Q11(i,j)))/pow(L + L*rgrid[i]*Qrr(i,j),2) + 
  Qrr.diff(d01,i,j)*(-((pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))*
          Qr1(i,j)*((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
            pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
               pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)))/
        (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qrr(i,j),2))) \
- (2*Qr1.diff(d01,i,j)*pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))*
        Qr1(i,j))/pow(L + L*rgrid[i]*Qrr(i,j),2)) - 
  (Qtt.diff(d01,i,j)*rgrid[i]*
     (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))) + 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
  (Q11.diff(d10,i,j)*((pow(rgrid[i],3)*
          (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j)))/2. - 
       (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        (3 + 2*rgrid[i]*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q11(i,j)*(2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)))))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        (-1 + rgrid[i])*rgrid[i]*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(2 + rgrid[i]*Q11(i,j))*
        (3 + 2*rgrid[i]*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q11(i,j)*(2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j))) - 
       (rgrid[i]*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                  pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))\
) + 3*pow(A1,2)*(24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],3)*
                   (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
               24*(3*pow(A1,2) + 
                  2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))\
)*rgrid[i]*Qtt(i,j)) + rgrid[i]*Q11(i,j)*
             (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                  pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) \
+ 3*pow(A1,2)*(48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))\
) + 48*(3*pow(A1,2) + 2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                     (3.*A1)))*rgrid[i]*Qtt(i,j))))/2.))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[2][2]+=
(-2*pow(Qr1.diff(d01,i,j),2)*pow(rgrid[i],2))/pow(L,2) + 
  (pow(Q11.diff(d10,i,j),2)*pow(rgrid[i],2))/
   pow(1 + rgrid[i]*Q11(i,j),2) + 
  (Qr1.diff(d01,i,j)*(-4 + (-8*pow(rgrid[i],2) + 
          2*(4 + pow(mu,2))*pow(rgrid[i],5) - 
          2*pow(mu,2)*pow(rgrid[i],6))*pow(Qr1(i,j),2)))/L - 
  (2*Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*Qrr(i,j))/
   (L*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (2*Q11.diff(d02,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (2*Q22.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],2)*
     (6 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 6*rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),3)) - 
  (2*Qr1.diff(d10,i,j)*Qrr.diff(d01,i,j)*pow(rgrid[i],2))/
   (L + L*rgrid[i]*Qrr(i,j)) + 
  Q11.diff(d10,i,j)*((-2*Q11.diff(d01,i,j)*L*pow(rgrid[i],3)*
        Qr1(i,j))/pow(L + L*rgrid[i]*Q11(i,j),2) + 
     (-1 + rgrid[i]*Qrr(i,j))/pow(1 + rgrid[i]*Q11(i,j),2)) + 
  Qrr.diff(d01,i,j)*((2*Qr1.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
     (pow(rgrid[i],3)*Qr1(i,j)*
        ((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
          pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j)))) - 
  (Qtt.diff(d01,i,j)*rgrid[i]*
     (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],4)*pow(Qr1(i,j),2))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  a0.diff(d10,i,j)*((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        pow(rgrid[i],2))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((rgrid[i]*(1 + rgrid[i]*Q11(i,j))*
          (8*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
             exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
             (1 + rgrid[i]*Q22(i,j)) + 
            (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
             pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
               pow(mu,2)*pow(rgrid[i],3),2)*
             (3 + 2*rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
             pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)) - 
            (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
             (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                     pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                  3*pow(A1,2)*
                   (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
               rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                     pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                  3*pow(A1,2)*
                   (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) + 
               rgrid[i]*Qrr(i,j)*
                (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                     pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                  3*pow(A1,2)*
                   (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                  48*(3*pow(A1,2) + 
                     2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)))*rgrid[i]*Qtt(i,j))) + 
            2*(2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             (8 + 5*rgrid[i]*Qtt(i,j) + 
               2*rgrid[i]*Q11(i,j)*(1 + rgrid[i]*Qrr(i,j))*
                (2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)) + 
               rgrid[i]*Qrr(i,j)*(7 + 4*rgrid[i]*Qtt(i,j)) + 
               rgrid[i]*Q22(i,j)*
                (5 + 2*rgrid[i]*Qtt(i,j) + 
                  rgrid[i]*Qrr(i,j)*(4 + rgrid[i]*Qtt(i,j)))) - 
            (1 + rgrid[i]*Q22(i,j))*
             (2*(2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (12*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (12 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                  rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ rgrid[i]*Qrr(i,j)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     24*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j))) + 
               rgrid[i]*Q11(i,j)*
                (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                  rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ rgrid[i]*Qrr(i,j)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     48*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j))))))/2. - 
       rgrid[i]*(2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
           exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*pow(1 + rgrid[i]*Q11(i,j),2)*
           (1 + rgrid[i]*Q22(i,j)) + 
          ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
             pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
               pow(mu,2)*pow(rgrid[i],3),2)*
             (2 + 3*rgrid[i]*Q11(i,j) + 
               pow(rgrid[i],2)*pow(Q11(i,j),2))*
             (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
             (1 + rgrid[i]*Qtt(i,j)))/2. - 
          ((1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
             (2*(2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (12*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (12 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                  rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ rgrid[i]*Qrr(i,j)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     24*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j))) + 
               rgrid[i]*Q11(i,j)*
                (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                  rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ rgrid[i]*Qrr(i,j)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     48*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j)))))/2. + 
          (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*
           (pow(rgrid[i],2)*pow(Q11(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
              (2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)) + 
             2*(1 + rgrid[i]*Qrr(i,j))*
              (3 + 2*rgrid[i]*Qtt(i,j) + 
                rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
             rgrid[i]*Q11(i,j)*
              (8 + 5*rgrid[i]*Qtt(i,j) + 
                rgrid[i]*Qrr(i,j)*(7 + 4*rgrid[i]*Qtt(i,j)) + 
                rgrid[i]*Q22(i,j)*
                 (5 + 2*rgrid[i]*Qtt(i,j) + 
                   rgrid[i]*Qrr(i,j)*(4 + rgrid[i]*Qtt(i,j))))))))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[2][3]+=
(-2*Q11.diff(d11,i,j)*rgrid[i])/L + 
  (2*Q11.diff(d02,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) - 
  (2*pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*(1 + rgrid[i]*Q11(i,j))) - 
  (4*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/L + 
  ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(2 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   rgrid[i] + Q11.diff(d10,i,j)*
   ((2*Q11.diff(d01,i,j)*pow(rgrid[i],2))/(L + L*rgrid[i]*Q11(i,j)) - 
     (-1 + rgrid[i])*rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)) - 
  (2*Q11.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j) - rgrid[i]*Qrr(i,j)))/
   (L + L*rgrid[i]*Q11(i,j)) + 
  (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  Qrr.diff(d01,i,j)*((2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
        (1 + rgrid[i]*Q11(i,j)))/(pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
     (pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*
        ((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
          3*pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j)))) - 
  (Qtt.diff(d01,i,j)*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*a0.diff(d10,i,j)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (1 + rgrid[i]*Q11(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[2][4]+=
(Q11.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q22(i,j),2) - 
  ((2 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(rgrid[i],2)*pow(1 + rgrid[i]*Q22(i,j),2)) - 
  (2*pow(Q22.diff(d01,i,j),2)*L*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q22(i,j),3)) - 
  (4*Q22.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),3))
;
elem[2][5]+=
(4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
     (1 + rgrid[i]*Q11(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
     (1 + rgrid[i]*Q11(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[2][6]+=
(2*pow(a0.diff(d10,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*B1*pow(mu,3)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (2*pow(a0.diff(d01,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) - 2*(1 + rgrid[i]*Qrr(i,j))))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))) + 
  a0.diff(d10,i,j)*((4*B1*pow(mu,3)*a0(i,j)*
        exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
        (1 + rgrid[i]*Q11(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
     (4*a0.diff(d01,i,j)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (2*mu*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*(1 + rgrid[i]*Q11(i,j))*
     ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4) - 
       24*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[3][0]+=
(-2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (Q22.diff(d10,i,j)*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Qrr(i,j)))/
   (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*
     (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (-4*pow(mu,2)*pow(a0(i,j),2)*exp(B1*mu*h(i,j)*rgrid[i])*
      pow(rgrid[i],4)*(1 + rgrid[i]*Q22(i,j)) + 
     (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
        2*pow(mu,2)*pow(rgrid[i],4))*(2 + rgrid[i]*Q22(i,j))*
      (1 + rgrid[i]*Qrr(i,j)))/
   (2.*pow(rgrid[i],2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  a0.diff(d10,i,j)*((-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)))
;
elem[3][1]+=
(2*Q22.diff(d02,i,j)*rgrid[i])/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (2*pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],2))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))) + 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (Q22.diff(d10,i,j)*((pow(rgrid[i],3)*
          (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j)))/2. - 
       (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        (3 + 2*rgrid[i]*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q11(i,j)*(2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)))))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        (-1 + rgrid[i])*rgrid[i]*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(2 + rgrid[i]*Q22(i,j))*
        (3 + 2*rgrid[i]*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q11(i,j)*(2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j))) - 
       (rgrid[i]*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                  pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))\
) + 3*pow(A1,2)*(24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],3)*
                   (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
               24*(3*pow(A1,2) + 
                  2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))\
)*rgrid[i]*Qtt(i,j)) + rgrid[i]*Q22(i,j)*
             (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                  pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) \
+ 3*pow(A1,2)*(48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))\
) + 48*(3*pow(A1,2) + 2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                     (3.*A1)))*rgrid[i]*Qtt(i,j))))/2.))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[3][2]+=
(Q22.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q11(i,j),2) - 
  (2*Q22.diff(d02,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (2*pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Q22(i,j))) - 
  ((2 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(rgrid[i],2)*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[3][3]+=
(-2*Q22.diff(d11,i,j)*rgrid[i])/L - 
  (2*Q22.diff(d01,i,j)*(2 + rgrid[i]*Q22(i,j)))/
   (L + L*rgrid[i]*Q22(i,j)) + 
  (2*Q22.diff(d02,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) - 
  (2*pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*(1 + rgrid[i]*Q22(i,j))) + 
  ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(2 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   rgrid[i] + Q22.diff(d10,i,j)*
   ((2*Q22.diff(d01,i,j)*pow(rgrid[i],2))/(L + L*rgrid[i]*Q22(i,j)) - 
     (-1 + rgrid[i])*rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*a0.diff(d10,i,j)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (1 + rgrid[i]*Q22(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[3][4]+=
(pow(Q22.diff(d10,i,j),2)*pow(rgrid[i],2))/
   pow(1 + rgrid[i]*Q22(i,j),2) + 
  (2*Q22.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],2)*
     (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Q22(i,j),2)) + 
  Q22.diff(d10,i,j)*((-2*Q22.diff(d01,i,j)*L*pow(rgrid[i],3)*
        Qr1(i,j))/pow(L + L*rgrid[i]*Q22(i,j),2) + 
     (-1 + rgrid[i]*Qrr(i,j))/pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  a0.diff(d10,i,j)*((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        pow(rgrid[i],2))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((rgrid[i]*(1 + rgrid[i]*Q22(i,j))*
          (8*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
             exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
             (1 + rgrid[i]*Q22(i,j)) + 
            (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
             pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
               pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
             (3 + 2*rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
             (1 + rgrid[i]*Qtt(i,j)) - 
            (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
             (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                     pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                  3*pow(A1,2)*
                   (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
               rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                     pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                  3*pow(A1,2)*
                   (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) + 
               rgrid[i]*Qrr(i,j)*
                (2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                     pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                  3*pow(A1,2)*
                   (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],3)*
                      (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                  48*(3*pow(A1,2) + 
                     2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)))*rgrid[i]*Qtt(i,j))) + 
            2*(2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             (8 + 5*rgrid[i]*Qtt(i,j) + 
               2*rgrid[i]*Q22(i,j)*(1 + rgrid[i]*Qrr(i,j))*
                (2 + rgrid[i]*Qtt(i,j)) + 
               rgrid[i]*Qrr(i,j)*(7 + 4*rgrid[i]*Qtt(i,j)) + 
               rgrid[i]*Q11(i,j)*
                (5 + 2*rgrid[i]*Q22(i,j)*(1 + rgrid[i]*Qrr(i,j)) + 
                  2*rgrid[i]*Qtt(i,j) + 
                  rgrid[i]*Qrr(i,j)*(4 + rgrid[i]*Qtt(i,j)))) - 
            (1 + rgrid[i]*Q11(i,j))*
             (2*(2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (12*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (12 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                  rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ rgrid[i]*Qrr(i,j)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     24*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j))) + 
               rgrid[i]*Q22(i,j)*
                (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                  rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ rgrid[i]*Qrr(i,j)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     48*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j))))))/2. - 
       rgrid[i]*(2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
           exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
           pow(1 + rgrid[i]*Q22(i,j),2) + 
          ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
             pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
               pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
             (2 + 3*rgrid[i]*Q22(i,j) + 
               pow(rgrid[i],2)*pow(Q22(i,j),2))*pow(Qr1(i,j),2)*
             (1 + rgrid[i]*Qtt(i,j)))/2. - 
          ((1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
             (2*(2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (12*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (12 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                  rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ rgrid[i]*Qrr(i,j)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     24*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j))) + 
               rgrid[i]*Q22(i,j)*
                (2*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (24*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (24 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i])))) + 
                  rgrid[i]*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))))*Qtt(i,j) \
+ rgrid[i]*Qrr(i,j)*(2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      (48*exp(A1*mu*h(i,j)*rgrid[i]) + 
                        pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     3*pow(A1,2)*
                      (48 + exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],3)*
                         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))) + 
                     48*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*rgrid[i]*Qtt(i,j)))))/2. + 
          (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*
           (pow(rgrid[i],2)*pow(Q22(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
              (2 + rgrid[i]*Qtt(i,j)) + 
             2*(1 + rgrid[i]*Qrr(i,j))*(3 + 2*rgrid[i]*Qtt(i,j)) + 
             rgrid[i]*Q22(i,j)*
              (8 + 5*rgrid[i]*Qtt(i,j) + 
                rgrid[i]*Qrr(i,j)*(7 + 4*rgrid[i]*Qtt(i,j))) + 
             rgrid[i]*Q11(i,j)*
              (pow(rgrid[i],2)*pow(Q22(i,j),2)*
                 (1 + rgrid[i]*Qrr(i,j)) + 
                2*(1 + rgrid[i]*Qrr(i,j))*(2 + rgrid[i]*Qtt(i,j)) + 
                rgrid[i]*Q22(i,j)*
                 (5 + 2*rgrid[i]*Qtt(i,j) + 
                   rgrid[i]*Qrr(i,j)*(4 + rgrid[i]*Qtt(i,j))))))))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Q22(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[3][5]+=
(4*a0.diff(d10,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
     (1 + rgrid[i]*Q22(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
     (1 + rgrid[i]*Q22(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[3][6]+=
(2*pow(a0.diff(d10,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*B1*pow(mu,3)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (2*pow(a0.diff(d01,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*
     (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  a0.diff(d10,i,j)*((4*B1*pow(mu,3)*a0(i,j)*
        exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
        (1 + rgrid[i]*Q22(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
     (4*a0.diff(d01,i,j)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (2*mu*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*(1 + rgrid[i]*Q22(i,j))*
     ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4) - 
       24*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[4][0]+=
-(Qr1.diff(d10,i,j)*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
        2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Qrr(i,j)))/
   (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (8*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (8*a0.diff(d01,i,j)*a0.diff(d10,i,j)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(L + L*rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (8*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (2*pow(Qtt.diff(d01,i,j),2)*L*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Qtt(i,j),3)) - 
  (Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
     (-12 - 3*(4 + pow(mu,2))*pow(rgrid[i],3) + 
       7*pow(mu,2)*pow(rgrid[i],4) + 
       rgrid[i]*(-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          5*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j)))/
   ((-1 + rgrid[i])*rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),3)) - 
  (Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
     ((-1 + rgrid[i])*pow(rgrid[i],2)*
        (32 + 32*rgrid[i] + 32*pow(rgrid[i],2) - 
          4*(-4 + pow(mu,2))*pow(rgrid[i],3) - 
          4*(-4 + pow(mu,2))*pow(rgrid[i],4) - 
          4*(-4 + pow(mu,2))*pow(rgrid[i],5) - 
          pow(mu,2)*(12 + pow(mu,2))*pow(rgrid[i],6) + 
          2*pow(mu,4)*pow(rgrid[i],7))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*
        (8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
          2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
       (-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          3*pow(mu,2)*pow(rgrid[i],4))*(-1 + rgrid[i]*Qtt(i,j))))/
   (pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4),2)*(L + L*rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),3)) + 
  Qtt.diff(d10,i,j)*((-2*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*
        (1 + rgrid[i]*Qrr(i,j)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(1 + rgrid[i]*Qtt(i,j),3)) + 
     (rgrid[i]*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          2*pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),3)))
;
elem[4][1]+=
(2*Qr1.diff(d02,i,j)*rgrid[i])/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (4*h.diff(d01,i,j)*pow(mu,2)*h(i,j)*rgrid[i])/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (4*h.diff(d01,i,j)*h.diff(d10,i,j)*pow(mu,2)*pow(rgrid[i],2))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),3)) - 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Q22(i,j),2)) + 
  Q11.diff(d10,i,j)*(-((Q11.diff(d01,i,j)*pow(rgrid[i],2))/
        (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3))) \
+ (rgrid[i]*Qr1(i,j))/pow(1 + rgrid[i]*Q11(i,j),2)) + 
  Q22.diff(d10,i,j)*((Q22.diff(d01,i,j)*pow(rgrid[i],2))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(1 + rgrid[i]*Q22(i,j),2)) + 
     (rgrid[i]*Qr1(i,j))/pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (Q22.diff(d01,i,j)*(2 - rgrid[i]*Q22(i,j) - 
       (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 4*rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Q22(i,j),2)) - 
  (pow(Qrr.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Qrr(i,j),2)) - 
  (Qr1.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*
     (1 - rgrid[i]*Q11(i,j) + 4*rgrid[i]*Qrr(i,j) + 
       2*pow(rgrid[i],2)*pow(Qrr(i,j),2)))/
   ((L + L*rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  Q11.diff(d01,i,j)*((2*Qr1.diff(d01,i,j)*pow(rgrid[i],2))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) \
+ (6 + (-4*pow(rgrid[i],2) + (4 + pow(mu,2))*pow(rgrid[i],5) - 
           pow(mu,2)*pow(rgrid[i],6))*pow(Qr1(i,j),2) + 
        Q11(i,j)*(rgrid[i] + (-4*pow(rgrid[i],3) + 
              (4 + pow(mu,2))*pow(rgrid[i],6) - 
              pow(mu,2)*pow(rgrid[i],7))*pow(Qr1(i,j),2)) + 
        4*rgrid[i]*Qrr(i,j))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3))) + 
  Qrr.diff(d10,i,j)*((Qr1.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(1 + rgrid[i]*Qrr(i,j),2) + 
     (Qrr.diff(d01,i,j)*pow(rgrid[i],2))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(1 + rgrid[i]*Qrr(i,j),2)) - 
     (pow(rgrid[i],3)*Qr1(i,j)*
        ((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
          pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)))/
      (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
     (Qr1.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      pow(L + L*rgrid[i]*Qrr(i,j),2)) + 
  Qrr.diff(d01,i,j)*((-2 - pow(rgrid[i],5)*
         (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*(1 + rgrid[i]*Q11(i,j))*
         pow(Qr1(i,j),2) + pow(rgrid[i],4)*
         pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
           pow(mu,2)*pow(rgrid[i],4),2)*(1 + rgrid[i]*Q11(i,j))*
         pow(Qr1(i,j),4))/
      (2.*L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(1 + rgrid[i]*Qrr(i,j),2)) + 
     (Qr1.diff(d01,i,j)*pow(rgrid[i],4)*pow(Qr1(i,j),2))/
      pow(L + L*rgrid[i]*Qrr(i,j),2)) - 
  (8*a0.diff(d01,i,j)*a0.diff(d10,i,j)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2))/
   (L*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (8*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2))/
   (L*(-1 + rgrid[i])*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (8*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
  Qtt.diff(d10,i,j)*((Qtt.diff(d01,i,j)*pow(rgrid[i],2))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(1 + rgrid[i]*Qtt(i,j),2)) - 
     (rgrid[i]*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          2*pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j))/
      (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
  (Qtt.diff(d01,i,j)*((-1 + rgrid[i])*pow(rgrid[i],2)*
        (32 + 32*rgrid[i] + 32*pow(rgrid[i],2) - 
          4*(-4 + pow(mu,2))*pow(rgrid[i],3) - 
          4*(-4 + pow(mu,2))*pow(rgrid[i],4) - 
          4*(-4 + pow(mu,2))*pow(rgrid[i],5) - 
          pow(mu,2)*(12 + pow(mu,2))*pow(rgrid[i],6) + 
          2*pow(mu,4)*pow(rgrid[i],7))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*(8 + 4*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],3) - 
          2*pow(mu,2)*pow(rgrid[i],4) + 
          2*rgrid[i]*(8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
             2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
          rgrid[i]*(-4 - 2*(4 + pow(mu,2))*pow(rgrid[i],3) + 
             3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j))))/
   (2.*pow(-1 + rgrid[i],2)*pow(4 + 4*rgrid[i] + 
       4*pow(rgrid[i],2) - pow(mu,2)*pow(rgrid[i],3),2)*
     (L + L*rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  Qr1.diff(d10,i,j)*(-((Qrr.diff(d01,i,j)*L*pow(rgrid[i],3)*
          Qr1(i,j))/pow(L + L*rgrid[i]*Qrr(i,j),2)) + 
     (-32 + 20*pow(rgrid[i],3) + 5*pow(mu,2)*pow(rgrid[i],3) - 
        4*pow(mu,2)*pow(rgrid[i],4) - 48*rgrid[i]*Qrr(i,j) + 
        24*pow(rgrid[i],4)*Qrr(i,j) + 
        6*pow(mu,2)*pow(rgrid[i],4)*Qrr(i,j) - 
        4*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j) - 
        24*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
        12*pow(rgrid[i],5)*pow(Qrr(i,j),2) + 
        3*pow(mu,2)*pow(rgrid[i],5)*pow(Qrr(i,j),2) - 
        2*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2) - 
        24*rgrid[i]*Qtt(i,j) + 24*pow(rgrid[i],4)*Qtt(i,j) + 
        6*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
        6*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
        32*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
        32*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) + 
        8*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
        8*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) - 
        16*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) + 
        16*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) + 
        4*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) - 
        4*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*Qtt(i,j) + 
        rgrid[i]*Q22(i,j)*(-24 + 12*pow(rgrid[i],3) + 
           3*pow(mu,2)*pow(rgrid[i],3) - 
           2*pow(mu,2)*pow(rgrid[i],4) - 
           4*(4*rgrid[i] - (4 + pow(mu,2))*pow(rgrid[i],4) + 
              pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j) + 
           2*rgrid[i]*Qrr(i,j)*
            (-16 + (4 + pow(mu,2))*pow(rgrid[i],3) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j)) + 
           pow(rgrid[i],2)*pow(Qrr(i,j),2)*
            (-16 + (4 + pow(mu,2))*pow(rgrid[i],3) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j))) + 
        rgrid[i]*Q11(i,j)*(-24 + 12*pow(rgrid[i],3) + 
           3*pow(mu,2)*pow(rgrid[i],3) - 
           2*pow(mu,2)*pow(rgrid[i],4) - 16*rgrid[i]*Qtt(i,j) + 
           16*pow(rgrid[i],4)*Qtt(i,j) + 
           4*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
           4*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) + 
           2*rgrid[i]*Qrr(i,j)*
            (-16 + (4 + pow(mu,2))*pow(rgrid[i],3) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j)) + 
           pow(rgrid[i],2)*pow(Qrr(i,j),2)*
            (-16 + (4 + pow(mu,2))*pow(rgrid[i],3) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j)) + 
           rgrid[i]*Q22(i,j)*(-16 + 4*pow(rgrid[i],3) + 
              pow(mu,2)*pow(rgrid[i],3) + 
              2*rgrid[i]*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                 2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
              pow(rgrid[i],2)*
               (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                 2*pow(mu,2)*pow(rgrid[i],4))*pow(Qrr(i,j),2) + 
              (-8*rgrid[i] + 2*(4 + pow(mu,2))*pow(rgrid[i],4) - 
                 2*pow(mu,2)*pow(rgrid[i],5))*Qtt(i,j))))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
        (1 + rgrid[i]*Qtt(i,j)))) + 
  (Qr1(i,j)*(-48 + 72*pow(rgrid[i],3) + 
       18*pow(mu,2)*pow(rgrid[i],3) - 
       16*pow(mu,2)*pow(rgrid[i],4) + 
       16*pow(rgrid[i],2)*pow(Qr1(i,j),2) - 
       32*pow(rgrid[i],5)*pow(Qr1(i,j),2) - 
       8*pow(mu,2)*pow(rgrid[i],5)*pow(Qr1(i,j),2) + 
       8*pow(mu,2)*pow(rgrid[i],6)*pow(Qr1(i,j),2) + 
       16*pow(rgrid[i],8)*pow(Qr1(i,j),2) + 
       8*pow(mu,2)*pow(rgrid[i],8)*pow(Qr1(i,j),2) + 
       pow(mu,4)*pow(rgrid[i],8)*pow(Qr1(i,j),2) - 
       8*pow(mu,2)*pow(rgrid[i],9)*pow(Qr1(i,j),2) - 
       2*pow(mu,4)*pow(rgrid[i],9)*pow(Qr1(i,j),2) + 
       pow(mu,4)*pow(rgrid[i],10)*pow(Qr1(i,j),2) - 
       96*rgrid[i]*Qrr(i,j) + 120*pow(rgrid[i],4)*Qrr(i,j) + 
       30*pow(mu,2)*pow(rgrid[i],4)*Qrr(i,j) - 
       24*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j) - 
       48*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
       60*pow(rgrid[i],5)*pow(Qrr(i,j),2) + 
       15*pow(mu,2)*pow(rgrid[i],5)*pow(Qrr(i,j),2) - 
       12*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2) - 
       72*rgrid[i]*Qtt(i,j) + 168*pow(rgrid[i],4)*Qtt(i,j) + 
       42*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
       46*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) + 
       32*pow(rgrid[i],3)*pow(Qr1(i,j),2)*Qtt(i,j) - 
       64*pow(rgrid[i],6)*pow(Qr1(i,j),2)*Qtt(i,j) - 
       16*pow(mu,2)*pow(rgrid[i],6)*pow(Qr1(i,j),2)*Qtt(i,j) + 
       16*pow(mu,2)*pow(rgrid[i],7)*pow(Qr1(i,j),2)*Qtt(i,j) + 
       32*pow(rgrid[i],9)*pow(Qr1(i,j),2)*Qtt(i,j) + 
       16*pow(mu,2)*pow(rgrid[i],9)*pow(Qr1(i,j),2)*Qtt(i,j) + 
       2*pow(mu,4)*pow(rgrid[i],9)*pow(Qr1(i,j),2)*Qtt(i,j) - 
       16*pow(mu,2)*pow(rgrid[i],10)*pow(Qr1(i,j),2)*Qtt(i,j) - 
       4*pow(mu,4)*pow(rgrid[i],10)*pow(Qr1(i,j),2)*Qtt(i,j) + 
       2*pow(mu,4)*pow(rgrid[i],11)*pow(Qr1(i,j),2)*Qtt(i,j) - 
       144*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
       288*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) + 
       72*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
       76*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) - 
       72*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) + 
       144*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) + 
       36*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) - 
       38*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*Qtt(i,j) - 
       32*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
       92*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
       23*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
       28*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
       16*pow(rgrid[i],4)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
       32*pow(rgrid[i],7)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
       8*pow(mu,2)*pow(rgrid[i],7)*pow(Qr1(i,j),2)*
        pow(Qtt(i,j),2) + 8*pow(mu,2)*pow(rgrid[i],8)*
        pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
       16*pow(rgrid[i],10)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
       8*pow(mu,2)*pow(rgrid[i],10)*pow(Qr1(i,j),2)*
        pow(Qtt(i,j),2) + pow(mu,4)*pow(rgrid[i],10)*
        pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
       8*pow(mu,2)*pow(rgrid[i],11)*pow(Qr1(i,j),2)*
        pow(Qtt(i,j),2) - 2*pow(mu,4)*pow(rgrid[i],11)*
        pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
       pow(mu,4)*pow(rgrid[i],12)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
       64*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
       160*pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) + 
       40*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) - 
       48*pow(mu,2)*pow(rgrid[i],7)*Qrr(i,j)*pow(Qtt(i,j),2) - 
       32*pow(rgrid[i],4)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
       80*pow(rgrid[i],7)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
       20*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*
        pow(Qtt(i,j),2) - 24*pow(mu,2)*pow(rgrid[i],8)*
        pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
       pow(rgrid[i],2)*pow(Q22(i,j),2)*
        (-32 + 32*pow(rgrid[i],3) + 8*pow(mu,2)*pow(rgrid[i],3) - 
          4*pow(mu,2)*pow(rgrid[i],4) - 40*rgrid[i]*Qtt(i,j) + 
          88*pow(rgrid[i],4)*Qtt(i,j) + 
          22*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
          22*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
          16*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
          52*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
          13*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
          16*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
          pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
           pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
             pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
           pow(1 + rgrid[i]*Qtt(i,j),2) + 
          2*rgrid[i]*Qrr(i,j)*
           (-32 + 5*(4 + pow(mu,2))*pow(rgrid[i],3) - 
             2*rgrid[i]*(20 - 8*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                7*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j) - 
             2*pow(rgrid[i],2)*
              (8 - 5*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                6*pow(mu,2)*pow(rgrid[i],4))*pow(Qtt(i,j),2)) + 
          pow(rgrid[i],2)*pow(Qrr(i,j),2)*
           (-32 + 5*(4 + pow(mu,2))*pow(rgrid[i],3) - 
             2*rgrid[i]*(20 - 8*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                7*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j) - 
             2*pow(rgrid[i],2)*
              (8 - 5*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                6*pow(mu,2)*pow(rgrid[i],4))*pow(Qtt(i,j),2))) + 
       2*rgrid[i]*Q22(i,j)*(-36 + 48*pow(rgrid[i],3) + 
          12*pow(mu,2)*pow(rgrid[i],3) - 
          9*pow(mu,2)*pow(rgrid[i],4) - 48*rgrid[i]*Qtt(i,j) + 
          120*pow(rgrid[i],4)*Qtt(i,j) + 
          30*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
          32*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
          20*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
          68*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
          17*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
          21*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
          pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
           pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
             pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
           pow(1 + rgrid[i]*Qtt(i,j),2) - 
          2*rgrid[i]*Qrr(i,j)*
           (36 - 9*(4 + pow(mu,2))*pow(rgrid[i],3) + 
             5*pow(mu,2)*pow(rgrid[i],4) + 
             24*rgrid[i]*(2 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j) + 
             pow(rgrid[i],2)*
              (20 - 14*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                17*pow(mu,2)*pow(rgrid[i],4))*pow(Qtt(i,j),2)) + 
          pow(rgrid[i],2)*pow(Qrr(i,j),2)*
           (-36 + 9*(4 + pow(mu,2))*pow(rgrid[i],3) - 
             5*pow(mu,2)*pow(rgrid[i],4) - 
             24*rgrid[i]*(2 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j) + 
             (-20*pow(rgrid[i],2) + 
                14*(4 + pow(mu,2))*pow(rgrid[i],5) - 
                17*pow(mu,2)*pow(rgrid[i],6))*pow(Qtt(i,j),2))) + 
       pow(rgrid[i],2)*pow(Q11(i,j),2)*
        (-32 + 32*pow(rgrid[i],3) + 8*pow(mu,2)*pow(rgrid[i],3) - 
          4*pow(mu,2)*pow(rgrid[i],4) - 64*rgrid[i]*Qrr(i,j) + 
          40*pow(rgrid[i],4)*Qrr(i,j) + 
          10*pow(mu,2)*pow(rgrid[i],4)*Qrr(i,j) - 
          32*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
          20*pow(rgrid[i],5)*pow(Qrr(i,j),2) + 
          5*pow(mu,2)*pow(rgrid[i],5)*pow(Qrr(i,j),2) - 
          40*rgrid[i]*Qtt(i,j) + 88*pow(rgrid[i],4)*Qtt(i,j) + 
          22*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
          22*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
          80*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          128*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) + 
          32*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
          28*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) - 
          40*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) + 
          64*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) + 
          16*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) - 
          14*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*Qtt(i,j) - 
          16*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
          52*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
          13*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
          16*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) - 
          32*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          80*pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          20*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) - 
          24*pow(mu,2)*pow(rgrid[i],7)*Qrr(i,j)*pow(Qtt(i,j),2) - 
          16*pow(rgrid[i],4)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          40*pow(rgrid[i],7)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          10*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) - 
          12*pow(mu,2)*pow(rgrid[i],8)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) + 
          pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
           pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
             pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
           pow(1 + rgrid[i]*Qtt(i,j),2) + 
          pow(rgrid[i],2)*pow(Q22(i,j),2)*
           (-16 - 8*pow(rgrid[i],3) - 2*pow(mu,2)*pow(rgrid[i],3) + 
             8*pow(mu,2)*pow(rgrid[i],4) - 32*rgrid[i]*Qrr(i,j) - 
             40*pow(rgrid[i],4)*Qrr(i,j) - 
             10*pow(mu,2)*pow(rgrid[i],4)*Qrr(i,j) + 
             24*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j) - 
             16*pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
             20*pow(rgrid[i],5)*pow(Qrr(i,j),2) - 
             5*pow(mu,2)*pow(rgrid[i],5)*pow(Qrr(i,j),2) + 
             12*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2) - 
             8*rgrid[i]*Qtt(i,j) + 8*pow(rgrid[i],4)*Qtt(i,j) + 
             2*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) + 
             2*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
             16*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
             32*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
             8*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) + 
             20*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) - 
             8*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) - 
             16*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) - 
             4*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) + 
             10*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*
              Qtt(i,j) + 12*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
             3*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
             4*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
             pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
              pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
                pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
              pow(1 + rgrid[i]*Qtt(i,j),2)) + 
          2*rgrid[i]*Q22(i,j)*(-20 + 8*pow(rgrid[i],3) + 
             2*pow(mu,2)*pow(rgrid[i],3) + 
             3*pow(mu,2)*pow(rgrid[i],4) - 16*rgrid[i]*Qtt(i,j) + 
             40*pow(rgrid[i],4)*Qtt(i,j) + 
             10*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
             8*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
             4*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             28*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
             7*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
             9*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
             pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
              pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
                pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
              pow(1 + rgrid[i]*Qtt(i,j),2) - 
             2*rgrid[i]*Qrr(i,j)*
              (20 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                7*pow(mu,2)*pow(rgrid[i],4) - 
                4*rgrid[i]*(-4 + (4 + pow(mu,2))*pow(rgrid[i],3))*
                 Qtt(i,j) + pow(rgrid[i],2)*
                 (4 - 4*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                   5*pow(mu,2)*pow(rgrid[i],4))*pow(Qtt(i,j),2)) - 
             pow(rgrid[i],2)*pow(Qrr(i,j),2)*
              (20 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                7*pow(mu,2)*pow(rgrid[i],4) - 
                4*rgrid[i]*(-4 + (4 + pow(mu,2))*pow(rgrid[i],3))*
                 Qtt(i,j) + pow(rgrid[i],2)*
                 (4 - 4*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                   5*pow(mu,2)*pow(rgrid[i],4))*pow(Qtt(i,j),2)))) + 
       2*rgrid[i]*Q11(i,j)*(-36 + 48*pow(rgrid[i],3) + 
          12*pow(mu,2)*pow(rgrid[i],3) - 
          9*pow(mu,2)*pow(rgrid[i],4) + 
          16*pow(rgrid[i],2)*pow(Qr1(i,j),2) - 
          32*pow(rgrid[i],5)*pow(Qr1(i,j),2) - 
          8*pow(mu,2)*pow(rgrid[i],5)*pow(Qr1(i,j),2) + 
          8*pow(mu,2)*pow(rgrid[i],6)*pow(Qr1(i,j),2) + 
          16*pow(rgrid[i],8)*pow(Qr1(i,j),2) + 
          8*pow(mu,2)*pow(rgrid[i],8)*pow(Qr1(i,j),2) + 
          pow(mu,4)*pow(rgrid[i],8)*pow(Qr1(i,j),2) - 
          8*pow(mu,2)*pow(rgrid[i],9)*pow(Qr1(i,j),2) - 
          2*pow(mu,4)*pow(rgrid[i],9)*pow(Qr1(i,j),2) + 
          pow(mu,4)*pow(rgrid[i],10)*pow(Qr1(i,j),2) - 
          72*rgrid[i]*Qrr(i,j) + 72*pow(rgrid[i],4)*Qrr(i,j) + 
          18*pow(mu,2)*pow(rgrid[i],4)*Qrr(i,j) - 
          10*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j) - 
          36*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
          36*pow(rgrid[i],5)*pow(Qrr(i,j),2) + 
          9*pow(mu,2)*pow(rgrid[i],5)*pow(Qrr(i,j),2) - 
          5*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2) - 
          48*rgrid[i]*Qtt(i,j) + 120*pow(rgrid[i],4)*Qtt(i,j) + 
          30*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
          32*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) + 
          32*pow(rgrid[i],3)*pow(Qr1(i,j),2)*Qtt(i,j) - 
          64*pow(rgrid[i],6)*pow(Qr1(i,j),2)*Qtt(i,j) - 
          16*pow(mu,2)*pow(rgrid[i],6)*pow(Qr1(i,j),2)*Qtt(i,j) + 
          16*pow(mu,2)*pow(rgrid[i],7)*pow(Qr1(i,j),2)*Qtt(i,j) + 
          32*pow(rgrid[i],9)*pow(Qr1(i,j),2)*Qtt(i,j) + 
          16*pow(mu,2)*pow(rgrid[i],9)*pow(Qr1(i,j),2)*Qtt(i,j) + 
          2*pow(mu,4)*pow(rgrid[i],9)*pow(Qr1(i,j),2)*Qtt(i,j) - 
          16*pow(mu,2)*pow(rgrid[i],10)*pow(Qr1(i,j),2)*Qtt(i,j) - 
          4*pow(mu,4)*pow(rgrid[i],10)*pow(Qr1(i,j),2)*Qtt(i,j) + 
          2*pow(mu,4)*pow(rgrid[i],11)*pow(Qr1(i,j),2)*Qtt(i,j) - 
          96*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          192*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) + 
          48*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j)*Qtt(i,j) - 
          48*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*Qtt(i,j) - 
          48*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) + 
          96*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) + 
          24*pow(mu,2)*pow(rgrid[i],6)*pow(Qrr(i,j),2)*Qtt(i,j) - 
          24*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*Qtt(i,j) - 
          20*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
          68*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
          17*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
          21*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
          16*pow(rgrid[i],4)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
          32*pow(rgrid[i],7)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
          8*pow(mu,2)*pow(rgrid[i],7)*pow(Qr1(i,j),2)*
           pow(Qtt(i,j),2) + 
          8*pow(mu,2)*pow(rgrid[i],8)*pow(Qr1(i,j),2)*
           pow(Qtt(i,j),2) + 
          16*pow(rgrid[i],10)*pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
          8*pow(mu,2)*pow(rgrid[i],10)*pow(Qr1(i,j),2)*
           pow(Qtt(i,j),2) + 
          pow(mu,4)*pow(rgrid[i],10)*pow(Qr1(i,j),2)*
           pow(Qtt(i,j),2) - 
          8*pow(mu,2)*pow(rgrid[i],11)*pow(Qr1(i,j),2)*
           pow(Qtt(i,j),2) - 
          2*pow(mu,4)*pow(rgrid[i],11)*pow(Qr1(i,j),2)*
           pow(Qtt(i,j),2) + 
          pow(mu,4)*pow(rgrid[i],12)*pow(Qr1(i,j),2)*
           pow(Qtt(i,j),2) - 
          40*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          112*pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          28*pow(mu,2)*pow(rgrid[i],6)*Qrr(i,j)*pow(Qtt(i,j),2) - 
          34*pow(mu,2)*pow(rgrid[i],7)*Qrr(i,j)*pow(Qtt(i,j),2) - 
          20*pow(rgrid[i],4)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          56*pow(rgrid[i],7)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
          14*pow(mu,2)*pow(rgrid[i],7)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) - 
          17*pow(mu,2)*pow(rgrid[i],8)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) + 
          pow(rgrid[i],2)*pow(Q22(i,j),2)*
           (-20 + 8*pow(rgrid[i],3) + 2*pow(mu,2)*pow(rgrid[i],3) + 
             3*pow(mu,2)*pow(rgrid[i],4) - 16*rgrid[i]*Qtt(i,j) + 
             40*pow(rgrid[i],4)*Qtt(i,j) + 
             10*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
             8*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
             4*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             28*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
             7*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
             9*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
             pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
              pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
                pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
              pow(1 + rgrid[i]*Qtt(i,j),2) - 
             2*rgrid[i]*Qrr(i,j)*
              (20 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                7*pow(mu,2)*pow(rgrid[i],4) - 
                4*rgrid[i]*(-4 + (4 + pow(mu,2))*pow(rgrid[i],3))*
                 Qtt(i,j) + pow(rgrid[i],2)*
                 (4 - 4*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                   5*pow(mu,2)*pow(rgrid[i],4))*pow(Qtt(i,j),2)) - 
             pow(rgrid[i],2)*pow(Qrr(i,j),2)*
              (20 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                7*pow(mu,2)*pow(rgrid[i],4) - 
                4*rgrid[i]*(-4 + (4 + pow(mu,2))*pow(rgrid[i],3))*
                 Qtt(i,j) + pow(rgrid[i],2)*
                 (4 - 4*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                   5*pow(mu,2)*pow(rgrid[i],4))*pow(Qtt(i,j),2))) + 
          2*rgrid[i]*Q22(i,j)*(-24 + 24*pow(rgrid[i],3) + 
             6*pow(mu,2)*pow(rgrid[i],3) - 
             2*pow(mu,2)*pow(rgrid[i],4) - 24*rgrid[i]*Qtt(i,j) + 
             72*pow(rgrid[i],4)*Qtt(i,j) + 
             18*pow(mu,2)*pow(rgrid[i],4)*Qtt(i,j) - 
             18*pow(mu,2)*pow(rgrid[i],5)*Qtt(i,j) - 
             8*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             44*pow(rgrid[i],5)*pow(Qtt(i,j),2) + 
             11*pow(mu,2)*pow(rgrid[i],5)*pow(Qtt(i,j),2) - 
             14*pow(mu,2)*pow(rgrid[i],6)*pow(Qtt(i,j),2) + 
             pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
              pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
                pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
              pow(1 + rgrid[i]*Qtt(i,j),2) + 
             2*rgrid[i]*Qrr(i,j)*
              (-24 + 3*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                2*pow(mu,2)*pow(rgrid[i],4) - 
                2*rgrid[i]*(12 - 6*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                   5*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j) + 
                (-8*pow(rgrid[i],2) + 
                   8*(4 + pow(mu,2))*pow(rgrid[i],5) - 
                   10*pow(mu,2)*pow(rgrid[i],6))*pow(Qtt(i,j),2)) + 
             pow(rgrid[i],2)*pow(Qrr(i,j),2)*
              (-24 + 3*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                2*pow(mu,2)*pow(rgrid[i],4) - 
                2*rgrid[i]*(12 - 6*(4 + pow(mu,2))*pow(rgrid[i],3) + 
                   5*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j) + 
                (-8*pow(rgrid[i],2) + 
                   8*(4 + pow(mu,2))*pow(rgrid[i],5) - 
                   10*pow(mu,2)*pow(rgrid[i],6))*pow(Qtt(i,j),2))))))/
   (2.*(-1 + rgrid[i])*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
     pow(1 + rgrid[i]*Qrr(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2))
;
elem[4][2]+=
(Qrr.diff(d01,i,j)*Qrr.diff(d10,i,j)*pow(rgrid[i],2))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j))) - 
  (pow(Qrr.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (Qr1.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q11(i,j),2) - 
  (4*h.diff(d01,i,j)*pow(mu,2)*h(i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (4*h.diff(d01,i,j)*h.diff(d10,i,j)*pow(mu,2)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (2*Qr1.diff(d02,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) - 
  (Q22.diff(d01,i,j)*Q22.diff(d10,i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Q22(i,j),2)) - 
  (3*pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),4)) + 
  (2*Qr1.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q11(i,j),2)) + 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  ((12 - 6*(4 + pow(mu,2))*pow(rgrid[i],3) + 
       7*pow(mu,2)*pow(rgrid[i],4) + 
       rgrid[i]*(4 - 4*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          5*pow(mu,2)*pow(rgrid[i],4))*Q11(i,j))*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3)) + 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(L + L*rgrid[i]*Q22(i,j),2)) + 
  (Q22.diff(d01,i,j)*rgrid[i]*(Q22(i,j) - 2*Qrr(i,j))*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Q22(i,j),2)) + 
  Q11.diff(d10,i,j)*((-4*Qrr.diff(d01,i,j)*pow(rgrid[i],2))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3)) + 
     (3*Q11.diff(d01,i,j)*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),4)) - 
     (2*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
      pow(1 + rgrid[i]*Q11(i,j),3)) + 
  Q11.diff(d01,i,j)*((-4*Qr1.diff(d01,i,j)*L*pow(rgrid[i],2)*
        (1 + rgrid[i]*Qrr(i,j)))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),3)) \
+ ((-11 + 2*(-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2) + 
          2*rgrid[i]*Q11(i,j)*
           (-1 + (-1 + rgrid[i])*pow(rgrid[i],2)*
              (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)) - 
          6*rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),4))) + 
  Qrr.diff(d01,i,j)*((4*Qr1.diff(d01,i,j)*pow(rgrid[i],2))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) \
+ (4*Q11.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),3)) \
+ (6 + 7*rgrid[i]*Qrr(i,j) - (-1 + rgrid[i])*pow(rgrid[i],2)*
         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
         (1 + rgrid[i]*Qrr(i,j)) + 
        rgrid[i]*Q11(i,j)*(2 + 3*rgrid[i]*Qrr(i,j) - 
           (-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
            (1 + rgrid[i]*Qrr(i,j))))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3)*
        (1 + rgrid[i]*Qrr(i,j)))) - 
  (Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (8*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (8*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
  (8*a0.diff(d01,i,j)*a0.diff(d10,i,j)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     (L + L*rgrid[i]*Qtt(i,j))) + 
  (Qtt.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j))*
     ((-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          2*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j) + 
       (4 + 2*(4 + pow(mu,2))*pow(rgrid[i],3) - 
          3*pow(mu,2)*pow(rgrid[i],4))*Qtt(i,j)))/
   (L*pow(-1 + rgrid[i],2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Qtt(i,j),2))
;
elem[4][3]+=
(-2*Qr1.diff(d11,i,j)*rgrid[i])/L + 
  (2*Qr1.diff(d02,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) + 
  (pow(Qrr.diff(d01,i,j),2)*pow(rgrid[i],2))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),3)) + 
  (Q11.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q11(i,j),2) - 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (Q22.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q22(i,j),2) - 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Q22(i,j),2)) - 
  (2*Q11.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (2*Q22.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  Qrr.diff(d10,i,j)*((pow(rgrid[i],2)*
        ((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
          3*pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) + 
     (Qr1.diff(d01,i,j)*pow(rgrid[i],2))/(L + L*rgrid[i]*Qrr(i,j))) + 
  Qr1.diff(d10,i,j)*(-3*(-1 + rgrid[i])*rgrid[i]*
      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j) + 
     (Qrr.diff(d01,i,j)*pow(rgrid[i],2))/(L + L*rgrid[i]*Qrr(i,j))) + 
  Qrr.diff(d01,i,j)*((-2*Q11.diff(d01,i,j)*pow(rgrid[i],2))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) \
- (2*Qr1.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) - 
     (rgrid[i]*Qr1(i,j)*(-8 + 20*pow(rgrid[i],3) + 
          5*pow(mu,2)*pow(rgrid[i],3) - 
          6*pow(mu,2)*pow(rgrid[i],4) + 
          2*pow(rgrid[i],2)*
           pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2) + 
          pow(rgrid[i],3)*Q11(i,j)*
           ((12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
             2*pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),2)) - 
          8*rgrid[i]*Qrr(i,j) + 8*pow(rgrid[i],4)*Qrr(i,j) + 
          2*pow(mu,2)*pow(rgrid[i],4)*Qrr(i,j) - 
          2*pow(mu,2)*pow(rgrid[i],5)*Qrr(i,j)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j)))) + 
  Qr1.diff(d01,i,j)*((2*(-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2))/L + 
     (-(pow(rgrid[i],3)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
           (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
        pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
         pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
           pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
         pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j)) - 
        (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (6 + 7*rgrid[i]*Qrr(i,j) + 
           2*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
           rgrid[i]*Q11(i,j)*(4 + 3*rgrid[i]*Qrr(i,j))))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j)))) - 
  (Qtt.diff(d10,i,j)*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Qrr(i,j)))/
   (2.*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (Qtt.diff(d01,i,j)*rgrid[i]*
     (-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       2*pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (8*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Qtt(i,j),2)) - 
  (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
      pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
        pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
      pow(1 + rgrid[i]*Q22(i,j),2)*pow(Qr1(i,j),2)*
      (4 + 3*rgrid[i]*Qrr(i,j))*pow(1 + rgrid[i]*Qtt(i,j),2) + 
     2*(-1 + rgrid[i])*pow(rgrid[i],2)*
      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
      pow(1 + rgrid[i]*Q22(i,j),2)*pow(Qr1(i,j),2)*
      (2*(8 - 5*(4 + pow(mu,2))*pow(rgrid[i],3) + 
           6*pow(mu,2)*pow(rgrid[i],4)) + 
        rgrid[i]*(12 - 9*(4 + pow(mu,2))*pow(rgrid[i],3) + 
           11*pow(mu,2)*pow(rgrid[i],4))*Qrr(i,j))*
      pow(1 + rgrid[i]*Qtt(i,j),2) + 
     2*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))*
      (6 + pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         pow(Qr1(i,j),2) + 6*rgrid[i]*Qrr(i,j) + 9*rgrid[i]*Qtt(i,j) + 
        2*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         pow(Qr1(i,j),2)*Qtt(i,j) + 
        9*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
        4*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
        pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
         pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
        4*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
        pow(rgrid[i],2)*pow(Q22(i,j),2)*
         (pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
            pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) + 
           (1 + rgrid[i]*Qrr(i,j))*
            (4 + 5*rgrid[i]*Qtt(i,j) + 
              2*pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
        rgrid[i]*Q22(i,j)*(2*pow(rgrid[i],5)*
            (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qr1(i,j),2)*
            pow(1 + rgrid[i]*Qtt(i,j),2) + 
           (1 + rgrid[i]*Qrr(i,j))*
            (9 + 12*rgrid[i]*Qtt(i,j) + 
              5*pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
        pow(rgrid[i],2)*pow(Q11(i,j),2)*
         (4 + pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
            pow(Qr1(i,j),2) + 4*rgrid[i]*Qrr(i,j) + 
           5*rgrid[i]*Qtt(i,j) + 
           2*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
            pow(Qr1(i,j),2)*Qtt(i,j) + 
           5*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
           2*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
           pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
            pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
           2*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
           pow(rgrid[i],2)*pow(Q22(i,j),2)*
            (pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
               pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) + 
              (1 + rgrid[i]*Qrr(i,j))*(2 + rgrid[i]*Qtt(i,j))) + 
           rgrid[i]*Q22(i,j)*(2*pow(rgrid[i],5)*
               (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*pow(Qr1(i,j),2)*
               pow(1 + rgrid[i]*Qtt(i,j),2) + 
              (1 + rgrid[i]*Qrr(i,j))*
               (5 + 4*rgrid[i]*Qtt(i,j) + 
                 pow(rgrid[i],2)*pow(Qtt(i,j),2)))) + 
        rgrid[i]*Q11(i,j)*(9 + 
           2*pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
            pow(Qr1(i,j),2) + 9*rgrid[i]*Qrr(i,j) + 
           12*rgrid[i]*Qtt(i,j) + 
           4*pow(rgrid[i],6)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
            pow(Qr1(i,j),2)*Qtt(i,j) + 
           12*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
           5*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
           2*pow(rgrid[i],7)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
            pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
           5*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
           4*rgrid[i]*Q22(i,j)*
            (pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
               pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) + 
              (1 + rgrid[i]*Qrr(i,j))*
               (3 + 3*rgrid[i]*Qtt(i,j) + 
                 pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
           pow(rgrid[i],2)*pow(Q22(i,j),2)*
            (2*pow(rgrid[i],5)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
               pow(Qr1(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) + 
              (1 + rgrid[i]*Qrr(i,j))*
               (5 + 4*rgrid[i]*Qtt(i,j) + 
                 pow(rgrid[i],2)*pow(Qtt(i,j),2))))) - 
     pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
      (-48 - 12*pow(mu,2) + 24*pow(mu,2)*rgrid[i] - 
        24*rgrid[i]*Qrr(i,j) - 6*pow(mu,2)*rgrid[i]*Qrr(i,j) + 
        20*pow(mu,2)*pow(rgrid[i],2)*Qrr(i,j) + 
        12*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
        3*pow(mu,2)*pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
        48*rgrid[i]*Qtt(i,j) - 12*pow(mu,2)*rgrid[i]*Qtt(i,j) + 
        28*pow(mu,2)*pow(rgrid[i],2)*Qtt(i,j) + 
        48*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
        12*pow(mu,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
        72*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) + 
        18*pow(mu,2)*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) - 
        20*pow(mu,2)*pow(rgrid[i],4)*pow(Qrr(i,j),2)*Qtt(i,j) - 
        12*pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
        3*pow(mu,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
        8*pow(mu,2)*pow(rgrid[i],3)*pow(Qtt(i,j),2) + 
        48*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
        12*pow(mu,2)*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) - 
        12*pow(mu,2)*pow(rgrid[i],4)*Qrr(i,j)*pow(Qtt(i,j),2) + 
        48*pow(rgrid[i],4)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
        12*pow(mu,2)*pow(rgrid[i],4)*pow(Qrr(i,j),2)*
         pow(Qtt(i,j),2) - 16*pow(mu,2)*pow(rgrid[i],5)*
         pow(Qrr(i,j),2)*pow(Qtt(i,j),2) - 
        rgrid[i]*Q22(i,j)*(72 + 18*pow(mu,2) - 32*pow(mu,2)*rgrid[i] + 
           4*(24 + pow(mu,2)*(6 - 11*rgrid[i]))*rgrid[i]*Qtt(i,j) + 
           (36 + pow(mu,2)*(9 - 16*rgrid[i]))*pow(rgrid[i],2)*
            pow(Qtt(i,j),2) + 
           2*rgrid[i]*Qrr(i,j)*
            (9*(4 + pow(mu,2)*(1 - 2*rgrid[i])) + 
              2*(12 + pow(mu,2)*(3 - 8*rgrid[i]))*rgrid[i]*Qtt(i,j) - 
              2*pow(mu,2)*pow(rgrid[i],3)*pow(Qtt(i,j),2)) + 
           pow(rgrid[i],2)*pow(Qrr(i,j),2)*
            (12 + pow(mu,2)*(3 - 8*rgrid[i]) + 
              2*rgrid[i]*(-12 + pow(mu,2)*(-3 + 2*rgrid[i]))*Qtt(i,j) + 
              2*pow(rgrid[i],2)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*
               pow(Qtt(i,j),2))) + 
        rgrid[i]*Q11(i,j)*(-72 - 18*pow(mu,2) + 32*pow(mu,2)*rgrid[i] - 
           96*rgrid[i]*Qtt(i,j) - 24*pow(mu,2)*rgrid[i]*Qtt(i,j) + 
           44*pow(mu,2)*pow(rgrid[i],2)*Qtt(i,j) - 
           36*pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
           9*pow(mu,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
           16*pow(mu,2)*pow(rgrid[i],3)*pow(Qtt(i,j),2) + 
           pow(rgrid[i],2)*pow(Qrr(i,j),2)*
            (-12 + pow(mu,2)*(-3 + 8*rgrid[i]) + 
              2*(12 + pow(mu,2)*(3 - 2*rgrid[i]))*rgrid[i]*Qtt(i,j) + 
              2*(12 + pow(mu,2)*(3 - 4*rgrid[i]))*pow(rgrid[i],2)*
               pow(Qtt(i,j),2)) + 
           2*rgrid[i]*Qrr(i,j)*
            (9*(-4 + pow(mu,2)*(-1 + 2*rgrid[i])) + 
              2*rgrid[i]*(-12 + pow(mu,2)*(-3 + 8*rgrid[i]))*Qtt(i,j) + 
              2*pow(mu,2)*pow(rgrid[i],3)*pow(Qtt(i,j),2)) + 
           rgrid[i]*Q22(i,j)*(8*(-12 + pow(mu,2)*(-3 + 5*rgrid[i])) + 
              12*rgrid[i]*(-12 + pow(mu,2)*(-3 + 5*rgrid[i]))*Qtt(i,j) + 
              3*pow(rgrid[i],2)*(-20 + pow(mu,2)*(-5 + 8*rgrid[i]))*
               pow(Qtt(i,j),2) + 
              pow(rgrid[i],2)*pow(Qrr(i,j),2)*
               (-36 + pow(mu,2)*(-9 + 16*rgrid[i]) + 
                 6*rgrid[i]*(-4 + pow(mu,2)*(-1 + 2*rgrid[i]))*Qtt(i,j)) \
+ 2*rgrid[i]*Qrr(i,j)*(-60 + pow(mu,2)*(-15 + 26*rgrid[i]) + 
                 2*rgrid[i]*(-36 + pow(mu,2)*(-9 + 16*rgrid[i]))*
                  Qtt(i,j) + 2*pow(rgrid[i],2)*
                  (-12 + pow(mu,2)*(-3 + 5*rgrid[i]))*pow(Qtt(i,j),2))))\
))/(2.*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j))*pow(1 + rgrid[i]*Qtt(i,j),2))
;
elem[4][4]+=
(Qr1.diff(d10,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   pow(1 + rgrid[i]*Q22(i,j),2) + 
  (2*pow(Q22.diff(d01,i,j),2)*L*pow(rgrid[i],3)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Q22(i,j),3)) + 
  ((12 - 6*(4 + pow(mu,2))*pow(rgrid[i],3) + 
       7*pow(mu,2)*pow(rgrid[i],4) + 
       rgrid[i]*(4 - 4*(4 + pow(mu,2))*pow(rgrid[i],3) + 
          5*pow(mu,2)*pow(rgrid[i],4))*Q22(i,j))*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q22(i,j),3)) + 
  (Q22.diff(d01,i,j)*(-1 + rgrid[i]*Q22(i,j) + 
       2*(-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) - 4*rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Q22(i,j),3)) + 
  Q22.diff(d10,i,j)*((-2*Q22.diff(d01,i,j)*pow(rgrid[i],2)*
        (1 + rgrid[i]*Qrr(i,j)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(1 + rgrid[i]*Q22(i,j),3)) - 
     (2*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
      pow(1 + rgrid[i]*Q22(i,j),3))
;
elem[4][5]+=
(-8*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
    (1 + rgrid[i]*Qrr(i,j)))/
  (L*(-1 + rgrid[i])*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
      pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
    (1 + rgrid[i]*Qtt(i,j)))
;
elem[4][6]+=
(4*h.diff(d01,i,j)*pow(mu,2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (8*a0.diff(d01,i,j)*a0.diff(d10,i,j)*B1*pow(mu,3)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (8*a0.diff(d01,i,j)*B1*pow(mu,3)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (8*pow(a0.diff(d01,i,j),2)*B1*pow(mu,3)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[5][0]+=
-a0(i,j)/(2.*(-1 + rgrid[i])*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (a0.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (2.*L*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  Qtt.diff(d10,i,j)*((a0(i,j)*pow(rgrid[i],2))/
      (2.*(-1 + rgrid[i])*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
     (2*a0.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(2 + 2*rgrid[i]*Qtt(i,j),2) - 
     (a0.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qtt(i,j),2))) + 
  a0.diff(d10,i,j)*(-1/(2.*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
     (Qtt.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qtt(i,j),2))) + 
  Qtt.diff(d01,i,j)*(-(a0(i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*L*(-1 + rgrid[i])*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
     (a0.diff(d01,i,j)*pow(rgrid[i],2)*
        (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(L + L*rgrid[i]*Qtt(i,j),2)))
;
elem[5][1]+=
(2*a0.diff(d02,i,j)*rgrid[i])/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (a0.diff(d01,i,j)*Q11.diff(d01,i,j)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (a0.diff(d01,i,j)*Q22.diff(d01,i,j)*pow(rgrid[i],2))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))) - 
  a0(i,j)/(2.*(-1 + rgrid[i])*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  a0.diff(d01,i,j)*((2*h.diff(d01,i,j)*B1*mu*pow(rgrid[i],2))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
     (rgrid[i]*Qr1(i,j))/(2.*L*pow(1 + rgrid[i]*Qrr(i,j),2))) + 
  Qrr.diff(d10,i,j)*((a0(i,j)*pow(rgrid[i],2))/
      (2.*(-1 + rgrid[i])*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
     (2*a0.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(2 + 2*rgrid[i]*Qrr(i,j),2) - 
     (a0.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qrr(i,j),2))) + 
  a0.diff(d10,i,j)*(-1/(2.*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
     (Qrr.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qrr(i,j),2))) + 
  Qrr.diff(d01,i,j)*(-(a0(i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*L*(-1 + rgrid[i])*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
     (a0.diff(d01,i,j)*pow(rgrid[i],4)*pow(Qr1(i,j),2))/
      (2.*pow(L + L*rgrid[i]*Qrr(i,j),2))) - 
  (a0.diff(d01,i,j)*Qtt.diff(d01,i,j)*pow(rgrid[i],2))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[5][2]+=
a0(i,j)/(2.*(-1 + rgrid[i])*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (a0.diff(d01,i,j)*Qrr.diff(d01,i,j)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  Q11.diff(d10,i,j)*(-(a0(i,j)*pow(rgrid[i],2))/
      (2.*(-1 + rgrid[i])*pow(1 + rgrid[i]*Q11(i,j),2)) - 
     (2*a0.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(2 + 2*rgrid[i]*Q11(i,j),2) + 
     (a0.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q11(i,j),2))) + 
  a0.diff(d10,i,j)*(1/(2.*pow(1 + rgrid[i]*Q11(i,j),2)) + 
     (Q11.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q11(i,j),2))) - 
  (2*a0.diff(d02,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) - 
  (a0.diff(d01,i,j)*Q22.diff(d01,i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Q22(i,j))) + 
  a0.diff(d01,i,j)*(-(rgrid[i]*Qr1(i,j))/
      (2.*L*pow(1 + rgrid[i]*Q11(i,j),2)) - 
     (2*h.diff(d01,i,j)*B1*mu*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
      ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2))) \
+ Q11.diff(d01,i,j)*((a0(i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*L*(-1 + rgrid[i])*pow(1 + rgrid[i]*Q11(i,j),2)) + 
     (a0.diff(d01,i,j)*pow(rgrid[i],2)*
        (4 - (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 4*rgrid[i]*Qrr(i,j)))/
      (2.*pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),3))) + 
  (a0.diff(d01,i,j)*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[5][3]+=
(-2*a0.diff(d11,i,j)*rgrid[i])/L + 
  (2*a0.diff(d01,i,j)*Qr1.diff(d01,i,j)*pow(rgrid[i],2))/pow(L,2) - 
  (a0.diff(d01,i,j)*h.diff(d10,i,j)*B1*mu*pow(rgrid[i],2))/L + 
  (h.diff(d01,i,j)*B1*mu*a0(i,j)*pow(rgrid[i],2))/(L - L*rgrid[i]) - 
  (a0.diff(d01,i,j)*Q11.diff(d10,i,j)*pow(rgrid[i],2))/
   (2.*(L + L*rgrid[i]*Q11(i,j))) - 
  (a0.diff(d01,i,j)*Q22.diff(d10,i,j)*pow(rgrid[i],2))/
   (2.*(L + L*rgrid[i]*Q22(i,j))) + 
  (2*a0.diff(d02,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) + 
  Q11.diff(d01,i,j)*(-(a0(i,j)*pow(rgrid[i],2))/
      (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))) + 
     (a0.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Q11(i,j)))) + 
  Q22.diff(d01,i,j)*(-(a0(i,j)*pow(rgrid[i],2))/
      (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q22(i,j))) + 
     (a0.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Q22(i,j)))) + 
  (a0.diff(d01,i,j)*Qrr.diff(d10,i,j)*pow(rgrid[i],2))/
   (2*L + 2*L*rgrid[i]*Qrr(i,j)) + 
  Qrr.diff(d01,i,j)*((a0(i,j)*pow(rgrid[i],2))/
      (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j))) - 
     (a0.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Qrr(i,j)))) + 
  (a0.diff(d01,i,j)*Qtt.diff(d10,i,j)*pow(rgrid[i],2))/
   (2*L + 2*L*rgrid[i]*Qtt(i,j)) + 
  Qtt.diff(d01,i,j)*((a0(i,j)*pow(rgrid[i],2))/
      (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Qtt(i,j))) - 
     (a0.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Qtt(i,j)))) + 
  a0.diff(d10,i,j)*(-((h.diff(d01,i,j)*B1*mu*pow(rgrid[i],2))/L) - 
     (Q11.diff(d01,i,j)*pow(rgrid[i],2))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) - 
     (Q22.diff(d01,i,j)*pow(rgrid[i],2))/
      (2.*(L + L*rgrid[i]*Q22(i,j))) + 
     (Qrr.diff(d01,i,j)*pow(rgrid[i],2))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j)) + 
     (Qtt.diff(d01,i,j)*pow(rgrid[i],2))/(2*L + 2*L*rgrid[i]*Qtt(i,j))) \
+ a0.diff(d01,i,j)*((2*h.diff(d01,i,j)*B1*mu*pow(rgrid[i],3)*
        Qr1(i,j))/pow(L,2) + 
     (2 - 6*rgrid[i] + 3*rgrid[i]*Q22(i,j) - 7*pow(rgrid[i],2)*Q22(i,j) + 
        rgrid[i]*Qrr(i,j) - 5*pow(rgrid[i],2)*Qrr(i,j) + 
        2*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) - 
        6*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j) + rgrid[i]*Qtt(i,j) - 
        5*pow(rgrid[i],2)*Qtt(i,j) + 
        2*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
        6*pow(rgrid[i],3)*Q22(i,j)*Qtt(i,j) - 
        4*pow(rgrid[i],3)*Qrr(i,j)*Qtt(i,j) + 
        pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
        5*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
        2*B1*mu*h(i,j)*(-1 + rgrid[i])*rgrid[i]*(1 + rgrid[i]*Q11(i,j))*
         (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
         (1 + rgrid[i]*Qtt(i,j)) + 
        rgrid[i]*Q11(i,j)*(3 - 7*rgrid[i] + 2*rgrid[i]*Qtt(i,j) - 
           6*pow(rgrid[i],2)*Qtt(i,j) + 
           rgrid[i]*Qrr(i,j)*(2 - 6*rgrid[i] + 
              (rgrid[i] - 5*pow(rgrid[i],2))*Qtt(i,j)) + 
           rgrid[i]*Q22(i,j)*(4 - 8*rgrid[i] + 
              (3 - 7*rgrid[i])*rgrid[i]*Qtt(i,j) + 
              rgrid[i]*Qrr(i,j)*
               (3 - 7*rgrid[i] + 2*(1 - 3*rgrid[i])*rgrid[i]*Qtt(i,j)))))/
      (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))
;
elem[5][4]+=
a0(i,j)/(2.*(-1 + rgrid[i])*pow(1 + rgrid[i]*Q22(i,j),2)) - 
  (a0.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (2.*L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  Q22.diff(d10,i,j)*(-(a0(i,j)*pow(rgrid[i],2))/
      (2.*(-1 + rgrid[i])*pow(1 + rgrid[i]*Q22(i,j),2)) - 
     (2*a0.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(2 + 2*rgrid[i]*Q22(i,j),2) + 
     (a0.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q22(i,j),2))) + 
  a0.diff(d10,i,j)*(1/(2.*pow(1 + rgrid[i]*Q22(i,j),2)) + 
     (Q22.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q22(i,j),2))) + 
  Q22.diff(d01,i,j)*((a0(i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*L*(-1 + rgrid[i])*pow(1 + rgrid[i]*Q22(i,j),2)) - 
     (a0.diff(d01,i,j)*pow(rgrid[i],2)*
        (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(L + L*rgrid[i]*Q22(i,j),2)))
;
elem[5][5]+=
(h.diff(d10,i,j)*B1*mu*rgrid[i])/(-1 + rgrid[i]) + 
  (Qr1.diff(d01,i,j)*rgrid[i])/(L - L*rgrid[i]) + 
  (Q11.diff(d10,i,j)*rgrid[i])/
   (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))) + 
  (Q22.diff(d10,i,j)*rgrid[i])/
   (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q22(i,j))) + 
  (h.diff(d01,i,j)*B1*mu*pow(rgrid[i],2)*Qr1(i,j))/(L - L*rgrid[i]) - 
  (Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
   (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))) - 
  (Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
   (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q22(i,j))) - 
  (Qrr.diff(d10,i,j)*rgrid[i])/
   (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j))) + 
  (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
   (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j))) - 
  (Qtt.diff(d10,i,j)*rgrid[i])/
   (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Qtt(i,j))) + 
  (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
   (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Qtt(i,j))) - 
  (-Q22(i,j) + Qrr(i,j) + Qtt(i,j) + 2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
     pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
     2*B1*mu*h(i,j)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
      (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) - 
     Q11(i,j)*(1 - pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
        rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j))))/
   (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[5][6]+=
a0.diff(d10,i,j)*B1*mu + (B1*mu*a0(i,j))/(-1 + rgrid[i]) - 
  (a0.diff(d01,i,j)*B1*mu*rgrid[i]*Qr1(i,j))/L
;
elem[6][0]+=
-((pow(a0.diff(d10,i,j),2)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
       (-1 + rgrid[i])*pow(rgrid[i],2))/
     ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
         pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
  (-2*B1*mu*pow(a0(i,j),2)*exp(B1*mu*h(i,j)*rgrid[i])*
      pow(rgrid[i],3) + h(i,j)*
      (4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
        pow(mu,2)*pow(rgrid[i],4)))/
   (2.*(-1 + rgrid[i])*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (h.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (2.*L*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (2*a0.diff(d01,i,j)*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],3)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (pow(a0.diff(d01,i,j),2)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  a0.diff(d10,i,j)*((-2*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        pow(rgrid[i],2))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
     (2*a0.diff(d01,i,j)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
  Qtt.diff(d10,i,j)*((-2*h(i,j)*rgrid[i])/
      pow(2 + 2*rgrid[i]*Qtt(i,j),2) - 
     (2*h.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(2 + 2*rgrid[i]*Qtt(i,j),2) + 
     (h.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qtt(i,j),2))) + 
  h.diff(d10,i,j)*(1/(2.*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
     (Qtt.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qtt(i,j),2))) + 
  Qtt.diff(d01,i,j)*((L*h(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qtt(i,j),2)) - 
     (h.diff(d01,i,j)*pow(rgrid[i],2)*
        (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(L + L*rgrid[i]*Qtt(i,j),2)))
;
elem[6][1]+=
(2*h.diff(d02,i,j)*rgrid[i])/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (h.diff(d01,i,j)*Q11.diff(d01,i,j)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (h.diff(d01,i,j)*Q22.diff(d01,i,j)*pow(rgrid[i],2))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))) + 
  (h.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (2.*L*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     (-((2 + 3*pow(A1,2))*mu*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*h(i,j)*
          (-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       48*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)))*
        pow(1 + rgrid[i]*Qrr(i,j),2)))/
   (2.*(2 + 3*pow(A1,2))*mu*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  Qrr.diff(d10,i,j)*((2*h(i,j)*rgrid[i])/
      pow(2 + 2*rgrid[i]*Qrr(i,j),2) + 
     (2*h.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(2 + 2*rgrid[i]*Qrr(i,j),2) - 
     (h.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qrr(i,j),2))) + 
  h.diff(d10,i,j)*(-1/(2.*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
     (Qrr.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qrr(i,j),2))) + 
  Qrr.diff(d01,i,j)*(-(L*h(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Qrr(i,j),2)) + 
     (h.diff(d01,i,j)*pow(rgrid[i],4)*pow(Qr1(i,j),2))/
      (2.*pow(L + L*rgrid[i]*Qrr(i,j),2))) + 
  (2*pow(a0.diff(d01,i,j),2)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (h.diff(d01,i,j)*Qtt.diff(d01,i,j)*pow(rgrid[i],2))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[6][2]+=
h(i,j)/(2.*rgrid[i]*pow(1 + rgrid[i]*Q11(i,j),2)) - 
  (h.diff(d01,i,j)*Qrr.diff(d01,i,j)*pow(rgrid[i],2))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) - 
  (h.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (2.*L*pow(1 + rgrid[i]*Q11(i,j),2)) + 
  Q11.diff(d10,i,j)*((-2*h(i,j)*rgrid[i])/
      pow(2 + 2*rgrid[i]*Q11(i,j),2) - 
     (2*h.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(2 + 2*rgrid[i]*Q11(i,j),2) + 
     (h.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q11(i,j),2))) + 
  h.diff(d10,i,j)*(1/(2.*pow(1 + rgrid[i]*Q11(i,j),2)) + 
     (Q11.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q11(i,j),2))) - 
  (2*h.diff(d02,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4))*pow(L + L*rgrid[i]*Q11(i,j),2)) - 
  (h.diff(d01,i,j)*Q22.diff(d01,i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Q22(i,j))) + 
  Q11.diff(d01,i,j)*((L*h(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
     (h.diff(d01,i,j)*pow(rgrid[i],2)*
        (4 - (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 4*rgrid[i]*Qrr(i,j)))/
      (2.*pow(L,2)*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
          pow(mu,2)*pow(rgrid[i],4))*pow(1 + rgrid[i]*Q11(i,j),3))) - 
  (2*pow(a0.diff(d01,i,j),2)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (h.diff(d01,i,j)*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j)))
;
elem[6][3]+=
(-2*h.diff(d11,i,j)*rgrid[i])/L + 
  (2*h.diff(d01,i,j)*Qr1.diff(d01,i,j)*pow(rgrid[i],2))/pow(L,2) - 
  (h.diff(d01,i,j)*Q11.diff(d10,i,j)*pow(rgrid[i],2))/
   (2.*(L + L*rgrid[i]*Q11(i,j))) - 
  (h.diff(d01,i,j)*Q22.diff(d10,i,j)*pow(rgrid[i],2))/
   (2.*(L + L*rgrid[i]*Q22(i,j))) + 
  (2*h.diff(d02,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2) + 
  Q11.diff(d01,i,j)*(-(h(i,j)*rgrid[i])/(2.*(L + L*rgrid[i]*Q11(i,j))) + 
     (h.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Q11(i,j)))) + 
  Q22.diff(d01,i,j)*(-(h(i,j)*rgrid[i])/(2.*(L + L*rgrid[i]*Q22(i,j))) + 
     (h.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Q22(i,j)))) + 
  (h.diff(d01,i,j)*Qrr.diff(d10,i,j)*pow(rgrid[i],2))/
   (2*L + 2*L*rgrid[i]*Qrr(i,j)) + 
  Qrr.diff(d01,i,j)*(-((h.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
        (pow(L,2)*(1 + rgrid[i]*Qrr(i,j)))) + 
     (h(i,j)*rgrid[i])/(2*L + 2*L*rgrid[i]*Qrr(i,j))) - 
  (2*a0.diff(d01,i,j)*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (2*a0.diff(d01,i,j)*a0.diff(d10,i,j)*B1*mu*
     exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (2*pow(a0.diff(d01,i,j),2)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],3)*Qr1(i,j))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (h.diff(d01,i,j)*Qtt.diff(d10,i,j)*pow(rgrid[i],2))/
   (2.*(L + L*rgrid[i]*Qtt(i,j))) + 
  Qtt.diff(d01,i,j)*((h.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Qtt(i,j))) - 
     (h(i,j)*rgrid[i])/(2.*(L + L*rgrid[i]*Qtt(i,j)))) + 
  h.diff(d10,i,j)*(-(Q11.diff(d01,i,j)*pow(rgrid[i],2))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) - 
     (Q22.diff(d01,i,j)*pow(rgrid[i],2))/
      (2.*(L + L*rgrid[i]*Q22(i,j))) + 
     (Qrr.diff(d01,i,j)*pow(rgrid[i],2))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j)) - 
     (Qtt.diff(d01,i,j)*pow(rgrid[i],2))/(2.*(L + L*rgrid[i]*Qtt(i,j)))) \
- (h.diff(d01,i,j)*(pow(rgrid[i],3)*
        (-12 + pow(mu,2)*(-3 + 4*rgrid[i]))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j)) + 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (2 + rgrid[i]*Qrr(i,j) + 3*rgrid[i]*Qtt(i,j) + 
            2*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
            rgrid[i]*Q22(i,j)*
             (3 + 4*rgrid[i]*Qtt(i,j) + 
               rgrid[i]*Qrr(i,j)*(2 + 3*rgrid[i]*Qtt(i,j))) + 
            rgrid[i]*Q11(i,j)*(3 + 4*rgrid[i]*Qtt(i,j) + 
               rgrid[i]*Qrr(i,j)*(2 + 3*rgrid[i]*Qtt(i,j)) + 
               rgrid[i]*Q22(i,j)*
                (4 + 5*rgrid[i]*Qtt(i,j) + 
                  rgrid[i]*Qrr(i,j)*(3 + 4*rgrid[i]*Qtt(i,j))))))/2.))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[6][4]+=
h(i,j)/(2.*rgrid[i]*pow(1 + rgrid[i]*Q22(i,j),2)) - 
  (h.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
   (2.*L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  Q22.diff(d10,i,j)*((-2*h(i,j)*rgrid[i])/
      pow(2 + 2*rgrid[i]*Q22(i,j),2) - 
     (2*h.diff(d10,i,j)*pow(rgrid[i],2))/
      pow(2 + 2*rgrid[i]*Q22(i,j),2) + 
     (h.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q22(i,j),2))) + 
  h.diff(d10,i,j)*(1/(2.*pow(1 + rgrid[i]*Q22(i,j),2)) + 
     (Q22.diff(d01,i,j)*L*pow(rgrid[i],3)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q22(i,j),2))) + 
  Q22.diff(d01,i,j)*((L*h(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*pow(L + L*rgrid[i]*Q22(i,j),2)) - 
     (h.diff(d01,i,j)*pow(rgrid[i],2)*
        (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
      (2.*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(L + L*rgrid[i]*Q22(i,j),2)))
;
elem[6][5]+=
(2*a0.diff(d10,i,j)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i])/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (2*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i])/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (2*a0.diff(d01,i,j)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     Qr1(i,j))/(L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[6][6]+=
-(Qr1.diff(d01,i,j)/L) + Q11.diff(d10,i,j)/(2 + 2*rgrid[i]*Q11(i,j)) + 
  Q22.diff(d10,i,j)/(2 + 2*rgrid[i]*Q22(i,j)) - 
  (Q11.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/(2.*(L + L*rgrid[i]*Q11(i,j))) - 
  (Q22.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/(2.*(L + L*rgrid[i]*Q22(i,j))) - 
  Qrr.diff(d10,i,j)/(2 + 2*rgrid[i]*Qrr(i,j)) + 
  (Qrr.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/(2*L + 2*L*rgrid[i]*Qrr(i,j)) + 
  (pow(a0.diff(d10,i,j),2)*pow(B1,2)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (2*a0.diff(d01,i,j)*pow(B1,2)*pow(mu,2)*a0(i,j)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],3)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (pow(a0.diff(d01,i,j),2)*pow(B1,2)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (2 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 2*rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  Qtt.diff(d10,i,j)/(2 + 2*rgrid[i]*Qtt(i,j)) - 
  (Qtt.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/(2.*(L + L*rgrid[i]*Qtt(i,j))) + 
  a0.diff(d10,i,j)*((2*pow(B1,2)*pow(mu,2)*a0(i,j)*
        exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
     (2*a0.diff(d01,i,j)*pow(B1,2)*pow(mu,2)*
        exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],3)*
        Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  ((-4*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
           exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4) + 
          2*(12*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                  (3.*A1))) + 
             ((2 + 3*pow(A1,2))*mu*
                exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*h(i,j)*
                pow(rgrid[i],4)*(-12 + pow(mu,2)*(-3 + 4*rgrid[i]))\
)/2. + 12*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                  (3.*A1)))*rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))/
      (A1*(2 + 3*pow(A1,2))*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))) + 
     (2*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        (B1*(2 + 3*A1*B1)*pow(mu,2)*pow(a0(i,j),2)*
           exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4) + 
          exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           (2*mu*h(i,j)*pow(rgrid[i],4)*
              (-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
             3*A1*(8*exp(A1*mu*h(i,j)*rgrid[i]) + 
                pow(rgrid[i],3)*
                 (-12 + pow(mu,2)*(-3 + 4*rgrid[i])) + 
                8*exp(A1*mu*h(i,j)*rgrid[i])*rgrid[i]*Qrr(i,j)))*
           (1 + rgrid[i]*Qtt(i,j))))/
      (A1*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))) - 
     (3*(4 + 5*rgrid[i]*Qrr(i,j) + 3*rgrid[i]*Qtt(i,j) + 
          4*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(3 + 2*rgrid[i]*Qtt(i,j) + 
             rgrid[i]*Qrr(i,j)*(4 + 3*rgrid[i]*Qtt(i,j))) + 
          rgrid[i]*Q11(i,j)*(3 + 2*rgrid[i]*Qtt(i,j) + 
             rgrid[i]*Qrr(i,j)*(4 + 3*rgrid[i]*Qtt(i,j)) + 
             rgrid[i]*Q22(i,j)*
              (2 + rgrid[i]*Qtt(i,j) + 
                rgrid[i]*Qrr(i,j)*(3 + 2*rgrid[i]*Qtt(i,j))))))/
      ((1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (1 + rgrid[i]*Qrr(i,j))))/
   (6.*pow(rgrid[i],2)*(1 + rgrid[i]*Qtt(i,j)))
;
}
}
