
if(rgrid.IsRightBoundary(i)){
elem[0]=
Qtt(i,j)-Qtt(i+1,j)
;
elem[1]=
Qrr(i,j)-Qrr(i+1,j)
;
elem[2]=
Q11(i,j)-Q11(i+1,j)
;
elem[3]=
Qr1(i,j)-Qr1(i+1,j)
;
elem[4]=
Q22(i,j)-Q22(i+1,j)
;
elem[5]=
a0(i,j)-a0(i+1,j)
;
elem[6]=
h(i,j)-h(i+1,j)
;
}
else if(rgrid.IsLeftBoundary(i)){
elem[0]=
Qtt.diff(d10,i-1,j)-Qtt.diff(d10,i,j)
;
elem[1]=
Qrr.diff(d10,i-1,j)-Qrr.diff(d10,i,j)
;
elem[2]=
Q11.diff(d10,i-1,j)-Q11.diff(d10,i,j)
;
elem[3]=
Qr1.diff(d10,i-1,j)-Qr1.diff(d10,i,j)
;
elem[4]=
Q22.diff(d10,i-1,j)-Q22.diff(d10,i,j)
;
elem[5]=
a0.diff(d10,i-1,j)-a0.diff(d10,i,j)
;
elem[6]=
h.diff(d10,i-1,j)-h.diff(d10,i,j)
;
}
else{
elem[0]=
Qtt.diff(d20,i,j) - (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*rgrid[i])/
   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3)) - 
  (2*Qtt.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))) + 
  a0.diff(d10,i,j)*((-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        rgrid[i])/
      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3)) + 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],2)*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3)))) - 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     rgrid[i]*(1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))) + 
  (2*Qtt.diff(d02,i,j)*(1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (pow(Qtt.diff(d10,i,j),2)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j)) - 
  (2*pow(Qtt.diff(d01,i,j),2)*rgrid[i]*
     (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (2*Qtt.diff(d01,i,j)*Qr1(i,j)*(2 + rgrid[i]*Qtt(i,j)))/
   (L + L*rgrid[i]*Qtt(i,j)) + 
  Qtt.diff(d10,i,j)*((2*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qtt(i,j)) + 
     (-(pow(-1 + rgrid[i],2)*rgrid[i]*
            pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
            (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
            (1 + rgrid[i]*Qtt(i,j)))/2. + 
        (((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
                3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))/2.)*(1 + rgrid[i]*Q11(i,j))*
         (1 + rgrid[i]*Q22(i,j))*
         (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
        (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (3*Qrr(i,j) + Qtt(i,j) + 2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           Q22(i,j)*(-1 + rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
           Q11(i,j)*(-1 + rgrid[i]*Q22(i,j)*
               (-2 + rgrid[i]*Qrr(i,j) - rgrid[i]*Qtt(i,j)) + 
              rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j)))))/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (2*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        pow(rgrid[i],3)*pow(((-1 + rgrid[i])*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))/2.,2)*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(Qrr(i,j) - Qtt(i,j))*
        (1 + rgrid[i]*Qtt(i,j)) + 
       ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
          pow(-1 + rgrid[i],3)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),3)*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
          (2 + 3*rgrid[i]*Qtt(i,j) + pow(rgrid[i],2)*pow(Qtt(i,j),2)))/
        4. - ((-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (72*pow(A1,2) + 48*
             exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
            16*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.) + 
            24*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.) + 
            72*pow(A1,2)*rgrid[i]*Q22(i,j) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             rgrid[i]*Q22(i,j) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j) - 6*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j) + 12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j) + 
            18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j) + 
            72*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             rgrid[i]*Qrr(i,j) + 
            16*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
            24*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
            72*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + 
            12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j) + 
            18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j) + 
            144*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             rgrid[i]*Qtt(i,j) - 
            8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Qtt(i,j) - 12*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Qtt(i,j) + 20*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
            30*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
            144*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
            8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*Qtt(i,j) - 
            12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*Qtt(i,j) + 
            12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qtt(i,j) + 
            18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qtt(i,j) + 
            144*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
            22*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*Qtt(i,j) + 
            33*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*Qtt(i,j) + 
            144*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*
             Qtt(i,j) + 96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
               (3.*A1))*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
            14*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j)*
             Qtt(i,j) + 21*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j)*
             Qtt(i,j) + 72*pow(A1,2)*pow(rgrid[i],2)*
             pow(Qtt(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             pow(Qtt(i,j),2) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             pow(Qtt(i,j),2) + 
            6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qtt(i,j),2) + 
            9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qtt(i,j),2) + 
            72*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*pow(Qtt(i,j),2) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*pow(Qtt(i,j),2) + 
            2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qtt(i,j),2) + 
            3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qtt(i,j),2) + 
            72*pow(A1,2)*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
            8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            72*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
             exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
             (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
            rgrid[i]*Q11(i,j)*
             (72*pow(A1,2) + 
               48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) - 
               4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                  ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) - 
               6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],2)*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                  ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
               12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.) + 
               18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                rgrid[i]*(((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.) + 
               72*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
               48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                rgrid[i]*Qrr(i,j) + 
               12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
               18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],2)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
               144*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
               96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                rgrid[i]*Qtt(i,j) - 
               8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                  ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
                Qtt(i,j) - 12*pow(A1,2)*
                exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                  ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
                Qtt(i,j) + 12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],2)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
               18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],2)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
               144*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
               96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
               14*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*Qtt(i,j) \
+ 21*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*Qtt(i,j) \
+ 72*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
               48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
               4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                  ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
                pow(Qtt(i,j),2) - 
               6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],4)*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                  ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
                pow(Qtt(i,j),2) + 
               2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qtt(i,j),2) \
+ 3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qtt(i,j),2) \
+ 72*pow(A1,2)*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
               48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
               4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*
                pow(Qtt(i,j),2) + 
               6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                pow(rgrid[i],4)*
                (((-1 + rgrid[i])*
                     (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*
                pow(Qtt(i,j),2) + 
               rgrid[i]*Q22(i,j)*
                (72*pow(A1,2) + 
                  48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                     (3.*A1)) - 
                  4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],2)*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                     ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) \
- 6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                     ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) \
+ 8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                  12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   rgrid[i]*(((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                  2*rgrid[i]*(72*pow(A1,2) + 
                     48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)) - 
                     2*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
                     (2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*Qtt(i,j) - 
                  pow(rgrid[i],2)*
                   (-24*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1))) + 
                     2*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
                     (2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*
                   pow(Qtt(i,j),2) + 
                  rgrid[i]*Qrr(i,j)*
                   (72*pow(A1,2) + 
                     48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)) + 
                     4*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))/2.) + 
                     3*rgrid[i]*
                      (48*pow(A1,2) + 
                        32*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                        rgrid[i])/(3.*A1)) + 
                        (2 + 3*pow(A1,2))*
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*Qtt(i,j) + 
                     24*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*pow(rgrid[i],2)*pow(Qtt(i,j),2)))))\
)/2. - ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
          pow(-1 + rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          (-6 + pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) - 
            6*rgrid[i]*Qrr(i,j) - 8*rgrid[i]*Qtt(i,j) + 
            2*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             Qtt(i,j) - 7*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
            2*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
            pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             pow(Qtt(i,j),2) - 
            2*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
            rgrid[i]*Q22(i,j)*
             (-(rgrid[i]*Qrr(i,j)*pow(2 + rgrid[i]*Qtt(i,j),2)) + 
               (1 + rgrid[i]*Qtt(i,j))*
                (-4 + pow(rgrid[i],3)*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                   pow(Qr1(i,j),2) + 
                  rgrid[i]*(-1 + 
                     pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*Qtt(i,j))) + 
            rgrid[i]*Q11(i,j)*(-(rgrid[i]*Qrr(i,j)*
                  pow(2 + rgrid[i]*Qtt(i,j),2)) + 
               (1 + rgrid[i]*Qtt(i,j))*
                (-4 + pow(rgrid[i],3)*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                   pow(Qr1(i,j),2) + 
                  rgrid[i]*(-1 + 
                     pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*Qtt(i,j)) + 
               rgrid[i]*Q22(i,j)*
                (-(rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
                  (1 + rgrid[i]*Qtt(i,j))*
                   (-2 + pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                         (-4 - 8*rgrid[i] + 
                         3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))))))/2.))/
   ((2 + 3*pow(A1,2))*pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))
;
elem[1]=
Qrr.diff(d20,i,j) - (2*Qrr.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  (2*pow(Qr1.diff(d01,i,j),2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   pow(L,2) + 2*pow(h.diff(d10,i,j),2)*pow(mu,2)*rgrid[i]*
   (1 + rgrid[i]*Qrr(i,j)) + (pow(Q11.diff(d10,i,j),2)*rgrid[i]*
     (1 + rgrid[i]*Qrr(i,j)))/(2.*pow(1 + rgrid[i]*Q11(i,j),2)) + 
  (pow(Q22.diff(d10,i,j),2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   (2.*pow(1 + rgrid[i]*Q22(i,j),2)) - 
  (4*h.diff(d01,i,j)*pow(mu,2)*h(i,j)*rgrid[i]*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/L - 
  2*Qr1.diff(d10,i,j)*(-1 + rgrid[i])*
   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
   Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)) + 
  (2*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(rgrid[i],3)*
     pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j)))/pow(L,2) + 
  (pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],3)*pow(Qr1(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j)))/(2.*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],3)*pow(Qr1(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j)))/(2.*pow(L + L*rgrid[i]*Q22(i,j),2)) + 
  (2*Qr1.diff(d01,i,j)*(2 + 
       (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2) + 
       Q11(i,j)*(rgrid[i] + (-1 + rgrid[i])*pow(rgrid[i],3)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)))*
     (1 + rgrid[i]*Qrr(i,j)))/(L*rgrid[i]*(1 + rgrid[i]*Q11(i,j))) + 
  (Q22.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*(Q22(i,j) - 2*Qrr(i,j))*
     (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (2*Qrr.diff(d02,i,j)*(1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (3*pow(Qrr.diff(d10,i,j),2)*rgrid[i])/(2 + 2*rgrid[i]*Qrr(i,j)) - 
  (pow(Qrr.diff(d01,i,j),2)*rgrid[i]*
     (2 + (3*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + 2*rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  h.diff(d10,i,j)*(4*pow(mu,2)*h(i,j)*(1 + rgrid[i]*Qrr(i,j)) - 
     (4*h.diff(d01,i,j)*pow(mu,2)*pow(rgrid[i],2)*Qr1(i,j)*
        (1 + rgrid[i]*Qrr(i,j)))/L) - 
  (Qrr.diff(d01,i,j)*Qr1(i,j)*
     (4 + rgrid[i]*Qrr(i,j) + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
        (1 + rgrid[i]*Qrr(i,j))))/(L + L*rgrid[i]*Qrr(i,j)) + 
  Q11.diff(d10,i,j)*((-2*Qr1.diff(d01,i,j)*rgrid[i]*
        (1 + rgrid[i]*Qrr(i,j)))/(L + L*rgrid[i]*Q11(i,j)) - 
     (Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
        (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q11(i,j),2)) - 
     ((Q11(i,j) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
      pow(1 + rgrid[i]*Q11(i,j),2)) + 
  Q11.diff(d01,i,j)*((2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
        (1 + rgrid[i]*Qrr(i,j)))/(pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
     (rgrid[i]*Qr1(i,j)*(Q11(i,j) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
      (L*pow(1 + rgrid[i]*Q11(i,j),2))) + 
  Q22.diff(d10,i,j)*(-((Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
          (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),2))) - 
     ((Q22(i,j) - 2*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
      pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (pow(Qtt.diff(d10,i,j),2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   (2.*pow(1 + rgrid[i]*Qtt(i,j),2)) - 
  (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     rgrid[i]*(-1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. - rgrid[i]*Qrr(i,j))*
     (1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],3)*pow(Qr1(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j)))/(2.*pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
  (2*Qtt.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
     (rgrid[i]*(((-1 + rgrid[i])*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))/2.)*(Qrr(i,j) - Qtt(i,j)) + 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(-2*Qrr(i,j) + Qtt(i,j)))/2.))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  a0.diff(d10,i,j)*((-4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],2)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  Qtt.diff(d10,i,j)*(-((Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)*
          (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
     (2*(1 + rgrid[i]*Qrr(i,j))*
        (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(2*Qrr(i,j) - Qtt(i,j)))/2. \
+ rgrid[i]*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*(-Qrr(i,j) + Qtt(i,j))))/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
  Qrr.diff(d10,i,j)*((3*Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qrr(i,j)) + 
     ((pow(-1 + rgrid[i],2)*rgrid[i]*
           pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
           (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/2. + 
        (((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
                3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))/2.)*(1 + rgrid[i]*Q11(i,j))*
         (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
         (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
        (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (6*Qrr(i,j) + 3*rgrid[i]*pow(Qrr(i,j),2) - Qtt(i,j) + 
           4*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           2*pow(rgrid[i],2)*pow(Qrr(i,j),2)*Qtt(i,j) + 
           Q22(i,j)*(-1 - 2*rgrid[i]*Qtt(i,j) + 
              2*rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
              pow(rgrid[i],2)*pow(Qrr(i,j),2)*(2 + rgrid[i]*Qtt(i,j))) \
+ Q11(i,j)*(-1 - 2*rgrid[i]*Qtt(i,j) + 
              rgrid[i]*Q22(i,j)*
               (-2 + 2*rgrid[i]*Qrr(i,j) + 
                 pow(rgrid[i],2)*pow(Qrr(i,j),2) - 3*rgrid[i]*Qtt(i,j)\
) + 2*rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j)) + 
              pow(rgrid[i],2)*pow(Qrr(i,j),2)*(2 + rgrid[i]*Qtt(i,j))))\
)/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j)))) + 
  (2*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     ((-2 - 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        pow(rgrid[i],3)*pow(((-1 + rgrid[i])*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))/2.,2)*
        pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
        pow(1 + rgrid[i]*Qrr(i,j),2)*(Qrr(i,j) - Qtt(i,j))*
        (1 + rgrid[i]*Qtt(i,j)) - 
       ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
          pow(-1 + rgrid[i],3)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),3)*
          pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
          pow(Qr1(i,j),2)*(4 + 7*rgrid[i]*Qrr(i,j) + 
            3*pow(rgrid[i],2)*pow(Qrr(i,j),2))*
          pow(1 + rgrid[i]*Qtt(i,j),2))/4. - 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (72*pow(A1,2) + 48*
             exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
            16*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.) + 
            24*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.) + 
            72*pow(A1,2)*rgrid[i]*Q22(i,j) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             rgrid[i]*Q22(i,j) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j) - 6*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j) + 12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j) + 
            18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j) + 
            144*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             rgrid[i]*Qrr(i,j) - 
            8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Qrr(i,j) - 12*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Qrr(i,j) + 20*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
            30*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
            144*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) - 
            8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*Qrr(i,j) - 
            12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*Qrr(i,j) + 
            12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j) + 
            18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j) + 
            72*pow(A1,2)*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             pow(Qrr(i,j),2) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             pow(Qrr(i,j),2) + 
            6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2) + 
            9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2) + 
            72*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*pow(Qrr(i,j),2) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*pow(Qrr(i,j),2) + 
            2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qrr(i,j),2) + 
            3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qrr(i,j),2) + 144*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             rgrid[i]*Qtt(i,j) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Qtt(i,j) - 6*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Qtt(i,j) + 32*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
            48*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
            144*pow(A1,2)*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*Qtt(i,j) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*Qtt(i,j) + 
            24*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qtt(i,j) + 
            36*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qtt(i,j) + 
            288*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
            192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
            8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Qrr(i,j)*Qtt(i,j) - 
            12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Qrr(i,j)*Qtt(i,j) + 
            46*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*Qtt(i,j) + 
            69*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*Qtt(i,j) + 
            288*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*
             Qtt(i,j) + 192*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                 rgrid[i])/(3.*A1))*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*
             Qtt(i,j) - 8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
            12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
            30*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j)*
             Qtt(i,j) + 45*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j)*
             Qtt(i,j) + 144*pow(A1,2)*pow(rgrid[i],3)*
             pow(Qrr(i,j),2)*Qtt(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             pow(Qrr(i,j),2)*Qtt(i,j) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             pow(Qrr(i,j),2)*Qtt(i,j) + 
            18*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2)*
             Qtt(i,j) + 27*pow(A1,2)*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2)*
             Qtt(i,j) + 144*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*
             pow(Qrr(i,j),2)*Qtt(i,j) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) - 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) - 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],6)*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
               ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
             Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) + 
            10*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qrr(i,j),2)*Qtt(i,j) + 
            15*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qrr(i,j),2)*Qtt(i,j) + 
            72*pow(A1,2)*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
            12*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qtt(i,j),2) + 
            18*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qtt(i,j),2) + 
            72*pow(A1,2)*pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) + 
            8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qtt(i,j),2) + 
            12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qtt(i,j),2) + 
            144*pow(A1,2)*pow(rgrid[i],3)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
            18*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            27*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            144*pow(A1,2)*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            96*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) + 
            10*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            15*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            72*pow(A1,2)*pow(rgrid[i],4)*pow(Qrr(i,j),2)*
             pow(Qtt(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
            8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2)*
             pow(Qtt(i,j),2) + 
            12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2)*
             pow(Qtt(i,j),2) + 
            72*pow(A1,2)*pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*
             pow(Qtt(i,j),2) + 
            48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*
             pow(Qtt(i,j),2) + 
            4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
            6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
            2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
             exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
             pow(rgrid[i],4)*(1 + rgrid[i]*Q11(i,j))*
             (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
             (1 + rgrid[i]*Qtt(i,j)) + 
            rgrid[i]*Q11(i,j)*
             (pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                (72*pow(A1,2) + 
                  48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                     (3.*A1)) - 
                  4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],2)*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                     ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) \
- 6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                     ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) \
+ 2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                  3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   rgrid[i]*(((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                  rgrid[i]*(48*
                      (3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                       rgrid[i])/(3.*A1))) - 
                     2*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
                     5*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*Qtt(i,j) + 
                  2*pow(rgrid[i],2)*
                   (36*pow(A1,2) + 
                     24*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)) + 
                     (2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*
                   pow(Qtt(i,j),2)) + 
               2*(36*pow(A1,2) + 
                  24*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                     (3.*A1)) - 
                  2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   pow(rgrid[i],2)*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                     ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) \
- 3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
                     ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) \
+ 6*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                  9*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                   rgrid[i]*(((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                  rgrid[i]*(72*pow(A1,2) + 
                     48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)) - 
                     (2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
                     6*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*Qtt(i,j) + 
                  2*pow(rgrid[i],2)*
                   (6*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1))) + 
                     (2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*
                   pow(Qtt(i,j),2)) + 
               rgrid[i]*Qrr(i,j)*
                (2*(72*pow(A1,2) + 
                     48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)) - 
                     2*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
                     3*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)) + 
                  rgrid[i]*(15*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                     4*(72*pow(A1,2) + 
                        48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                        rgrid[i])/(3.*A1)) - 
                        (2 + 3*pow(A1,2))*
                        exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],2)*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.)))*Qtt(i,j) + 
                  pow(rgrid[i],2)*
                   (48*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1))) + 
                     5*(2 + 3*pow(A1,2))*
                      exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*
                   pow(Qtt(i,j),2)) + 
               rgrid[i]*Q22(i,j)*
                (pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                   (72*pow(A1,2) + 
                     48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)) - 
                     4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/
                         2.) - 
                     6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/
                         2.) - 
                     2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) - 
                     3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                     rgrid[i]*
                      (48*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                       rgrid[i])/(3.*A1))) - 
                        2*(2 + 3*pow(A1,2))*
                        exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],2)*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
                        (2 + 3*pow(A1,2))*
                        exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                        (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*Qtt(i,j) + 
                     24*(3*pow(A1,2) + 
                        2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)))*pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
                  2*(36*pow(A1,2) + 
                     24*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                        (3.*A1)) - 
                     2*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/
                         2.) - 
                     3*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      pow(rgrid[i],2)*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/
                         2.) + 
                     4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                     6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                      rgrid[i]*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                     rgrid[i]*
                      (72*pow(A1,2) + 
                        48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                        rgrid[i])/(3.*A1)) - 
                        (2 + 3*pow(A1,2))*
                        exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                        pow(rgrid[i],2)*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
                        4*(2 + 3*pow(A1,2))*
                        exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                        (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (36*pow(A1,2) + 
                        24*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                        rgrid[i])/(3.*A1)) + 
                        (2 + 3*pow(A1,2))*
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.))*
                      pow(Qtt(i,j),2)) + 
                  rgrid[i]*Qrr(i,j)*
                   (2*(72*pow(A1,2) + 
                        48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1)) - 
                        2*(2 + 3*pow(A1,2))*
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],2)*
                         (-4 - 8*rgrid[i] + 
                         3*pow(mu,2)*pow(rgrid[i],2) + 
                         ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
                        (2 + 3*pow(A1,2))*
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))/2.)) + 
                     rgrid[i]*
                      (7*(2 + 3*pow(A1,2))*
                        exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                        (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.) + 
                        4*(72*pow(A1,2) + 
                         48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
                        rgrid[i])/(3.*A1)) - 
                         (2 + 3*pow(A1,2))*
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                         pow(rgrid[i],2)*
                         (-4 - 8*rgrid[i] + 
                         3*pow(mu,2)*pow(rgrid[i],2) + 
                         ((-1 + rgrid[i])*
                        (-8 + 6*pow(mu,2)*rgrid[i]))/2.)))*Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (48*(3*pow(A1,2) + 
                         2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                         (3.*A1))) + 
                        (2 + 3*pow(A1,2))*
                         exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))/2.))*
                      pow(Qtt(i,j),2))))))/2. - 
       ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
          pow(-1 + rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          (-12 - 24*rgrid[i]*Q22(i,j) - 
            9*pow(rgrid[i],2)*pow(Q22(i,j),2) + 
            6*pow(rgrid[i],3)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) + 
            12*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2) + 
            6*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2) - 20*rgrid[i]*Qrr(i,j) - 
            46*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) - 
            16*pow(rgrid[i],3)*pow(Q22(i,j),2)*Qrr(i,j) + 
            12*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             Qrr(i,j) + 24*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2)*Qrr(i,j) + 
            12*pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2)*Qrr(i,j) + 
            pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
            10*pow(rgrid[i],3)*Q22(i,j)*pow(Qrr(i,j),2) + 
            6*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             pow(Qrr(i,j),2) + 
            12*pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2)*pow(Qrr(i,j),2) + 
            6*pow(rgrid[i],7)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2)*pow(Qrr(i,j),2) + 
            6*pow(rgrid[i],3)*pow(Qrr(i,j),3) + 
            6*pow(rgrid[i],4)*Q22(i,j)*pow(Qrr(i,j),3) + 
            4*pow(rgrid[i],5)*pow(Q22(i,j),2)*pow(Qrr(i,j),3) - 
            24*rgrid[i]*Qtt(i,j) - 
            48*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) - 
            18*pow(rgrid[i],3)*pow(Q22(i,j),2)*Qtt(i,j) + 
            12*pow(rgrid[i],4)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             Qtt(i,j) + 24*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2)*Qtt(i,j) + 
            12*pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2)*Qtt(i,j) - 
            46*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) - 
            104*pow(rgrid[i],3)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) - 
            38*pow(rgrid[i],4)*pow(Q22(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
            24*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             Qrr(i,j)*Qtt(i,j) + 
            48*pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2)*Qrr(i,j)*Qtt(i,j) + 
            24*pow(rgrid[i],7)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2)*Qrr(i,j)*Qtt(i,j) - 
            10*pow(rgrid[i],3)*pow(Qrr(i,j),2)*Qtt(i,j) - 
            44*pow(rgrid[i],4)*Q22(i,j)*pow(Qrr(i,j),2)*Qtt(i,j) - 
            12*pow(rgrid[i],5)*pow(Q22(i,j),2)*pow(Qrr(i,j),2)*
             Qtt(i,j) + 12*pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             pow(Qrr(i,j),2)*Qtt(i,j) + 
            24*pow(rgrid[i],7)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*Qtt(i,j) + 
            12*pow(rgrid[i],8)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*Qtt(i,j) + 
            6*pow(rgrid[i],4)*pow(Qrr(i,j),3)*Qtt(i,j) + 
            2*pow(rgrid[i],6)*pow(Q22(i,j),2)*pow(Qrr(i,j),3)*
             Qtt(i,j) - 9*pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
            18*pow(rgrid[i],3)*Q22(i,j)*pow(Qtt(i,j),2) - 
            6*pow(rgrid[i],4)*pow(Q22(i,j),2)*pow(Qtt(i,j),2) + 
            6*pow(rgrid[i],5)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             pow(Qtt(i,j),2) + 
            12*pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
            6*pow(rgrid[i],7)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2)*pow(Qtt(i,j),2) - 
            16*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) - 
            38*pow(rgrid[i],4)*Q22(i,j)*Qrr(i,j)*pow(Qtt(i,j),2) - 
            12*pow(rgrid[i],5)*pow(Q22(i,j),2)*Qrr(i,j)*
             pow(Qtt(i,j),2) + 
            12*pow(rgrid[i],6)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             Qrr(i,j)*pow(Qtt(i,j),2) + 
            24*pow(rgrid[i],7)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2)*Qrr(i,j)*pow(Qtt(i,j),2) + 
            12*pow(rgrid[i],8)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2)*Qrr(i,j)*pow(Qtt(i,j),2) - 
            12*pow(rgrid[i],5)*Q22(i,j)*pow(Qrr(i,j),2)*
             pow(Qtt(i,j),2) - 
            pow(rgrid[i],6)*pow(Q22(i,j),2)*pow(Qrr(i,j),2)*
             pow(Qtt(i,j),2) + 
            6*pow(rgrid[i],7)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
             pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
            12*pow(rgrid[i],8)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*Q22(i,j)*
             pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
            6*pow(rgrid[i],9)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Q22(i,j),2)*
             pow(Qr1(i,j),2)*pow(Qrr(i,j),2)*pow(Qtt(i,j),2) + 
            4*pow(rgrid[i],5)*pow(Qrr(i,j),3)*pow(Qtt(i,j),2) + 
            2*pow(rgrid[i],6)*Q22(i,j)*pow(Qrr(i,j),3)*
             pow(Qtt(i,j),2) + 
            2*pow(rgrid[i],7)*pow(Q22(i,j),2)*pow(Qrr(i,j),3)*
             pow(Qtt(i,j),2) - 
            4*pow(mu,2)*pow(h(i,j),2)*pow(rgrid[i],2)*
             pow(1 + rgrid[i]*Q11(i,j),2)*
             pow(1 + rgrid[i]*Q22(i,j),2)*
             pow(1 + rgrid[i]*Qrr(i,j),2)*pow(1 + rgrid[i]*Qtt(i,j),2) \
+ 2*rgrid[i]*Q11(i,j)*(pow(rgrid[i],3)*pow(Qrr(i,j),3)*
                (3 + pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
               pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                (-5 + 6*pow(rgrid[i],3)*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                   pow(Qr1(i,j),2) + 
                  2*rgrid[i]*
                   (-11 + 6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*Qtt(i,j) + 
                  6*pow(rgrid[i],2)*
                   (-1 + pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
               3*(-4 + 2*pow(rgrid[i],3)*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                   pow(Qr1(i,j),2) + 
                  4*rgrid[i]*
                   (-2 + pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*Qtt(i,j) + 
                  pow(rgrid[i],2)*
                   (-3 + 2*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
               rgrid[i]*Qrr(i,j)*
                (-23 + 12*pow(rgrid[i],3)*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                   pow(Qr1(i,j),2) + 
                  4*rgrid[i]*
                   (-13 + 6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*Qtt(i,j) + 
                  pow(rgrid[i],2)*
                   (-19 + 12*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
               2*rgrid[i]*Q22(i,j)*
                (-2*pow(rgrid[i],4)*pow(Qrr(i,j),3)*Qtt(i,j)*
                   (3 + rgrid[i]*Qtt(i,j)) + 
                  pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                   (-11 + 6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     2*rgrid[i]*
                      (-17 + 
                        6*pow(rgrid[i],3)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*Qtt(i,j) + 
                     6*pow(rgrid[i],2)*
                      (-2 + pow(rgrid[i],3)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
                  3*(-4 + 2*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     4*rgrid[i]*
                      (-2 + pow(rgrid[i],3)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (-3 + 2*pow(rgrid[i],3)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
                  2*rgrid[i]*Qrr(i,j)*
                   (-13 + 6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     rgrid[i]*
                      (-29 + 
                        12*pow(rgrid[i],3)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (-11 + 
                        6*pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2))) + 
               pow(rgrid[i],2)*pow(Q22(i,j),2)*
                (-(pow(rgrid[i],3)*pow(Qrr(i,j),3)*
                     (-1 + 4*rgrid[i]*Qtt(i,j) + 
                       pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
                  3*(-3 + 2*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     (-6*rgrid[i] + 
                        4*pow(rgrid[i],4)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*Qtt(i,j) + 
                     2*pow(rgrid[i],2)*
                      (-1 + pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
                  rgrid[i]*Qrr(i,j)*
                   (-19 + 12*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     4*rgrid[i]*
                      (-11 + 
                        6*pow(rgrid[i],3)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*Qtt(i,j) + 
                     3*pow(rgrid[i],2)*
                      (-5 + 4*pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
                  pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                   (-6 + 6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     12*rgrid[i]*
                      (-2 + pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (-7 + 6*pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2)))) + 
            pow(rgrid[i],2)*pow(Q11(i,j),2)*
             (2*pow(rgrid[i],3)*pow(Qrr(i,j),3)*
                (2 + rgrid[i]*Qtt(i,j) + 
                  pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
               3*(-3 + 2*pow(rgrid[i],3)*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                   pow(Qr1(i,j),2) + 
                  (-6*rgrid[i] + 
                     4*pow(rgrid[i],4)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*Qtt(i,j) + 
                  2*pow(rgrid[i],2)*
                   (-1 + pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
               2*rgrid[i]*Qrr(i,j)*
                (-8 + 6*pow(rgrid[i],3)*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                   pow(Qr1(i,j),2) + 
                  rgrid[i]*(-19 + 
                     12*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*Qtt(i,j) + 
                  6*pow(rgrid[i],2)*
                   (-1 + pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
               pow(rgrid[i],3)*pow(Qrr(i,j),2)*
                (6*pow(rgrid[i],2)*
                   (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                   pow(Qr1(i,j),2) + 
                  12*(-1 + pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*Qtt(i,j) + 
                  rgrid[i]*(-1 + 
                     6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
               2*rgrid[i]*Q22(i,j)*
                (-(pow(rgrid[i],3)*pow(Qrr(i,j),3)*
                     (-1 + 4*rgrid[i]*Qtt(i,j) + 
                       pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
                  3*(-3 + 2*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     (-6*rgrid[i] + 
                        4*pow(rgrid[i],4)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*Qtt(i,j) + 
                     2*pow(rgrid[i],2)*
                      (-1 + pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
                  rgrid[i]*Qrr(i,j)*
                   (-19 + 12*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     4*rgrid[i]*
                      (-11 + 
                        6*pow(rgrid[i],3)*
                        (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                        pow(Qr1(i,j),2))*Qtt(i,j) + 
                     3*pow(rgrid[i],2)*
                      (-5 + 4*pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
                  pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                   (-6 + 6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     12*rgrid[i]*
                      (-2 + pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (-7 + 6*pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2))) + 
               pow(rgrid[i],2)*pow(Q22(i,j),2)*
                (-2*pow(rgrid[i],3)*pow(Qrr(i,j),3)*
                   (-1 + rgrid[i]*Qtt(i,j)) + 
                  3*(-2 + 2*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     4*rgrid[i]*
                      (-1 + pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*Qtt(i,j) + 
                     pow(rgrid[i],2)*
                      (-1 + 2*pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
                  2*rgrid[i]*Qrr(i,j)*
                   (-6 + 6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     3*rgrid[i]*
                      (-5 + 4*pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*Qtt(i,j) + 
                     (-4*pow(rgrid[i],2) + 
                        6*pow(rgrid[i],5)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2)) + 
                  pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                   (-1 + 6*pow(rgrid[i],3)*
                      (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))/2.)*
                      pow(Qr1(i,j),2) + 
                     2*rgrid[i]*
                      (-7 + 6*pow(rgrid[i],3)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*Qtt(i,j) + 
                     (-2*pow(rgrid[i],2) + 
                        6*pow(rgrid[i],5)*
                         (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))/2.)*
                         pow(Qr1(i,j),2))*pow(Qtt(i,j),2))))))/4.))/
   ((2 + 3*pow(A1,2))*pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Q22(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),2))
;
elem[2]=
Q11.diff(d20,i,j) - (pow(Q11.diff(d10,i,j),2)*rgrid[i])/
   (1 + rgrid[i]*Q11(i,j)) - (2*pow(Qr1.diff(d01,i,j),2)*rgrid[i]*
     (1 + rgrid[i]*Q11(i,j)))/pow(L,2) - 
  (2*Q11.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L - 
  (4*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
     (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2))/2.))/
   (L*rgrid[i]) - (2*Q11.diff(d01,i,j)*Qr1(i,j)*
     (1 + rgrid[i]*Q11(i,j) - rgrid[i]*Qrr(i,j)))/(L + L*rgrid[i]*Q11(i,j)) \
+ (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*rgrid[i]*
     (1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))\
) + (pow(Q22.diff(d01,i,j),2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q22(i,j),2)) + 
  (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/(L*pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (2*Q11.diff(d02,i,j)*(1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (pow(Q11.diff(d01,i,j),2)*rgrid[i]*
     (3 + (-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2) + 3*rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)) + 
  (pow(Qrr.diff(d01,i,j),2)*rgrid[i])/
   (pow(L,2)*(-1 + rgrid[i])*
      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3)) + 
     pow(L,2)*(-1 + rgrid[i])*rgrid[i]*
      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
        pow(mu,2)*pow(rgrid[i],3))*Qrr(i,j)) + 
  Qrr.diff(d01,i,j)*((2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
        (1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
      (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
     (2*rgrid[i]*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
        (-((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
                3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (4 + 4*rgrid[i] + 4*pow(rgrid[i],2) - 
             pow(mu,2)*pow(rgrid[i],3))/2. + 
          (pow(-1 + rgrid[i],2)*rgrid[i]*
             pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2))/2.))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) - 
     (2*Qr1.diff(d10,i,j)*rgrid[i]*(1 + rgrid[i]*Q11(i,j)))/
      (L + L*rgrid[i]*Qrr(i,j))) + 
  (2*Qtt.diff(d01,i,j)*((-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3)) - 
       rgrid[i]*(((-1 + rgrid[i])*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))/2.))*(1 + rgrid[i]*Q11(i,j))*
     Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*rgrid[i]*(1 + rgrid[i]*Q11(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     rgrid[i]*(-1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. - rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))) + 
  (pow(Qtt.diff(d01,i,j),2)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
  a0.diff(d10,i,j)*((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        rgrid[i]*(1 + rgrid[i]*Q11(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     (2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
        pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j)) + 
       ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
          pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          (2 + 3*rgrid[i]*Q11(i,j) + 
            pow(rgrid[i],2)*pow(Q11(i,j),2))*(1 + rgrid[i]*Q22(i,j))*
          pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/2. - 
       (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (72*pow(A1,2) + 48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
               rgrid[i])/(3.*A1)) + 
          8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.) + 
          12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.) + 
          72*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*Qrr(i,j) + 
          4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
          6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
          72*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*Qtt(i,j) + 
          4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
          6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
          72*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          rgrid[i]*Q11(i,j)*(72*pow(A1,2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
             4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                 2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.) + 
             6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                 2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.) + 
             rgrid[i]*(72*pow(A1,2) + 
                48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                   (3.*A1)) + 
                (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 rgrid[i]*(((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.))*Qtt(i,j) + 
             rgrid[i]*Qrr(i,j)*
              (72*pow(A1,2) + 
                48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                   (3.*A1)) + 
                (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 rgrid[i]*(((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.) + 
                24*(3*pow(A1,2) + 
                   2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                      (3.*A1)))*rgrid[i]*Qtt(i,j)))) + 
       (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        (pow(rgrid[i],2)*pow(Q11(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
           (2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)) + 
          2*(1 + rgrid[i]*Qrr(i,j))*
           (3 + 2*rgrid[i]*Qtt(i,j) + 
             rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
          rgrid[i]*Q11(i,j)*(8 + 5*rgrid[i]*Qtt(i,j) + 
             rgrid[i]*Qrr(i,j)*(7 + 4*rgrid[i]*Qtt(i,j)) + 
             rgrid[i]*Q22(i,j)*
              (5 + 2*rgrid[i]*Qtt(i,j) + 
                rgrid[i]*Qrr(i,j)*(4 + rgrid[i]*Qtt(i,j)))))))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  Q11.diff(d10,i,j)*((2*Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Q11(i,j)) + 
     (-(pow(-1 + rgrid[i],2)*rgrid[i]*
            pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
            (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
            (1 + rgrid[i]*Qtt(i,j)))/2. + 
        (((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
                3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))/2.)*(1 + rgrid[i]*Q11(i,j))*
         (1 + rgrid[i]*Q22(i,j))*
         (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
        (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (3*Qrr(i,j) - Qtt(i,j) + 2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           Q22(i,j)*(-1 - 2*rgrid[i]*Qtt(i,j) + 
              rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
           Q11(i,j)*(1 - pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
              rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)))\
))/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))))
;
elem[3]=
Q22.diff(d20,i,j) - (pow(Q22.diff(d10,i,j),2)*rgrid[i])/
   (1 + rgrid[i]*Q22(i,j)) - (2*Q22.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L - 
  (2*Q22.diff(d01,i,j)*(2 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   (L + L*rgrid[i]*Q22(i,j)) + 
  (2*Q22.diff(d02,i,j)*(1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (2*pow(Q22.diff(d01,i,j),2)*rgrid[i]*
     (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))) + 
  (2*pow(a0.diff(d10,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*rgrid[i]*(1 + rgrid[i]*Q22(i,j)))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     rgrid[i]*(1 + rgrid[i]*Q22(i,j))*
     (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  a0.diff(d10,i,j)*((4*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        rgrid[i]*(1 + rgrid[i]*Q22(i,j)))/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
     (4*a0.diff(d01,i,j)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
     (2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
        exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],4)*
        (1 + rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Q22(i,j),2) + 
       ((2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
          pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
          (2 + 3*rgrid[i]*Q22(i,j) + 
            pow(rgrid[i],2)*pow(Q22(i,j),2))*pow(Qr1(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j)))/2. - 
       (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (72*pow(A1,2) + 48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*
               rgrid[i])/(3.*A1)) + 
          8*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.) + 
          12*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.) + 
          72*pow(A1,2)*rgrid[i]*Qrr(i,j) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*Qrr(i,j) + 
          4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
          6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
          72*pow(A1,2)*rgrid[i]*Qtt(i,j) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           rgrid[i]*Qtt(i,j) + 
          4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],2)*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
          6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
          72*pow(A1,2)*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(72*pow(A1,2) + 
             48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
             4*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                 2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.) + 
             6*pow(A1,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*rgrid[i]*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                 2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.) + 
             rgrid[i]*(72*pow(A1,2) + 
                48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                   (3.*A1)) + 
                (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 rgrid[i]*(((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.))*Qtt(i,j) + 
             rgrid[i]*Qrr(i,j)*
              (72*pow(A1,2) + 
                48*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                   (3.*A1)) + 
                (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
                 rgrid[i]*(((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.) + 
                24*(3*pow(A1,2) + 
                   2*exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                      (3.*A1)))*rgrid[i]*Qtt(i,j)))) + 
       (2 + 3*pow(A1,2))*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        (pow(rgrid[i],2)*pow(Q22(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
           (2 + rgrid[i]*Qtt(i,j)) + 
          2*(1 + rgrid[i]*Qrr(i,j))*(3 + 2*rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q22(i,j)*(8 + 5*rgrid[i]*Qtt(i,j) + 
             rgrid[i]*Qrr(i,j)*(7 + 4*rgrid[i]*Qtt(i,j))) + 
          rgrid[i]*Q11(i,j)*(pow(rgrid[i],2)*pow(Q22(i,j),2)*
              (1 + rgrid[i]*Qrr(i,j)) + 
             2*(1 + rgrid[i]*Qrr(i,j))*(2 + rgrid[i]*Qtt(i,j)) + 
             rgrid[i]*Q22(i,j)*
              (5 + 2*rgrid[i]*Qtt(i,j) + 
                rgrid[i]*Qrr(i,j)*(4 + rgrid[i]*Qtt(i,j)))))))/
   ((2 + 3*pow(A1,2))*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  Q22.diff(d10,i,j)*((2*Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Q22(i,j)) + 
     (-(pow(-1 + rgrid[i],2)*rgrid[i]*
            pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
            (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
            (1 + rgrid[i]*Qtt(i,j)))/2. + 
        (((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
                3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))/2.)*(1 + rgrid[i]*Q11(i,j))*
         (1 + rgrid[i]*Q22(i,j))*
         (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
        (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*
         (3*Qrr(i,j) - Qtt(i,j) + 2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
           Q22(i,j)*(1 + rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
           Q11(i,j)*(-1 - rgrid[i]*(2 + rgrid[i]*Q22(i,j))*Qtt(i,j) + 
              rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Q22(i,j) + rgrid[i]*Qtt(i,j)))\
))/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))))
;
elem[4]=
Qr1.diff(d20,i,j) - (2*Qr1.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  (pow(Qrr.diff(d01,i,j),2)*pow(rgrid[i],2)*Qr1(i,j))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (4*h.diff(d01,i,j)*pow(mu,2)*h(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (4*h.diff(d01,i,j)*h.diff(d10,i,j)*pow(mu,2)*rgrid[i]*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],2)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3)) - 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(rgrid[i],2)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (pow(Q22.diff(d01,i,j),2)*pow(rgrid[i],2)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Q22(i,j),2)) - 
  (Q22.diff(d01,i,j)*(Q22(i,j) + 
       2*(((-1 + rgrid[i])*rgrid[i]*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. - Qrr(i,j)))*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (2*Qr1.diff(d02,i,j)*(1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  Q11.diff(d10,i,j)*(-((Q11.diff(d01,i,j)*rgrid[i]*
          (1 + rgrid[i]*Qrr(i,j)))/
        (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3))) \
+ (Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/pow(1 + rgrid[i]*Q11(i,j),2)) + 
  Q22.diff(d10,i,j)*((Q22.diff(d01,i,j)*rgrid[i]*
        (1 + rgrid[i]*Qrr(i,j)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(1 + rgrid[i]*Q22(i,j),2)) + 
     (Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/pow(1 + rgrid[i]*Q22(i,j),2)) + 
  Q11.diff(d01,i,j)*((2*Qr1.diff(d01,i,j)*rgrid[i]*
        (1 + rgrid[i]*Qrr(i,j)))/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)) \
- ((-4 + (-1 + rgrid[i])*pow(rgrid[i],2)*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2) + 
          rgrid[i]*Q11(i,j)*(-1 + 
             (-1 + rgrid[i])*pow(rgrid[i],2)*
              (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)) - 
          2*rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
      (L*(-1 + rgrid[i])*rgrid[i]*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),3))) + 
  Qrr.diff(d10,i,j)*(-((Qr1.diff(d10,i,j)*rgrid[i])/
        (1 + rgrid[i]*Qrr(i,j))) - 
     (Qrr.diff(d01,i,j)*rgrid[i])/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j))) + 
     (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qrr(i,j)) + 
     (-(rgrid[i]*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qr1(i,j)) + 
        (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
           pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),3))/2.)/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3)) + 
        (-1 + rgrid[i])*rgrid[i]*
         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*Qrr(i,j))) + 
  (2*Qr1.diff(d01,i,j)*Qr1(i,j)*
     (-(rgrid[i]*(((-1 + rgrid[i])*
               (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))/2.)*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j)))/2. - 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (6 + 7*rgrid[i]*Qrr(i,j) + 
            2*pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
            rgrid[i]*Q11(i,j)*(4 + 3*rgrid[i]*Qrr(i,j))))/2.))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  Qrr.diff(d01,i,j)*((2*Q11.diff(d10,i,j)*rgrid[i])/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) - 
     (2*Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2)) \
- (2*Qr1.diff(d01,i,j)*rgrid[i]*
        (2 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + 2*rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j))) + 
     (Qr1.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (L + L*rgrid[i]*Qrr(i,j)) + 
     (-4 + pow(rgrid[i],3)*(((-1 + rgrid[i])*
              (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) - 
        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
           pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),4))/2. + 
        pow(rgrid[i],5)*pow(Q11(i,j),2)*pow(Qr1(i,j),2)*
         (((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
                3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))/2. - 
           (pow(-1 + rgrid[i],2)*rgrid[i]*
              pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2))/2.) - 
        5*rgrid[i]*Qrr(i,j) + (-1 + rgrid[i])*pow(rgrid[i],2)*
         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
         (1 + rgrid[i]*Qrr(i,j)) + 
        rgrid[i]*Q11(i,j)*(-2 + 
           2*pow(rgrid[i],3)*
            (((-1 + rgrid[i])*
                 (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) - 
           pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
            pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),4) - 
           3*rgrid[i]*Qrr(i,j) + 
           (-1 + rgrid[i])*pow(rgrid[i],2)*
            (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
              pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
            (1 + rgrid[i]*Qrr(i,j))))/
      (L*(-1 + rgrid[i])*rgrid[i]*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
        (1 + rgrid[i]*Qrr(i,j)))) - 
  (8*a0.diff(d01,i,j)*a0.diff(d10,i,j)*pow(mu,2)*
     exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (8*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (8*pow(a0.diff(d01,i,j),2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],2)*Qr1(i,j)*
     (1 + rgrid[i]*Qrr(i,j)))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(L + L*rgrid[i]*Qtt(i,j),2)) + 
  (2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
     (-(pow(-1 + rgrid[i],2)*rgrid[i]*
           pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
           pow(Qr1(i,j),2))/2. + 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2))\
)/2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))/2.)*
             (1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 2*Qrr(i,j) - 
            Qtt(i,j)))/2. + rgrid[i]*
        (((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
               3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))/2.)*(-Qrr(i,j) + Qtt(i,j))))/
   (pow(-1 + rgrid[i],2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(L + L*rgrid[i]*Q11(i,j))*
     pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  Qtt.diff(d10,i,j)*((Qtt.diff(d01,i,j)*rgrid[i]*
        (1 + rgrid[i]*Qrr(i,j)))/
      (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(1 + rgrid[i]*Qtt(i,j),2)) + 
     (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3)) - 
          rgrid[i]*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.))*Qr1(i,j)*
        (1 + rgrid[i]*Qrr(i,j)))/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Qtt(i,j),2))) + 
  (Qr1.diff(d10,i,j)*((-3*pow(-1 + rgrid[i],2)*rgrid[i]*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/2. + 
       (((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
               3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))/2.)*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (4 + rgrid[i]*Qrr(i,j) + 3*rgrid[i]*Qtt(i,j)) - 
       (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        (4*Qrr(i,j) + 3*rgrid[i]*pow(Qrr(i,j),2) - Qtt(i,j) + 
          2*rgrid[i]*Qrr(i,j)*Qtt(i,j) + 
          2*pow(rgrid[i],2)*pow(Qrr(i,j),2)*Qtt(i,j) + 
          Q22(i,j)*(-1 + 2*rgrid[i]*Qrr(i,j) - 2*rgrid[i]*Qtt(i,j) + 
             pow(rgrid[i],2)*pow(Qrr(i,j),2)*(2 + rgrid[i]*Qtt(i,j))) \
+ Q11(i,j)*(-1 + 2*rgrid[i]*Qrr(i,j) - 2*rgrid[i]*Qtt(i,j) + 
             pow(rgrid[i],2)*pow(Qrr(i,j),2)*(2 + rgrid[i]*Qtt(i,j)) + 
             rgrid[i]*Q22(i,j)*
              (-2 + pow(rgrid[i],2)*pow(Qrr(i,j),2) - 
                3*rgrid[i]*Qtt(i,j) - 2*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j)\
)))))/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))\
) - (Qr1(i,j)*((pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
          pow(Qr1(i,j),2)*(4 + 3*rgrid[i]*Qrr(i,j))*
          pow(1 + rgrid[i]*Qtt(i,j),2))/2. + 
       (-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))*
        (6 + 2*pow(rgrid[i],3)*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) + 
          6*rgrid[i]*Qrr(i,j) + 9*rgrid[i]*Qtt(i,j) + 
          4*pow(rgrid[i],4)*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
           Qtt(i,j) + 9*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          4*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
          2*pow(rgrid[i],5)*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
           pow(Qtt(i,j),2) + 
          4*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
          pow(rgrid[i],2)*pow(Q22(i,j),2)*
           (4 + 2*pow(rgrid[i],3)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2))\
)/2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) + 
             rgrid[i]*(5 + 4*pow(rgrid[i],3)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*
                 pow(Qr1(i,j),2))*Qtt(i,j) + 
             2*(pow(rgrid[i],2) + 
                pow(rgrid[i],5)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*
                 pow(Qr1(i,j),2))*pow(Qtt(i,j),2) + 
             rgrid[i]*Qrr(i,j)*
              (4 + 5*rgrid[i]*Qtt(i,j) + 
                2*pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
          rgrid[i]*Q22(i,j)*(9 + 
             4*pow(rgrid[i],3)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2))\
)/2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) + 
             4*rgrid[i]*(3 + 
                2*pow(rgrid[i],3)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*
                 pow(Qr1(i,j),2))*Qtt(i,j) + 
             pow(rgrid[i],2)*
              (5 + 4*pow(rgrid[i],3)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*
                 pow(Qr1(i,j),2))*pow(Qtt(i,j),2) + 
             rgrid[i]*Qrr(i,j)*
              (9 + 12*rgrid[i]*Qtt(i,j) + 
                5*pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
          pow(rgrid[i],2)*pow(Q11(i,j),2)*
           (4 + 2*pow(rgrid[i],3)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2))\
)/2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) + 
             4*rgrid[i]*Qrr(i,j) + 5*rgrid[i]*Qtt(i,j) + 
             4*pow(rgrid[i],4)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2))\
)/2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
              Qtt(i,j) + 5*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
             2*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             2*pow(rgrid[i],5)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2))\
)/2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
              pow(Qtt(i,j),2) + 
             2*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             pow(rgrid[i],2)*pow(Q22(i,j),2)*
              (2 + 2*pow(rgrid[i],3)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*
                 pow(Qr1(i,j),2) + 
                (rgrid[i] + 4*pow(rgrid[i],4)*
                    (((-1 + rgrid[i])*
                       (-4 - 8*rgrid[i] + 
                       3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                    pow(Qr1(i,j),2))*Qtt(i,j) + 
                2*pow(rgrid[i],5)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*
                 pow(Qr1(i,j),2)*pow(Qtt(i,j),2) + 
                rgrid[i]*Qrr(i,j)*(2 + rgrid[i]*Qtt(i,j))) + 
             rgrid[i]*Q22(i,j)*
              (5 + 4*pow(rgrid[i],3)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*
                 pow(Qr1(i,j),2) + 
                4*(rgrid[i] + 
                   2*pow(rgrid[i],4)*
                    (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                    pow(Qr1(i,j),2))*Qtt(i,j) + 
                (pow(rgrid[i],2) + 
                   4*pow(rgrid[i],5)*
                    (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                    pow(Qr1(i,j),2))*pow(Qtt(i,j),2) + 
                rgrid[i]*Qrr(i,j)*
                 (5 + 4*rgrid[i]*Qtt(i,j) + 
                   pow(rgrid[i],2)*pow(Qtt(i,j),2)))) + 
          rgrid[i]*Q11(i,j)*(9 + 
             4*pow(rgrid[i],3)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                 2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) + 
             9*rgrid[i]*Qrr(i,j) + 12*rgrid[i]*Qtt(i,j) + 
             8*pow(rgrid[i],4)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                 2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
              Qtt(i,j) + 12*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
             5*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
             4*pow(rgrid[i],5)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                 2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2)*
              pow(Qtt(i,j),2) + 
             5*pow(rgrid[i],3)*Qrr(i,j)*pow(Qtt(i,j),2) + 
             4*rgrid[i]*Q22(i,j)*
              (3 + 2*pow(rgrid[i],3)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*
                 pow(Qr1(i,j),2) + 
                rgrid[i]*(3 + 
                   4*pow(rgrid[i],3)*
                    (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                    pow(Qr1(i,j),2))*Qtt(i,j) + 
                (pow(rgrid[i],2) + 
                   2*pow(rgrid[i],5)*
                    (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                    pow(Qr1(i,j),2))*pow(Qtt(i,j),2) + 
                rgrid[i]*Qrr(i,j)*
                 (3 + 3*rgrid[i]*Qtt(i,j) + 
                   pow(rgrid[i],2)*pow(Qtt(i,j),2))) + 
             pow(rgrid[i],2)*pow(Q22(i,j),2)*
              (5 + 4*pow(rgrid[i],3)*
                 (((-1 + rgrid[i])*
                      (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qr1(i,j),2) \
+ 4*(rgrid[i] + 2*pow(rgrid[i],4)*
                    (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                    pow(Qr1(i,j),2))*Qtt(i,j) + 
                (pow(rgrid[i],2) + 
                   4*pow(rgrid[i],5)*
                    (((-1 + rgrid[i])*
                        (-4 - 8*rgrid[i] + 
                        3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))/2.)*
                    pow(Qr1(i,j),2))*pow(Qtt(i,j),2) + 
                rgrid[i]*Qrr(i,j)*
                 (5 + 4*rgrid[i]*Qtt(i,j) + 
                   pow(rgrid[i],2)*pow(Qtt(i,j),2))))) - 
       pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (2*(-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
             ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.) + 
          3*rgrid[i]*(-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
             ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*Qrr(i,j) - 
          4*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j) + 
          pow(rgrid[i],2)*(-4 - 8*rgrid[i] + 
             3*pow(mu,2)*pow(rgrid[i],2) + 
             ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
           pow(Qrr(i,j),2) - 
          3*rgrid[i]*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2) + 
          3*rgrid[i]*(-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
             ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*Qtt(i,j) - 
          2*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qtt(i,j) + 
          4*pow(rgrid[i],2)*(-4 - 8*rgrid[i] + 
             3*pow(mu,2)*pow(rgrid[i],2) + 
             ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*Qrr(i,j)*
           Qtt(i,j) - 12*rgrid[i]*
           (((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*Qtt(i,j) + 
          pow(rgrid[i],3)*(-4 - 8*rgrid[i] + 
             3*pow(mu,2)*pow(rgrid[i],2) + 
             ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
           pow(Qrr(i,j),2)*Qtt(i,j) - 
          8*pow(rgrid[i],2)*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2)*
           Qtt(i,j) + pow(rgrid[i],2)*
           (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2) + 
             ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
           pow(Qtt(i,j),2) - 
          rgrid[i]*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qtt(i,j),2) + 
          pow(rgrid[i],3)*(-4 - 8*rgrid[i] + 
             3*pow(mu,2)*pow(rgrid[i],2) + 
             ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*Qrr(i,j)*
           pow(Qtt(i,j),2) - 
          6*pow(rgrid[i],2)*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*Qrr(i,j)*
           pow(Qtt(i,j),2) - 
          4*pow(rgrid[i],3)*(((-1 + rgrid[i])*
                (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))/2.)*pow(Qrr(i,j),2)*
           pow(Qtt(i,j),2) + 
          Q22(i,j)*(rgrid[i]*(-4 - 8*rgrid[i] + 
                3*pow(mu,2)*pow(rgrid[i],2) + 
                ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
              (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))*
              (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) - 
             (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/
                 2. + (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.)*
              (-2 - 2*rgrid[i]*Qtt(i,j) - 
                pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                2*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j)*
                 (2 + rgrid[i]*Qtt(i,j)) + 
                pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                 (1 + 4*rgrid[i]*Qtt(i,j) + 
                   2*pow(rgrid[i],2)*pow(Qtt(i,j),2)))) + 
          Q11(i,j)*(rgrid[i]*(-4 - 8*rgrid[i] + 
                3*pow(mu,2)*pow(rgrid[i],2) + 
                ((-1 + rgrid[i])*(-8 + 6*pow(mu,2)*rgrid[i]))/2.)*
              (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
              (1 + rgrid[i]*Qtt(i,j))*
              (2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)) + 
             (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. \
+ (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))/
                 2.)*(2 + 2*rgrid[i]*Qtt(i,j) + 
                pow(rgrid[i],2)*pow(Qtt(i,j),2) - 
                2*pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j)*
                 (2 + rgrid[i]*Qtt(i,j)) - 
                pow(rgrid[i],2)*pow(Qrr(i,j),2)*
                 (1 + 4*rgrid[i]*Qtt(i,j) + 
                   2*pow(rgrid[i],2)*pow(Qtt(i,j),2)) + 
                rgrid[i]*Q22(i,j)*
                 (4 + pow(rgrid[i],2)*pow(Qrr(i,j),2) + 
                   6*rgrid[i]*Qtt(i,j) + 
                   3*pow(rgrid[i],2)*pow(Qtt(i,j),2) + 
                   2*rgrid[i]*Qrr(i,j)*
                    (2 + 2*rgrid[i]*Qtt(i,j) + 
                      pow(rgrid[i],2)*pow(Qtt(i,j),2))))))))/
   ((-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
     (1 + rgrid[i]*Qrr(i,j))*pow(1 + rgrid[i]*Qtt(i,j),2))
;
elem[5]=
a0.diff(d20,i,j) - (a0.diff(d01,i,j)*Qr1.diff(d10,i,j)*rgrid[i])/L + 
  (h.diff(d10,i,j)*B1*mu*a0(i,j)*rgrid[i])/(-1 + rgrid[i]) - 
  (2*a0.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  (h.diff(d01,i,j)*B1*mu*a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
   (L - L*rgrid[i]) + Qr1.diff(d01,i,j)*
   (-((a0.diff(d10,i,j)*rgrid[i])/L) + 
     (a0(i,j)*rgrid[i])/(L - L*rgrid[i]) + 
     (2*a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2)) + 
  Q11.diff(d10,i,j)*((a0(i,j)*rgrid[i])/
      (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))) + 
     (a0.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q11(i,j)) - 
     (a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j)))) + 
  Q22.diff(d10,i,j)*((a0(i,j)*rgrid[i])/
      (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q22(i,j))) + 
     (a0.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q22(i,j)) - 
     (a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q22(i,j)))) + 
  (2*a0.diff(d02,i,j)*(1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  Q11.diff(d01,i,j)*(-(a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))) - 
     (a0.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) + 
     (a0.diff(d01,i,j)*rgrid[i]*
        (-1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. - rgrid[i]*Qrr(i,j)))/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2))) \
+ Q22.diff(d01,i,j)*(-(a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
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
        (1 + rgrid[i]*Q22(i,j)))) + 
  Qrr.diff(d10,i,j)*(-(a0(i,j)*rgrid[i])/
      (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j))) - 
     (a0.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qrr(i,j)) + 
     (a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j))) + 
  Qrr.diff(d01,i,j)*((a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
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
      (2*L + 2*L*rgrid[i]*Qrr(i,j))) + 
  Qtt.diff(d10,i,j)*(-(a0(i,j)*rgrid[i])/
      (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Qtt(i,j))) - 
     (a0.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qtt(i,j)) + 
     (a0.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qtt(i,j))) + 
  Qtt.diff(d01,i,j)*((a0(i,j)*pow(rgrid[i],2)*Qr1(i,j))/
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
      (2*L + 2*L*rgrid[i]*Qtt(i,j))) + 
  (a0(i,j)*(Q22(i,j) - Qrr(i,j) - Qtt(i,j) - 
       2*rgrid[i]*Qrr(i,j)*Qtt(i,j) - 
       pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j)*Qtt(i,j) + 
       2*B1*mu*h(i,j)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
       Q11(i,j)*(1 - pow(rgrid[i],2)*Qrr(i,j)*Qtt(i,j) + 
          rgrid[i]*Q22(i,j)*(2 + rgrid[i]*Qrr(i,j) + rgrid[i]*Qtt(i,j)))))/
   (2.*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  a0.diff(d01,i,j)*(-((h.diff(d10,i,j)*B1*mu*pow(rgrid[i],2)*
          Qr1(i,j))/L) + (2*h.diff(d01,i,j)*B1*mu*rgrid[i]*
        (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
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
          2*B1*mu*h(i,j)*(-1 + rgrid[i])*rgrid[i]*(1 + rgrid[i]*Q11(i,j))*
           (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
           (1 + rgrid[i]*Qtt(i,j)) + 
          rgrid[i]*Q11(i,j)*(3 - 7*rgrid[i] + 2*rgrid[i]*Qtt(i,j) - 
             6*pow(rgrid[i],2)*Qtt(i,j) + 
             rgrid[i]*Qrr(i,j)*
              (2 - 6*rgrid[i] + (rgrid[i] - 5*pow(rgrid[i],2))*Qtt(i,j)) \
+ rgrid[i]*Q22(i,j)*(4 - 8*rgrid[i] + (3 - 7*rgrid[i])*rgrid[i]*Qtt(i,j) + 
                rgrid[i]*Qrr(i,j)*
                 (3 - 7*rgrid[i] + 2*(1 - 3*rgrid[i])*rgrid[i]*Qtt(i,j)))))\
)/(2.*L*(-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))) + 
  a0.diff(d10,i,j)*(h.diff(d10,i,j)*B1*mu*rgrid[i] - 
     (h.diff(d01,i,j)*B1*mu*pow(rgrid[i],2)*Qr1(i,j))/L + 
     (4 - Q22(i,j) + 5*rgrid[i]*Q22(i,j) + Qrr(i,j) + 3*rgrid[i]*Qrr(i,j) + 
        4*pow(rgrid[i],2)*Q22(i,j)*Qrr(i,j) + Qtt(i,j) + 
        3*rgrid[i]*Qtt(i,j) + 4*pow(rgrid[i],2)*Q22(i,j)*Qtt(i,j) + 
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
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))
;
elem[6]=
h.diff(d20,i,j) - (h.diff(d01,i,j)*Qr1.diff(d10,i,j)*rgrid[i])/L - 
  (2*h.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  Qr1.diff(d01,i,j)*(-(h(i,j)/L) - (h.diff(d10,i,j)*rgrid[i])/L + 
     (2*h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/pow(L,2)) + 
  Q11.diff(d10,i,j)*(h(i,j)/(2 + 2*rgrid[i]*Q11(i,j)) + 
     (h.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q11(i,j)) - 
     (h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j)))) + 
  Q22.diff(d10,i,j)*(h(i,j)/(2 + 2*rgrid[i]*Q22(i,j)) + 
     (h.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Q22(i,j)) - 
     (h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q22(i,j)))) + 
  (2*h.diff(d02,i,j)*(1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  Q11.diff(d01,i,j)*(-(h(i,j)*rgrid[i]*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) - 
     (h.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q11(i,j))) + 
     (h.diff(d01,i,j)*rgrid[i]*
        (-1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. - rgrid[i]*Qrr(i,j)))/
      ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*pow(L + L*rgrid[i]*Q11(i,j),2))) \
+ Q22.diff(d01,i,j)*(-(h(i,j)*rgrid[i]*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Q22(i,j))) - 
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
        (1 + rgrid[i]*Q22(i,j)))) + 
  Qrr.diff(d10,i,j)*(-(h(i,j)/(2 + 2*rgrid[i]*Qrr(i,j))) - 
     (h.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qrr(i,j)) + 
     (h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j))) + 
  Qrr.diff(d01,i,j)*((h.diff(d01,i,j)*rgrid[i]*
        (1 - ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j))) + 
     (h(i,j)*rgrid[i]*Qr1(i,j))/(2*L + 2*L*rgrid[i]*Qrr(i,j)) + 
     (h.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2*L + 2*L*rgrid[i]*Qrr(i,j))) + 
  (pow(a0.diff(d10,i,j),2)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*rgrid[i])/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (2*a0.diff(d01,i,j)*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*Qr1(i,j))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (2*pow(a0.diff(d01,i,j),2)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*rgrid[i]*
     (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  a0.diff(d10,i,j)*((2*B1*mu*a0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
        rgrid[i])/
      ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
     (2*a0.diff(d01,i,j)*B1*mu*exp(B1*mu*h(i,j)*rgrid[i])*
        (-1 + rgrid[i])*pow(rgrid[i],2)*Qr1(i,j))/
      (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
  Qtt.diff(d10,i,j)*(h(i,j)/(2 + 2*rgrid[i]*Qtt(i,j)) + 
     (h.diff(d10,i,j)*rgrid[i])/(2 + 2*rgrid[i]*Qtt(i,j)) - 
     (h.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Qtt(i,j)))) + 
  Qtt.diff(d01,i,j)*((h.diff(d01,i,j)*rgrid[i]*
        (1 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
             pow(Qr1(i,j),2))/2. + rgrid[i]*Qrr(i,j)))/
      (pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qtt(i,j))) - 
     (h(i,j)*rgrid[i]*Qr1(i,j))/(2.*(L + L*rgrid[i]*Qtt(i,j))) - 
     (h.diff(d10,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
      (2.*(L + L*rgrid[i]*Qtt(i,j)))) + 
  (h.diff(d10,i,j)*(2*(((-1 + rgrid[i])*
             (-4 - 8*rgrid[i] + 3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))/2.)*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j)) + 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (-Qrr(i,j) + Qtt(i,j) + 
            Q22(i,j)*(1 + rgrid[i]*(2 + rgrid[i]*Qrr(i,j))*Qtt(i,j)) + 
            Q11(i,j)*(1 + rgrid[i]*(2 + rgrid[i]*Qrr(i,j))*Qtt(i,j) + 
               rgrid[i]*Q22(i,j)*
                (2 + 3*rgrid[i]*Qtt(i,j) + 
                  rgrid[i]*Qrr(i,j)*(1 + 2*rgrid[i]*Qtt(i,j))))))/2.))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))\
) + ((2*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*
        ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
           exp(((2 + 3*A1*B1)*mu*h(i,j)*rgrid[i])/(3.*A1))*
           pow(rgrid[i],4) + 
          2*(12*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                  (3.*A1))) + 
             (2 + 3*pow(A1,2))*mu*
              exp((2*mu*h(i,j)*rgrid[i])/(3.*A1))*h(i,j)*
              pow(rgrid[i],2)*
              (((-1 + rgrid[i])*
                   (-4 - 8*rgrid[i] + 
                     3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))/2.) + 
             12*A1*(-1 + exp(((2 + 3*pow(A1,2))*mu*h(i,j)*rgrid[i])/
                  (3.*A1)))*rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))/
      ((2 + 3*pow(A1,2))*mu*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))) - 
     (h(i,j)*rgrid[i]*(4 + 5*rgrid[i]*Qrr(i,j) + 3*rgrid[i]*Qtt(i,j) + 
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
   (2.*pow(rgrid[i],3)*(1 + rgrid[i]*Qtt(i,j))) - 
  (h.diff(d01,i,j)*Qr1(i,j)*(2*rgrid[i]*
        (((-1 + rgrid[i])*(-4 - 8*rgrid[i] + 
               3*pow(mu,2)*pow(rgrid[i],2)))/2. + 
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))/2.)*(1 + rgrid[i]*Q11(i,j))*
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
}
