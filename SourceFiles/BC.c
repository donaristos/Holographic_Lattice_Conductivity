
elem[0]=
-Qrr(i,j) + Qtt(i,j)
;
elem[1]=
(-8*pow(Qtt.diff(d01,i,j),2) + 8*Qtt.diff(d02,i,j) - 
    12*pow(Qtt.diff(d01,i,j),2)*pow(A1,2) + 
    12*Qtt.diff(d02,i,j)*pow(A1,2) + 96*pow(L,2) - 
    24*Qrr.diff(d10,i,j)*pow(L,2) - 24*Qtt.diff(d10,i,j)*pow(L,2) + 
    144*pow(A1,2)*pow(L,2) - 
    36*Qrr.diff(d10,i,j)*pow(A1,2)*pow(L,2) - 
    36*Qtt.diff(d10,i,j)*pow(A1,2)*pow(L,2) + 
    8*pow(L,2)*pow(mu,2) + 
    2*Qrr.diff(d10,i,j)*pow(L,2)*pow(mu,2) + 
    2*Qtt.diff(d10,i,j)*pow(L,2)*pow(mu,2) + 
    12*pow(A1,2)*pow(L,2)*pow(mu,2) + 
    3*Qrr.diff(d10,i,j)*pow(A1,2)*pow(L,2)*pow(mu,2) + 
    3*Qtt.diff(d10,i,j)*pow(A1,2)*pow(L,2)*pow(mu,2) - 
    144*pow(A1,2)*pow(L,2)*exp((-2*mu*h(i,j))/(3.*A1)) - 
    96*pow(L,2)*exp(A1*mu*h(i,j)) - 
    8*pow(Qtt.diff(d01,i,j),2)*Q22(i,j) + 
    8*Qtt.diff(d02,i,j)*Q22(i,j) - 
    12*pow(Qtt.diff(d01,i,j),2)*pow(A1,2)*Q22(i,j) + 
    12*Qtt.diff(d02,i,j)*pow(A1,2)*Q22(i,j) + 48*pow(L,2)*Q22(i,j) - 
    24*Qrr.diff(d10,i,j)*pow(L,2)*Q22(i,j) - 
    24*Qtt.diff(d10,i,j)*pow(L,2)*Q22(i,j) + 
    72*pow(A1,2)*pow(L,2)*Q22(i,j) - 
    36*Qrr.diff(d10,i,j)*pow(A1,2)*pow(L,2)*Q22(i,j) - 
    36*Qtt.diff(d10,i,j)*pow(A1,2)*pow(L,2)*Q22(i,j) + 
    12*pow(L,2)*pow(mu,2)*Q22(i,j) + 
    2*Qrr.diff(d10,i,j)*pow(L,2)*pow(mu,2)*Q22(i,j) + 
    2*Qtt.diff(d10,i,j)*pow(L,2)*pow(mu,2)*Q22(i,j) + 
    18*pow(A1,2)*pow(L,2)*pow(mu,2)*Q22(i,j) + 
    3*Qrr.diff(d10,i,j)*pow(A1,2)*pow(L,2)*pow(mu,2)*Q22(i,j) + 
    3*Qtt.diff(d10,i,j)*pow(A1,2)*pow(L,2)*pow(mu,2)*Q22(i,j) - 
    144*pow(A1,2)*pow(L,2)*exp((-2*mu*h(i,j))/(3.*A1))*Q22(i,j) - 
    96*pow(L,2)*exp(A1*mu*h(i,j))*Q22(i,j) - 
    4*(2 + 3*pow(A1,2))*pow(L,2)*pow(mu,2)*pow(a0(i,j),2)*
     exp(B1*mu*h(i,j))*(1 + Q11(i,j))*(1 + Q22(i,j)) + 
    8*Qtt.diff(d02,i,j)*Qtt(i,j) + 
    12*Qtt.diff(d02,i,j)*pow(A1,2)*Qtt(i,j) + 144*pow(L,2)*Qtt(i,j) + 
    216*pow(A1,2)*pow(L,2)*Qtt(i,j) + 
    4*pow(L,2)*pow(mu,2)*Qtt(i,j) + 
    6*pow(A1,2)*pow(L,2)*pow(mu,2)*Qtt(i,j) - 
    288*pow(A1,2)*pow(L,2)*exp((-2*mu*h(i,j))/(3.*A1))*Qtt(i,j) - 
    192*pow(L,2)*exp(A1*mu*h(i,j))*Qtt(i,j) + 
    8*Qtt.diff(d02,i,j)*Q22(i,j)*Qtt(i,j) + 
    12*Qtt.diff(d02,i,j)*pow(A1,2)*Q22(i,j)*Qtt(i,j) + 
    48*pow(L,2)*Q22(i,j)*Qtt(i,j) + 
    72*pow(A1,2)*pow(L,2)*Q22(i,j)*Qtt(i,j) + 
    12*pow(L,2)*pow(mu,2)*Q22(i,j)*Qtt(i,j) + 
    18*pow(A1,2)*pow(L,2)*pow(mu,2)*Q22(i,j)*Qtt(i,j) - 
    288*pow(A1,2)*pow(L,2)*exp((-2*mu*h(i,j))/(3.*A1))*Q22(i,j)*
     Qtt(i,j) - 192*pow(L,2)*exp(A1*mu*h(i,j))*Q22(i,j)*Qtt(i,j) + 
    96*pow(L,2)*pow(Qtt(i,j),2) + 
    144*pow(A1,2)*pow(L,2)*pow(Qtt(i,j),2) - 
    8*pow(L,2)*pow(mu,2)*pow(Qtt(i,j),2) - 
    12*pow(A1,2)*pow(L,2)*pow(mu,2)*pow(Qtt(i,j),2) - 
    144*pow(A1,2)*pow(L,2)*exp((-2*mu*h(i,j))/(3.*A1))*
     pow(Qtt(i,j),2) - 96*pow(L,2)*exp(A1*mu*h(i,j))*
     pow(Qtt(i,j),2) + 48*pow(L,2)*Q22(i,j)*pow(Qtt(i,j),2) + 
    72*pow(A1,2)*pow(L,2)*Q22(i,j)*pow(Qtt(i,j),2) - 
    4*pow(L,2)*pow(mu,2)*Q22(i,j)*pow(Qtt(i,j),2) - 
    6*pow(A1,2)*pow(L,2)*pow(mu,2)*Q22(i,j)*pow(Qtt(i,j),2) - 
    144*pow(A1,2)*pow(L,2)*exp((-2*mu*h(i,j))/(3.*A1))*Q22(i,j)*
     pow(Qtt(i,j),2) - 96*pow(L,2)*exp(A1*mu*h(i,j))*Q22(i,j)*
     pow(Qtt(i,j),2) + pow(L,2)*Q11(i,j)*
     (48 - 24*Qrr.diff(d10,i,j) - 24*Qtt.diff(d10,i,j) + 
       72*pow(A1,2) - 36*Qrr.diff(d10,i,j)*pow(A1,2) - 
       36*Qtt.diff(d10,i,j)*pow(A1,2) + 12*pow(mu,2) + 
       2*Qrr.diff(d10,i,j)*pow(mu,2) + 
       2*Qtt.diff(d10,i,j)*pow(mu,2) + 18*pow(A1,2)*pow(mu,2) + 
       3*Qrr.diff(d10,i,j)*pow(A1,2)*pow(mu,2) + 
       3*Qtt.diff(d10,i,j)*pow(A1,2)*pow(mu,2) - 
       144*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1)) - 96*exp(A1*mu*h(i,j)) + 
       48*Qtt(i,j) + 72*pow(A1,2)*Qtt(i,j) + 12*pow(mu,2)*Qtt(i,j) + 
       18*pow(A1,2)*pow(mu,2)*Qtt(i,j) - 
       288*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1))*Qtt(i,j) - 
       192*exp(A1*mu*h(i,j))*Qtt(i,j) + 48*pow(Qtt(i,j),2) + 
       72*pow(A1,2)*pow(Qtt(i,j),2) - 4*pow(mu,2)*pow(Qtt(i,j),2) - 
       6*pow(A1,2)*pow(mu,2)*pow(Qtt(i,j),2) - 
       144*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1))*pow(Qtt(i,j),2) - 
       96*exp(A1*mu*h(i,j))*pow(Qtt(i,j),2) + 
       Q22(i,j)*(-24*Qtt.diff(d10,i,j) - 
          36*Qtt.diff(d10,i,j)*pow(A1,2) + 16*pow(mu,2) + 
          2*Qtt.diff(d10,i,j)*pow(mu,2) + 24*pow(A1,2)*pow(mu,2) + 
          3*Qtt.diff(d10,i,j)*pow(A1,2)*pow(mu,2) + 
          Qrr.diff(d10,i,j)*(2 + 3*pow(A1,2))*(-12 + pow(mu,2)) - 
          144*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1)) - 
          96*exp(A1*mu*h(i,j)) + 
          (pow(A1,2)*(-72 + 30*pow(mu,2) - 
                288*exp((-2*mu*h(i,j))/(3.*A1))) - 
             4*(12 - 5*pow(mu,2) + 48*exp(A1*mu*h(i,j))))*Qtt(i,j) + 
          (-144*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1)) - 
             96*exp(A1*mu*h(i,j)))*pow(Qtt(i,j),2))))/
  (2.*(2 + 3*pow(A1,2))*pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*
    (1 + Q22(i,j)))
;
elem[2]=
Q11.diff(d10,i,j) + (2*pow(Qtt.diff(d01,i,j),2))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Qtt(i,j))) - 
  (2*Qtt.diff(d01,i,j)*(1 + Q11(i,j))*Qr1(i,j))/(L*(1 + Qtt(i,j))) + 
  (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))) - 
  (3*pow(Q11.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
  (2*Q11.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (pow(Q22.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
   ((-12 + pow(mu,2))*pow(L + L*Q22(i,j),2)) + 
  (exp((-2*mu*h(i,j))/(3.*A1))*
     (2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
        exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j))*(1 + Q11(i,j)) - 
       (1 + Qtt(i,j))*(72*pow(A1,2) - 48*exp((2*mu*h(i,j))/(3.*A1)) - 
          72*pow(A1,2)*exp((2*mu*h(i,j))/(3.*A1)) + 
          4*pow(mu,2)*exp((2*mu*h(i,j))/(3.*A1)) + 
          6*pow(A1,2)*pow(mu,2)*exp((2*mu*h(i,j))/(3.*A1)) + 
          48*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)) + 
          72*pow(A1,2)*Qtt(i,j) + 
          48*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j))*Qtt(i,j) + 
          Q11(i,j)*(3*pow(A1,2)*
              (24 + (-12 + pow(mu,2))*exp((2*mu*h(i,j))/(3.*A1))) + 
             2*exp((2*mu*h(i,j))/(3.*A1))*
              (-12 + pow(mu,2) + 24*exp(A1*mu*h(i,j))) + 
             24*(3*pow(A1,2) + 2*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*
              Qtt(i,j)))))/
   ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[3]=
Q22.diff(d10,i,j) + (2*Q22.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (2*pow(Q22.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))) + 
  (exp((-2*mu*h(i,j))/(3.*A1))*
     (2*(2 + 3*pow(A1,2))*pow(mu,2)*pow(a0(i,j),2)*
        exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j))*(1 + Q22(i,j)) - 
       (1 + Qtt(i,j))*(72*pow(A1,2) - 48*exp((2*mu*h(i,j))/(3.*A1)) - 
          72*pow(A1,2)*exp((2*mu*h(i,j))/(3.*A1)) + 
          4*pow(mu,2)*exp((2*mu*h(i,j))/(3.*A1)) + 
          6*pow(A1,2)*pow(mu,2)*exp((2*mu*h(i,j))/(3.*A1)) + 
          48*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)) + 
          72*pow(A1,2)*Qtt(i,j) + 
          48*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j))*Qtt(i,j) + 
          Q22(i,j)*(3*pow(A1,2)*
              (24 + (-12 + pow(mu,2))*exp((2*mu*h(i,j))/(3.*A1))) + 
             2*exp((2*mu*h(i,j))/(3.*A1))*
              (-12 + pow(mu,2) + 24*exp(A1*mu*h(i,j))) + 
             24*(3*pow(A1,2) + 2*exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*
              Qtt(i,j)))))/
   ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
elem[4]=
(4*Qr1.diff(d10,i,j) - (16*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*
       exp(B1*mu*h(i,j)))/(L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))) - 
    (2*Qr1.diff(d01,i,j)*Qr1(i,j))/L - 
    (Qtt.diff(d10,i,j)*((-2*Qtt.diff(d01,i,j))/
          (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + Qr1(i,j)))/
     (1 + Qtt(i,j)) - (Qrr.diff(d10,i,j)*
       ((2*Qtt.diff(d01,i,j))/(L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
         Qr1(i,j)))/(1 + Qtt(i,j)) - 
    (2*Q11.diff(d01,i,j)*Q11.diff(d10,i,j)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) + 
    (4*Qr1.diff(d02,i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
    (8*h.diff(d01,i,j)*h.diff(d10,i,j)*pow(mu,2)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
    (8*h.diff(d01,i,j)*pow(mu,2)*h(i,j)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
    (2*Q22.diff(d01,i,j)*Q22.diff(d10,i,j)*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
    (2*pow(Q11.diff(d01,i,j),2)*Qr1(i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) - 
    (8*pow(h.diff(d01,i,j),2)*pow(mu,2)*Qr1(i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
    (2*pow(Q22.diff(d01,i,j),2)*Qr1(i,j)*(1 + Qtt(i,j)))/
     ((-12 + pow(mu,2))*(1 + Q11(i,j))*pow(L + L*Q22(i,j),2)) - 
    (2*Q22.diff(d01,i,j)*(Q22(i,j) - 2*Qtt(i,j))*(1 + Qtt(i,j)))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
    (Qtt.diff(d01,i,j)*((-12 + pow(mu,2))*(1 + Q11(i,j))*
          pow(Qr1(i,j),2) + 
         2*(-Qrr.diff(d10,i,j) + Qtt.diff(d10,i,j) + Qtt(i,j))))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) + 
    (Qtt.diff(d01,i,j)*(-8*Qr1.diff(d01,i,j) - 8*L - 
         4*Q11.diff(d01,i,j)*Qr1(i,j) - 12*L*pow(Qr1(i,j),2) + 
         L*pow(mu,2)*pow(Qr1(i,j),2) + 
         L*(-12 + pow(mu,2))*pow(Q11(i,j),2)*pow(Qr1(i,j),2) - 
         8*Qr1.diff(d01,i,j)*Qtt(i,j) - 10*L*Qtt(i,j) - 
         4*Q11.diff(d01,i,j)*Qr1(i,j)*Qtt(i,j) + 
         4*Q11.diff(d10,i,j)*L*(1 + Qtt(i,j)) + 
         2*Q11(i,j)*(-4*Qr1.diff(d01,i,j) - 2*L + 
            L*(-12 + pow(mu,2))*pow(Qr1(i,j),2) - 
            4*Qr1.diff(d01,i,j)*Qtt(i,j) - 3*L*Qtt(i,j))))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Qtt(i,j))) \
+ (2*Q11.diff(d01,i,j)*(1 + Qtt(i,j))*
       ((2*Qr1.diff(d01,i,j) + L)*Q11(i,j) + 
         2*(Qr1.diff(d01,i,j) + 2*L + L*Qtt(i,j))))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)) + 
    (2*Qr1(i,j)*(-24 + 6*pow(mu,2) + 12*Qtt(i,j) + 
         3*pow(mu,2)*Qtt(i,j) + 24*pow(Qtt(i,j),2) - 
         2*pow(mu,2)*pow(Qtt(i,j),2) + 
         Q22(i,j)*(-36 + 7*pow(mu,2) + (-12 + 5*pow(mu,2))*Qtt(i,j) - 
            (-12 + pow(mu,2))*pow(Qtt(i,j),2)) + 
         Q11(i,j)*(-36 + 7*pow(mu,2) + (-12 + 5*pow(mu,2))*Qtt(i,j) - 
            (-12 + pow(mu,2))*pow(Qtt(i,j),2) + 
            Q22(i,j)*(8*(-6 + pow(mu,2)) + 
               (-36 + 7*pow(mu,2))*Qtt(i,j)))))/
     ((-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j))))/2.
;
elem[5]=
(4*a0.diff(d10,i,j) - (2*Qr1.diff(d01,i,j)*a0(i,j))/L + 
    2*h.diff(d10,i,j)*B1*mu*a0(i,j) + 
    (Q11.diff(d10,i,j)*a0(i,j))/(1 + Q11(i,j)) + 
    (Q22.diff(d10,i,j)*a0(i,j))/(1 + Q22(i,j)) - 
    (2*h.diff(d01,i,j)*B1*mu*a0(i,j)*Qr1(i,j))/L - 
    (Qrr.diff(d10,i,j)*a0(i,j))/(1 + Qtt(i,j)) - 
    (Qtt.diff(d10,i,j)*a0(i,j))/(1 + Qtt(i,j)) + 
    (4*a0.diff(d02,i,j)*(1 + Qtt(i,j)))/
     (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
    (Qtt.diff(d01,i,j)*((-2*a0.diff(d01,i,j))/
          ((-12 + pow(mu,2))*(1 + Q11(i,j))) + 
         (L*a0(i,j)*Qr1(i,j))/(1 + Qtt(i,j))))/pow(L,2) + 
    (Qtt.diff(d01,i,j)*((2*a0.diff(d01,i,j))/
          ((-12 + pow(mu,2))*(1 + Q11(i,j))) + 
         (L*a0(i,j)*Qr1(i,j))/(1 + Qtt(i,j))))/pow(L,2) - 
    (Q11.diff(d01,i,j)*(L*(-12 + pow(mu,2))*a0(i,j)*(1 + Q11(i,j))*
          Qr1(i,j) + 2*a0.diff(d01,i,j)*(1 + Qtt(i,j))))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
    (Q22.diff(d01,i,j)*(-(L*a0(i,j)*Qr1(i,j)) + 
         (2*a0.diff(d01,i,j)*(1 + Qtt(i,j)))/
          ((-12 + pow(mu,2))*(1 + Q11(i,j)))))/
     (pow(L,2)*(1 + Q22(i,j))) + 
    (4*a0.diff(d01,i,j)*(-(L*Qr1(i,j)) + 
         (h.diff(d01,i,j)*B1*mu*(1 + Qtt(i,j)))/
          ((-12 + pow(mu,2))*(1 + Q11(i,j)))))/pow(L,2) + 
    (a0(i,j)*(Q22(i,j) + Q11(i,j)*(1 + 2*Q22(i,j) - Qtt(i,j)) - 
         2*Qtt(i,j) - Q22(i,j)*Qtt(i,j) + 
         2*B1*mu*h(i,j)*(1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j))))/
     ((1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j))))/2.
;
elem[6]=
h.diff(d10,i,j) + (2*h.diff(d01,i,j)*Qtt.diff(d01,i,j))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (h.diff(d01,i,j)*Qr1(i,j))/L - 
  (h.diff(d01,i,j)*Q11.diff(d01,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
  (2*h.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (h.diff(d01,i,j)*Q22.diff(d01,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))) + 
  (exp((-2*mu*h(i,j))/(3.*A1))*
     ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
        exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j)) + 
       2*(1 + Qtt(i,j))*(((2 + 3*pow(A1,2))*mu*(-12 + pow(mu,2))*
             exp((2*mu*h(i,j))/(3.*A1))*h(i,j))/2. + 
          12*A1*(-1 + exp(((2*mu)/(3.*A1) + A1*mu)*h(i,j)))*(1 + Qtt(i,j))))\
)/((2 + 3*pow(A1,2))*mu*(-12 + pow(mu,2))*(1 + Qtt(i,j)))
;
