
elem[0]=
(complex(0,1)*a0.diff(d01,i,j)*mu*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2))/(L*(1 + rgrid[i]*Q11(i,j))) + 
  (complex(0,1)*a0.diff(d01,i,j)*mu*w*Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-1 + rgrid[i])*pow(rgrid[i],2))/(L*(1 + rgrid[i]*Q11(i,j))) - 
  (h11.diff(d01,i,j)*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3)))/
   (4.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
  (A0.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(-1 + rgrid[i],2)*
     pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (2.*L*(1 + rgrid[i]*Qrr(i,j))) + 
  (A0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (2.*L*(-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (A0.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(-1 + rgrid[i],2)*
     pow(rgrid[i],4)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
  (complex(0,0.5)*ht1.diff(d01,i,j)*w*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))/
   (L*(1 + rgrid[i]*Q11(i,j))) - 
  (htt.diff(d20,i,j)*pow(-1 + rgrid[i],2)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j)))/
   (8.*(1 + rgrid[i]*Qrr(i,j))) + 
  (htt.diff(d11,i,j)*pow(-1 + rgrid[i],2)*rgrid[i]*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*(1 + rgrid[i]*Qtt(i,j)))/
   (4.*L*(1 + rgrid[i]*Qrr(i,j))) - 
  (htt.diff(d02,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))))*(1 + rgrid[i]*Qtt(i,j)))/
   (8.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. - 
       (L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Qtt(i,j)))/2. - 
            ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                 3*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))*
               (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
   (4.*L*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
     (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. + 
       (L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Qtt(i,j)))/2. - 
            ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                 3*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))*
               (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
   (4.*pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
  (Dh(i,j)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*(-1 + rgrid[i])*
     pow(rgrid[i],5)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((-2*pow(a0.diff(d01,i,j),2)*(2 + 3*pow(A1,2))*B1*pow(mu,2)*
          exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + B1*mu*h(i,j)*rgrid[i])*
          (-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) - 
       ((1 + rgrid[i]*Q11(i,j))*
          ((2 + 3*pow(A1,2))*B1*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + B1*mu*h(i,j)*rgrid[i])*
             pow(-(L*(-(mu*a0(i,j)) - 
                    a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
               a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j),2) \
+ (24*A1*pow(L,2)*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)*rgrid[i]))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],4)))/pow(rgrid[i],2)))/
   (4.*(2 + 3*pow(A1,2))*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) - 
  (complex(0,0.25)*w*ht1(i,j)*(-1 + rgrid[i])*pow(rgrid[i],8)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],7) + ((1 + rgrid[i]*Q22(i,j))*
          ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)))/
   (L*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (htt.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((-2*Qrr.diff(d01,i,j)*Qr1(i,j))/
             ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (4*L*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qrr(i,j)))/2.))/
             (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)))*
          (1 + rgrid[i]*Qtt(i,j)))/(2.*L*pow(rgrid[i],2)) + 
       (2*(1 + rgrid[i]*Qrr(i,j))*
          ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/L + 
            ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               ((2*Qr1.diff(d01,i,j)*rgrid[i])/L - 
                 (complex(0,24)*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                 (pow(rgrid[i],2)*
                    (-((-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))/pow(rgrid[i],3)) + 
                      (Q11.diff(d01,i,j)*Qr1(i,j))/L))/
                  (1 + rgrid[i]*Q11(i,j)) + 
                 (pow(rgrid[i],2)*
                    (-((-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j))/pow(rgrid[i],3)) + 
                      (Q22.diff(d01,i,j)*Qr1(i,j))/L))/
                  (1 + rgrid[i]*Q22(i,j)))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],2)) - 
            (2*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qtt(i,j)))/2. - 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (htt.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
        ((-2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
           pow(rgrid[i],3) + 
          (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
           pow(rgrid[i],3) - 
          (Q11.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
             (1 + rgrid[i]*Qtt(i,j)))/
           (pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))) - 
          (Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
             (1 + rgrid[i]*Qtt(i,j)))/
           (pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j)))) + 
       (2*pow(rgrid[i],4)*(1 + rgrid[i]*Qrr(i,j))*
          (-((Q22.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                 (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
               pow(rgrid[i],7)) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Q11(i,j))*
                    ((-2*Qrr.diff(d01,i,j))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*L*(1 + rgrid[i]*Q11(i,j))*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                       pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                  (2.*pow(rgrid[i],4)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (-((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                      (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)\
))/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
          (1 + rgrid[i]*Q22(i,j))) + 
       rgrid[i]*Qr1(i,j)*((-2*L*
             (((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Qrr(i,j)))/2. + 
               ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                    3*(-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                  (1 + rgrid[i]*Qrr(i,j)))/2.)*(1 + rgrid[i]*Qtt(i,j)))/
           ((-1 + rgrid[i])*pow(rgrid[i],5)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (2*(1 + rgrid[i]*Qrr(i,j))*
             (((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (-4*Qr1.diff(d01,i,j)*rgrid[i] + 
                    (complex(0,24)*L*w*pow(rgrid[i],2))/
                     ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                    (L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/
                     (rgrid[i]*(1 + rgrid[i]*Q11(i,j))) + 
                    (L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/
                     (rgrid[i]*(1 + rgrid[i]*Q22(i,j))))*
                  (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
               (2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
           ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))))))/
   (16.*pow(L,2)*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h11(i,j)*pow(rgrid[i],10)*
     ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        (2.*pow(rgrid[i],8)) + 
       ((1 + rgrid[i]*Q11(i,j))*
          (-(Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*
                (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],6)) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (-(Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                     (-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                     (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],4)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((pow(Qtt.diff(d01,i,j),2)*
                        pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2))/
                       (4.*pow(rgrid[i],2)) + 
                      ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-1 + rgrid[i],2) - 
                        (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/pow(rgrid[i],2)))/
   (4.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) - 
  (h22(i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((complex(0,12)*L*w*pow(rgrid[i],2)*
          (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. + 
            (L*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qtt(i,j)))/2. - 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))) + 
       (2*pow(rgrid[i],8)*((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
               pow(-1 + rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],8)) - 
            (L*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qtt(i,j))*
               ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. - 
                 (4*L*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                  (pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 (L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3))\
)/pow(rgrid[i],8) + ((-1 + rgrid[i])*(1 + rgrid[i]*Q11(i,j))*
               (-(Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                     (-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  (2.*pow(rgrid[i],6)) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    (-(Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],4)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((pow(Qtt.diff(d01,i,j),2)*
                        pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2))/
                        (4.*pow(rgrid[i],2)) + 
                        ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2) - 
                        (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
          (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/
   (8.*pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
  (htt(i,j)*((12*(3*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
            2*exp(A1*mu*h(i,j)*rgrid[i]))*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))/
        ((2 + 3*pow(A1,2))*pow(rgrid[i],2)) + 
       (9*w*(complex(0,1)*pow(mu,2) + 2*(complex(0,-6) + w))*
          pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j)))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)*
          (1 + rgrid[i]*Qrr(i,j))) - 
       (complex(0,6)*w*pow(-1 + rgrid[i],2)*rgrid[i]*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (pow(rgrid[i],6)*((pow(Qtt.diff(d01,i,j),2)*
               (-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j)))/
             (2.*pow(rgrid[i],4)) + 
            ((1 + rgrid[i]*Q11(i,j))*
               pow((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. - 
                 (L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3),
                2))/pow(rgrid[i],2)))/
        (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j))) + 
       (complex(0,1.5)*w*pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          (((1 + rgrid[i]*Q11(i,j))*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                      (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    ((2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (4*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          ((2*Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
               (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2))/
             ((-1 + rgrid[i])*pow(rgrid[i],8)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) - 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                    (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],4) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                       pow(rgrid[i],2) + 
                      (Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        pow(Qr1(i,j),2))/2. + 
                      (2*Qtt.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       pow(rgrid[i],3) + 
                      (pow(L,2)*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6) \
- L*rgrid[i]*Qr1(i,j)*((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) + 
                         (Q11.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4))\
))/pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  ((Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],4) - 
                    (2*((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/
                        (2.*pow(rgrid[i],2)) + 
                        (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  ((2*L*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) + 
                        (Q22.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) - 
                    ((1 + rgrid[i]*Q22(i,j))*
                       ((4*Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],2) - 
                         (8*L*
                        ((Qtt.diff(d11,i,j)*(-1 + rgrid[i])*
                       rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                        Qtt.diff(d01,i,j)*
                        (-((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                        (-12 + pow(mu,2)*pow(rgrid[i],4) + 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2.))*
                        (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*Qtt.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*Qrr.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                         ((-1 + rgrid[i])*pow(rgrid[i],4)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2)) + 
                 L*((-2*L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                     ((-1 + rgrid[i])*pow(rgrid[i],8)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                    ((1 + rgrid[i]*Q22(i,j))*
                       ((4*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                         (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                         (4*(1 + rgrid[i]*Qrr(i,j))*
                         ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        (2.*rgrid[i]) + 
                         (Qr1.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2) \
- (L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (6 - 
                        2*Qtt.diff(d10,i,j)*pow(rgrid[i],2) + 
                        Qtt.diff(d20,i,j)*pow(rgrid[i],3) + 
                        2*rgrid[i]*Qtt(i,j)))/2. + 
                        rgrid[i]*
                        (-(((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/rgrid[i]) + 
                        (3*(-4 + pow(mu,2)*pow(rgrid[i],4) + 
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i])))/
                         pow(rgrid[i],4)))/
                         ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],4)))/
        (4.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
          (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2))))/4.
;
elem[1]=
complex(0,-1)*h22.diff(d10,i,j)*w*(-1 + rgrid[i]) + 
  complex(0,1)*mu*w*Dh(i,j)*rgrid[i]*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
  (complex(0,0.5)*h11.diff(d01,i,j)*w*rgrid[i]*Qr1(i,j))/L + 
  (complex(0,0.5)*h22.diff(d01,i,j)*w*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))/
   L + (Ax.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],4)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
  (Ax.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],5)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
  (muj*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],4)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (complex(0,-6)*L*w*pow(rgrid[i],2) + 
       Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*rgrid[i]*
        (-1 + pow(rgrid[i],3)))*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] + (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (2.*pow(L,2)*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
  (Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],4)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (complex(0,-6)*L*w*pow(rgrid[i],2) + 
       Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*rgrid[i]*
        (-1 + pow(rgrid[i],3)))*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] + (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (2.*pow(L,2)*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
     (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
  (ht1.diff(d20,i,j)*pow(-1 + rgrid[i],2)*rgrid[i]*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*(1 + rgrid[i]*Qtt(i,j)))/
   (8.*(1 + rgrid[i]*Qrr(i,j))) + 
  (ht1.diff(d11,i,j)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (pow(rgrid[i],2)/(1 + rgrid[i]*Q11(i,j)) + 
       ((-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2))/
        (1 + rgrid[i]*Qrr(i,j)))*(1 + rgrid[i]*Qtt(i,j)))/
   (4.*L*pow(rgrid[i],2)) - (ht1.diff(d02,i,j)*pow(-1 + rgrid[i],2)*
     pow(rgrid[i],3)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))))*(1 + rgrid[i]*Qtt(i,j)))/
   (8.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
  (complex(0,0.25)*w*h11(i,j)*pow(rgrid[i],8)*
     (((1 + rgrid[i]*Q11(i,j))*
          ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) + 
            Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
       ((1 + rgrid[i]*Q22(i,j))*
          (((-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/
                  rgrid[i])*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
            rgrid[i]*Qr1(i,j)*
             ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                  (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
               (2*(1 + rgrid[i]*Qrr(i,j))*
                  (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],3)) \
+ (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))))/
                ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))))))/pow(rgrid[i],2)))/
   (L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (ht1.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((pow(rgrid[i],2)*pow(Qr1(i,j),2)*
          ((3*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
             pow(rgrid[i],3) - 
            (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],3) - 
            (Q11.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/
             (pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))) + 
            (Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/
             (pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j)))))/L + 
       (2*pow(rgrid[i],4)*(1 + rgrid[i]*Qrr(i,j))*
          ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],7) - 
            ((1 + rgrid[i]*Q22(i,j))*
               ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((-3*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                       (2.*pow(rgrid[i],3)) + 
                      (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)\
))/(L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
          (1 + rgrid[i]*Q22(i,j))) + 
       rgrid[i]*Qr1(i,j)*((2*(((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Qrr(i,j)))/2. + 
               ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                    3*(-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                  (1 + rgrid[i]*Qrr(i,j)))/2.)*(1 + rgrid[i]*Qtt(i,j)))/
           ((-1 + rgrid[i])*pow(rgrid[i],5)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (2*(1 + rgrid[i]*Qrr(i,j))*
             (((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  ((2*Qr1.diff(d01,i,j)*rgrid[i])/L - 
                    (complex(0,24)*w*pow(rgrid[i],2))/
                     ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))/
                     (rgrid[i]*(1 + rgrid[i]*Q11(i,j))) - 
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j))/
                     (rgrid[i]*(1 + rgrid[i]*Q22(i,j))))*
                  (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) - 
               (3*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
           ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))))))/
   (16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (ht1.diff(d01,i,j)*((Q11.diff(d01,i,j)*(-1 + rgrid[i])*
          pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
          (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(1 + rgrid[i]*Q11(i,j),2)) + 
       (pow(rgrid[i],2)*((-3*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. + 
            ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-4*Qr1.diff(d01,i,j)*rgrid[i] + 
                 (complex(0,12)*L*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
                 (L*(-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2))/
                  (2.*(1 + rgrid[i]*Qrr(i,j))) + 
                 (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                    pow(rgrid[i],4)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    pow(Qr1(i,j),3))/(2.*(1 + rgrid[i]*Qrr(i,j))) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((Q22.diff(d01,i,j)*rgrid[i])/
                       (1 + rgrid[i]*Q22(i,j))) + 
                    (Qrr.diff(d01,i,j)*rgrid[i])/
                     (1 + rgrid[i]*Qrr(i,j))))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],2))))/(1 + rgrid[i]*Q11(i,j)) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],5)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*
          ((2*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],4) + 
            pow(rgrid[i],2)*pow(Qr1(i,j),2)*
             ((-3*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                pow(rgrid[i],3) + 
               (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],3) - 
               (Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                  (1 + rgrid[i]*Qtt(i,j)))/
                (pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j)))) + 
            rgrid[i]*Qr1(i,j)*
             ((-2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                  (1 + rgrid[i]*Qtt(i,j)))/
                ((-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) + 
               (2*(1 + rgrid[i]*Qrr(i,j))*
                  (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-4*Qr1.diff(d01,i,j)*rgrid[i] + 
                        (complex(0,24)*L*w*pow(rgrid[i],2))/
                        ((-12 + pow(mu,2))*
                       (-1 + pow(rgrid[i],3))) + 
                        (L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/
                        (rgrid[i]*(1 + rgrid[i]*Q22(i,j))))*
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
                    (3*L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
                ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))))))/
        (4.*pow(1 + rgrid[i]*Qrr(i,j),2))))/(4.*pow(L,2)) + 
  (w*h22(i,j)*((24*w*pow(rgrid[i],2))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))) - 
       complex(0,1)*(4 + ((1 - rgrid[i])*pow(rgrid[i],8)*
             (((1 + rgrid[i]*Q11(i,j))*
                  (-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                    Q22.diff(d01,i,j)*Qr1(i,j))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],6) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  (rgrid[i]*Qr1(i,j)*
                     ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                       (2*(1 + rgrid[i]*Qrr(i,j))*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       (-(L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],5)) + 
                         (2*(1 + rgrid[i]*Q11(i,j))*
                        ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/pow(rgrid[i],2)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                pow(rgrid[i],2)))/
           (L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
             (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))))))/4. + 
  (ht1(i,j)*((9*w*(complex(0,1)*pow(mu,2) + 2*(complex(0,-6) + w))*
          pow(-1 + rgrid[i],2)*pow(rgrid[i],5)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*
          (1 + rgrid[i]*Qtt(i,j)))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)*
          (1 + rgrid[i]*Qrr(i,j))) - 
       (complex(0,6)*w*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*
          (1 + rgrid[i]*Qtt(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],8)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (-((Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                 (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                 (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6)) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                    (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (((-4*a0.diff(d01,i,j)*L*mu*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                       (-1 + rgrid[i]) - 
                        (3*Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2.)*
                        (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) + 
                      ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Q11.diff(d01,i,j)*Qr1.diff(d01,i,j) - 
                        (2*Qr1.diff(d02,i,j)*
                        (1 + rgrid[i]*Q11(i,j)))/rgrid[i])*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/(2.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))) - 
       ((-1 + rgrid[i])*pow(rgrid[i],13)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),3)*
          ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],8)) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                    (-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  (2.*pow(rgrid[i],6)) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    (-(Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],4)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                        (-(pow(Qtt.diff(d01,i,j),2)*
                        pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2))/
                        (4.*pow(rgrid[i],2)) + 
                        ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2) + 
                        (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j))) + 
       (complex(0,1.5)*w*pow(-1 + rgrid[i],2)*pow(rgrid[i],12)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          ((-2*Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],9)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((3*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       pow(rgrid[i],3) - 
                      (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                       pow(rgrid[i],3) + 
                      (L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*Qr1(i,j)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       (2.*pow(rgrid[i],4)) - 
                      (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
                       (2.*rgrid[i])))/pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*Qr1(i,j)*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                      (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((3*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (3*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],3)))/
        (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)) + 
       ((-1 + rgrid[i])*pow(rgrid[i],12)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
          (((1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j))*
               ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) - 
                 (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                    (-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
                 (Q22.diff(d01,i,j)*L*
                    (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4))\
)/pow(rgrid[i],6) + ((1 + rgrid[i]*Q22(i,j))*
               (((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j))*
                    ((Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) - 
                      L*((2*Qtt.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (2*Qrr.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))))))/
                  (2.*pow(rgrid[i],4)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)*
                        (Q11.diff(d01,i,j)*Qr1.diff(d01,i,j) - 
                        (2*Qr1.diff(d02,i,j)*
                       (1 + rgrid[i]*Q11(i,j)))/rgrid[i])*
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                       (4.*pow(rgrid[i],4)) - 
                      (Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j))*
                        ((((-7*Qr1.diff(d01,i,j)*Qtt.diff(d01,\
i,j)*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                       4*L*(-(a0.diff(d01,i,j)*mu*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                       (-1 + rgrid[i])) + 
                       ((Qtt.diff(d11,i,j)*(-1 + rgrid[i])*
                      rgrid[i]*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))/2. + 
                       Qtt.diff(d01,i,j)*
                       (-((-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))/2. + 
                       (-12 + pow(mu,2)*pow(rgrid[i],4) + 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2.))/
                       pow(rgrid[i],2)))*(1 + rgrid[i]*Q11(i,j)))/
                        pow(rgrid[i],2) + 
                        L*((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) + 
                        (Q11.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4)\
)))/(2.*pow(rgrid[i],2))))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/(2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j))) + 
       rgrid[i]*Qr1(i,j)*(((-1 + rgrid[i])*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             ((72*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
                (2 + 3*pow(A1,2)) + 
               (48*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) + 
               ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))/
                (L*(1 + rgrid[i]*Qrr(i,j))) - 
               (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q11(i,j)))/
                (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
               (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q22(i,j)))/
                (2.*L*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
               (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
                  (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                (L*pow(1 + rgrid[i]*Qrr(i,j),2)))*
             (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
          (2*pow(rgrid[i],2)*
             ((pow(Qtt.diff(d01,i,j),2)*pow(-1 + rgrid[i],2)*
                  pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3),2))/
                (4.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
               ((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  pow(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.,2))/
                (2.*pow(rgrid[i],4)*(1 + rgrid[i]*Qrr(i,j)))))/
           ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
          (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
             pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3),2)*
             ((2*Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                  (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2))/
                ((-1 + rgrid[i])*pow(rgrid[i],8)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) - 
               (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                  ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],4) + 
                    ((1 + rgrid[i]*Q22(i,j))*
                       ((Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                         (4*(pow(a0.diff(d01,i,j),2)*
                       pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2) + 
                        (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]))*
                        (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (pow(L,2)*
                        (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6))\
)/pow(rgrid[i],2)))/
                ((-1 + rgrid[i])*pow(rgrid[i],4)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) - 
               (L*pow(1 + rgrid[i]*Q11(i,j),2)*
                  ((2*L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                     ((-1 + rgrid[i])*pow(rgrid[i],8)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                    ((1 + rgrid[i]*Q22(i,j))*
                       ((-4*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                         (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                         (2*L*exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                         (Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i] - 
                         (5*Qr1.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2) \
+ (2*L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (6 - 
                        2*Qtt.diff(d10,i,j)*pow(rgrid[i],2) + 
                        Qtt.diff(d20,i,j)*pow(rgrid[i],3) + 
                        2*rgrid[i]*Qtt(i,j)))/2. + 
                        rgrid[i]*
                        (-(((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/rgrid[i]) + 
                        (3*(-4 + pow(mu,2)*pow(rgrid[i],4) + 
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i])))/
                         pow(rgrid[i],4)))/
                         ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2)))/pow(rgrid[i],4)))/
           (4.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
             (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)))))/4.
;
elem[2]=
complex(0,-0.5)*h11.diff(d01,i,j)*w - 
  complex(0,0.5)*h22.diff(d01,i,j)*w*(-1 + rgrid[i]) + 
  complex(0,1)*h.diff(d01,i,j)*mu*w*Dh(i,j)*pow(rgrid[i],2) + 
  (Ax.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (2.*(1 + rgrid[i]*Qrr(i,j))) + 
  (muj*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (complex(0,-6)*L*w*pow(rgrid[i],2) + 
       Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*rgrid[i]*
        (-1 + pow(rgrid[i],3)))*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (2.*L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (complex(0,-6)*L*w*pow(rgrid[i],2) + 
       Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*rgrid[i]*
        (-1 + pow(rgrid[i],3)))*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (2.*L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (Ax.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],3)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (2.*L*(1 + rgrid[i]*Qrr(i,j))) - 
  (ht1.diff(d20,i,j)*L*pow(-1 + rgrid[i],2)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j)))/
   (8.*(1 + rgrid[i]*Qrr(i,j))) + 
  (ht1.diff(d11,i,j)*pow(-1 + rgrid[i],2)*rgrid[i]*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*(1 + rgrid[i]*Qtt(i,j)))/
   (4.*(1 + rgrid[i]*Qrr(i,j))) - 
  (ht1.diff(d02,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j)))/(8.*L*(1 + rgrid[i]*Qrr(i,j))) - 
  complex(0,0.25)*w*h11(i,j)*((2*Q22.diff(d01,i,j)*rgrid[i])/
      (1 + rgrid[i]*Q22(i,j)) + 
     (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) - 
     (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))) - 
  complex(0,0.25)*w*h22(i,j)*(1 - rgrid[i])*
   ((Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) + 
     (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))) + 
  (ht1.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((-2*Qrr.diff(d01,i,j)*Qr1(i,j))/
             ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (4*L*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qrr(i,j)))/2.))/
             (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)))*
          (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
       (2*(1 + rgrid[i]*Qrr(i,j))*
          ((3*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. + 
            ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (2*Qr1.diff(d01,i,j)*rgrid[i] - 
                 (complex(0,24)*L*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                 (pow(rgrid[i],2)*
                    ((L*(-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                      Q11.diff(d01,i,j)*Qr1(i,j)))/
                  (1 + rgrid[i]*Q11(i,j)) + 
                 (pow(rgrid[i],2)*
                    (-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j)))/
                  (1 + rgrid[i]*Q22(i,j)))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],2)) - 
            (3*L*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qtt(i,j)))/2. - 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (ht1.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((2*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],4) + (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
          ((-3*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
             pow(rgrid[i],3) + 
            (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],3) + 
            (Q11.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/
             (pow(rgrid[i],3)*(1 + rgrid[i]*Q11(i,j))) - 
            (Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/
             (pow(rgrid[i],3)*(1 + rgrid[i]*Q22(i,j)))))/L + 
       rgrid[i]*Qr1(i,j)*((-2*
             (((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Qrr(i,j)))/2. + 
               ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                    3*(-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                  (1 + rgrid[i]*Qrr(i,j)))/2.)*(1 + rgrid[i]*Qtt(i,j)))/
           ((-1 + rgrid[i])*pow(rgrid[i],5)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (2*(1 + rgrid[i]*Qrr(i,j))*
             (((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  ((-4*Qr1.diff(d01,i,j)*rgrid[i])/L + 
                    (complex(0,24)*w*pow(rgrid[i],2))/
                     ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))/
                     (rgrid[i]*(1 + rgrid[i]*Q11(i,j))) + 
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j))/
                     (rgrid[i]*(1 + rgrid[i]*Q22(i,j))))*
                  (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
               (3*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
           ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))))))/
   (16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (ht1(i,j)*((9*L*w*(complex(0,1)*pow(mu,2) + 2*(complex(0,-6) + w))*
          pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j)))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)*
          (1 + rgrid[i]*Qrr(i,j))) - 
       (complex(0,6)*L*w*pow(-1 + rgrid[i],2)*rgrid[i]*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((72*pow(A1,2)*L*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
             (2 + 3*pow(A1,2)) + 
            (48*L*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) + 
            (Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*(-1 + rgrid[i])*
               pow(rgrid[i],5)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
             (2.*L*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
            ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))/(1 + rgrid[i]*Qrr(i,j)) \
- (Qr1.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
             (L*(1 + rgrid[i]*Qrr(i,j))) + 
            (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],5)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 Q11.diff(d01,i,j)*Qr1(i,j)))/
             (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
            (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],5)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j)))/
             (2.*L*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
            (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
               (((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qrr(i,j)))/2.))/
             pow(1 + rgrid[i]*Qrr(i,j),2))*(1 + rgrid[i]*Qtt(i,j)))/
        (2.*pow(rgrid[i],2)) + 
       (pow(rgrid[i],6)*((pow(Qtt.diff(d01,i,j),2)*
               (-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j)))/
             (2.*pow(rgrid[i],4)) + 
            ((1 + rgrid[i]*Q11(i,j))*
               pow((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. - 
                 (L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3),
                2))/pow(rgrid[i],2)))/
        (L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j))) + 
       (complex(0,1.5)*w*pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          (((1 + rgrid[i]*Q11(i,j))*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) - 
                 rgrid[i]*Qr1(i,j)*
                  ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((-3*Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    ((2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (3*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          ((2*Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
               (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2))/
             ((-1 + rgrid[i])*pow(rgrid[i],8)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (-((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                      (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],4)) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    (-((Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                        pow(rgrid[i],2)) - 
                      (Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        pow(Qr1(i,j),2))/2. + 
                      (4*(pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2) - 
                        (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]))*
                        (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) - 
                      (pow(L,2)*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6) \
+ L*rgrid[i]*Qr1(i,j)*((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) + 
                         (Q11.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4))\
))/pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               ((-2*(-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*
                    ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. - 
                      (L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                       ((Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                        pow(rgrid[i],2) - 
                         (4*(pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-1 + rgrid[i],2) + 
                        (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]))*
                        (1 + rgrid[i]*Qrr(i,j)))/
                         ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) - 
                      rgrid[i]*Qr1(i,j)*
                       ((2*((7*Qr1.diff(d01,i,j)*Qtt.diff(d01,i\
,j)*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. - 
                        4*L*
                        (-(a0.diff(d01,i,j)*mu*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                       (-1 + rgrid[i])) + 
                        ((Qtt.diff(d11,i,j)*(-1 + rgrid[i])*
                       rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                       Qtt.diff(d01,i,j)*
                       (-((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                       (-12 + pow(mu,2)*pow(rgrid[i],4) + 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2.))/
                        pow(rgrid[i],2)))*(1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         L*((2*Qtt.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*Qrr.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                         ((-1 + rgrid[i])*pow(rgrid[i],4)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))) - 
                      L*((-4*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                         (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                         (2*L*exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                         (Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i] - 
                         (5*Qr1.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2) \
+ (2*L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (6 - 
                        2*Qtt.diff(d10,i,j)*pow(rgrid[i],2) + 
                        Qtt.diff(d20,i,j)*pow(rgrid[i],3) + 
                        2*rgrid[i]*Qtt(i,j)))/2. + 
                        rgrid[i]*
                        (-(((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/rgrid[i]) + 
                        (3*(-4 + pow(mu,2)*pow(rgrid[i],4) + 
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i])))/
                         pow(rgrid[i],4)))/
                         ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],4)))/
        (4.*L*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2))))/4.
;
elem[3]=
htt.diff(d20,i,j)/2. - 2*Dh.diff(d10,i,j)*mu*rgrid[i]*
   (h(i,j) + h.diff(d10,i,j)*rgrid[i]) + 
  complex(0,1)*ht1.diff(d10,i,j)*w*rgrid[i]*Qr1(i,j) - 
  (h11.diff(d20,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(Qr1(i,j),2))/(4.*(1 + rgrid[i]*Qrr(i,j))) + 
  (h11.diff(d11,i,j)*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/(2.*L*(1 + rgrid[i]*Qrr(i,j))) \
- (h22.diff(d11,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/(2.*L*(1 + rgrid[i]*Qrr(i,j))) \
- (h11.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (4.*pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d02,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (4.*pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d20,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (4*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/(4.*(1 + rgrid[i]*Qrr(i,j))) + 
  (htt.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (2*Qr1.diff(d01,i,j)*pow(1 + rgrid[i]*Q11(i,j),2)*
        pow(Qr1(i,j),3) + (2*
          ((2*Qrr.diff(d01,i,j))/
             ((-1 + rgrid[i])*rgrid[i]*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) - 
            (2*L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q11(i,j))*Qr1(i,j))/pow(rgrid[i],2) + 
            Q11.diff(d01,i,j)*rgrid[i]*pow(Qr1(i,j),2))*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((2*Qrr.diff(d01,i,j)*rgrid[i]*pow(Qr1(i,j),2))/
             ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) - 
            L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
               rgrid[i]*Q11(i,j))*pow(Qr1(i,j),3) + 
            Q11.diff(d01,i,j)*pow(rgrid[i],3)*pow(Qr1(i,j),4) - 
            (4*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
               (1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            rgrid[i]*Qr1(i,j)*
             ((4*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                ((-1 + rgrid[i])*rgrid[i]*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) + 
               (4*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                  pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3),2)))))/pow(rgrid[i],2))\
)/(8.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
  (complex(0,2)*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],4)*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-2*L*(-(mu*a0(i,j)) - 
               a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,2)*w*Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],4)*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-2*L*(-(mu*a0(i,j)) - 
               a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],6)*(((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))))*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(1 - rgrid[i])*
     pow(rgrid[i],4)*((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
        (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
          a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)) + 
       (2*(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (L*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) - 
  (A0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],4)*
     (-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
        (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
          a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)) + 
       (2*(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (L*(-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (htt.diff(d10,i,j)*((complex(0,24)*w*pow(rgrid[i],2))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
       (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
        (L*(1 + rgrid[i]*Qrr(i,j))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i])*
          pow(Qr1(i,j),2))/(2.*L*(1 + rgrid[i]*Qrr(i,j))) - 
       (Q11.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),3))/
        (2.*L*(1 + rgrid[i]*Qrr(i,j))) - 
       (2*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Qrr(i,j)))/2. + 
            ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                 3*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))*
               (1 + rgrid[i]*Qrr(i,j)))/2.))/
        ((-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) + 
       (4*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Qtt(i,j)))/2. - 
            ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                 3*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))*
               (1 + rgrid[i]*Qtt(i,j)))/2.))/
        ((-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))/4. + 
  (w*ht1(i,j)*((complex(0,2)*Qrr.diff(d01,i,j)*rgrid[i])/
        (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
       (complex(0,1)*Q11.diff(d01,i,j)*pow(rgrid[i],3)*
          pow(Qr1(i,j),2))/(L*(1 + rgrid[i]*Q11(i,j))) + 
       rgrid[i]*Qr1(i,j)*((complex(0,2)*Qr1.diff(d01,i,j)*rgrid[i])/L - 
          (12*w*pow(rgrid[i],2))/
           ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
          (complex(0,2)*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
               rgrid[i]*Q11(i,j)))/(rgrid[i]*(1 + rgrid[i]*Q11(i,j))) + 
          (complex(0,4)*(((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Qtt(i,j)))/2. - 
               ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                    3*(-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                  (1 + rgrid[i]*Qtt(i,j)))/2.))/
           ((-1 + rgrid[i])*rgrid[i]*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))))))/2. \
+ Dh(i,j)*(-2*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
     (complex(0,12)*mu*w*pow(rgrid[i],3)*
        (h(i,j) + h.diff(d10,i,j)*rgrid[i]))/
      ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
     (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],7)*
        ((-4*pow(a0.diff(d01,i,j),2)*(2 + 3*pow(A1,2))*B1*
             pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + 
               B1*mu*h(i,j)*rgrid[i])*pow(1 + rgrid[i]*Qrr(i,j),2))/
           (pow(rgrid[i],4)*
             pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3),2)) - 
          (pow(1 + rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2)*
             ((2 + 3*pow(A1,2))*B1*
                exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + 
                  B1*mu*h(i,j)*rgrid[i])*
                pow(-(L*(-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                  a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                   Qr1(i,j),2) - 
               (24*A1*pow(L,2)*
                  (-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)*rgrid[i]))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],4)))/pow(rgrid[i],2) + 
          (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
             (-((2 + 3*pow(A1,2))*B1*
                  exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + 
                    B1*mu*h(i,j)*rgrid[i])*
                  (-(pow(L,2)*
                       pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2)) + 
                    2*a0.diff(d01,i,j)*L*mu*
                     (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                     (-1 + rgrid[i])*rgrid[i]*Qr1(i,j) + 
                    2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                     pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                     pow(Qr1(i,j),2))) + 
               (24*A1*pow(L,2)*
                  (-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)*rgrid[i]))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],4)))/
           ((-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3)))))/
      (2.*(2 + 3*pow(A1,2))*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))) + 
  (htt(i,j)*((complex(0,-18)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          pow(rgrid[i],4))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)) + 
       (complex(0,12)*w*rgrid[i])/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
       (exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],6)*
          ((pow(1 + rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2)*
               pow(-(L*(-(mu*a0(i,j)) - 
                      a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                 a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                  Qr1(i,j),2))/pow(rgrid[i],2) + 
            (2*(1 + rgrid[i]*Q11(i,j))*
               (-(pow(L,2)*
                    pow(-(mu*a0(i,j)) - 
                      a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2)) + 
                 2*a0.diff(d01,i,j)*L*mu*
                  (-(mu*a0(i,j)) - 
                    a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                  (-1 + rgrid[i])*rgrid[i]*Qr1(i,j) + 
                 2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                  pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                  pow(Qr1(i,j),2))*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
               pow(1 + rgrid[i]*Qrr(i,j),2))/
             (pow(rgrid[i],4)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2))))/
        (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j))) - 
       (complex(0,3)*w*pow(rgrid[i],6)*
          ((Qrr.diff(d01,i,j)*Qr1(i,j)*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],2) + 
            ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/
                  rgrid[i])*pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/2. \
+ (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),3)*
               (1 + rgrid[i]*Qtt(i,j)))/2. + 
            L*((2*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                  (1 + rgrid[i]*Qtt(i,j)))/
                ((-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) - 
               (4*(1 + rgrid[i]*Qrr(i,j))*
                  (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/
                ((-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))))))/
        (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/2. + 
  (h11.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((8*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            Q11.diff(d01,i,j)*Qr1(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2))/
        (L*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)) + 
       (4*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          ((2*((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (-(L*(-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                     (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
                     (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3)) + 
                 (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    pow(Qr1(i,j),3)*(1 + rgrid[i]*Qtt(i,j)))/2. + 
                 (4*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3) + 
                 rgrid[i]*Qr1(i,j)*
                  ((2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                     pow(rgrid[i],3) - 
                    (2*Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],3))))/pow(rgrid[i],2)))/
        (L*pow(-1 + rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qtt(i,j))) + 
       (pow(1 + rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2)*
          ((-2*Qrr.diff(d01,i,j)*Qr1(i,j))/
             (L*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (4*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qrr(i,j)))/2.))/
             (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               ((2*Qr1.diff(d01,i,j)*rgrid[i])/L - 
                 (complex(0,24)*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                 (pow(rgrid[i],2)*
                    (-((-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))/pow(rgrid[i],3)) + 
                      (Q22.diff(d01,i,j)*Qr1(i,j))/L))/
                  (1 + rgrid[i]*Q22(i,j)) + 
                 (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
                  (L*(1 + rgrid[i]*Qtt(i,j))) - 
                 (2*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)))/
   (16.*(1 + rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h11.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((4*((4*Qrr.diff(d01,i,j))/
             ((-1 + rgrid[i])*rgrid[i]*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) - 
            Q11.diff(d01,i,j)*rgrid[i]*pow(Qr1(i,j),2))*
          pow(1 + rgrid[i]*Qrr(i,j),2))/
        (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))) + 
       (2*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
          (rgrid[i]*Qr1(i,j)*
             ((6*Qrr.diff(d01,i,j))/
                ((-1 + rgrid[i])*rgrid[i]*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) + 
               (L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q11(i,j))*Qr1(i,j))/pow(rgrid[i],2) - 
               Q11.diff(d01,i,j)*rgrid[i]*pow(Qr1(i,j),2)) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               (-4*Qr1.diff(d01,i,j)*rgrid[i] + 
                 (complex(0,24)*L*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
                 (Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
                  (1 + rgrid[i]*Q22(i,j)) - 
                 (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
                  (1 + rgrid[i]*Qtt(i,j))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
        ((4*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
             (1 + rgrid[i]*Qrr(i,j)))/
           ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          pow(rgrid[i],2)*pow(Qr1(i,j),2)*
           ((2*Qrr.diff(d01,i,j))/
              ((-1 + rgrid[i])*rgrid[i]*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))) - 
             (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
              ((-1 + rgrid[i])*rgrid[i]*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))) \
- (2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
              ((-1 + rgrid[i])*rgrid[i]*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) \
+ rgrid[i]*Qr1(i,j)*((-4*L*(((-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                     (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                  ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                     (1 + rgrid[i]*Qrr(i,j)))/2.))/
              (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3),2)) + 
             (2*(1 + rgrid[i]*Qrr(i,j))*
                (-4*Qr1.diff(d01,i,j)*rgrid[i] + 
                  (complex(0,24)*L*w*pow(rgrid[i],2))/
                   ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                  (L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/
                   (rgrid[i]*(1 + rgrid[i]*Q22(i,j))) + 
                  (2*L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                         (1 + rgrid[i]*Qtt(i,j)))/2.))/
                   ((-1 + rgrid[i])*rgrid[i]*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                     (1 + rgrid[i]*Qtt(i,j)))))/
              ((-1 + rgrid[i])*pow(rgrid[i],2)*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3)))))))/
   (16.*pow(L,2)*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  h22.diff(d10,i,j)*(2 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2))/(2.*(1 + rgrid[i]*Qrr(i,j))) + 
     (complex(0,3)*w*(-1 + rgrid[i])*pow(rgrid[i],4)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
          (4*(1 + rgrid[i]*Qrr(i,j)))/
           ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3)))))/
      ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
        (1 + rgrid[i]*Qrr(i,j))) - 
     (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        ((4*(1 + rgrid[i]*Q22(i,j))*
             (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
               Q11.diff(d01,i,j)*Qr1(i,j))*
             pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
           ((-1 + rgrid[i])*pow(rgrid[i],8)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
             ((2*(-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                    Q22.diff(d01,i,j)*Qr1(i,j))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],4) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  ((2*Qtt.diff(d01,i,j)*Qr1(i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],2) - 
                    (3*L*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3)) \
+ (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       pow(Qr1(i,j),3)*(1 + rgrid[i]*Qtt(i,j)))/2. + 
                    ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       ((4*Qr1.diff(d01,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)))/
                pow(rgrid[i],2)))/
           ((-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (pow(1 + rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2)*
             (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
               rgrid[i]*Qr1(i,j)*
                (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                  (2*(1 + rgrid[i]*Qrr(i,j))*
                     ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],3)) \
+ (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  ((2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                       (1 + rgrid[i]*Qtt(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],5)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                         (L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                pow(rgrid[i],2)))/pow(rgrid[i],2)))/
      (8.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))) + 
  (h11(i,j)*((-2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
          pow(rgrid[i],2))/
        (pow(L,2)*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) \
+ (18*w*(complex(0,1)*pow(mu,2) + 2*(complex(0,-6) + w))*
          (-1 + rgrid[i])*pow(rgrid[i],6)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)*
          (1 + rgrid[i]*Qrr(i,j))) - 
       (complex(0,12)*w*(-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (pow(rgrid[i],8)*((2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j)*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],6)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((8*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                    exp(B1*mu*h(i,j)*rgrid[i])*
                    pow(1 + rgrid[i]*Qrr(i,j),2))/
                  (pow(rgrid[i],4)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) + 
                 ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    ((-4*pow(Qrr.diff(d01,i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (pow(L,2)*
                       pow(-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2))/
                       pow(rgrid[i],4) - 
                      (2*Q11.diff(d01,i,j)*L*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*pow(Qr1(i,j),3))/
                       rgrid[i] + 
                      pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],2)*
                       pow(Qr1(i,j),4))*(1 + rgrid[i]*Qtt(i,j)))/
                  (2.*pow(rgrid[i],2)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                       pow(rgrid[i],2) + 
                      (2*Qrr.diff(d02,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                       pow(rgrid[i],3)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)\
))/(pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
       2*(1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
        ((36*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
           (2 + 3*pow(A1,2)) + 
          (24*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) + 
          (Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*(-1 + rgrid[i])*
             pow(rgrid[i],5)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
           (2.*pow(L,2)*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
          (pow(Qr1.diff(d01,i,j),2)*(-1 + rgrid[i])*
             pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3)))/
           (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
          ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3)))/
           (L*(1 + rgrid[i]*Qrr(i,j))) - 
          (Qr1.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
           (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
          (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],5)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
               Q22.diff(d01,i,j)*Qr1(i,j)))/
           (2.*pow(L,2)*(1 + rgrid[i]*Q22(i,j))*
             (1 + rgrid[i]*Qrr(i,j))) - 
          (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
             (((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Qrr(i,j)))/2. + 
               ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                    3*(-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                  (1 + rgrid[i]*Qrr(i,j)))/2.))/
           (L*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
          (pow(rgrid[i],4)*(rgrid[i]*
                (2*a0.diff(d01,i,j)*L*mu*exp(B1*mu*h(i,j)*rgrid[i])*
                   (-(mu*a0(i,j)) - 
                     a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                   (-1 + rgrid[i]) + 
                  (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                     (-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2.)*Qr1(i,j) + 
               pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                exp(B1*mu*h(i,j)*rgrid[i])*pow(-1 + rgrid[i],2)*
                pow(rgrid[i],2)*pow(Qr1(i,j),2) + 
               (2*pow(L,2)*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) + 
               L*(L*exp(B1*mu*h(i,j)*rgrid[i])*
                   pow(-(mu*a0(i,j)) - 
                     a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                  (Qr1.diff(d01,i,j)*
                     (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2)))\
)/(pow(L,2)*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))) + 
       (complex(0,3)*w*(-1 + rgrid[i])*pow(rgrid[i],12)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((4*(1 + rgrid[i]*Q22(i,j))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 Q11.diff(d01,i,j)*Qr1(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],8)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((2*((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    (-(L*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
                       (1 + rgrid[i]*Qtt(i,j)))/
                       (2.*pow(rgrid[i],3)) + 
                      (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       pow(Qr1(i,j),3)*(1 + rgrid[i]*Qtt(i,j)))/2. \
+ (4*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3) + 
                      rgrid[i]*Qr1(i,j)*
                       ((2*Qtt.diff(d01,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],3) - 
                        (2*Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3))))/
                  pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2)*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                      (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q22(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],3))))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],10)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
          (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
               (-((Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                      (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                    pow(rgrid[i],6)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((Q11.diff(d01,i,j)*Q22.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                        (2*Q11.diff(d02,i,j)*
                        (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))))) + 
            rgrid[i]*Qr1(i,j)*
             ((L*((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                    (Q11.diff(d01,i,j)*
                       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],4) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  (-(L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q11.diff(d01,i,j)*
                       (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) - 
                        (3*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        (2*L*(-1 + rgrid[i])*
                        (-Q11.diff(d01,i,j) + 
                       Q11.diff(d11,i,j)*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                pow(rgrid[i],2)) + 
            L*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],10)) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  ((2*L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                       (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                       (1 + rgrid[i]*Qtt(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],8)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) - 
                        (L*
                       (6 - 
                       2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],2) - 
                         (L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                pow(rgrid[i],2))))/
        (2.*pow(L,2)*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))))/4. + 
  (h22.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((-8*L*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
          ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
            (2*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) - 
       (complex(0,48)*L*w*rgrid[i]*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
          ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
            (2*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-12 + pow(mu,2))*(-1 + rgrid[i])*
          (1 + rgrid[i] + pow(rgrid[i],2))*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (2*pow(rgrid[i],6)*((2*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j)*
               ((-4*L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) + 
                 3*Q11.diff(d01,i,j)*Qr1(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],7)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((Q22.diff(d01,i,j)*pow(Qr1(i,j),2)*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],3) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((-3*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*pow(Qr1(i,j),3)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       (2.*pow(rgrid[i],2)) + 
                      (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        pow(Qr1(i,j),4)*(1 + rgrid[i]*Qtt(i,j)))/2. \
- (4*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
                        (1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                      ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
                        ((8*Qr1.diff(d01,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (4*L*
                       (((-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                      (-2 + 
                      Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] + 
                      pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                       ((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                        pow(rgrid[i],3) - 
                        (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],3))))/pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2)*
               ((-2*L*(1 + rgrid[i]*Q22(i,j))*
                    (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],6) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                    ((1 + rgrid[i]*Q22(i,j))*
                       ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((4*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],2)))/
        ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/
   (16.*pow(L,2)*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h22(i,j)*((complex(0,-18)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          pow(rgrid[i],6)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
            (4*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        (pow(-12 + pow(mu,2),2)*
          pow(1 + rgrid[i] + pow(rgrid[i],2),2)*(1 + rgrid[i]*Qrr(i,j))\
) + (complex(0,12)*w*(-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
            (4*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
          (1 + rgrid[i]*Qrr(i,j))) - 
       ((-1 + rgrid[i])*pow(rgrid[i],10)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((4*(1 + rgrid[i]*Q22(i,j))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 Q11.diff(d01,i,j)*Qr1(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],8)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((2*(-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((2*Qtt.diff(d01,i,j)*Qr1(i,j)*
                        (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],2) - 
                      (3*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       (2.*pow(rgrid[i],3)) + 
                      (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        pow(Qr1(i,j),3)*(1 + rgrid[i]*Qtt(i,j)))/2. \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((4*Qr1.diff(d01,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*L*
                       (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)))/
                  pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2)*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                      (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j))) + 
       (complex(0,24)*w*pow(rgrid[i],2)*
          (2 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
               pow(Qr1(i,j),2))/(2.*(1 + rgrid[i]*Qrr(i,j))) - 
            (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               ((4*(1 + rgrid[i]*Q22(i,j))*
                    (-((L*(-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                      Q11.diff(d01,i,j)*Qr1(i,j))*
                    pow(1 + rgrid[i]*Qrr(i,j),2)*
                    (1 + rgrid[i]*Qtt(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],8)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    ((2*(-((L*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                       Q22.diff(d01,i,j)*Qr1(i,j))*
                        (1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                      ((1 + rgrid[i]*Q22(i,j))*
                        ((2*Qtt.diff(d01,i,j)*Qr1(i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],2) - 
                        (3*L*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
                       (1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       pow(Qr1(i,j),3)*(1 + rgrid[i]*Qtt(i,j)))/2. \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((4*Qr1.diff(d01,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (4*L*
                       (((-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                      (-2 + 
                      Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)))/
                       pow(rgrid[i],2)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],4)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (pow(1 + rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2)*
                    (-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                      rgrid[i]*Qr1(i,j)*
                       (-((Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Q22(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q22(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],3))))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                      ((1 + rgrid[i]*Q22(i,j))*
                        ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                       pow(rgrid[i],2)))/pow(rgrid[i],2)))/
             (8.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (1 - rgrid[i])*((-2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
             pow(rgrid[i],2))/
           (pow(L,2)*(-1 + rgrid[i])*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             pow(1 + rgrid[i]*Q11(i,j),2)) + 
          (pow(rgrid[i],8)*((2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j\
)*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                ((-1 + rgrid[i])*pow(rgrid[i],6)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  ((8*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(1 + rgrid[i]*Qrr(i,j),2))/
                     (pow(rgrid[i],4)*
                       pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                    ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       ((-4*pow(Qrr.diff(d01,i,j),2))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) + 
                        (pow(L,2)*
                       pow(-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j),2)*pow(Qr1(i,j),2))/
                        pow(rgrid[i],4) - 
                        (2*Q11.diff(d01,i,j)*L*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*pow(Qr1(i,j),3))/rgrid[i] \
+ pow(Q11.diff(d01,i,j),2)*pow(rgrid[i],2)*pow(Qr1(i,j),4))*
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                         (2*Qrr.diff(d02,i,j)*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                pow(rgrid[i],2)))/
           (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
             (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
          2*(1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2)*
           ((36*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
              (2 + 3*pow(A1,2)) + 
             (24*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) + 
             (Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*(-1 + rgrid[i])*
                pow(rgrid[i],5)*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
              (2.*pow(L,2)*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
             (pow(Qr1.diff(d01,i,j),2)*(-1 + rgrid[i])*
                pow(rgrid[i],4)*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3)))/
              (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
                (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3)))/
              (L*(1 + rgrid[i]*Qrr(i,j))) - 
             (Qr1.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
              (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
             (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],5)*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*
                ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                  Q22.diff(d01,i,j)*Qr1(i,j)))/
              (2.*pow(L,2)*(1 + rgrid[i]*Q22(i,j))*
                (1 + rgrid[i]*Qrr(i,j))) - 
             (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
                (((-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                     (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                  ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                     (1 + rgrid[i]*Qrr(i,j)))/2.))/
              (L*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
             (pow(rgrid[i],4)*
                (rgrid[i]*(2*a0.diff(d01,i,j)*L*mu*
                      exp(B1*mu*h(i,j)*rgrid[i])*
                      (-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                      (-1 + rgrid[i]) + 
                     (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2.)*Qr1(i,j) + 
                  pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                   exp(B1*mu*h(i,j)*rgrid[i])*pow(-1 + rgrid[i],2)*
                   pow(rgrid[i],2)*pow(Qr1(i,j),2) + 
                  (2*pow(L,2)*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                   ((-1 + rgrid[i])*pow(rgrid[i],2)*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                  L*(L*exp(B1*mu*h(i,j)*rgrid[i])*
                      pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                     (Qr1.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                         (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2)))\
)/(pow(L,2)*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))) + 
          ((-1 + rgrid[i])*pow(rgrid[i],10)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
             (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                         (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))\
)/pow(rgrid[i],6)) + (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                         ((Q11.diff(d01,i,j)*Q22.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                        (2*Q11.diff(d02,i,j)*
                        (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                         (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))) + 
               rgrid[i]*Qr1(i,j)*
                ((L*((Q22.diff(d01,i,j)*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                       (Q11.diff(d01,i,j)*
                        (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                   pow(rgrid[i],4) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     (-(L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q11.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                         ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) - 
                         (3*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                         (2*L*(-1 + rgrid[i])*
                        (-Q11.diff(d01,i,j) + 
                        Q11.diff(d11,i,j)*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                         (Q11.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2)) + 
               L*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                         rgrid[i]*Q11(i,j))*
                       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                         rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],10)) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     ((2*L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                         (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],8)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                       (2*(1 + rgrid[i]*Qrr(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) - 
                        (L*(6 - 
                       2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                         pow(rgrid[i],2) - 
                         (L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2))))/
           (2.*pow(L,2)*(1 + rgrid[i]*Q22(i,j))*
             pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j))))))/4.
;
elem[4]=
htt.diff(d11,i,j)/2. + complex(0,0.5)*ht1.diff(d10,i,j)*L*w - 
  Dh.diff(d10,i,j)*h.diff(d01,i,j)*mu*pow(rgrid[i],2) - 
  Dh.diff(d01,i,j)*mu*rgrid[i]*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) + 
  complex(0,0.5)*ht1.diff(d01,i,j)*w*rgrid[i]*Qr1(i,j) - 
  (h11.diff(d20,i,j)*L*(-1 + rgrid[i])*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (4.*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d20,i,j)*L*pow(-1 + rgrid[i],2)*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (4.*(1 + rgrid[i]*Qrr(i,j))) - 
  (h11.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(Qr1(i,j),3))/(4.*L*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d02,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     pow(Qr1(i,j),3))/(4.*L*(1 + rgrid[i]*Qrr(i,j))) + 
  h22.diff(d11,i,j)*(1 - rgrid[i])*
   (-0.5 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2))/(2.*(1 + rgrid[i]*Qrr(i,j)))) + 
  h11.diff(d11,i,j)*(0.5 + ((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2))/(2.*(1 + rgrid[i]*Qrr(i,j)))) + 
  (complex(0,2)*L*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*
     (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
     pow(rgrid[i],2))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,2)*L*w*Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
     pow(rgrid[i],2))/
   ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(1 - rgrid[i])*
     pow(rgrid[i],4)*(((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   ((1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) - 
  (A0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],4)*
     (-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],4)*(-(L*(-(mu*a0(i,j)) - 
            a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (L*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (htt.diff(d10,i,j)*pow(rgrid[i],4)*
     ((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],3) - 
       (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3) + 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i])*
          Qr1(i,j)*(1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) - 
       (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i])))/
   (4.*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (htt.diff(d01,i,j)*((complex(0,12)*w*pow(rgrid[i],2))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
       (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - rgrid[i]*Q11(i,j))/
        (rgrid[i]*(1 + rgrid[i]*Q11(i,j))) + 
       (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
        (L*(1 + rgrid[i]*Qrr(i,j))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i])*
          pow(Qr1(i,j),2))/(2.*L*(1 + rgrid[i]*Qrr(i,j))) + 
       (Q11.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),3))/
        (2.*L*(1 + rgrid[i]*Qrr(i,j))) + 
       (2*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Qtt(i,j)))/2. - 
            ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                 3*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))*
               (1 + rgrid[i]*Qtt(i,j)))/2.))/
        ((-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))/4. - 
  (htt(i,j)*pow(rgrid[i],4)*((exp(B1*mu*h(i,j)*rgrid[i])*
          (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
               pow(-(L*(-(mu*a0(i,j)) - 
                      a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                 a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                  Qr1(i,j),2))/rgrid[i] - 
            (2*a0.diff(d01,i,j)*mu*
               (-2*L*(-(mu*a0(i,j)) - 
                    a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) - 
                 a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                  Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
             (pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/L + 
       (complex(0,3)*w*pow(rgrid[i],2)*
          (-((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
               pow(rgrid[i],3)) + 
            (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],3) + 
            ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/
                  rgrid[i])*Qr1(i,j)*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*rgrid[i]) + (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
               (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i])))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3)))))/
   (2.*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  Dh(i,j)*(-(h.diff(d01,i,j)*mu*rgrid[i]) - 
     (complex(0,6)*h.diff(d01,i,j)*mu*w*pow(rgrid[i],4))/
      ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
     (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
        ((2*a0.diff(d01,i,j)*(2 + 3*pow(A1,2))*B1*mu*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + B1*mu*h(i,j)*rgrid[i])*
             (-2*L*(-(mu*a0(i,j)) - 
                  a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) - 
               a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))*
             (1 + rgrid[i]*Qrr(i,j)))/
           (pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) - 
          ((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
             ((2 + 3*pow(A1,2))*B1*
                exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + 
                  B1*mu*h(i,j)*rgrid[i])*
                pow(-(L*(-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                  a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                   Qr1(i,j),2) - 
               (24*A1*pow(L,2)*
                  (-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)*rgrid[i]))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],4)))/rgrid[i]))/
      (2.*(2 + 3*pow(A1,2))*L*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j)))) + 
  (complex(0,1)*w*ht1(i,j)*pow(rgrid[i],4)*
     (-(L*(-12 + pow(mu,2))*(-1 + rgrid[i])*(-1 + pow(rgrid[i],3))*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*
           (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
             rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        (2.*pow(rgrid[i],5)) + 
       ((1 + rgrid[i]*Q11(i,j))*
          (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (complex(0,6)*L*w*pow(rgrid[i],2) + 
                 Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*rgrid[i]*
                  (-1 + pow(rgrid[i],3)))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],2)) + 
            (-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
             ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. + 
               (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3))))/
        pow(rgrid[i],2)))/
   ((-12 + pow(mu,2))*(-1 + rgrid[i])*(-1 + pow(rgrid[i],3))*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (h11.diff(d10,i,j)*((2*Q22.diff(d01,i,j)*rgrid[i])/
        (1 + rgrid[i]*Q22(i,j)) - 
       (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) - 
       (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
            rgrid[i]*Q11(i,j))*Qr1(i,j))/(2.*(1 + rgrid[i]*Qrr(i,j))) + 
       (Q11.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2))/
        (2.*(1 + rgrid[i]*Qrr(i,j))) + 
       (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j)) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
          Qr1(i,j)*((-2*Qrr.diff(d01,i,j)*Qr1(i,j))/
             ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (4*L*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qrr(i,j)))/2.))/
             (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               (2*Qr1.diff(d01,i,j)*rgrid[i] - 
                 (complex(0,24)*L*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                 (pow(rgrid[i],2)*
                    (-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j)))/
                  (1 + rgrid[i]*Q22(i,j)) + 
                 (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
                  (1 + rgrid[i]*Qtt(i,j)) - 
                 (2*L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        (4.*pow(1 + rgrid[i]*Qrr(i,j),2))))/4. + 
  (h11.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((-4*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
            rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2))/
        (pow(-1 + rgrid[i],2)*pow(rgrid[i],5)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))) + 
       (2*(1 + rgrid[i]*Qrr(i,j))*
          ((rgrid[i]*Qr1(i,j)*
               ((4*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j))/pow(rgrid[i],2) - 
                 Q11.diff(d01,i,j)*rgrid[i]*pow(Qr1(i,j),2)))/L + 
            (2*((complex(0,12)*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                 (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q22(i,j))/
                  (rgrid[i]*(1 + rgrid[i]*Q22(i,j))))*
               (1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       ((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          ((4*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
               (1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
               ((2*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Q22(i,j))) - 
                 (2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))))/L + 
            rgrid[i]*Qr1(i,j)*
             ((-4*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                  pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3),2)) + 
               (2*(1 + rgrid[i]*Qrr(i,j))*
                  ((-4*Qr1.diff(d01,i,j)*rgrid[i])/L + 
                    (complex(0,24)*w*pow(rgrid[i],2))/
                     ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j))/
                     (rgrid[i]*(1 + rgrid[i]*Q22(i,j))) + 
                    (2*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                     ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))))/
                ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))))))/rgrid[i]))/
   (16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h22.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((8*L*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (complex(0,48)*L*w*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-12 + pow(mu,2))*(-1 + rgrid[i])*rgrid[i]*
          (1 + rgrid[i] + pow(rgrid[i],2))*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (2*(1 - rgrid[i])*pow(rgrid[i],4)*
          (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((2*(1 + rgrid[i]*Qrr(i,j))*
                    ((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       pow(rgrid[i],3) + 
                      (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                       pow(rgrid[i],3)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  ((-3*L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) + 
                    ((1 + rgrid[i]*Q11(i,j))*
                       ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],2)))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))))/(16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h22.diff(d01,i,j)*(2 - (2*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          pow(Qr1(i,j),2))/(1 + rgrid[i]*Qrr(i,j)) - 
       (complex(0,6)*w*(-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (2*(1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) - 
            (2*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((-2*L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],11)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            ((1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (3*(1 + rgrid[i]*Q22(i,j))*pow(Qr1(i,j),2)*
                  (-((L*(-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                    Q11.diff(d01,i,j)*Qr1(i,j)) + 
                 (2*L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],5)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*Qr1(i,j)*
               ((-2*L*(1 + rgrid[i]*Q22(i,j))*
                    (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],6) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                    ((1 + rgrid[i]*Q22(i,j))*
                       ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((4*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],3)))/
        (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))))/4. + 
  (h11(i,j)*((18*L*w*(complex(0,1)*pow(mu,2) + 2*(complex(0,-6) + w))*
          (-1 + rgrid[i])*pow(rgrid[i],5)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          Qr1(i,j))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)*
          (1 + rgrid[i]*Qrr(i,j))) - 
       (complex(0,12)*L*w*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          Qr1(i,j))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (complex(0,3)*w*(-1 + rgrid[i])*pow(rgrid[i],10)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Qrr(i,j))*
               (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
                    (-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j)))/rgrid[i] + 
                 (4*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((2*(1 + rgrid[i]*Qrr(i,j))*
                    ((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       pow(rgrid[i],3) - 
                      (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                       pow(rgrid[i],3)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],3))))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                    ((1 + rgrid[i]*Q11(i,j))*
                       ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],2)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j))) + 
       rgrid[i]*Qr1(i,j)*(((-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             pow(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
               Q11.diff(d01,i,j)*Qr1(i,j),2))/
           (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
          (2*(1 + rgrid[i]*Q11(i,j))*
             ((36*pow(A1,2)*L*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
                (2 + 3*pow(A1,2)) + 
               (24*L*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) + 
               (Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                  (-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
                (2.*L*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
               (pow(Qr1.diff(d01,i,j),2)*(-1 + rgrid[i])*
                  pow(rgrid[i],4)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))/
                (L*(1 + rgrid[i]*Qrr(i,j))) + 
               ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))/
                (1 + rgrid[i]*Qrr(i,j)) - 
               (Qr1.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
                (L*(1 + rgrid[i]*Qrr(i,j))) + 
               (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                    Q22.diff(d01,i,j)*Qr1(i,j)))/
                (2.*L*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
               (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
                  (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                pow(1 + rgrid[i]*Qrr(i,j),2) - 
               (pow(rgrid[i],4)*
                  (rgrid[i]*(2*a0.diff(d01,i,j)*L*mu*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                        (-1 + rgrid[i]) + 
                       (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2.)*Qr1(i,j) + 
                    pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                     exp(B1*mu*h(i,j)*rgrid[i])*pow(-1 + rgrid[i],2)*
                     pow(rgrid[i],2)*pow(Qr1(i,j),2) + 
                    (2*pow(L,2)*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                    L*(L*exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                       (Qr1.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2))\
))/(L*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/pow(rgrid[i],2) \
+ ((-1 + rgrid[i])*pow(rgrid[i],8)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))\
)/pow(rgrid[i],6)) + (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((Q11.diff(d01,i,j)*Q22.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                        (2*Q11.diff(d02,i,j)*
                        (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))) + 
               rgrid[i]*Qr1(i,j)*
                ((L*((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                       (Q11.diff(d01,i,j)*
                        (-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                   pow(rgrid[i],4) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     (-(L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q11.diff(d01,i,j)*
                       (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) - 
                        (3*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        (2*L*(-1 + rgrid[i])*
                        (-Q11.diff(d01,i,j) + 
                       Q11.diff(d11,i,j)*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2)) + 
               L*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],10)) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     ((2*L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],8)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                       (2*(1 + rgrid[i]*Qrr(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) - 
                        (L*
                       (6 - 
                       2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],2) - 
                         (L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2))))/
           (2.*L*(1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
             (1 + rgrid[i]*Qtt(i,j))))))/4. + 
  (h22(i,j)*((complex(0,-18)*L*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          pow(rgrid[i],5)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          Qr1(i,j))/
        (pow(-12 + pow(mu,2),2)*
          pow(1 + rgrid[i] + pow(rgrid[i],2),2)*(1 + rgrid[i]*Qrr(i,j))\
) + (complex(0,12)*L*w*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          Qr1(i,j))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
          (1 + rgrid[i]*Qrr(i,j))) - 
       ((-1 + rgrid[i])*pow(rgrid[i],8)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((2*(1 + rgrid[i]*Qrr(i,j))*
                    ((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       pow(rgrid[i],3) + 
                      (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                       pow(rgrid[i],3)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  ((-3*L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) + 
                    ((1 + rgrid[i]*Q11(i,j))*
                       ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],2)))/
        (2.*(1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j))) + 
       (1 - rgrid[i])*rgrid[i]*Qr1(i,j)*
        (((-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             pow(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
               Q11.diff(d01,i,j)*Qr1(i,j),2))/
           (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
          (2*(1 + rgrid[i]*Q11(i,j))*
             ((36*pow(A1,2)*L*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
                (2 + 3*pow(A1,2)) + 
               (24*L*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) + 
               (Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                  (-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
                (2.*L*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
               (pow(Qr1.diff(d01,i,j),2)*(-1 + rgrid[i])*
                  pow(rgrid[i],4)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))/
                (L*(1 + rgrid[i]*Qrr(i,j))) + 
               ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))/
                (1 + rgrid[i]*Qrr(i,j)) - 
               (Qr1.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
                (L*(1 + rgrid[i]*Qrr(i,j))) + 
               (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                    Q22.diff(d01,i,j)*Qr1(i,j)))/
                (2.*L*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
               (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
                  (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                pow(1 + rgrid[i]*Qrr(i,j),2) - 
               (pow(rgrid[i],4)*
                  (rgrid[i]*(2*a0.diff(d01,i,j)*L*mu*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                        (-1 + rgrid[i]) + 
                       (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2.)*Qr1(i,j) + 
                    pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                     exp(B1*mu*h(i,j)*rgrid[i])*pow(-1 + rgrid[i],2)*
                     pow(rgrid[i],2)*pow(Qr1(i,j),2) + 
                    (2*pow(L,2)*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                    L*(L*exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                       (Qr1.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2))\
))/(L*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/pow(rgrid[i],2) \
+ ((-1 + rgrid[i])*pow(rgrid[i],8)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))\
)/pow(rgrid[i],6)) + (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((Q11.diff(d01,i,j)*Q22.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                        (2*Q11.diff(d02,i,j)*
                        (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))) + 
               rgrid[i]*Qr1(i,j)*
                ((L*((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                       (Q11.diff(d01,i,j)*
                        (-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                   pow(rgrid[i],4) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     (-(L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q11.diff(d01,i,j)*
                       (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) - 
                        (3*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        (2*L*(-1 + rgrid[i])*
                        (-Q11.diff(d01,i,j) + 
                       Q11.diff(d11,i,j)*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2)) + 
               L*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],10)) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     ((2*L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],8)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                       (2*(1 + rgrid[i]*Qrr(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) - 
                        (L*
                       (6 - 
                       2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],2) - 
                         (L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2))))/
           (2.*L*(1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
             (1 + rgrid[i]*Qtt(i,j)))) + 
       (complex(0,1.5)*w*pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          ((8*L*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],3)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 - rgrid[i])*pow(rgrid[i],4)*
               (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
                    (-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],5) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],3) + 
                        (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],3)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                      pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                       (-((Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                      rgrid[i]*Qr1(i,j)*
                       ((-3*L*
                        (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) + 
                         ((1 + rgrid[i]*Q11(i,j))*
                         ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                         pow(rgrid[i],2))))/pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          pow(1 + rgrid[i]*Qrr(i,j),2))))/4.
;
elem[5]=
htt.diff(d02,i,j)/2. + complex(0,1)*ht1.diff(d01,i,j)*L*w - 
  2*Dh.diff(d01,i,j)*h.diff(d01,i,j)*mu*pow(rgrid[i],2) - 
  (h11.diff(d20,i,j)*pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j)))/
   (4.*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d20,i,j)*pow(L,2)*pow(-1 + rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j)))/
   (4.*(1 + rgrid[i]*Qrr(i,j))) + 
  (h11.diff(d11,i,j)*L*(-1 + rgrid[i])*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (2.*(1 + rgrid[i]*Qrr(i,j))) - 
  (h22.diff(d11,i,j)*L*pow(-1 + rgrid[i],2)*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
   (2.*(1 + rgrid[i]*Qrr(i,j))) + 
  (htt.diff(d10,i,j)*L*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
            rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
       (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i] - 
       Q11.diff(d01,i,j)*Qr1(i,j)))/(8.*(1 + rgrid[i]*Qrr(i,j))) + 
  h11.diff(d02,i,j)*(0.5 - ((-1 + rgrid[i])*pow(rgrid[i],2)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        pow(Qr1(i,j),2))/(4.*(1 + rgrid[i]*Qrr(i,j)))) + 
  (h22.diff(d02,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/(4.*(1 + rgrid[i]*Qrr(i,j))) - 
  (complex(0,2)*a0.diff(d01,i,j)*L*mu*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,2)*a0.diff(d01,i,j)*L*mu*w*Ax(i,j)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2))/
   ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d10,i,j)*L*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q11(i,j))*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   ((1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (L*A0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
     (1 + rgrid[i]*Q11(i,j))*(L*
        (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(1 - rgrid[i])*
     pow(rgrid[i],4)*(((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   ((1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) - 
  complex(0,0.5)*L*w*ht1(i,j)*((Q11.diff(d01,i,j)*rgrid[i])/
      (1 + rgrid[i]*Q11(i,j)) - 
     (2*Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))) + 
  (htt.diff(d01,i,j)*(-((Q11.diff(d01,i,j)*rgrid[i])/
          (1 + rgrid[i]*Q11(i,j))) - 
       (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
            rgrid[i]*Q11(i,j))*Qr1(i,j))/(2.*(1 + rgrid[i]*Qrr(i,j))) + 
       (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          Qr1(i,j))/(1 + rgrid[i]*Qrr(i,j)) + 
       (Q11.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2))/
        (2.*(1 + rgrid[i]*Qrr(i,j))) + 
       (2*Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))))/4. + 
  (htt(i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((complex(0,3)*L*w*pow(rgrid[i],2)*
          ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i] - 
            Q11.diff(d01,i,j)*Qr1(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (2*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
          (-(((1 + rgrid[i]*Q11(i,j))*
                 pow(-(L*(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                   a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                    Qr1(i,j),2))/pow(rgrid[i],2)) + 
            (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*(-1 + rgrid[i])*
               (1 + rgrid[i]*Qrr(i,j)))/
             (pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))/
   (4.*(1 + rgrid[i]*Qrr(i,j))) + 
  (Dh(i,j)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
     ((2*pow(a0.diff(d01,i,j),2)*(2 + 3*pow(A1,2))*B1*pow(mu,2)*
          exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + B1*mu*h(i,j)*rgrid[i])*
          (-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) - 
       ((1 + rgrid[i]*Q11(i,j))*
          ((2 + 3*pow(A1,2))*B1*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + B1*mu*h(i,j)*rgrid[i])*
             pow(-(L*(-(mu*a0(i,j)) - 
                    a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
               a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j),2) \
- (24*A1*pow(L,2)*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)*rgrid[i]))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],4)))/pow(rgrid[i],2)))/
   (2.*(2 + 3*pow(A1,2))*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) \
+ (h11.diff(d10,i,j)*L*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((2*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((-2*Qrr.diff(d01,i,j)*Qr1(i,j))/
             ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (4*L*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qrr(i,j)))/2.))/
             (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               (2*Qr1.diff(d01,i,j)*rgrid[i] - 
                 (complex(0,24)*L*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                 (pow(rgrid[i],2)*
                    (-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j)))/
                  (1 + rgrid[i]*Q22(i,j)) + 
                 (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
                  (1 + rgrid[i]*Qtt(i,j)) - 
                 (2*L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)))/
   (16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h11.diff(d01,i,j)*(-((Q11.diff(d01,i,j)*rgrid[i])/
          (1 + rgrid[i]*Q11(i,j))) + 
       (3*Q22.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) + 
       (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) + 
       (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
            rgrid[i]*Q11(i,j))*Qr1(i,j))/(2.*(1 + rgrid[i]*Qrr(i,j))) - 
       (Q11.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2))/
        (2.*(1 + rgrid[i]*Qrr(i,j))) + 
       (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j)) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
          ((4*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
               (1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            pow(rgrid[i],2)*pow(Qr1(i,j),2)*
             ((2*Qrr.diff(d01,i,j))/
                ((-1 + rgrid[i])*rgrid[i]*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) - 
               (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                ((-1 + rgrid[i])*rgrid[i]*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Q22(i,j))) - 
               (2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                ((-1 + rgrid[i])*rgrid[i]*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Qtt(i,j)))) + 
            rgrid[i]*Qr1(i,j)*
             ((-4*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                  pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3),2)) + 
               (2*(1 + rgrid[i]*Qrr(i,j))*
                  (-4*Qr1.diff(d01,i,j)*rgrid[i] + 
                    (complex(0,24)*L*w*pow(rgrid[i],2))/
                     ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                    (L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/
                     (rgrid[i]*(1 + rgrid[i]*Q22(i,j))) + 
                    (2*L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                     ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))))/
                ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))))))/
        (4.*pow(1 + rgrid[i]*Qrr(i,j),2))))/4. + 
  (h22.diff(d10,i,j)*L*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((8*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (complex(0,48)*L*w*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-12 + pow(mu,2))*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (1 + rgrid[i] + pow(rgrid[i],2))*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (2*(1 - rgrid[i])*pow(rgrid[i],4)*
          (((1 + rgrid[i]*Q11(i,j))*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((-3*L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    ((2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))))/(16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h11(i,j)*((18*pow(L,2)*w*
          (complex(0,1)*pow(mu,2) + 2*(complex(0,-6) + w))*
          (-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j)))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)*
          (1 + rgrid[i]*Qrr(i,j))) - 
       (complex(0,12)*pow(L,2)*w*(-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          pow(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            Q11.diff(d01,i,j)*Qr1(i,j),2))/
        (2.*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
       (2*(1 + rgrid[i]*Q11(i,j))*
          ((36*pow(A1,2)*pow(L,2)*
               exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/(2 + 3*pow(A1,2)) \
+ (24*pow(L,2)*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) + 
            (Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*(-1 + rgrid[i])*
               pow(rgrid[i],5)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
             (2.*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
            (pow(Qr1.diff(d01,i,j),2)*(-1 + rgrid[i])*
               pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))/(1 + rgrid[i]*Qrr(i,j)) \
+ (L*(-1 + rgrid[i])*pow(rgrid[i],2)*
               (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))/(1 + rgrid[i]*Qrr(i,j)) \
- (Qr1.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
             (1 + rgrid[i]*Qrr(i,j)) + 
            (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],5)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j)))/
             (2.*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
            (Qr1.diff(d01,i,j)*L*pow(rgrid[i],2)*
               (((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qrr(i,j)))/2.))/
             pow(1 + rgrid[i]*Qrr(i,j),2) - 
            (pow(rgrid[i],4)*
               (rgrid[i]*(2*a0.diff(d01,i,j)*L*mu*
                     exp(B1*mu*h(i,j)*rgrid[i])*
                     (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                     (-1 + rgrid[i]) + 
                    (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2.)*Qr1(i,j) + 
                 pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                  exp(B1*mu*h(i,j)*rgrid[i])*pow(-1 + rgrid[i],2)*
                  pow(rgrid[i],2)*pow(Qr1(i,j),2) + 
                 (2*pow(L,2)*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 L*(L*exp(B1*mu*h(i,j)*rgrid[i])*
                     pow(-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                    (Qr1.diff(d01,i,j)*
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2))\
))/((1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/pow(rgrid[i],2) - 
       (complex(0,3)*L*w*(-1 + rgrid[i])*pow(rgrid[i],10)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Q11(i,j))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) - 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],3))))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    ((-2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (-((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]) + 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],8)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
               (-((Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                      (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                    pow(rgrid[i],6)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((Q11.diff(d01,i,j)*Q22.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                        (2*Q11.diff(d02,i,j)*
                        (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))))) + 
            rgrid[i]*Qr1(i,j)*
             ((L*((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                    (Q11.diff(d01,i,j)*
                       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],4) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  (-(L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q11.diff(d01,i,j)*
                       (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) - 
                        (3*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        (2*L*(-1 + rgrid[i])*
                        (-Q11.diff(d01,i,j) + 
                       Q11.diff(d11,i,j)*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                pow(rgrid[i],2)) + 
            L*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],10)) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  ((2*L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                       (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                       (1 + rgrid[i]*Qtt(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],8)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) - 
                        (L*
                       (6 - 
                       2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],2) - 
                         (L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                pow(rgrid[i],2))))/
        (2.*(1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j)))))/4. + 
  h22.diff(d01,i,j)*(-(L*(-1 + rgrid[i])*rgrid[i]*
         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
           pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)\
)/(2.*(1 + rgrid[i]*Qrr(i,j))) - 
     (complex(0,3)*L*w*(-1 + rgrid[i])*pow(rgrid[i],3)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j))/
      ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
        (1 + rgrid[i]*Qrr(i,j))) + 
     ((1 - rgrid[i])*(-1 + rgrid[i])*pow(rgrid[i],10)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        ((2*Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
             pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
           ((-1 + rgrid[i])*pow(rgrid[i],9)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
             (-((Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  ((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                     pow(rgrid[i],3) + 
                    (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],3) + 
                    (3*L*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*Qr1(i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],4)) \
- (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
                     (2.*rgrid[i])))/pow(rgrid[i],2)))/
           ((-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (pow(1 + rgrid[i]*Q11(i,j),2)*
             ((2*L*(1 + rgrid[i]*Q22(i,j))*
                  (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],6) - 
               pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                  (2*(1 + rgrid[i]*Qrr(i,j))*
                     ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],3)) \
+ (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))) + 
               rgrid[i]*Qr1(i,j)*
                ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                     (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     ((-2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                       (2*(1 + rgrid[i]*Qrr(i,j))*
                         ((-4*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] + 
                         (L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2))))/pow(rgrid[i],4)))/
      (8.*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))) + 
  (h22(i,j)*((complex(0,-18)*pow(L,2)*
          (-12 + pow(mu,2) - complex(0,2)*w)*w*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j)))/
        (pow(-12 + pow(mu,2),2)*
          pow(1 + rgrid[i] + pow(rgrid[i],2),2)*(1 + rgrid[i]*Qrr(i,j))\
) + (complex(0,12)*pow(L,2)*w*(-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j)))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (L*(-1 + rgrid[i])*pow(rgrid[i],8)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Q11(i,j))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((3*L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) + 
                 rgrid[i]*Qr1(i,j)*
                  ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) - 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    ((-2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((-3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] + 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (2.*(1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j))) + 
       (complex(0,1.5)*L*w*pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          ((8*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 - rgrid[i])*pow(rgrid[i],4)*
               (((1 + rgrid[i]*Q11(i,j))*
                    (-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],6) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((-3*L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) + 
                      rgrid[i]*Qr1(i,j)*
                       (-((Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (3*Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],3))))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                      ((1 + rgrid[i]*Q11(i,j))*
                        ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((3*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                       pow(rgrid[i],2)))/pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          pow(1 + rgrid[i]*Qrr(i,j),2)) + 
       (1 - rgrid[i])*(((-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             pow(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
               Q11.diff(d01,i,j)*Qr1(i,j),2))/
           (2.*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
          (2*(1 + rgrid[i]*Q11(i,j))*
             ((36*pow(A1,2)*pow(L,2)*
                  exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
                (2 + 3*pow(A1,2)) + 
               (24*pow(L,2)*exp(A1*mu*h(i,j)*rgrid[i]))/
                (2 + 3*pow(A1,2)) + 
               (Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*(-1 + rgrid[i])*
                  pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
                (2.*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
               (pow(Qr1.diff(d01,i,j),2)*(-1 + rgrid[i])*
                  pow(rgrid[i],4)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))/
                (1 + rgrid[i]*Qrr(i,j)) + 
               (L*(-1 + rgrid[i])*pow(rgrid[i],2)*
                  (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))/
                (1 + rgrid[i]*Qrr(i,j)) - 
               (Qr1.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
                (1 + rgrid[i]*Qrr(i,j)) + 
               (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                    Q22.diff(d01,i,j)*Qr1(i,j)))/
                (2.*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
               (Qr1.diff(d01,i,j)*L*pow(rgrid[i],2)*
                  (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                pow(1 + rgrid[i]*Qrr(i,j),2) - 
               (pow(rgrid[i],4)*
                  (rgrid[i]*(2*a0.diff(d01,i,j)*L*mu*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        (-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                        (-1 + rgrid[i]) + 
                       (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2.)*Qr1(i,j) + 
                    pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                     exp(B1*mu*h(i,j)*rgrid[i])*pow(-1 + rgrid[i],2)*
                     pow(rgrid[i],2)*pow(Qr1(i,j),2) + 
                    (2*pow(L,2)*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                    L*(L*exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                       (Qr1.diff(d01,i,j)*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2)))\
)/((1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/pow(rgrid[i],2) + 
          ((-1 + rgrid[i])*pow(rgrid[i],8)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                         (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))\
)/pow(rgrid[i],6)) + (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                         ((Q11.diff(d01,i,j)*Q22.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                        (2*Q11.diff(d02,i,j)*
                        (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                         (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))) + 
               rgrid[i]*Qr1(i,j)*
                ((L*((Q22.diff(d01,i,j)*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                       (Q11.diff(d01,i,j)*
                        (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                   pow(rgrid[i],4) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     (-(L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q11.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                         ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) - 
                         (3*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                         (2*L*(-1 + rgrid[i])*
                        (-Q11.diff(d01,i,j) + 
                        Q11.diff(d11,i,j)*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                         (Q11.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2)) + 
               L*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                         rgrid[i]*Q11(i,j))*
                       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                         rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],10)) + 
                  ((1 + rgrid[i]*Q22(i,j))*
                     ((2*L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                         (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],8)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                       (2*(1 + rgrid[i]*Qrr(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) - 
                        (L*(6 - 
                       2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                        (Q11.diff(d01,i,j)*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                         pow(rgrid[i],2) - 
                         (L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2))))/
           (2.*(1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
             (1 + rgrid[i]*Qtt(i,j))))))/4.
;
elem[6]=
(complex(0,0.5)*Q22.diff(d01,i,j)*w*ht1(i,j)*rgrid[i])/
   (L*(1 + rgrid[i]*Q11(i,j))) + 
  (h11.diff(d20,i,j)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j)))/
   (4.*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d20,i,j)*pow(-1 + rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j)))/
   (4.*(1 + rgrid[i]*Qrr(i,j))) - 
  (h11.diff(d11,i,j)*(-1 + rgrid[i])*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   (2.*L*(1 + rgrid[i]*Qrr(i,j))) - 
  (h22.diff(d11,i,j)*pow(-1 + rgrid[i],2)*rgrid[i]*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j))/
   (2.*L*(1 + rgrid[i]*Qrr(i,j))) + 
  (htt.diff(d10,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
            rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
       Q22.diff(d01,i,j)*Qr1(i,j)))/(8.*L*(1 + rgrid[i]*Qrr(i,j))) + 
  (htt.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
            Q22.diff(d01,i,j)*Qr1(i,j)))/rgrid[i] + 
       (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (8.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
  (h11.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (4.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d02,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (4.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
  (complex(0,2)*a0.diff(d01,i,j)*mu*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,2)*a0.diff(d01,i,j)*mu*w*Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Q22(i,j))*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (A0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
     (1 + rgrid[i]*Q22(i,j))*(L*
        (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],4)*(1 + rgrid[i]*Q22(i,j))*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (htt(i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((complex(0,3)*L*w*pow(rgrid[i],2)*
          ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
            Q22.diff(d01,i,j)*Qr1(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
       (2*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
          (1 + rgrid[i]*Q22(i,j))*
          (((1 + rgrid[i]*Q11(i,j))*
               pow(-(L*(-(mu*a0(i,j)) - 
                      a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                 a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                  Qr1(i,j),2))/pow(rgrid[i],2) + 
            (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*(-1 + rgrid[i])*
               (1 + rgrid[i]*Qrr(i,j)))/
             (pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))))/(4.*pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) \
- (Dh(i,j)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
     (1 + rgrid[i]*Q22(i,j))*((2*pow(a0.diff(d01,i,j),2)*
          (2 + 3*pow(A1,2))*B1*pow(mu,2)*
          exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + B1*mu*h(i,j)*rgrid[i])*
          (-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((2 + 3*pow(A1,2))*B1*
             exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + B1*mu*h(i,j)*rgrid[i])*
             pow(-(L*(-(mu*a0(i,j)) - 
                    a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
               a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j),2) \
- (24*A1*pow(L,2)*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)*rgrid[i]))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],4)))/pow(rgrid[i],2)))/
   (2.*(2 + 3*pow(A1,2))*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (h11.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((2*((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
            Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
        (L*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       ((1 + rgrid[i]*Q22(i,j))*
          ((2*Qrr.diff(d01,i,j)*Qr1(i,j))/
             (L*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) - 
            (4*(((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qrr(i,j)))/2.))/
             (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               ((-2*Qr1.diff(d01,i,j)*rgrid[i])/L + 
                 (complex(0,24)*w*pow(rgrid[i],2))/
                  ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                 (pow(rgrid[i],2)*
                    ((-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))/pow(rgrid[i],3) - 
                      (Q11.diff(d01,i,j)*Qr1(i,j))/L))/
                  (1 + rgrid[i]*Q11(i,j)) - 
                 (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
                  (L*(1 + rgrid[i]*Qtt(i,j))) + 
                 (2*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)))/
   (16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h11.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((2*(1 + rgrid[i]*Qrr(i,j))*
          (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j)))/rgrid[i] + 
            (6*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],3)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
       ((1 + rgrid[i]*Q22(i,j))*
          (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
             ((-2*Qrr.diff(d01,i,j))/
                ((-1 + rgrid[i])*rgrid[i]*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) + 
               (2*Q11.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                ((-1 + rgrid[i])*rgrid[i]*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Q11(i,j))) + 
               (2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                ((-1 + rgrid[i])*rgrid[i]*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Qtt(i,j)))) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               (-2*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)) - 
                 (2*Q11.diff(d01,i,j)*rgrid[i]*
                    (1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    pow(1 + rgrid[i]*Q11(i,j),2)) + 
                 (pow(rgrid[i],2)*
                    ((2*Qrr.diff(d01,i,j))/
                       ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))))/
                  (1 + rgrid[i]*Q11(i,j))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            rgrid[i]*Qr1(i,j)*
             ((4*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                    ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                (pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
                  pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3),2)) + 
               (2*(1 + rgrid[i]*Qrr(i,j))*
                  (4*Qr1.diff(d01,i,j)*rgrid[i] - 
                    (complex(0,24)*L*w*pow(rgrid[i],2))/
                     ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
                    (L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/
                     (rgrid[i]*(1 + rgrid[i]*Q11(i,j))) - 
                    (2*L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                     ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))))/
                ((-1 + rgrid[i])*pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))))))/pow(rgrid[i],2)))/
   (16.*pow(L,2)*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h22.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((8*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (complex(0,48)*w*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
        ((-12 + pow(mu,2))*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (1 + rgrid[i] + pow(rgrid[i],2))*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (2*pow(rgrid[i],4)*((3*(1 + rgrid[i]*Q11(i,j))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) - 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    ((-2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (-((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]) + 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))))/(16.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h11(i,j)*((complex(0,-18)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          (-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j)))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (complex(0,12)*w*(-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))) - 
       ((-1 + rgrid[i])*pow(rgrid[i],6)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Q11(i,j))*
               pow(-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j),2))/pow(rgrid[i],2) + 
            (4*pow(Q22.diff(d01,i,j),2)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        (2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (2*(1 + rgrid[i]*Q22(i,j))*
          (-((36*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
                 24*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2))) + 
            (pow(rgrid[i],6)*
               ((4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                    exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
                    (1 + rgrid[i]*Qrr(i,j)))/
                  (pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    (exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-(L*
                       (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                        a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
                        rgrid[i]*Qr1(i,j),2) + 
                      (2*pow(L,2)*pow(w,2)*
                        (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/
             (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))))/pow(rgrid[i],2) - 
       (complex(0,3)*w*(-1 + rgrid[i])*pow(rgrid[i],10)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Q11(i,j))*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                      (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],3))))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    ((2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],10)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((-4*Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*
               pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],8)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((2*Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                    (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],4) - 
                 (L*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    ((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                      (Q11.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                    Qr1(i,j)*(1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                 (Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*
                    (-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
                  (2.*pow(rgrid[i],2)) + 
                 ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    ((4*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (pow(L,2)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],6) + 
                      (8*Q22.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],3)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))))*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2))))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],4) + 
                    ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       ((-2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q22.diff(d02,i,j)*
                        (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],3)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))))*
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2))) \
+ L*((-2*L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*
                       (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                       (1 + rgrid[i]*Qtt(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],8)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       (-(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],2) - 
                       (L*(6 - 
                       2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q22.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q22(i,j)))/pow(rgrid[i],4) + 
                       (Q22.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                       rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],2)) + 
                        (L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  ((L*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q22.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)))*
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) - 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-4*Q22.diff(d01,i,j)*Qr1.diff(d01,i,j\
) + (4*L*(-Q22.diff(d01,i,j) + Q22.diff(d11,i,j)*rgrid[i]))/
                        pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
                         L*((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) + 
                         (Q22.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))))/
             pow(rgrid[i],4)))/
        (2.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))))/4. + 
  (h22.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((-8*L*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) - 
       (complex(0,48)*L*w*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j)*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-12 + pow(mu,2))*(-1 + rgrid[i])*rgrid[i]*
          (1 + rgrid[i] + pow(rgrid[i],2))*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (2*pow(rgrid[i],6)*((-2*Q11.diff(d01,i,j)*
               (1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
               (1 + rgrid[i]*Qtt(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],9)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       pow(rgrid[i],3) + 
                      (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                       pow(rgrid[i],3) - 
                      (L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*Qr1(i,j)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       (2.*pow(rgrid[i],4)) + 
                      (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))/
                       (2.*rgrid[i])))/pow(rgrid[i],2)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               ((-2*L*(1 + rgrid[i]*Q22(i,j))*
                    (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],6) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q22(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (3*Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  ((-3*L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) + 
                    ((1 + rgrid[i]*Q22(i,j))*
                       ((2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((2*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] - 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],4)))/
        ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j)))))/
   (16.*pow(L,2)*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h22(i,j)*((complex(0,-18)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          pow(rgrid[i],4)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j)))/
        (pow(-12 + pow(mu,2),2)*
          pow(1 + rgrid[i] + pow(rgrid[i],2),2)*(1 + rgrid[i]*Qrr(i,j))\
) + (complex(0,12)*w*(-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j)))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],8)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((3*(1 + rgrid[i]*Q11(i,j))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) - 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    ((-2*L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (-((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]) + 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (2.*L*(1 + rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)*
          (1 + rgrid[i]*Qtt(i,j))) + 
       (complex(0,1.5)*w*pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          ((8*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*pow(rgrid[i],4)*
               ((3*(1 + rgrid[i]*Q11(i,j))*
                    ((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],6) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) - 
                      rgrid[i]*Qr1(i,j)*
                       (-((Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Q11(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5)) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q11(i,j)))/
                        (2.*pow(rgrid[i],3)) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],3))))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                      ((1 + rgrid[i]*Q11(i,j))*
                        ((-2*L*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                       ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        (-((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]) + 
                        (L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                       pow(rgrid[i],2)))/pow(rgrid[i],2)))/
             (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          pow(1 + rgrid[i]*Qrr(i,j),2)) + 
       (1 - rgrid[i])*(((-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             pow(-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
               Q22.diff(d01,i,j)*Qr1(i,j),2))/
           (2.*pow(L,2)*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) \
+ ((1 + rgrid[i]*Q22(i,j))*((72*pow(A1,2)*
                   exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
                  48*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) - 
               (2*pow(rgrid[i],4)*
                  (exp(B1*mu*h(i,j)*rgrid[i])*
                     pow(-(L*
                        (-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                        Qr1(i,j),2) + 
                    (2*pow(L,2)*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))*
                  (1 + rgrid[i]*Qtt(i,j)))))/pow(rgrid[i],2) + 
          ((-1 + rgrid[i])*pow(rgrid[i],8)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  (-((Q22.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                         (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j))\
)/pow(rgrid[i],6)) + (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                         ((Q11.diff(d01,i,j)*Q22.diff(d01,i,j))/
                        pow(rgrid[i],2) + 
                        (2*Q22.diff(d02,i,j)*
                        (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3))*
                         (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))) + 
               rgrid[i]*Qr1(i,j)*
                ((L*((Q22.diff(d01,i,j)*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                       (Q11.diff(d01,i,j)*
                        (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                   pow(rgrid[i],4) + 
                  ((1 + rgrid[i]*Q11(i,j))*
                     (-(L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q22.diff(d01,i,j)*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                         ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) - 
                         (2*Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                         (2*L*(-1 + rgrid[i])*
                        (-Q22.diff(d01,i,j) + 
                        Q22.diff(d11,i,j)*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                         (Q22.diff(d01,i,j)*L*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2)) + 
               L*(-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                         rgrid[i]*Q11(i,j))*
                       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                         rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],10)) + 
                  ((1 + rgrid[i]*Q11(i,j))*
                     ((2*L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.)*
                         (1 + rgrid[i]*Qtt(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],8)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                       (2*(1 + rgrid[i]*Qrr(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],2) - 
                        (L*(6 - 
                       2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q22.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q22(i,j)))/pow(rgrid[i],4) + 
                        (Q22.diff(d01,i,j)*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                         pow(rgrid[i],2) - 
                         (L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                         ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                   pow(rgrid[i],2))))/
           (2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
             pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j))))))/4.
;
