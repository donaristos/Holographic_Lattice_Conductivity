
elem[0]=
A0.diff(d20,i,j)*(1 - rgrid[i]) - 
  (2*A0.diff(d10,i,j)*(-12 + pow(mu,2) + 
       (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2)))/
   ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))) - 
  (complex(0,6)*w*A0(i,j)*rgrid[i]*
     (24 - 2*pow(mu,2) + 2*(-12 + pow(mu,2))*rgrid[i] + 
       2*(-12 + pow(mu,2))*pow(rgrid[i],2) + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],3)))/
   (pow(-12 + pow(mu,2),2)*(-1 + rgrid[i])*
     pow(1 + rgrid[i] + pow(rgrid[i],2),2)) + 
  (2*A0.diff(d11,i,j)*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))/L + 
  (h22.diff(d01,i,j)*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/pow(L,2) \
+ (htt.diff(d10,i,j)*(-(mu*a0(i,j)) - 
       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]) + 
       (a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))/L))/2. + 
  h22.diff(d10,i,j)*(1 - rgrid[i])*
   (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]) + 
     (a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))/L) + 
  Dh.diff(d10,i,j)*B1*rgrid[i]*
   (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]) + 
     (a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))/L) + 
  h22(i,j)*(mu*a0(i,j) + a0.diff(d10,i,j)*mu*(-1 + rgrid[i]) - 
     (a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))/L + 
     (complex(0,6)*w*pow(rgrid[i],2)*
        (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
          a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
      (L*(-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))) + 
     ((1 - rgrid[i])*(rgrid[i]*
           (-2*L*(-(a0.diff(d01,i,j)*mu) - 
                a0.diff(d11,i,j)*mu*(-1 + rgrid[i])) - 
             a0.diff(d01,i,j)*Qr1.diff(d01,i,j)*mu*(-1 + rgrid[i])*
              rgrid[i])*Qr1(i,j) - 
          a0.diff(d02,i,j)*mu*(-1 + rgrid[i])*pow(rgrid[i],2)*
           pow(Qr1(i,j),2) + 
          L*(L*(-2*a0.diff(d10,i,j)*mu - 
                a0.diff(d20,i,j)*mu*(-1 + rgrid[i])) + 
             a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
              (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))))/pow(L,2)) + 
  (2*a0.diff(d01,i,j)*h11.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (htt.diff(d01,i,j)*pow(rgrid[i],2)*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
  (Dh.diff(d01,i,j)*B1*pow(rgrid[i],3)*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (pow(L,2)*(1 + rgrid[i]*Q11(i,j))) - 
  (A0.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     ((1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
  (htt(i,j)*((complex(0,12)*L*w*pow(rgrid[i],2)*
          (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       2*(rgrid[i]*(-2*L*(-(a0.diff(d01,i,j)*mu) - 
                a0.diff(d11,i,j)*mu*(-1 + rgrid[i])) - 
             a0.diff(d01,i,j)*Qr1.diff(d01,i,j)*mu*(-1 + rgrid[i])*
              rgrid[i])*Qr1(i,j) - 
          a0.diff(d02,i,j)*mu*(-1 + rgrid[i])*pow(rgrid[i],2)*
           pow(Qr1(i,j),2) + 
          L*(L*(-2*a0.diff(d10,i,j)*mu - 
                a0.diff(d20,i,j)*mu*(-1 + rgrid[i])) + 
             a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
              (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))) + 
       (2*a0.diff(d01,i,j)*Q11.diff(d01,i,j)*mu*rgrid[i]*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) \
+ (pow(rgrid[i],2)*((-2*a0.diff(d01,i,j)*Qrr.diff(d01,i,j)*mu)/
             (rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               (-2*a0.diff(d02,i,j)*mu*(-1 + rgrid[i]) - 
                 a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
                  (2*h.diff(d01,i,j)*B1*mu*rgrid[i] + 
                    (Q22.diff(d01,i,j)*rgrid[i])/
                     (1 + rgrid[i]*Q22(i,j)) - 
                    (Qtt.diff(d01,i,j)*rgrid[i])/
                     (1 + rgrid[i]*Qtt(i,j)))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/(1 + rgrid[i]*Q11(i,j)))\
)/(4.*pow(L,2)) + (B1*Dh(i,j)*
     (rgrid[i]*(-4*L*(-(a0.diff(d01,i,j)*mu) - 
             a0.diff(d11,i,j)*mu*(-1 + rgrid[i]))*rgrid[i] - 
          a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
           (-2*L + 2*Qr1.diff(d01,i,j)*pow(rgrid[i],2) - 
             (complex(0,12)*L*w*pow(rgrid[i],3))/
              ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3)))))*Qr1(i,j) - 
       2*a0.diff(d02,i,j)*mu*(-1 + rgrid[i])*pow(rgrid[i],3)*
        pow(Qr1(i,j),2) + 2*L*
        (L*(-2*a0.diff(d10,i,j)*mu - 
             a0.diff(d20,i,j)*mu*(-1 + rgrid[i]))*rgrid[i] + 
          (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
           (L + (complex(0,6)*L*w*pow(rgrid[i],3))/
              ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3)))) + 
          a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
           (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))) + 
       (2*a0.diff(d01,i,j)*Q11.diff(d01,i,j)*mu*pow(rgrid[i],2)*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) \
+ (pow(rgrid[i],3)*((-2*a0.diff(d01,i,j)*Qrr.diff(d01,i,j)*mu)/
             (rgrid[i]*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Qrr(i,j))*
               (-2*a0.diff(d02,i,j)*mu*(-1 + rgrid[i]) - 
                 a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
                  (2*h.diff(d01,i,j)*B1*mu*rgrid[i] + 
                    (Q22.diff(d01,i,j)*rgrid[i])/
                     (1 + rgrid[i]*Q22(i,j)) - 
                    (Qtt.diff(d01,i,j)*rgrid[i])/
                     (1 + rgrid[i]*Qtt(i,j)))))/
             ((-1 + rgrid[i])*pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/(1 + rgrid[i]*Q11(i,j)))\
)/(2.*pow(L,2)) + (h11(i,j)*pow(rgrid[i],8)*
     (-((a0.diff(d01,i,j)*Q11.diff(d01,i,j)*mu*(-1 + rgrid[i])*
            (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
            (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) - 
       ((1 + rgrid[i]*Q11(i,j))*
          (-((a0.diff(d01,i,j)*Q22.diff(d01,i,j)*mu*(-1 + rgrid[i])*
                 (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
               pow(rgrid[i],5)) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (-((a0.diff(d01,i,j)*Qrr.diff(d01,i,j)*mu*
                      (-1 + rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                    pow(rgrid[i],3)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((a0.diff(d01,i,j)*Qtt.diff(d01,i,j)*mu*
                        pow(-1 + rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) + 
                      ((-1 + rgrid[i])*
                        (-(a0.diff(d02,i,j)*mu*(-1 + rgrid[i])) - 
                        a0.diff(d01,i,j)*h.diff(d01,i,j)*B1*
                        pow(mu,2)*(-1 + rgrid[i])*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/pow(rgrid[i],2)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d01,i,j)*(rgrid[i]*
        (-2*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*rgrid[i] + 
          4*L*(1 + (complex(0,6)*w*pow(rgrid[i],2))/
              ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2)))))*
        Qr1(i,j) + (2*pow(rgrid[i],8)*
          ((L*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
               (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],8) + 
            (Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],7) - 
            ((1 + rgrid[i]*Q11(i,j))*
               ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                       pow(rgrid[i],3) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
          (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/(2.*pow(L,2)) \
- (complex(0,0.25)*Ax.diff(d11,i,j)*(-1 + rgrid[i])*pow(rgrid[i],8)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
       ((1 + rgrid[i]*Q11(i,j))*
          (((-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
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
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 rgrid[i]*Qr1(i,j)*
                  ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],3) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
- (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j)*rgrid[i] - 
                        B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/pow(rgrid[i],2)))/
   (pow(L,2)*w*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (complex(0,0.25)*Ax.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],9)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
     (((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
       ((1 + rgrid[i]*Q11(i,j))*
          (((-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
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
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 rgrid[i]*Qr1(i,j)*
                  ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],3) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
- (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j)*rgrid[i] - 
                        B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/pow(rgrid[i],2)))/
   (pow(L,3)*w*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) + 
  (complex(0,0.25)*ht1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],8)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))*
     (((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
       ((1 + rgrid[i]*Q11(i,j))*
          (((-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
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
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 rgrid[i]*Qr1(i,j)*
                  ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],3) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
- (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j)*rgrid[i] - 
                        B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/pow(rgrid[i],2)))/
   (pow(L,3)*w*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))) - 
  (complex(0,0.125)*Ax.diff(d10,i,j)*(-1 + rgrid[i])*pow(rgrid[i],16)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (-((Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
            (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
          pow(rgrid[i],7)) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                    pow(rgrid[i],3)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)\
))/pow(rgrid[i],2))*(((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
       ((1 + rgrid[i]*Q11(i,j))*
          (((-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
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
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 rgrid[i]*Qr1(i,j)*
                  ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],3) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
- (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j)*rgrid[i] - 
                        B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/pow(rgrid[i],2)))/
   (pow(L,2)*w*pow(1 + rgrid[i]*Q11(i,j),3)*
     pow(1 + rgrid[i]*Q22(i,j),2)*pow(1 + rgrid[i]*Qrr(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,0.125)*ht1(i,j)*(-1 + rgrid[i])*pow(rgrid[i],16)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
          (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],7) + ((1 + rgrid[i]*Q11(i,j))*
          ((Q22.diff(d01,i,j)*
               (-(L*(-(mu*a0(i,j)) - 
                      a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                 a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                  Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((Qrr.diff(d01,i,j)*L*
                    (-(mu*a0(i,j)) - 
                      a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (-(Qtt.diff(d01,i,j)*L*
                       (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-(a0.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                       mu*(-1 + rgrid[i])*rgrid[i]) - 
                       L*(-(a0.diff(d01,i,j)*mu) - 
                       a0.diff(d11,i,j)*mu*(-1 + rgrid[i]) + 
                       h.diff(d01,i,j)*B1*mu*
                       (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                       rgrid[i]))*(1 + rgrid[i]*Qtt(i,j)))/
                       pow(rgrid[i],2)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 rgrid[i]*Qr1(i,j)*
                  ((a0.diff(d01,i,j)*Qrr.diff(d01,i,j)*mu*
                       (-1 + rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],3) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       (-(a0.diff(d01,i,j)*Qtt.diff(d01,i,j)*mu*
                        pow(-1 + rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ ((-1 + rgrid[i])*(-(a0.diff(d02,i,j)*mu*(-1 + rgrid[i])) - 
                        a0.diff(d01,i,j)*h.diff(d01,i,j)*B1*
                        pow(mu,2)*(-1 + rgrid[i])*rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))))/
             pow(rgrid[i],2)))/pow(rgrid[i],2))*
     (((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
       ((1 + rgrid[i]*Q11(i,j))*
          (((-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
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
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 rgrid[i]*Qr1(i,j)*
                  ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],3) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
- (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j)*rgrid[i] - 
                        B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/pow(rgrid[i],2)))/
   (pow(L,3)*w*pow(1 + rgrid[i]*Q11(i,j),3)*
     pow(1 + rgrid[i]*Q22(i,j),2)*pow(1 + rgrid[i]*Qrr(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (muj*(-1 + rgrid[i])*pow(rgrid[i],16)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((-6*L*pow(rgrid[i],2)*
          (-((Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                 (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
               pow(rgrid[i],7)) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],3)) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2))*
          (((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
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
                      rgrid[i]*Qr1(i,j)*
                       ((Qrr.diff(d01,i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
- (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) - 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (Qr1.diff(d01,i,j)*rgrid[i] - 
                       B1*L*mu*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (complex(0,1)*((Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
               pow(1 + rgrid[i]*Q22(i,j),2)*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*
               pow(1 + rgrid[i]*Qtt(i,j),2))/pow(rgrid[i],12) - 
            ((1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))*
               ((Qr1.diff(d01,i,j)*L*
                    ((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) - 
                      (Q11.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],3) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((8*Q11.diff(d01,i,j)*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],5)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (Qr1.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((-2*Qrr.diff(d01,i,j)*
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
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       ((L*
                      (Qr1.diff(d02,i,j)*rgrid[i] + 
                      h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*
                      mu*pow(rgrid[i],2))*
                      (-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) + 
                       (Q11.diff(d01,i,j)*
                       (pow(Qr1.diff(d01,i,j),2)*
                      pow(rgrid[i],2) - 
                       Qr1.diff(d01,i,j)*B1*L*mu*rgrid[i]*
                      (h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
                       Qr1.diff(d02,i,j)*pow(rgrid[i],2)*
                       Qr1(i,j)))/rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],2) + 
                        Qr1.diff(d01,i,j)*rgrid[i]*
                        ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) - 
                        (Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
                        rgrid[i] + 
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
)))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],8) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               ((Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                    (-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    pow(1 + rgrid[i]*Qrr(i,j),2)*
                    pow(1 + rgrid[i]*Qtt(i,j),2))/pow(rgrid[i],8) + 
                 ((1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j))*
                    ((8*Q22.diff(d01,i,j)*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],5)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*L*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) - 
                       (4*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                       Qr1(i,j))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (4*Q22.diff(d01,i,j)*L*
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
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                       (Qr1.diff(d01,i,j)*rgrid[i] - 
                       B1*L*mu*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i])) - 
                       (L*(Qr1.diff(d02,i,j)*rgrid[i] + 
                      h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*
                      mu*pow(rgrid[i],2))*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                       (Qr1.diff(d02,i,j)*rgrid[i] + 
                       2*h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*
                       mu*pow(rgrid[i],2))*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        Qr1.diff(d01,i,j)*L*rgrid[i]*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
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
)))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],6) + 
                 (pow(1 + rgrid[i]*Q22(i,j),2)*
                    ((Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                       ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) - 
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
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                       (2.*pow(rgrid[i],4)) + 
                      (16*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),3)*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                       (pow(-1 + rgrid[i],3)*pow(rgrid[i],6)*
                        pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),3)) + 
                      ((1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j))*
                        ((-2*Qrr.diff(d01,i,j)*
                       (Qr1.diff(d02,i,j)*rgrid[i] + 
                       2*h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*
                       mu*pow(rgrid[i],2))*Qr1(i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       ((-2*pow(Qr1.diff(d01,i,j),2)*Qrr.di\
ff(d01,i,j)*rgrid[i])/
                       ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (4*Qr1.diff(d02,i,j)*L*
                      (((-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 
                     4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))*
                      (-2 + 
                      Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))*
                      (1 + rgrid[i]*Qrr(i,j)))/2.))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                      4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) + 
                       Qr1.diff(d01,i,j)*B1*L*rgrid[i]*
                       ((2*Qrr.diff(d01,i,j)*mu*
                      (h(i,j) + h.diff(d10,i,j)*rgrid[i]))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (4*h.diff(d01,i,j)*mu*
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
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2))))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        Qr1.diff(d01,i,j)*L*rgrid[i]*
                        ((2*Qtt.diff(d01,i,j)*
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
                        pow(mu,2)*pow(rgrid[i],3))) - 
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
                       pow(rgrid[i],4) + 
                      (4*pow(1 + rgrid[i]*Qrr(i,j),2)*
                        (-((pow(-1 + rgrid[i],2)*
                        (Qr1.diff(d02,i,j)*rgrid[i] + 
                       h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                       pow(rgrid[i],2))*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)*
                        (-(Qr1.diff(d01,i,j)*rgrid[i]) + 
                       B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*pow(1 + rgrid[i]*Qtt(i,j),2))/pow(rgrid[i],4)) + 
                        rgrid[i]*Qr1(i,j)*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i])*
                        ((Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                        ((-1 + rgrid[i])*
                        (Qr1.diff(d02,i,j)*rgrid[i] + 
                       h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                       pow(rgrid[i],2))*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)) + 
                        (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*L*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
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
                        (2.*pow(rgrid[i],3)) + 
                        ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j))*
                        ((2*Qrr.diff(d01,i,j)*pow(L,2)*
                       pow(w,2))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (pow(Qr1.diff(d01,i,j),2)*Qtt.diff(d01\
,i,j)*(-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                        (Qr1.diff(d02,i,j)*L*
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
                        pow(rgrid[i],2) + 
                        Qr1.diff(d01,i,j)*B1*L*rgrid[i]*
                        (-(Qtt.diff(d01,i,j)*mu*(-1 + rgrid[i])*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*mu*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2)\
)))/pow(rgrid[i],2)))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2))))/
                  pow(rgrid[i],4)))/pow(rgrid[i],4)))/w))/
   (8.*pow(L,3)*pow(1 + rgrid[i]*Q11(i,j),3)*
     pow(1 + rgrid[i]*Q22(i,j),2)*pow(1 + rgrid[i]*Qrr(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (Ax(i,j)*(-1 + rgrid[i])*pow(rgrid[i],16)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((-6*L*pow(rgrid[i],2)*
          (-((Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
                 (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
               pow(rgrid[i],7)) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    (-((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],3)) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2))*
          (((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
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
                      rgrid[i]*Qr1(i,j)*
                       ((Qrr.diff(d01,i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
- (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) - 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (Qr1.diff(d01,i,j)*rgrid[i] - 
                       B1*L*mu*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (complex(0,1)*((Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
               pow(1 + rgrid[i]*Q22(i,j),2)*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*
               pow(1 + rgrid[i]*Qtt(i,j),2))/pow(rgrid[i],12) - 
            ((1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))*
               ((Qr1.diff(d01,i,j)*L*
                    ((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) - 
                      (Q11.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],3) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((8*Q11.diff(d01,i,j)*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],5)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (Qr1.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((-2*Qrr.diff(d01,i,j)*
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
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       ((L*
                      (Qr1.diff(d02,i,j)*rgrid[i] + 
                      h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*
                      mu*pow(rgrid[i],2))*
                      (-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) + 
                       (Q11.diff(d01,i,j)*
                       (pow(Qr1.diff(d01,i,j),2)*
                      pow(rgrid[i],2) - 
                       Qr1.diff(d01,i,j)*B1*L*mu*rgrid[i]*
                      (h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
                       Qr1.diff(d02,i,j)*pow(rgrid[i],2)*
                       Qr1(i,j)))/rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],2) + 
                        Qr1.diff(d01,i,j)*rgrid[i]*
                        ((Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) - 
                        (Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/
                        rgrid[i] + 
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
)))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],8) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               ((Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                    (-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    pow(1 + rgrid[i]*Qrr(i,j),2)*
                    pow(1 + rgrid[i]*Qtt(i,j),2))/pow(rgrid[i],8) + 
                 ((1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j))*
                    ((8*Q22.diff(d01,i,j)*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],5)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*L*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) - 
                       (4*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                       Qr1(i,j))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (4*Q22.diff(d01,i,j)*L*
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
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                       (Qr1.diff(d01,i,j)*rgrid[i] - 
                       B1*L*mu*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i])) - 
                       (L*(Qr1.diff(d02,i,j)*rgrid[i] + 
                      h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*
                      mu*pow(rgrid[i],2))*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                       (Qr1.diff(d02,i,j)*rgrid[i] + 
                       2*h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*
                       mu*pow(rgrid[i],2))*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        Qr1.diff(d01,i,j)*L*rgrid[i]*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
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
)))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],6) + 
                 (pow(1 + rgrid[i]*Q22(i,j),2)*
                    ((Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                       ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) - 
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
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                       (2.*pow(rgrid[i],4)) + 
                      (16*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),3)*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                       (pow(-1 + rgrid[i],3)*pow(rgrid[i],6)*
                        pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),3)) + 
                      ((1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j))*
                        ((-2*Qrr.diff(d01,i,j)*
                       (Qr1.diff(d02,i,j)*rgrid[i] + 
                       2*h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*
                       mu*pow(rgrid[i],2))*Qr1(i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       ((-2*pow(Qr1.diff(d01,i,j),2)*Qrr.di\
ff(d01,i,j)*rgrid[i])/
                       ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (4*Qr1.diff(d02,i,j)*L*
                      (((-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 
                     4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))*
                      (-2 + 
                      Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qrr(i,j)))/2. + 
                      ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))*
                      (1 + rgrid[i]*Qrr(i,j)))/2.))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                      4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) + 
                       Qr1.diff(d01,i,j)*B1*L*rgrid[i]*
                       ((2*Qrr.diff(d01,i,j)*mu*
                      (h(i,j) + h.diff(d10,i,j)*rgrid[i]))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (4*h.diff(d01,i,j)*mu*
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
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2))))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        Qr1.diff(d01,i,j)*L*rgrid[i]*
                        ((2*Qtt.diff(d01,i,j)*
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
                        pow(mu,2)*pow(rgrid[i],3))) - 
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
                       pow(rgrid[i],4) + 
                      (4*pow(1 + rgrid[i]*Qrr(i,j),2)*
                        (-((pow(-1 + rgrid[i],2)*
                        (Qr1.diff(d02,i,j)*rgrid[i] + 
                       h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                       pow(rgrid[i],2))*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)*
                        (-(Qr1.diff(d01,i,j)*rgrid[i]) + 
                       B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*pow(1 + rgrid[i]*Qtt(i,j),2))/pow(rgrid[i],4)) + 
                        rgrid[i]*Qr1(i,j)*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i])*
                        ((Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                        ((-1 + rgrid[i])*
                        (Qr1.diff(d02,i,j)*rgrid[i] + 
                       h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                       pow(rgrid[i],2))*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)) + 
                        (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*L*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
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
                        (2.*pow(rgrid[i],3)) + 
                        ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j))*
                        ((2*Qrr.diff(d01,i,j)*pow(L,2)*
                       pow(w,2))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (pow(Qr1.diff(d01,i,j),2)*Qtt.diff(d01\
,i,j)*(-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                        (Qr1.diff(d02,i,j)*L*
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
                        pow(rgrid[i],2) + 
                        Qr1.diff(d01,i,j)*B1*L*rgrid[i]*
                        (-(Qtt.diff(d01,i,j)*mu*(-1 + rgrid[i])*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*mu*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2)\
)))/pow(rgrid[i],2)))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2))))/
                  pow(rgrid[i],4)))/pow(rgrid[i],4)))/w))/
   (8.*pow(L,3)*pow(1 + rgrid[i]*Q11(i,j),3)*
     pow(1 + rgrid[i]*Q22(i,j),2)*pow(1 + rgrid[i]*Qrr(i,j),2)*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (Ax.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],14)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((-24*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*
          (((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
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
                      rgrid[i]*Qr1(i,j)*
                       ((Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],3) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
- (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) - 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j)*rgrid[i] - 
                       B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
        ((-12 + pow(mu,2))*(-1 + rgrid[i])*pow(rgrid[i],4)*
          (-1 + pow(rgrid[i],3))*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       (complex(0,2)*pow(rgrid[i],2)*
          ((Q11.diff(d01,i,j)*pow(1 + rgrid[i]*Q22(i,j),2)*Qr1(i,j)*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2)*
               pow(1 + rgrid[i]*Qtt(i,j),2))/pow(rgrid[i],12) + 
            ((1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))*
               ((L*(-((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4)) + 
                      (Q11.diff(d01,i,j)*
                        (-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                    Qr1(i,j)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((2*Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        pow(Qr1(i,j),2)*(1 + rgrid[i]*Qrr(i,j)))/
                       pow(rgrid[i],2) - 
                      (4*Qr1.diff(d01,i,j)*L*
                        (-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
                      rgrid[i]*Qr1(i,j)*
                       ((L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) - 
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
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((Q11.diff(d01,i,j)*
                       (Qr1.diff(d01,i,j)*rgrid[i] + 
                       B1*L*mu*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i])))/
                       rgrid[i] - 
                       (h.diff(d01,i,j)*B1*L*mu*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) - 
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
)))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],8) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               ((Q22.diff(d01,i,j)*Qr1(i,j)*
                    (-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    pow(1 + rgrid[i]*Qrr(i,j),2)*
                    pow(1 + rgrid[i]*Qtt(i,j),2))/pow(rgrid[i],8) + 
                 ((1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j))*
                    ((-4*Qr1.diff(d01,i,j)*L*
                        (-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
                      (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        pow(Qr1(i,j),2)*
                        ((-2*Qrr.diff(d01,i,j))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*h.diff(d01,i,j)*B1*mu*
                       (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] + 
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
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (2*(1 + rgrid[i]*Qrr(i,j))*
                        (-(Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/(2.*pow(rgrid[i],4)) - 
                        ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((Q22.diff(d01,i,j)*
                       (-3*Qr1.diff(d01,i,j)*rgrid[i] + 
                       B1*L*mu*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i])))/
                       rgrid[i] + 
                       (h.diff(d01,i,j)*B1*L*mu*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],2))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))))))/
                  pow(rgrid[i],6) + 
                 (pow(1 + rgrid[i]*Q22(i,j),2)*
                    ((16*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),3)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],8)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (Qrr.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
                        ((2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                        ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) - 
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
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                       (2.*pow(rgrid[i],4)) + 
                      (4*pow(1 + rgrid[i]*Qrr(i,j),2)*
                        ((2*Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i] + 
                        rgrid[i]*Qr1(i,j)*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j)*rgrid[i] - 
                       B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        rgrid[i]*Qr1(i,j)*
                        (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
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
))/(pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                      ((1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))*
                         (-4*h.diff(d01,i,j)*Qrr.diff(d01,i,j)*B1*
                        mu*pow(Qr1(i,j),2)*(1 + rgrid[i]*Qtt(i,j)) + 
                         (8*Qr1.diff(d01,i,j)*L*
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
                        ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                         rgrid[i]*Qr1(i,j)*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((-6*Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j\
))/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (2*Qrr.diff(d01,i,j)*B1*L*mu*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i]))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*h.diff(d01,i,j)*B1*L*mu*
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
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
                        pow(mu,2)*pow(rgrid[i],3))) - 
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
                         pow(mu,2)*pow(rgrid[i],3)))))))/
                       pow(rgrid[i],4)))/pow(rgrid[i],4)))/
             pow(rgrid[i],4)))/
        (w*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))/
   (16.*pow(L,3)*pow(1 + rgrid[i]*Q11(i,j),3)*
     pow(1 + rgrid[i]*Q22(i,j),2)*pow(1 + rgrid[i]*Qrr(i,j),2))
;
elem[1]=
Ax.diff(d20,i,j) - (2*Ax.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  (Ax.diff(d02,i,j)*pow(rgrid[i],2)*pow(Qr1(i,j),2))/pow(L,2) + 
  (ht1.diff(d01,i,j)*rgrid[i]*Qr1(i,j)*
     (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/pow(L,2) \
+ ht1.diff(d10,i,j)*(-(mu*a0(i,j)) - 
     a0.diff(d10,i,j)*mu*(-1 + rgrid[i]) + 
     (a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))/L) + 
  (complex(0,4)*A0.diff(d01,i,j)*w*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,4)*a0.diff(d01,i,j)*mu*w*h11(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,2)*a0.diff(d01,i,j)*mu*w*htt(i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,4)*a0.diff(d01,i,j)*B1*mu*w*Dh(i,j)*rgrid[i]*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))) + 
  (Ax.diff(d10,i,j)*((-2*Qr1.diff(d01,i,j)*rgrid[i])/L + 
       2*B1*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) + 
       (complex(0,24)*w*pow(rgrid[i],2))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
       (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - rgrid[i]*Q11(i,j))/
        (rgrid[i]*(1 + rgrid[i]*Q11(i,j))) + 
       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - rgrid[i]*Q22(i,j))/
        (rgrid[i]*(1 + rgrid[i]*Q22(i,j))) - 
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
       (rgrid[i]*Qr1(i,j)*(-2*h.diff(d01,i,j)*B1*mu*rgrid[i] + 
            (Q11.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q11(i,j)) - 
            (Q22.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) + 
            (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) - 
            (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))))/L + 
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
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))/2. + 
  (Ax.diff(d01,i,j)*(-2*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)) + 
       pow(rgrid[i],2)*pow(Qr1(i,j),2)*
        (2*h.diff(d01,i,j)*B1*mu*rgrid[i] - 
          (Q11.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q11(i,j)) + 
          (Q22.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) - 
          (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) + 
          (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))) + 
       rgrid[i]*Qr1(i,j)*(4*Qr1.diff(d01,i,j)*rgrid[i] + 
          L*(-2*B1*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
             (complex(0,24)*w*pow(rgrid[i],2))/
              ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
             (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q11(i,j))/(rgrid[i]*(1 + rgrid[i]*Q11(i,j))) - 
             (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q22(i,j))/(rgrid[i]*(1 + rgrid[i]*Q22(i,j))) + 
             (2*(((-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) \
- (2*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))\
))/(2.*pow(L,2)) + (ht1(i,j)*
     (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
        (-2*a0.diff(d02,i,j)*mu*(-1 + rgrid[i]) - 
          a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
           (2*h.diff(d01,i,j)*B1*mu*rgrid[i] - 
             (Q11.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q11(i,j)) + 
             (Q22.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) - 
             (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) + 
             (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j)))) + 
       L*(2*L*(-2*a0.diff(d10,i,j)*mu - 
             a0.diff(d20,i,j)*mu*(-1 + rgrid[i])) + 
          2*a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
           (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)) + 
          L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
           (2*B1*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) + 
             (complex(0,12)*w*pow(rgrid[i],2))/
              ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
             (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q11(i,j))/(rgrid[i]*(1 + rgrid[i]*Q11(i,j))) + 
             (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q22(i,j))/(rgrid[i]*(1 + rgrid[i]*Q22(i,j))) - 
             (2*(((-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) \
+ (2*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))))\
) + rgrid[i]*Qr1(i,j)*(L*(-4*(-(a0.diff(d01,i,j)*mu) - 
                a0.diff(d11,i,j)*mu*(-1 + rgrid[i])) - 
             2*h.diff(d01,i,j)*B1*mu*
              (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
              rgrid[i] + (Q11.diff(d01,i,j)*
                (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                rgrid[i])/(1 + rgrid[i]*Q11(i,j)) - 
             (Q22.diff(d01,i,j)*
                (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                rgrid[i])/(1 + rgrid[i]*Q22(i,j)) + 
             (Qrr.diff(d01,i,j)*
                (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) - 
             (Qtt.diff(d01,i,j)*
                (-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                rgrid[i])/(1 + rgrid[i]*Qtt(i,j))) - 
          a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*
           (2*Qr1.diff(d01,i,j)*rgrid[i] + 
             L*(-2*B1*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
                (complex(0,12)*w*pow(rgrid[i],2))/
                 ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
                (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j))/(rgrid[i]*(1 + rgrid[i]*Q11(i,j))) - 
                (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q22(i,j))/(rgrid[i]*(1 + rgrid[i]*Q22(i,j))) + 
                (2*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                     pow(mu,2)*pow(rgrid[i],3))*
                   (1 + rgrid[i]*Qrr(i,j))) - 
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
                     pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))\
))))))/(2.*pow(L,2)) + (muj*((complex(0,-36)*
          (-12 + pow(mu,2) - complex(0,2)*w)*w*pow(rgrid[i],4))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)) + 
       (complex(0,24)*w*rgrid[i])/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
       (complex(0,6)*w*pow(rgrid[i],10)*
          (((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q11(i,j))*
               (((-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
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
                      rgrid[i]*Qr1(i,j)*
                       (-((Qrr.diff(d01,i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3)) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (Qr1.diff(d01,i,j)*rgrid[i] - 
                       B1*L*mu*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) - 
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
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
       (pow(rgrid[i],8)*((Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((Qr1.diff(d01,i,j)*
                    (-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],3) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((8*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                        pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((-2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                       ((-1 + rgrid[i])*
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
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (rgrid[i]*Qr1(i,j)*
                        ((Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. + 
                        ((-1 + rgrid[i])*
                        (Qr1.diff(d02,i,j)*rgrid[i] + 
                       h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                       pow(rgrid[i],2))*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)) - 
                        L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j) + 
                       Qr1.diff(d11,i,j)*rgrid[i] + 
                       Qr1.diff(d01,i,j)*B1*mu*rgrid[i]*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2)\
)))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/2. + 
  (Ax(i,j)*((complex(0,-36)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          pow(rgrid[i],4))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)) + 
       (complex(0,24)*w*rgrid[i])/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
       (complex(0,6)*w*pow(rgrid[i],10)*
          (((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q11(i,j))*
               (((-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
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
                      rgrid[i]*Qr1(i,j)*
                       (-((Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3)) + 
                        (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]) \
+ (h.diff(d01,i,j)*B1*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j)*rgrid[i] - 
                       B1*L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) - 
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
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
       (pow(rgrid[i],8)*((Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((Qr1.diff(d01,i,j)*
                    (-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],3) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((8*pow(L,2)*pow(w,2)*
                        pow(1 + rgrid[i]*Qrr(i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((-2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                        ((-1 + rgrid[i])*
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
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                         (rgrid[i]*Qr1(i,j)*
                        ((Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2. + 
                        ((-1 + rgrid[i])*
                        (Qr1.diff(d02,i,j)*rgrid[i] + 
                        h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                        pow(rgrid[i],2))*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2)) - 
                         L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d01,i,j) + 
                        Qr1.diff(d11,i,j)*rgrid[i] + 
                        Qr1.diff(d01,i,j)*B1*mu*rgrid[i]*
                        (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
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
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/2.
;
elem[2]=
-((h11.diff(d01,i,j)*h.diff(d01,i,j)*mu*pow(rgrid[i],3))/
     (pow(L,2)*(1 + rgrid[i]*Q11(i,j)))) - 
  (complex(0,1)*h.diff(d01,i,j)*mu*w*ht1(i,j)*pow(rgrid[i],3))/
   (L*(1 + rgrid[i]*Q11(i,j))) + 
  (Dh.diff(d20,i,j)*(-1 + rgrid[i])*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3)))/(2.*(1 + rgrid[i]*Qrr(i,j))) - 
  (Dh.diff(d11,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/(L*(1 + rgrid[i]*Qrr(i,j))) \
+ (h22.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],3)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
     (L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
       h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/
   (2.*pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
  (htt.diff(d10,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (-(L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])) + 
       h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/
   (4.*L*(1 + rgrid[i]*Qrr(i,j))) + 
  (h22.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (-(L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])) + 
       h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/
   (2.*L*(1 + rgrid[i]*Qrr(i,j))) + 
  (Dh.diff(d02,i,j)*rgrid[i]*
     (pow(rgrid[i],2)/(1 + rgrid[i]*Q11(i,j)) + 
       ((-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2))/
        (2.*(1 + rgrid[i]*Qrr(i,j)))))/pow(L,2) + 
  (htt.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
            h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/rgrid[i] - 
       (2*h.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (4.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
  (complex(0,2)*a0.diff(d01,i,j)*B1*mu*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],4))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,2)*a0.diff(d01,i,j)*B1*mu*w*Ax(i,j)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],4))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d10,i,j)*B1*exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
     pow(rgrid[i],4)*(-(L*(-(mu*a0(i,j)) - 
            a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (B1*A0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],4)*
     (-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
     (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (A0.diff(d01,i,j)*B1*exp(B1*mu*h(i,j)*rgrid[i])*(1 - rgrid[i])*
     pow(rgrid[i],6)*(((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
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
     ((complex(0,-6)*L*w*pow(rgrid[i],2)*
          (L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
            h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (2*B1*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],4)*
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
+ (Dh.diff(d10,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     ((2*Qrr.diff(d01,i,j)*rgrid[i]*Qr1(i,j))/
        (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) - 
       (4*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Qrr(i,j)))/2. + 
            ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                 3*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))*
               (1 + rgrid[i]*Qrr(i,j)))/2.))/
        (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)) + 
       (2*(1 + rgrid[i]*Qrr(i,j))*
          (4 - (2*Qr1.diff(d01,i,j)*pow(rgrid[i],2))/L + 
            (complex(0,24)*w*pow(rgrid[i],3))/
             ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
            (pow(rgrid[i],3)*
               ((-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q11(i,j))/pow(rgrid[i],3) - 
                 (Q11.diff(d01,i,j)*Qr1(i,j))/L))/
             (1 + rgrid[i]*Q11(i,j)) + 
            (pow(rgrid[i],3)*
               ((-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q22(i,j))/pow(rgrid[i],3) - 
                 (Q22.diff(d01,i,j)*Qr1(i,j))/L))/
             (1 + rgrid[i]*Q22(i,j)) - 
            (Qtt.diff(d01,i,j)*pow(rgrid[i],3)*Qr1(i,j))/
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
             ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (8.*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (h11(i,j)*pow(rgrid[i],10)*
     ((h.diff(d01,i,j)*Q11.diff(d01,i,j)*mu*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],6) - ((1 + rgrid[i]*Q11(i,j))*
          ((h.diff(d01,i,j)*Q22.diff(d01,i,j)*mu*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],4) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((h.diff(d01,i,j)*Qrr.diff(d01,i,j)*mu*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (pow(a0.diff(d01,i,j),2)*B1*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2) + 
                      (h.diff(d01,i,j)*Qtt.diff(d01,i,j)*mu*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2. + 
                      (h.diff(d02,i,j)*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/pow(rgrid[i],2)))/
   (2.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (h22(i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*
     ((2*L*(-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
            (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
          (-(L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])) + 
            h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))) - 
       (2*h.diff(d01,i,j)*Q11.diff(d01,i,j)*mu*pow(rgrid[i],2)*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) \
+ (2*pow(rgrid[i],6)*((h.diff(d01,i,j)*Q22.diff(d01,i,j)*mu*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],4) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((h.diff(d01,i,j)*Qrr.diff(d01,i,j)*mu*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    (pow(a0.diff(d01,i,j),2)*B1*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2) + 
                      (h.diff(d01,i,j)*Qtt.diff(d01,i,j)*mu*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2. + 
                      (h.diff(d02,i,j)*mu*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/rgrid[i]))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2))\
)/((-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))*
          (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))))/(4.*pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) \
+ (Dh.diff(d01,i,j)*pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
     pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*
     (pow(rgrid[i],3)*pow(Qr1(i,j),2)*
        ((-2*Qrr.diff(d01,i,j))/
           ((-1 + rgrid[i])*rgrid[i]*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (2*Q11.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
           ((-1 + rgrid[i])*rgrid[i]*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
          (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
           ((-1 + rgrid[i])*rgrid[i]*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))) + 
          (2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
           ((-1 + rgrid[i])*rgrid[i]*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))) + 
       (4*pow(rgrid[i],7)*(1 + rgrid[i]*Qrr(i,j))*
          ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             pow(rgrid[i],7) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Q11(i,j))*
                    ((2*Qrr.diff(d01,i,j))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) - 
                      (2*L*(1 + rgrid[i]*Q11(i,j))*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                       pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                  (2.*pow(rgrid[i],4)) + 
                 (2*(1 + rgrid[i]*Qrr(i,j))*
                    ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/
                       (2.*pow(rgrid[i],3)) - 
                      (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)\
))/(pow(-1 + rgrid[i],2)*pow(-4 - 4*rgrid[i] - 
            4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3),2)*
          pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qtt(i,j))) + 
       rgrid[i]*Qr1(i,j)*((4*L*
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
           (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
             pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3),2)) + 
          (2*(1 + rgrid[i]*Qrr(i,j))*
             (4*Qr1.diff(d01,i,j)*pow(rgrid[i],2) + 
               L*(-4 - (complex(0,24)*w*pow(rgrid[i],3))/
                   ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
                  (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                     rgrid[i]*Q11(i,j))/(1 + rgrid[i]*Q11(i,j)) - 
                  (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                     rgrid[i]*Q22(i,j))/(1 + rgrid[i]*Q22(i,j)) - 
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
                   ((-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                     (1 + rgrid[i]*Qtt(i,j))))))/
           ((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))))))/
   (8.*pow(L,2)*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
  (Dh(i,j)*((complex(0,-18)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          (-1 + rgrid[i])*pow(rgrid[i],7)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (complex(0,12)*w*(-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       rgrid[i]*((16*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
             24*pow(A1,2)*exp(A1*mu*h(i,j)*rgrid[i]))/
           (2 + 3*pow(A1,2)) + 
          (pow(rgrid[i],6)*((2*pow(a0.diff(d01,i,j),2)*
                  pow(B1,2)*pow(mu,2)*exp(B1*mu*h(i,j)*rgrid[i])*
                  (-1 + rgrid[i])*(1 + rgrid[i]*Qrr(i,j)))/
                (pow(rgrid[i],2)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) + 
               ((1 + rgrid[i]*Q11(i,j))*
                  (pow(B1,2)*exp(B1*mu*h(i,j)*rgrid[i])*
                     pow(-(L*
                        (-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                        Qr1(i,j),2) + 
                    (4*pow(L,2)*pow(w,2)*(1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                pow(rgrid[i],2)))/
           (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
             (1 + rgrid[i]*Qtt(i,j)))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],10)*
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
        (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j))) + 
       (complex(0,3)*w*(-1 + rgrid[i])*pow(rgrid[i],12)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 Q11.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],3) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                        ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) - 
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
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                         (-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. + 
                         ((-1 + rgrid[i])*
                        (4*L - 
                        2*Qr1.diff(d01,i,j)*pow(rgrid[i],2))*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) \
+ (L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2 + 
                        Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                        3*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))))/2.
;
elem[3]=
h11.diff(d20,i,j) - (2*h11.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  (h11.diff(d02,i,j)*pow(rgrid[i],2)*pow(Qr1(i,j),2))/pow(L,2) + 
  (htt.diff(d10,i,j)*pow(rgrid[i],4)*
     (((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i] + 
            Q11.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
            Q22.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2)))/
   (4.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))) + 
  (h22.diff(d10,i,j)*(-1 + rgrid[i])*pow(rgrid[i],4)*
     (((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i] + 
            Q11.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
            Q22.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2)))/
   (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))) - 
  (htt.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (complex(0,2)*ht1.diff(d01,i,j)*w*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (4*Dh.diff(d01,i,j)*h.diff(d01,i,j)*mu*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (8*a0.diff(d01,i,j)*A0.diff(d01,i,j)*mu*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,8)*a0.diff(d01,i,j)*mu*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,8)*a0.diff(d01,i,j)*mu*w*Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (4*pow(a0.diff(d01,i,j),2)*B1*pow(mu,2)*Dh(i,j)*
     exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],3)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,2)*w*ht1(i,j)*pow(rgrid[i],6)*(1 + rgrid[i]*Qrr(i,j))*
     ((Q22.diff(d01,i,j)*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],5)) + 
       ((1 + rgrid[i]*Q22(i,j))*
          (-((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                 (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                   pow(mu,2)*pow(rgrid[i],3))*
                 (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],3))))/pow(rgrid[i],2)))/
   (L*pow(-1 + rgrid[i],2)*pow(-4 - 4*rgrid[i] - 
       4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3),2)*
     pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (h11.diff(d10,i,j)*((-2*Qr1.diff(d01,i,j)*rgrid[i])/L + 
       (complex(0,24)*w*pow(rgrid[i],2))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (pow(rgrid[i],2)*((-2 + 
               Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
               rgrid[i]*Q11(i,j))/pow(rgrid[i],3) - 
            (Q11.diff(d01,i,j)*Qr1(i,j))/L))/(1 + rgrid[i]*Q11(i,j)) + 
       (pow(rgrid[i],2)*((-2 + 
               Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
               rgrid[i]*Q22(i,j))/pow(rgrid[i],3) - 
            (Q22.diff(d01,i,j)*Qr1(i,j))/L))/(1 + rgrid[i]*Q22(i,j)) + 
       (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
        (L*(1 + rgrid[i]*Qrr(i,j))) - 
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
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) - 
       (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
        (L*(1 + rgrid[i]*Qtt(i,j))) + 
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
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))/2. - 
  (h22.diff(d01,i,j)*pow(rgrid[i],6)*
     (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          Qr1(i,j)*((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
            Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        (2.*pow(rgrid[i],3)) + 
       ((1 + rgrid[i]*Q22(i,j))*
          (-((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
               pow(rgrid[i],3)) - 
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
               (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i])))/pow(rgrid[i],2)))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (htt(i,j)*pow(rgrid[i],6)*((complex(0,1.5)*L*w*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
            Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],2) + ((1 + rgrid[i]*Q22(i,j))*
          ((-4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
               (-12 + pow(mu,2))*exp(B1*mu*h(i,j)*rgrid[i])*
               (-1 + rgrid[i])*(-1 + pow(rgrid[i],3))*
               (1 + rgrid[i]*Qrr(i,j)))/
             (pow(rgrid[i],2)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) - 
            complex(0,1.5)*L*w*(-1 + rgrid[i])*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*
             ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
               (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/
                rgrid[i] - Q11.diff(d01,i,j)*Qr1(i,j))*
             (1 + rgrid[i]*Qtt(i,j))))/pow(rgrid[i],2)))/
   (pow(L,2)*(-12 + pow(mu,2))*(-1 + rgrid[i])*(-1 + pow(rgrid[i],3))*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (htt.diff(d01,i,j)*pow(rgrid[i],8)*
     (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          pow(1 + rgrid[i]*Q11(i,j),2)*Qr1(i,j)*
          (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) - 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j)))/rgrid[i] + 
            Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        (2.*pow(rgrid[i],5)) + 
       (Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],7) + ((1 + rgrid[i]*Q11(i,j))*
          ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((-2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  pow(rgrid[i],3) + 
                 (L*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j)*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],4)) - 
                 (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i])))/
             pow(rgrid[i],2)))/pow(rgrid[i],2)))/
   (2.*pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (h11.diff(d01,i,j)*(-2*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)) + 
       pow(rgrid[i],2)*pow(Qr1(i,j),2)*
        ((Q11.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q11(i,j)) + 
          (Q22.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) - 
          (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) + 
          (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))) + 
       rgrid[i]*Qr1(i,j)*(4*Qr1.diff(d01,i,j)*rgrid[i] + 
          L*((complex(0,-24)*w*pow(rgrid[i],2))/
              ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
             (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q11(i,j))/(rgrid[i]*(1 + rgrid[i]*Q11(i,j))) - 
             (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q22(i,j))/(rgrid[i]*(1 + rgrid[i]*Q22(i,j))) + 
             (2*(((-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) \
- (2*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))\
))/(2.*pow(L,2)) + (h11(i,j)*
     ((complex(0,-72)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          pow(rgrid[i],4))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)) + 
       (complex(0,48)*w*rgrid[i])/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (4*(1 + rgrid[i]*Qrr(i,j))*
          ((72*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
               48*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) - 
            (2*pow(h.diff(d01,i,j),2)*pow(mu,2)*
               pow(rgrid[i],4))/(pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
            (pow(Qtt.diff(d01,i,j),2)*pow(rgrid[i],4))/
             (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
               pow(1 + rgrid[i]*Qtt(i,j),2)) + 
            (2*pow(rgrid[i],6)*
               ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                    (-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))/
                  (2.*pow(rgrid[i],2)) + 
                 (2*(pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2) - 
                      (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i]))*
                    (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) + 
                 (2*pow(L,2)*pow(w,2)*
                    pow(1 + rgrid[i]*Q11(i,j),2))/pow(rgrid[i],4)))/
             (pow(L,2)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) - 
       (complex(0,12)*w*pow(rgrid[i],10)*
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
          (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
       (pow(rgrid[i],8)*(((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (1 + rgrid[i]*Q11(i,j))*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*
               ((2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                  ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 (4*L*(((-1 + rgrid[i])*
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
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)))*
               (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],4)) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    ((8*pow(Qrr.diff(d01,i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (4*L*((L*
                      (-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                       (2*Qr1.diff(d01,i,j)*
                       (1 + rgrid[i]*Q11(i,j)))/rgrid[i])*
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
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
                 (2*Qrr.diff(d01,i,j)*rgrid[i]*pow(Qr1(i,j),2)*
                    ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) + 
                      (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       (2.*pow(rgrid[i],3))))/
                  ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (8*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
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
                  (pow(-1 + rgrid[i],2)*pow(rgrid[i],8)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) + 
                 rgrid[i]*Qr1(i,j)*
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
+ (2*(1 + rgrid[i]*Q11(i,j))*
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
                     pow(rgrid[i],2))))/pow(rgrid[i],2)))/
        (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
       (4*pow(rgrid[i],12)*((pow(-1 + rgrid[i],2)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)*
               pow(1 + rgrid[i]*Q22(i,j),2)*
               ((4*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (pow(L,2)*
                    pow(-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j),2))/pow(rgrid[i],6) - 
                 (2*Q11.diff(d01,i,j)*L*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j))/pow(rgrid[i],3) + 
                 pow(Q11.diff(d01,i,j),2)*pow(Qr1(i,j),2))*
               pow(1 + rgrid[i]*Qtt(i,j),2))/(4.*pow(rgrid[i],8)) - 
            ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
               pow(1 + rgrid[i]*Q22(i,j),2)*(1 + rgrid[i]*Qtt(i,j))*
               ((2*Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                  pow(rgrid[i],2) + 
                 (4*Qrr.diff(d02,i,j)*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],3) - 
                 (2*Qr1.diff(d01,i,j)*L*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
                 (pow(L,2)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (6 - 2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                      Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                      2*rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],6) - 
                 (Q11.diff(d01,i,j)*L*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/
                     (2.*pow(rgrid[i],2)) + 
                    (Q11.diff(d02,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3)) + 
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
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],6) \
- rgrid[i]*Qr1(i,j)*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (-6*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j) + 
                        (4*L*
                       (-Q11.diff(d01,i,j) + 
                       Q11.diff(d11,i,j)*rgrid[i]))/
                        pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                     (2.*pow(rgrid[i],2)) + 
                    L*((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
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
)))/(2.*pow(rgrid[i],8)) + (pow(1 + rgrid[i]*Q11(i,j),2)*
               ((pow(-1 + rgrid[i],2)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)*
                    pow(-((L*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j),2)*
                    pow(1 + rgrid[i]*Qtt(i,j),2))/
                  (4.*pow(rgrid[i],4)) + 
                 ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))*
                    (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                        ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/
                        (2.*pow(rgrid[i],2)) + 
                        (Q22.diff(d02,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3))) - 
                      L*(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       ((2*L*
                      (6 - 
                      2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
                      Q22.diff(d20,i,j)*pow(rgrid[i],3) + 
                      2*rgrid[i]*Q22(i,j)))/pow(rgrid[i],4) - 
                       (2*Q22.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                       rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
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
) + rgrid[i]*Qr1(i,j)*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-2*Q22.diff(d01,i,j)*Qr1.diff(d01,i,j\
) + (4*L*(-Q22.diff(d01,i,j) + Q22.diff(d11,i,j)*rgrid[i]))/
                       pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
                        L*((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
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
))))/(2.*pow(rgrid[i],4)) + 
                 (2*pow(1 + rgrid[i]*Q22(i,j),2)*
                    (-(pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)*
                        (pow(Qr1.diff(d01,i,j),2)*
                       pow(rgrid[i],2) + 
                        pow(L,2)*pow(mu,2)*
                       pow(h(i,j) + h.diff(d10,i,j)*rgrid[i],2) \
- L*(Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i]))*
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                       (2.*pow(rgrid[i],4)) + 
                      (pow(L,2)*
                        pow(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qtt(i,j)))/2.,2))/
                       pow(rgrid[i],6) + 
                      pow(rgrid[i],2)*pow(Qr1(i,j),2)*
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
                       pow(mu,2)*pow(rgrid[i],3)))/(2.*rgrid[i])\
)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) - 
                        (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                        pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)*
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                        (2.*pow(rgrid[i],2))) - 
                      rgrid[i]*Qr1(i,j)*
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (-4*a0.diff(d01,i,j)*L*mu*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                       (-1 + rgrid[i]) + 
                       (3*Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. - 
                       (4*L*
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
                       pow(mu,2)*pow(rgrid[i],3)))/2.)))/
                       pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
                        (pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)*
                        (Qr1.diff(d02,i,j)*rgrid[i] - 
                       2*h.diff(d01,i,j)*L*pow(mu,2)*rgrid[i]*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                        (2.*pow(rgrid[i],4)) + 
                        (Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
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
) + (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j))*
                        (2*L*exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) + 
                        (Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i] + 
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
- (2*L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                        pow(rgrid[i],4)))/(2.*pow(rgrid[i],2))))/
                  pow(rgrid[i],4)))/pow(rgrid[i],4)))/
        (pow(L,2)*pow(-1 + rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
          pow(1 + rgrid[i]*Qtt(i,j),2))))/4. + 
  (h22(i,j)*pow(rgrid[i],8)*((complex(0,12)*L*w*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j))*
          (((1 + rgrid[i]*Q22(i,j))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/
                  rgrid[i] + Q11.diff(d01,i,j)*Qr1(i,j)))/
             pow(rgrid[i],2) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2)))/
        ((-12 + pow(mu,2))*pow(rgrid[i],2)*
          (1 + rgrid[i] + pow(rgrid[i],2))) - 
       (2*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (((1 + rgrid[i]*Q22(i,j))*
               ((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                 (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/
                  rgrid[i] - Q11.diff(d01,i,j)*Qr1(i,j)))/
             pow(rgrid[i],2) + 
            ((1 + rgrid[i]*Q11(i,j))*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2)))/
        pow(rgrid[i],4) + ((-1 + rgrid[i])*pow(rgrid[i],4)*
          (-((pow(1 + rgrid[i]*Q11(i,j),2)*
                 pow(-((L*
                        (-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                   Q22.diff(d01,i,j)*Qr1(i,j),2)*
                 (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
               pow(rgrid[i],8)) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
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
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((2*L*
                       (6 - 
                       2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q22.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q22(i,j)))/pow(rgrid[i],4) - 
                       (2*Q22.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                       rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
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
                        (-2*Q22.diff(d01,i,j)*Qr1.diff(d01,i,j\
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
             pow(rgrid[i],6) + 
            (pow(1 + rgrid[i]*Q22(i,j),2)*
               ((pow(L,2)*pow(-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],10) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  ((Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                       (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                     pow(rgrid[i],6) + 
                    (2*(1 + rgrid[i]*Qrr(i,j))*
                       (-(Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/(2.*pow(rgrid[i],4)) \
+ ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (pow(Q11.diff(d01,i,j),2)/
                        pow(rgrid[i],2) - 
                        (2*Q11.diff(d02,i,j)*
                        (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 (pow(1 + rgrid[i]*Q11(i,j),2)*
                    ((-4*Qr1.diff(d01,i,j)*L*
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
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (2*(1 + rgrid[i]*Qrr(i,j))*
                        ((-2*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (pow(Qr1.diff(d01,i,j),2)*
                       pow(rgrid[i],2) - 
                       L*(Qr1.diff(d01,i,j) + 
                       Qr1.diff(d11,i,j)*rgrid[i]))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                        (2*Qr1.diff(d01,i,j)*L*
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
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],2)\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],4) + 
                 (L*(1 + rgrid[i]*Q11(i,j))*
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
                         pow(mu,2)*pow(rgrid[i],3)))))/
                  pow(rgrid[i],2) - 
                 rgrid[i]*Qr1(i,j)*
                  ((2*Q11.diff(d01,i,j)*L*
                       (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],8) + 
                    (pow(1 + rgrid[i]*Q11(i,j),2)*
                       ((2*Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],2) + 
                         ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((-4*Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j\
))/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (8*Qr1.diff(d02,i,j)*
                        (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
)/pow(rgrid[i],4) + ((1 + rgrid[i]*Q11(i,j))*
                       ((L*(-1 + rgrid[i])*
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
- (2*(1 + rgrid[i]*Qrr(i,j))*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-6*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j\
) + (4*L*(-Q11.diff(d01,i,j) + Q11.diff(d11,i,j)*rgrid[i]))/
                        pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
                         L*((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
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
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],4)))/
        ((1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/
   (4.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Q22(i,j),2))
;
elem[4]=
htt.diff(d20,i,j) - (2*htt.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  Dh.diff(d10,i,j)*rgrid[i]*(-(mu*
        (h(i,j) + h.diff(d10,i,j)*rgrid[i])) + 
     (h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j))/L) + 
  (h11.diff(d10,i,j)*pow(rgrid[i],4)*
     (((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i] + 
            Q11.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
            Q22.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2)))/
   (4.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))) - 
  (h22.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (h11.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (complex(0,2)*ht1.diff(d01,i,j)*w*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) + 
  (Dh.diff(d01,i,j)*pow(rgrid[i],3)*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
            h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/rgrid[i] + 
       (2*h.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
  (htt.diff(d02,i,j)*(2*pow(rgrid[i],2)*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j)))))/
   (2.*pow(L,2)) + (6*A0.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(-(L*(-(mu*a0(i,j)) - 
            a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (6*A0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
     (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(-12 + pow(mu,2))*(-1 + rgrid[i])*
     (1 + rgrid[i] + pow(rgrid[i],2))*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,4)*a0.diff(d01,i,j)*mu*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,4)*a0.diff(d01,i,j)*mu*w*Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (2*A0.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(1 - rgrid[i])*
     pow(rgrid[i],4)*((3*(1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (htt.diff(d10,i,j)*((-2*Qr1.diff(d01,i,j)*rgrid[i])/L + 
       (complex(0,48)*w*pow(rgrid[i],2))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (pow(rgrid[i],2)*((-2 + 
               Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
               rgrid[i]*Q11(i,j))/pow(rgrid[i],3) - 
            (Q11.diff(d01,i,j)*Qr1(i,j))/L))/(1 + rgrid[i]*Q11(i,j)) + 
       (pow(rgrid[i],2)*((-2 + 
               Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
               rgrid[i]*Q22(i,j))/pow(rgrid[i],3) - 
            (Q22.diff(d01,i,j)*Qr1(i,j))/L))/(1 + rgrid[i]*Q22(i,j)) + 
       (2*Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
        (L*(1 + rgrid[i]*Qrr(i,j))) - 
       (4*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) - 
       (4*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
        (L*(1 + rgrid[i]*Qtt(i,j))) + 
       (8*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
  (complex(0,1)*w*ht1(i,j)*pow(rgrid[i],8)*
     ((Q22.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],7) + ((1 + rgrid[i]*Q22(i,j))*
          ((2*Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j))*
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
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (h11.diff(d01,i,j)*pow(rgrid[i],8)*
     (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          pow(1 + rgrid[i]*Q11(i,j),2)*Qr1(i,j)*
          (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) - 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j)))/rgrid[i] + 
            Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        (2.*pow(rgrid[i],5)) + 
       (Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],7) + ((1 + rgrid[i]*Q11(i,j))*
          ((-3*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((2*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  pow(rgrid[i],3) + 
                 (L*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j)*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],4)) - 
                 (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i])))/
             pow(rgrid[i],2)))/pow(rgrid[i],2)))/
   (2.*pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  Dh(i,j)*(-(mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])) + 
     (h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j))/L + 
     (complex(0,6)*w*pow(rgrid[i],3)*
        (-(L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])) + 
          h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/
      (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
     (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
        ((2*pow(a0.diff(d01,i,j),2)*(2 + 3*pow(A1,2))*B1*
             pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + 
               B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
             (1 + rgrid[i]*Qrr(i,j)))/
           (pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (3*(1 + rgrid[i]*Q11(i,j))*
             ((2 + 3*pow(A1,2))*B1*
                exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + 
                  B1*mu*h(i,j)*rgrid[i])*
                pow(-(L*(-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
                  a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*
                   Qr1(i,j),2) + 
               (8*A1*pow(L,2)*
                  (-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)*rgrid[i]))*
                  (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                pow(rgrid[i],4)))/pow(rgrid[i],2)))/
      ((2 + 3*pow(A1,2))*pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qtt(i,j)))) + 
  (h22.diff(d10,i,j)*pow(rgrid[i],6)*
     (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
            Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        (2.*pow(rgrid[i],4)) + 
       ((1 + rgrid[i]*Q22(i,j))*
          (-(L*(-1 + rgrid[i])*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*
                (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                  rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],5)) + 
            rgrid[i]*Qr1(i,j)*
             (-((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
               (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))) + 
            (2*(1 + rgrid[i]*Q11(i,j))*
               ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
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
             pow(rgrid[i],2)))/pow(rgrid[i],2)))/
   (2.*L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (h22.diff(d01,i,j)*pow(rgrid[i],8)*
     ((Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
        pow(rgrid[i],7) - ((-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (((1 + rgrid[i]*Q22(i,j))*Qr1(i,j)*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 Q11.diff(d01,i,j)*Qr1(i,j)))/rgrid[i] + 
            (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],3)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))))*(1 + rgrid[i]*Qtt(i,j))\
)/(2.*pow(rgrid[i],4)) + (pow(1 + rgrid[i]*Q11(i,j),2)*Qr1(i,j)*
          ((L*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],5)) + 
            rgrid[i]*Qr1(i,j)*
             ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
               (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))) - 
            (2*(1 + rgrid[i]*Q22(i,j))*
               ((Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*rgrid[i]) + 
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
             pow(rgrid[i],2)))/pow(rgrid[i],3)))/
   (2.*pow(L,2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
     (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
  (htt.diff(d01,i,j)*(-4*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)) - 
       (2*Q11.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)) \
+ pow(rgrid[i],2)*pow(Qr1(i,j),2)*
        ((Q11.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q11(i,j)) + 
          (Q22.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) - 
          (2*Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) + 
          (4*Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))) + 
       (pow(rgrid[i],2)*((4*Qrr.diff(d01,i,j))/
             ((-1 + rgrid[i])*rgrid[i]*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*rgrid[i]*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))) \
+ (4*Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*rgrid[i]*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))))\
)/(1 + rgrid[i]*Q11(i,j)) + rgrid[i]*Qr1(i,j)*
        (6*Qr1.diff(d01,i,j)*rgrid[i] + 
          L*((complex(0,-48)*w*pow(rgrid[i],2))/
              ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
             (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q11(i,j))/(rgrid[i]*(1 + rgrid[i]*Q11(i,j))) - 
             (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q22(i,j))/(rgrid[i]*(1 + rgrid[i]*Q22(i,j))) + 
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
              ((-1 + rgrid[i])*rgrid[i]*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) \
- (8*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))\
))/(4.*pow(L,2)) + (htt(i,j)*
     ((complex(0,-72)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          pow(rgrid[i],4))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)) + 
       (complex(0,48)*w*rgrid[i])/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (2*((72*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
               48*exp(A1*mu*h(i,j)*rgrid[i]))/(2 + 3*pow(A1,2)) + 
            (Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*
               pow(rgrid[i],4))/
             (pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
               (1 + rgrid[i]*Q22(i,j))) + 
            (pow(rgrid[i],6)*
               (pow(Q22.diff(d01,i,j),2)/pow(rgrid[i],2) - 
                 (2*Q22.diff(d02,i,j)*(1 + rgrid[i]*Q22(i,j)))/
                  pow(rgrid[i],3) - 
                 (2*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                    pow(1 + rgrid[i]*Q22(i,j),2))/pow(rgrid[i],2)))/
             (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
               pow(1 + rgrid[i]*Q22(i,j),2)))*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       ((-1 + rgrid[i])*pow(rgrid[i],6)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Q11(i,j))*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*
               ((2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                  ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 (4*L*(((-1 + rgrid[i])*
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
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2))))/
             pow(rgrid[i],2) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((4*pow(Qrr.diff(d01,i,j),2))/
                  (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) + 
                 (2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                    pow(Qr1(i,j),2))/
                  ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (4*L*((L*(-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                      (2*Qr1.diff(d01,i,j)*
                       (1 + rgrid[i]*Q11(i,j)))/rgrid[i])*
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
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) + 
                 rgrid[i]*Qr1(i,j)*
                  ((4*Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                       (1 + rgrid[i]*Q11(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) - 
                    L*((2*Qrr.diff(d01,i,j)*
                        (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
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
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2))))))/
             pow(rgrid[i],2)))/
        (2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))) + 
       (complex(0,6)*w*pow(rgrid[i],2)*
          ((-2*Qr1.diff(d01,i,j)*rgrid[i])/L + 
            (pow(rgrid[i],2)*
               ((-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q11(i,j))/pow(rgrid[i],3) - 
                 (Q11.diff(d01,i,j)*Qr1(i,j))/L))/
             (1 + rgrid[i]*Q11(i,j)) + 
            (pow(rgrid[i],2)*
               ((-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q22(i,j))/pow(rgrid[i],3) - 
                 (Q22.diff(d01,i,j)*Qr1(i,j))/L))/
             (1 + rgrid[i]*Q22(i,j)) + 
            (2*Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
             (L*(1 + rgrid[i]*Qrr(i,j))) - 
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
             ((-1 + rgrid[i])*rgrid[i]*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) \
- (4*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
             (L*(1 + rgrid[i]*Qtt(i,j))) + 
            (8*(((-1 + rgrid[i])*
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
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))))\
)/((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (2*pow(rgrid[i],10)*(((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               pow(1 + rgrid[i]*Q22(i,j),2)*
               ((2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (pow(L,2)*
                    pow(-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j),2))/pow(rgrid[i],6) - 
                 (2*Q11.diff(d01,i,j)*L*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j))/pow(rgrid[i],3) + 
                 pow(Q11.diff(d01,i,j),2)*pow(Qr1(i,j),2))*
               (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],6)) - 
            ((-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
               (1 + rgrid[i]*Q22(i,j))*
               ((2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (pow(L,2)*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],6) + 
                 (4*Qrr.diff(d02,i,j)*(1 + rgrid[i]*Q22(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 (4*Qr1.diff(d01,i,j)*L*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j)))/
                  pow(rgrid[i],4) + 
                 (2*pow(L,2)*
                    (6 - 2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                      Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                      2*rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j)))/
                  pow(rgrid[i],6) - 
                 rgrid[i]*(((-6*Q11.diff(d01,i,j)*Qr1.diff(d0\
1,i,j) + (4*L*(-Q11.diff(d01,i,j) + 
                      Q11.diff(d11,i,j)*rgrid[i]))/
                       pow(rgrid[i],2))*(1 + rgrid[i]*Q22(i,j)))/
                     pow(rgrid[i],2) + 
                    L*((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                       (Q11.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],4)))*
                  Qr1(i,j) + 
                 pow(rgrid[i],2)*
                  ((Q11.diff(d01,i,j)*Q22.diff(d01,i,j))/
                     pow(rgrid[i],2) + 
                    (2*Q11.diff(d02,i,j)*(1 + rgrid[i]*Q22(i,j)))/
                     pow(rgrid[i],3))*pow(Qr1(i,j),2) - 
                 (2*Q11.diff(d01,i,j)*L*(1 + rgrid[i]*Q22(i,j))*
                    (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                  pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],6)) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               (((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    pow(-((L*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j),2)*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) - 
                 ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Q22(i,j))*
                    (2*rgrid[i]*
                       (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j) - 
                       (L*(-Q22.diff(d01,i,j) + 
                       Q22.diff(d11,i,j)*rgrid[i]))/
                       pow(rgrid[i],2))*Qr1(i,j) + 
                      Q22.diff(d02,i,j)*rgrid[i]*
                       pow(Qr1(i,j),2) + 
                      L*(-((Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],2)) + 
                        (L*(6 - 
                       2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q22.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q22(i,j)))/pow(rgrid[i],4) - 
                        (Q22.diff(d01,i,j)*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i]))*(1 + rgrid[i]*Qtt(i,j)))/
                  pow(rgrid[i],4) + 
                 (pow(1 + rgrid[i]*Q22(i,j),2)*
                    (4*pow(L,2)*exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) - 
                      ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (2*pow(Qr1.diff(d01,i,j),2)*
                       pow(rgrid[i],2) + 
                        pow(L,2)*pow(mu,2)*
                       pow(h(i,j) + h.diff(d10,i,j)*rgrid[i],2) \
- 2*L*(Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i]))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) + 
                      pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                       (4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-1 + rgrid[i],2) - 
                        pow(h.diff(d01,i,j),2)*pow(mu,2)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j))) - 
                      4*rgrid[i]*Qr1(i,j)*
                       (-2*a0.diff(d01,i,j)*L*mu*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        (-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                        (-1 + rgrid[i]) + 
                         ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d02,i,j)*rgrid[i] - 
                        h.diff(d01,i,j)*L*pow(mu,2)*rgrid[i]*
                        (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)))\
))/pow(rgrid[i],4)))/pow(rgrid[i],4)))/
        (pow(L,2)*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*pow(1 + rgrid[i]*Q11(i,j),2)*
          pow(1 + rgrid[i]*Q22(i,j),2)*(1 + rgrid[i]*Qtt(i,j)))))/4. + 
  (h11(i,j)*pow(rgrid[i],8)*((-2*(1 + rgrid[i]*Q22(i,j))*
          (-((Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*
                 (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
               pow(rgrid[i],6)) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((2*Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                    (1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],4) + 
                 ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    ((2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                      (pow(L,2)*
                       pow(-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j),2))/pow(rgrid[i],6))*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) - 
                 (Q11.diff(d01,i,j)*L*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j)*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],5) + 
                 (pow(Q11.diff(d01,i,j),2)*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
                    (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2))))/
             pow(rgrid[i],2)))/
        ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((2*pow(Q22.diff(d01,i,j),2)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) - 
            ((1 + rgrid[i]*Q22(i,j))*
               ((2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (4*Q22.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/
             pow(rgrid[i],2) + 
            (pow(1 + rgrid[i]*Q22(i,j),2)*
               ((4*Qrr.diff(d02,i,j))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 (4*Qr1.diff(d01,i,j)*L*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) - 
                 (complex(0,6)*pow(L,2)*w*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/
                  ((-12 + pow(mu,2))*rgrid[i]*
                    (-1 + pow(rgrid[i],3))) + 
                 (2*pow(L,2)*
                    (6 - 2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                      Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                      2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) - 
                 (2*Q11.diff(d01,i,j)*L*
                    (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/rgrid[i] \
- ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    ((4*pow(Qrr.diff(d01,i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                       pow(Qr1(i,j),2))/
                       ((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (4*pow(L,2)*
                       (-2 + 
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
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) - 
                      L*rgrid[i]*Qr1(i,j)*
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
                        pow(mu,2)*pow(rgrid[i],3),2)))))/
                  (2.*(1 + rgrid[i]*Qrr(i,j))) + 
                 (4*Qrr.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                  ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j))) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  ((2*Q11.diff(d02,i,j))/rgrid[i] + 
                    (Q11.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                     (1 + rgrid[i]*Qtt(i,j))) + 
                 (2*pow(L,2)*
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
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                  ((-1 + rgrid[i])*pow(rgrid[i],4)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (1 + rgrid[i]*Qtt(i,j))) - 
                 (16*pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j))*
                    ((pow(Qtt.diff(d01,i,j),2)*
                       pow(-1 + rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2))/
                       (4.*pow(rgrid[i],2)) + 
                      ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (3*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2) - 
                       (Qtt.diff(d02,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/rgrid[i])*
                       (1 + rgrid[i]*Qtt(i,j)))/
                       (2.*pow(rgrid[i],2)) - 
                      (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                        pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)*
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                       (4.*pow(rgrid[i],2))))/
                  (pow(-1 + rgrid[i],3)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),3)*
                    pow(1 + rgrid[i]*Qtt(i,j),2)) + 
                 rgrid[i]*Qr1(i,j)*
                  ((-4*L*(-Q11.diff(d01,i,j) + 
                        Q11.diff(d11,i,j)*rgrid[i]))/
                     pow(rgrid[i],2) - 
                    (Qtt.diff(d01,i,j)*L*
                       (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/
                     (pow(rgrid[i],2)*(1 + rgrid[i]*Qtt(i,j))) + 
                    (Q11.diff(d01,i,j)*
                       (6*Qr1.diff(d01,i,j)*rgrid[i] + 
                        (complex(0,6)*L*w*pow(rgrid[i],2))/
                        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) \
- (2*L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))))/rgrid[i])))/
             pow(rgrid[i],4)))/pow(rgrid[i],2) + 
       (pow(1 + rgrid[i]*Q11(i,j),2)*
          (pow(-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                     rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
              Q22.diff(d01,i,j)*Qr1(i,j),2) + 
            (2*pow(1 + rgrid[i]*Q22(i,j),2)*
               (2*pow(Qr1.diff(d01,i,j),2)*pow(rgrid[i],2) - 
                 2*L*(Qr1.diff(d01,i,j) + 
                    Qr1.diff(d11,i,j)*rgrid[i]) + 
                 2*Qr1.diff(d02,i,j)*pow(rgrid[i],2)*Qr1(i,j) + 
                 Qr1.diff(d01,i,j)*rgrid[i]*
                  (rgrid[i]*Qr1(i,j)*
                     (-((Qrr.diff(d01,i,j)*rgrid[i])/
                        (1 + rgrid[i]*Qrr(i,j))) + 
                       (Qtt.diff(d01,i,j)*rgrid[i])/
                        (1 + rgrid[i]*Qtt(i,j))) + 
                    L*((complex(0,6)*w*pow(rgrid[i],2))/
                        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) \
+ (2*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qrr(i,j))) - 
                       (2*(((-1 + rgrid[i])*
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
                        ((-1 + rgrid[i])*rgrid[i]*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                         (1 + rgrid[i]*Qtt(i,j)))))))/pow(rgrid[i],4) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  ((-2*Q22.diff(d02,i,j))/rgrid[i] + 
                    (Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                     (1 + rgrid[i]*Qrr(i,j)) - 
                    (Q22.diff(d01,i,j)*Qtt.diff(d01,i,j))/
                     (1 + rgrid[i]*Qtt(i,j))) + 
                 L*((-2*L*(6 - 
                        2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
                        Q22.diff(d20,i,j)*pow(rgrid[i],3) + 
                        2*rgrid[i]*Q22(i,j)))/pow(rgrid[i],4) + 
                    (2*Q22.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                     rgrid[i] + 
                    (L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*
                       ((complex(0,6)*w*pow(rgrid[i],2))/
                        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) \
+ (2*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qrr(i,j))) - 
                        (2*(((-1 + rgrid[i])*
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
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))))/pow(rgrid[i],3)) + 
                 rgrid[i]*Qr1(i,j)*
                  (L*((4*(-Q22.diff(d01,i,j) + 
                        Q22.diff(d11,i,j)*rgrid[i]))/
                        pow(rgrid[i],2) - 
                       (Qrr.diff(d01,i,j)*
                        (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/
                        (pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j))) + 
                       (Qtt.diff(d01,i,j)*
                        (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/
                        (pow(rgrid[i],2)*(1 + rgrid[i]*Qtt(i,j)))) + 
                    (Q22.diff(d01,i,j)*
                       (-2*Qr1.diff(d01,i,j)*rgrid[i] - 
                         (complex(0,6)*L*w*pow(rgrid[i],2))/
                        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) \
- (2*L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qrr(i,j)))/2. + 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qrr(i,j)))/2.))/
                        ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qrr(i,j))) + 
                         (2*L*
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
                         ((-1 + rgrid[i])*rgrid[i]*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                         (1 + rgrid[i]*Qtt(i,j)))))/rgrid[i])))/
             pow(rgrid[i],2)))/pow(rgrid[i],4)))/
   (4.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (h22(i,j)*((pow(rgrid[i],2)*
          (-((-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q11(i,j))/pow(rgrid[i],3)) + 
            (Q11.diff(d01,i,j)*Qr1(i,j))/L))/(1 + rgrid[i]*Q11(i,j)) + 
       (pow(rgrid[i],2)*(-((-2 + 
                 Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j))/pow(rgrid[i],3)) + 
            (Q22.diff(d01,i,j)*Qr1(i,j))/L))/(1 + rgrid[i]*Q22(i,j)) + 
       (4*pow(rgrid[i],2)*(-(Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j))/2. + 
            (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*rgrid[i]) + (L*
               (((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Qtt(i,j)))/2. - 
                 ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                      3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))*
                    (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3)))/
        (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
       (complex(0,12)*w*pow(rgrid[i],8)*
          (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
               (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],4)) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (-(L*(-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                     (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qtt(i,j)))/
                  (2.*pow(rgrid[i],5)) + 
                 rgrid[i]*Qr1(i,j)*
                  (-((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                    (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))) \
+ (2*(1 + rgrid[i]*Q11(i,j))*
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
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],3))\
)/pow(rgrid[i],2)))/pow(rgrid[i],2)))/
        (L*(-12 + pow(mu,2))*(-1 + rgrid[i])*
          (1 + rgrid[i] + pow(rgrid[i],2))*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
          (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
       (1 - rgrid[i])*((2*(1 + rgrid[i]*Qrr(i,j))*
             ((48*(3*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)) + 
                    2*exp(A1*mu*h(i,j)*rgrid[i])))/(2 + 3*pow(A1,2)) + 
               (Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*
                  pow(rgrid[i],4))/
                (pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
                  (1 + rgrid[i]*Q22(i,j))) + 
               (pow(rgrid[i],6)*
                  (pow(Q22.diff(d01,i,j),2)/pow(rgrid[i],2) - 
                    (2*Q22.diff(d02,i,j)*(1 + rgrid[i]*Q22(i,j)))/
                     pow(rgrid[i],3) - 
                    (2*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                       pow(1 + rgrid[i]*Q22(i,j),2))/pow(rgrid[i],2)\
))/(pow(L,2)*(1 + rgrid[i]*Q11(i,j))*pow(1 + rgrid[i]*Q22(i,j),2)) + 
               (2*pow(rgrid[i],2)*
                  (-4*pow(w,2) - 
                    (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2)*pow(rgrid[i],2))/
                     (pow(L,2)*(1 + rgrid[i]*Q11(i,j)))))/
                ((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))\
)))/((-1 + rgrid[i])*pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          (pow(rgrid[i],8)*(((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Q11(i,j))*
                  (-((L*(-2 + 
                        Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                    Q22.diff(d01,i,j)*Qr1(i,j))*
                  ((2*Qrr.diff(d01,i,j)*Qr1(i,j))/
                     ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) - 
                    (4*L*(((-1 + rgrid[i])*
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
                       pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)))*
                  (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],4)) + 
               ((1 + rgrid[i]*Q22(i,j))*
                  (((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       ((4*pow(Qrr.diff(d01,i,j),2))/
                        (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) + 
                        (4*L*
                        ((L*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) - 
                       (2*Qr1.diff(d01,i,j)*
                       (1 + rgrid[i]*Q11(i,j)))/rgrid[i])*
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
                       (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
                    (2*Qrr.diff(d01,i,j)*rgrid[i]*pow(Qr1(i,j),2)*
                       ((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) + 
                        (Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],3))\
))/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                    (8*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
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
                     (pow(-1 + rgrid[i],2)*pow(rgrid[i],8)*
                       pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)) + 
                    rgrid[i]*Qr1(i,j)*
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
+ (2*(1 + rgrid[i]*Q11(i,j))*
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
                        pow(rgrid[i],2))))/pow(rgrid[i],2)))/
           (pow(L,2)*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
             (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) + 
          (4*pow(rgrid[i],12)*
             ((pow(-1 + rgrid[i],2)*
                  pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3),2)*
                  pow(1 + rgrid[i]*Q22(i,j),2)*
                  ((2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                     ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                    (pow(L,2)*
                       pow(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j),2))/pow(rgrid[i],6) - 
                    (2*Q11.diff(d01,i,j)*L*
                       (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*Qr1(i,j))/pow(rgrid[i],3) \
+ pow(Q11.diff(d01,i,j),2)*pow(Qr1(i,j),2))*
                  pow(1 + rgrid[i]*Qtt(i,j),2))/(4.*pow(rgrid[i],8)) - 
               ((-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))*
                  (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
                  (1 + rgrid[i]*Qtt(i,j))*
                  ((Q22.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                       (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],4) + 
                    ((1 + rgrid[i]*Q22(i,j))*
                       (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        ((2*Qrr.diff(d02,i,j))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        L*((-2*Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) + 
                        (L*
                       (6 - 
                       2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) - 
                        (Q11.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i]))*(1 + rgrid[i]*Qtt(i,j)))/
                        pow(rgrid[i],2) + 
                         pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                        ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/
                        (2.*pow(rgrid[i],2)) + 
                        (Q11.diff(d02,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3)) + 
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
- rgrid[i]*Qr1(i,j)*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-6*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j\
) + (4*L*(-Q11.diff(d01,i,j) + Q11.diff(d11,i,j)*rgrid[i]))/
                        pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
                         L*((Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
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
)))/pow(rgrid[i],2)))/(2.*pow(rgrid[i],6)) + 
               (pow(1 + rgrid[i]*Q11(i,j),2)*
                  ((pow(-1 + rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)*
                       pow(-((L*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                        Q22.diff(d01,i,j)*Qr1(i,j),2)*
                       pow(1 + rgrid[i]*Qtt(i,j),2))/
                     (4.*pow(rgrid[i],4)) + 
                    ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                       (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qtt(i,j))*
                       (-(pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                        ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/
                        (2.*pow(rgrid[i],2)) + 
                        (Q22.diff(d02,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3))) - 
                         L*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                        ((2*L*
                       (6 - 
                       2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q22.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q22(i,j)))/pow(rgrid[i],4) - 
                       (2*Q22.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                       rgrid[i])*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
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
) + rgrid[i]*Qr1(i,j)*(((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-2*Q22.diff(d01,i,j)*Qr1.diff(d01,i,j\
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
)))/(2.*pow(rgrid[i],4)) + (2*pow(1 + rgrid[i]*Q22(i,j),2)*
                       (-(pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)*
                        (pow(Qr1.diff(d01,i,j),2)*
                        pow(rgrid[i],2) + 
                        pow(L,2)*pow(mu,2)*
                        pow(h(i,j) + h.diff(d10,i,j)*rgrid[i],2) \
- L*(Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i]))*
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                         (2.*pow(rgrid[i],4)) + 
                         (pow(L,2)*
                        pow(((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))*
                       (-2 + 
                       Qtt.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Qtt(i,j)))/2. - 
                        ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                       3*(-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.,2))/
                         pow(rgrid[i],6) + 
                         pow(rgrid[i],2)*pow(Qr1(i,j),2)*
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
                        (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],2) - 
                         (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                        pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                        4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)*
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                         (2.*pow(rgrid[i],2))) - 
                         rgrid[i]*Qr1(i,j)*
                         (((-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (-4*a0.diff(d01,i,j)*L*mu*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       (-(mu*a0(i,j)) - 
                       a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                       (-1 + rgrid[i]) + 
                        (3*Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))/2. - 
                        (4*L*
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
                       pow(mu,2)*pow(rgrid[i],3)))/2.)))/
                        pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
                        (2.*pow(rgrid[i],2)) + 
                         (pow(-1 + rgrid[i],2)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)*
                        (Qr1.diff(d02,i,j)*rgrid[i] - 
                        2*h.diff(d01,i,j)*L*pow(mu,2)*rgrid[i]*
                        (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                        pow(1 + rgrid[i]*Qtt(i,j),2))/
                        (2.*pow(rgrid[i],4)) + 
                         (Qtt.diff(d01,i,j)*L*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
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
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/pow(rgrid[i],4)) \
+ (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                         (1 + rgrid[i]*Qtt(i,j))*
                         (2*L*exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) + 
                         (Qtt.diff(d01,i,j)*(-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i] + 
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
- (2*L*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                         pow(rgrid[i],4)))/(2.*pow(rgrid[i],2))))/
                     pow(rgrid[i],4)))/pow(rgrid[i],4)))/
           (pow(L,2)*pow(-1 + rgrid[i],2)*
             pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3),2)*
             pow(1 + rgrid[i]*Q11(i,j),2)*pow(1 + rgrid[i]*Q22(i,j),2)*
             pow(1 + rgrid[i]*Qtt(i,j),2)))))/4.
;
elem[5]=
ht1.diff(d20,i,j) - (2*ht1.diff(d11,i,j)*rgrid[i]*Qr1(i,j))/L + 
  (ht1.diff(d02,i,j)*pow(rgrid[i],2)*pow(Qr1(i,j),2))/pow(L,2) + 
  (4*Ax.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],3)*
     Qr1(i,j)*(-(L*(-(mu*a0(i,j)) - 
            a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*Ax.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*muj*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (complex(0,6)*L*w*pow(rgrid[i],2) - 
       Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*rgrid[i]*
        (-1 + pow(rgrid[i],3)))*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(-1 + rgrid[i])*(-1 + pow(rgrid[i],3))*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (4*Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (complex(0,6)*L*w*pow(rgrid[i],2) - 
       Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*rgrid[i]*
        (-1 + pow(rgrid[i],3)))*
     (L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i])) + 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(-1 + rgrid[i])*(-1 + pow(rgrid[i],3))*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,4)*h11.diff(d01,i,j)*w*(1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(-1 + rgrid[i],2)*pow(-4 - 4*rgrid[i] - 
       4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3),2)*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,4)*h22.diff(d01,i,j)*w*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,8)*h.diff(d01,i,j)*mu*w*Dh(i,j)*pow(rgrid[i],2)*
     (1 + rgrid[i]*Qrr(i,j)))/
   (L*pow(-1 + rgrid[i],2)*pow(-4 - 4*rgrid[i] - 
       4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3),2)*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (complex(0,2)*w*h11(i,j)*(1 + rgrid[i]*Qrr(i,j))*
     ((2*Q22.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) + 
       (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) - 
       (Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))))/
   (L*pow(-1 + rgrid[i],2)*pow(-4 - 4*rgrid[i] - 
       4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3),2)*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,2)*w*h22(i,j)*pow(rgrid[i],4)*
     ((Qtt.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/pow(rgrid[i],3) + 
       (Qrr.diff(d01,i,j)*(1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],3)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*pow(1 + rgrid[i]*Qtt(i,j),2)) + 
  (ht1.diff(d10,i,j)*((-2*Qr1.diff(d01,i,j)*rgrid[i])/L + 
       (complex(0,24)*w*pow(rgrid[i],2))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (pow(rgrid[i],2)*(-((-2 + 
                 Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q11(i,j))/pow(rgrid[i],3)) + 
            (Q11.diff(d01,i,j)*Qr1(i,j))/L))/(1 + rgrid[i]*Q11(i,j)) + 
       (pow(rgrid[i],2)*((-2 + 
               Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
               rgrid[i]*Q22(i,j))/pow(rgrid[i],3) - 
            (Q22.diff(d01,i,j)*Qr1(i,j))/L))/(1 + rgrid[i]*Q22(i,j)) + 
       (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
        (L*(1 + rgrid[i]*Qrr(i,j))) - 
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
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) - 
       (3*Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))/
        (L*(1 + rgrid[i]*Qtt(i,j))) + 
       (6*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))/2. + 
  (ht1.diff(d01,i,j)*(-2*L*(Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)) + 
       pow(rgrid[i],2)*pow(Qr1(i,j),2)*
        (-((Q11.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q11(i,j))) + 
          (Q22.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Q22(i,j)) - 
          (Qrr.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qrr(i,j)) + 
          (3*Qtt.diff(d01,i,j)*rgrid[i])/(1 + rgrid[i]*Qtt(i,j))) + 
       rgrid[i]*Qr1(i,j)*(4*Qr1.diff(d01,i,j)*rgrid[i] + 
          L*((complex(0,-24)*w*pow(rgrid[i],2))/
              ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
             (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q11(i,j))/(rgrid[i]*(1 + rgrid[i]*Q11(i,j))) - 
             (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                rgrid[i]*Q22(i,j))/(rgrid[i]*(1 + rgrid[i]*Q22(i,j))) + 
             (2*(((-1 + rgrid[i])*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))) \
- (6*(((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
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
                  pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j)))))\
))/(2.*pow(L,2)) + (ht1(i,j)*
     ((complex(0,-36)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          pow(rgrid[i],4))/
        (pow(-12 + pow(mu,2),2)*pow(-1 + pow(rgrid[i],3),2)) + 
       (complex(0,24)*w*rgrid[i])/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (complex(0,6)*w*pow(rgrid[i],10)*
          (((1 + rgrid[i]*Q11(i,j))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
            ((1 + rgrid[i]*Q22(i,j))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                      (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
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
                         pow(mu,2)*pow(rgrid[i],3)))) - 
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
        (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))*
          (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))) - 
       (4*(1 + rgrid[i]*Qrr(i,j))*
          (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               ((-72*pow(A1,2)*L*
                    exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
                  (2 + 3*pow(A1,2)) - 
                 (48*L*exp(A1*mu*h(i,j)*rgrid[i]))/
                  (2 + 3*pow(A1,2)) + 
                 (2*pow(Qr1.diff(d01,i,j),2)*(-1 + rgrid[i])*
                    pow(rgrid[i],4)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))/
                  (L*(1 + rgrid[i]*Qrr(i,j))) + 
                 (L*pow(mu,2)*(-1 + rgrid[i])*pow(rgrid[i],2)*
                    pow(h(i,j) + h.diff(d10,i,j)*rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))/
                  (1 + rgrid[i]*Qrr(i,j)) - 
                 ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (Qr1.diff(d01,i,j) + 
                      Qr1.diff(d11,i,j)*rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))/
                  (1 + rgrid[i]*Qrr(i,j)) - 
                 (L*(-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    pow(-2 + 
                      Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j),2))/
                  (2.*pow(1 + rgrid[i]*Q22(i,j),2)*
                    (1 + rgrid[i]*Qrr(i,j))) - 
                 (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*
                    pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/
                  (2.*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) \
+ (L*(-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (6 - 2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
                      Q22.diff(d20,i,j)*pow(rgrid[i],3) + 
                      2*rgrid[i]*Q22(i,j)))/
                  ((1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
                 (Q22.diff(d01,i,j)*(-1 + rgrid[i])*
                    pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                  ((1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
                 ((-1 + rgrid[i])*pow(rgrid[i],8)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))*
                    (((1 + rgrid[i]*Q22(i,j))*
                       ((2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (pow(L,2)*
                       pow(-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j),2))/pow(rgrid[i],6) - 
                       (2*Q11.diff(d01,i,j)*L*
                       (-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j))/
                       pow(rgrid[i],3) + 
                       pow(Q11.diff(d01,i,j),2)*
                       pow(Qr1(i,j),2)))/pow(rgrid[i],2) + 
                      (2*Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*
                        (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                  (2.*L*pow(1 + rgrid[i]*Q11(i,j),2)*
                    (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
                 (Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
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
                  pow(1 + rgrid[i]*Qrr(i,j),2) - 
                 (L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
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
                        (1 + rgrid[i]*Qrr(i,j)))/2.))/
                  ((1 + rgrid[i]*Q22(i,j))*
                    pow(1 + rgrid[i]*Qrr(i,j),2)) + 
                 (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)*
                    pow(Qr1(i,j),2)*
                    ((-2*pow(Q22.diff(d01,i,j),2)*
                       (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                       pow(1 + rgrid[i]*Q22(i,j),2)*
                       (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      ((1 + rgrid[i]*Q22(i,j))*
                        ((-2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (4*Q22.diff(d02,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],3)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                       pow(rgrid[i],2)))/
                  (4.*L*pow(1 + rgrid[i]*Q22(i,j),2)*
                    pow(1 + rgrid[i]*Qrr(i,j),2)) + 
                 (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)*
                    ((-4*pow(Q22.diff(d01,i,j),2)*
                       pow(1 + rgrid[i]*Qrr(i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
                        pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) + 
                      (2*(1 + rgrid[i]*Q22(i,j))*
                       (1 + rgrid[i]*Qrr(i,j))*
                       ((2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j\
))/((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                       (pow(L,2)*
                       (-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],6) - 
                       L*rgrid[i]*
                       ((Q22.diff(d01,i,j)*
                      (-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                       (Q11.diff(d01,i,j)*
                      (-2 + 
                      Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*
                       Qr1(i,j) + 
                       Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*
                       pow(Qr1(i,j),2) + 
                       (4*Q22.diff(d02,i,j)*
                       (1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],3)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (pow(1 + rgrid[i]*Q22(i,j),2)*
                        ((-4*pow(Qrr.diff(d01,i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) - 
                        (2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                       pow(Qr1(i,j),2))/
                       ((-1 + rgrid[i])*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (2*
                       ((4*Qrr.diff(d02,i,j))/
                      ((-1 + rgrid[i])*rgrid[i]*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                       rgrid[i]*
                      (7*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j\
) - (4*L*(-Q11.diff(d01,i,j) + 
                     Q11.diff(d11,i,j)*rgrid[i]))/
                      pow(rgrid[i],2))*Qr1(i,j) + 
                       2*Q11.diff(d02,i,j)*rgrid[i]*
                      pow(Qr1(i,j),2) + 
                       L*((-5*Qr1.diff(d01,i,j)*
                      (-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) + 
                       (2*L*
                      (6 - 
                      2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                      Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                      2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) - 
                       (2*Q11.diff(d01,i,j)*
                      (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                       rgrid[i]))*(1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                        (8*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                       pow(1 + rgrid[i]*Qrr(i,j),2))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) - 
                        (4*pow(L,2)*
                       (-2 + 
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
                       (1 + rgrid[i]*Qrr(i,j)))/2.))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
                       pow(-4 - 4*rgrid[i] - 
                       4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3),2)) + 
                        L*rgrid[i]*Qr1(i,j)*
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
                       pow(mu,2)*pow(rgrid[i],3),2)))))/
                       pow(rgrid[i],4)))/
                  (4.*L*(1 + rgrid[i]*Q11(i,j))*
                    pow(1 + rgrid[i]*Q22(i,j),2)*
                    pow(1 + rgrid[i]*Qrr(i,j),2)) + 
                 (pow(-1 + rgrid[i],2)*pow(rgrid[i],9)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*
                    ((4*Q22.diff(d01,i,j)*L*
                        (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],6)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) - 
                      (pow(1 + rgrid[i]*Q22(i,j),2)*
                        ((2*Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j\
))/((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) - 
                        (4*(Qr1.diff(d02,i,j)*rgrid[i] - 
                       2*h.diff(d01,i,j)*L*pow(mu,2)*rgrid[i]*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                       (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))))/
                       pow(rgrid[i],4) + 
                      ((1 + rgrid[i]*Q22(i,j))*
                        ((2*
                       (3*Q22.diff(d01,i,j)*Qr1.diff(d01,i,j\
) - (4*L*(-Q22.diff(d01,i,j) + Q22.diff(d11,i,j)*rgrid[i]))/
                       pow(rgrid[i],2))*(1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        L*((2*Qrr.diff(d01,i,j)*
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
                        pow(mu,2)*pow(rgrid[i],3),2)))))/
                       pow(rgrid[i],2)))/
                  (4.*L*pow(1 + rgrid[i]*Q22(i,j),2)*
                    pow(1 + rgrid[i]*Qrr(i,j),2)))*
               (1 + rgrid[i]*Qtt(i,j)))/(2.*pow(rgrid[i],2)) + 
            ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               ((4*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                    exp(B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
                    (1 + rgrid[i]*Qrr(i,j)))/
                  (pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 ((1 + rgrid[i]*Q11(i,j))*
                    (rgrid[i]*
                       (-4*a0.diff(d01,i,j)*L*mu*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        (-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))*
                        (-1 + rgrid[i]) - 
                        (3*Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))/2.)*Qr1(i,j) - 
                      2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                       exp(B1*mu*h(i,j)*rgrid[i])*
                       pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                       pow(Qr1(i,j),2) + 
                      L*(-2*L*exp(B1*mu*h(i,j)*rgrid[i])*
                         pow(-(mu*a0(i,j)) - 
                        a0.diff(d10,i,j)*mu*(-1 + rgrid[i]),2) + 
                         (3*Qr1.diff(d01,i,j)*
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
)/pow(rgrid[i],2)))/
             (2.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j)))))/
        (L*pow(-1 + rgrid[i],2)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Qtt(i,j)))))/2.
;
elem[6]=
h22.diff(d20,i,j)*(1 - rgrid[i]) + 
  (2*h22.diff(d11,i,j)*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j))/L + 
  Dh.diff(d10,i,j)*rgrid[i]*(mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
     (h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j))/L) - 
  (h11.diff(d10,i,j)*pow(rgrid[i],4)*
     (((1 + rgrid[i]*Q22(i,j))*
          (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                   rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/rgrid[i] + 
            Q11.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2) + 
       ((1 + rgrid[i]*Q11(i,j))*
          ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
            Q22.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2)))/
   (4.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))) - 
  (h11.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))) - 
  (h22.diff(d02,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (2*(1 + rgrid[i]*Q11(i,j))*pow(Qr1(i,j),2) + 
       (2*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (2.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
  (Dh.diff(d01,i,j)*pow(rgrid[i],3)*
     (((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i])) + 
            h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/rgrid[i] + 
       (2*h.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*rgrid[i]*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (pow(L,2)*(1 + rgrid[i]*Q11(i,j))) + 
  (h11.diff(d01,i,j)*pow(rgrid[i],6)*
     ((pow(1 + rgrid[i]*Q11(i,j),2)*Qr1(i,j)*
          ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) + 
            (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j)))/rgrid[i] - 
            Q22.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],3) + 
       (2*Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j)))/
        ((-1 + rgrid[i])*pow(rgrid[i],5)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))) + 
       ((1 + rgrid[i]*Q11(i,j))*
          (((1 + rgrid[i]*Q22(i,j))*
               ((-4*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*rgrid[i]*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 (L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j))/pow(rgrid[i],2) + 
                 Q11.diff(d01,i,j)*rgrid[i]*pow(Qr1(i,j),2)))/
             pow(rgrid[i],2) - 
            (6*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],3)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)))/
   (4.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))) - 
  (htt(i,j)*(1 + rgrid[i]*Qrr(i,j))*
     ((-72*pow(A1,2)*exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1)))/
        (2 + 3*pow(A1,2)) - (48*exp(A1*mu*h(i,j)*rgrid[i]))/
        (2 + 3*pow(A1,2)) + (2*pow(Qr1.diff(d01,i,j),2)*
          (-1 + rgrid[i])*pow(rgrid[i],4)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))/
        (pow(L,2)*(1 + rgrid[i]*Qrr(i,j))) + 
       (pow(mu,2)*(-1 + rgrid[i])*pow(rgrid[i],2)*
          pow(h(i,j) + h.diff(d10,i,j)*rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))/(1 + rgrid[i]*Qrr(i,j)) - 
       (2*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j)*rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))/(L*(1 + rgrid[i]*Qrr(i,j))) - 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          pow(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
            rgrid[i]*Q22(i,j),2))/
        (2.*pow(1 + rgrid[i]*Q22(i,j),2)*(1 + rgrid[i]*Qrr(i,j))) - 
       (Qr1.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
            rgrid[i]*Q22(i,j)))/
        (L*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (6 - 2*Q22.diff(d10,i,j)*pow(rgrid[i],2) + 
            Q22.diff(d20,i,j)*pow(rgrid[i],3) + 2*rgrid[i]*Q22(i,j)))/
        ((1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
       (Q22.diff(d01,i,j)*(-1 + rgrid[i])*pow(rgrid[i],3)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
        (L*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) - 
       ((-1 + rgrid[i])*pow(rgrid[i],8)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          (((1 + rgrid[i]*Q22(i,j))*
               ((2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (pow(L,2)*
                    pow(-2 + 
                      Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j),2))/pow(rgrid[i],6) - 
                 (2*Q11.diff(d01,i,j)*L*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*Qr1(i,j))/pow(rgrid[i],3) + 
                 pow(Q11.diff(d01,i,j),2)*pow(Qr1(i,j),2)))/
             pow(rgrid[i],2) + 
            (2*Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*
               (1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3)))))/
        (2.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
          (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))) + 
       (2*Qr1.diff(d01,i,j)*pow(rgrid[i],2)*
          (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Qrr(i,j)))/2. + 
            ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                 3*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))*
               (1 + rgrid[i]*Qrr(i,j)))/2.))/
        (L*pow(1 + rgrid[i]*Qrr(i,j),2)) - 
       ((-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - rgrid[i]*Q22(i,j))*
          (((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               (-2 + Qrr.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Qrr(i,j)))/2. + 
            ((12 - pow(mu,2)*pow(rgrid[i],4) - 
                 3*(-1 + rgrid[i])*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3)))*
               (1 + rgrid[i]*Qrr(i,j)))/2.))/
        ((1 + rgrid[i]*Q22(i,j))*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*pow(Qr1(i,j),2)*
          ((-2*pow(Q22.diff(d01,i,j),2)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (4*pow(h.diff(d01,i,j),2)*pow(mu,2)*
               pow(1 + rgrid[i]*Q22(i,j),2)*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            ((1 + rgrid[i]*Q22(i,j))*
               ((-2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (4*Q22.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],2)\
))/(4.*pow(L,2)*pow(1 + rgrid[i]*Q22(i,j),2)*
          pow(1 + rgrid[i]*Qrr(i,j),2)) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],10)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*
          ((-4*pow(Q22.diff(d01,i,j),2)*
               pow(1 + rgrid[i]*Qrr(i,j),2))/
             (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)) + 
            (2*(1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (pow(L,2)*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q11(i,j))*
                    (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],6) - 
                 L*rgrid[i]*((Q22.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) + 
                    (Q11.diff(d01,i,j)*
                       (-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],4))*Qr1(i,j) \
+ Q11.diff(d01,i,j)*Q22.diff(d01,i,j)*pow(Qr1(i,j),2) + 
                 (4*Q22.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q22(i,j),2)*
               ((-4*pow(Qrr.diff(d01,i,j),2))/
                  (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) - 
                 (2*Q11.diff(d01,i,j)*Qrr.diff(d01,i,j)*
                    pow(Qr1(i,j),2))/
                  ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (4*((2*Qrr.diff(d02,i,j))/
                       ((-1 + rgrid[i])*rgrid[i]*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                      rgrid[i]*
                       (3*Q11.diff(d01,i,j)*Qr1.diff(d01,i,j\
) - (2*L*(-Q11.diff(d01,i,j) + Q11.diff(d11,i,j)*rgrid[i]))/
                       pow(rgrid[i],2))*Qr1(i,j) + 
                      Q11.diff(d02,i,j)*rgrid[i]*
                       pow(Qr1(i,j),2) + 
                      L*((-2*Qr1.diff(d01,i,j)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],2) + 
                        (L*
                       (6 - 
                       2*Q11.diff(d10,i,j)*pow(rgrid[i],2) + 
                       Q11.diff(d20,i,j)*pow(rgrid[i],3) + 
                       2*rgrid[i]*Q11(i,j)))/pow(rgrid[i],4) - 
                        (Q11.diff(d01,i,j)*
                       (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j)))/
                        rgrid[i]))*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (8*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                    pow(1 + rgrid[i]*Qrr(i,j),2))/
                  (pow(-1 + rgrid[i],2)*pow(rgrid[i],2)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) - 
                 (4*pow(L,2)*
                    (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
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
                        (1 + rgrid[i]*Qrr(i,j)))/2.))/
                  (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
                    pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3),2)) + 
                 L*rgrid[i]*Qr1(i,j)*
                  ((2*Qrr.diff(d01,i,j)*
                       (-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
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
                       pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2)))))/
             pow(rgrid[i],4)))/
        (4.*pow(L,2)*(1 + rgrid[i]*Q11(i,j))*
          pow(1 + rgrid[i]*Q22(i,j),2)*pow(1 + rgrid[i]*Qrr(i,j),2)) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],9)*
          pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3),2)*Qr1(i,j)*
          ((4*Q22.diff(d01,i,j)*L*
               (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                 rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
             ((-1 + rgrid[i])*pow(rgrid[i],6)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) - 
            (2*pow(1 + rgrid[i]*Q22(i,j),2)*
               ((2*Qr1.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                  ((-1 + rgrid[i])*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) - 
                 (4*(Qr1.diff(d02,i,j)*rgrid[i] - 
                      h.diff(d01,i,j)*L*pow(mu,2)*rgrid[i]*
                       (h(i,j) + h.diff(d10,i,j)*rgrid[i]))*
                    (1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/pow(rgrid[i],4) \
+ ((1 + rgrid[i]*Q22(i,j))*((8*
                    (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j) - 
                      (L*(-Q22.diff(d01,i,j) + 
                        Q22.diff(d11,i,j)*rgrid[i]))/
                       pow(rgrid[i],2))*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],2)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 L*((2*Qrr.diff(d01,i,j)*
                       (-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],4)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                    (4*Q22.diff(d01,i,j)*
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
                     (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                       pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3),2)))))/
             pow(rgrid[i],2)))/
        (4.*pow(L,2)*pow(1 + rgrid[i]*Q22(i,j),2)*
          pow(1 + rgrid[i]*Qrr(i,j),2))))/
   (2.*(-1 + rgrid[i])*pow(rgrid[i],2)*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + pow(mu,2)*pow(rgrid[i],3))\
) + h22.diff(d10,i,j)*(-2 - (complex(0,12)*w*pow(rgrid[i],2))/
      ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))) + 
     (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*
        ((6*(1 + rgrid[i]*Q11(i,j))*
             (-((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
               Q22.diff(d01,i,j)*Qr1(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
           ((-1 + rgrid[i])*pow(rgrid[i],4)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) + 
          ((1 + rgrid[i]*Q22(i,j))*
             ((-6*L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                    rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
                ((-1 + rgrid[i])*pow(rgrid[i],5)*
                  (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                    pow(mu,2)*pow(rgrid[i],3))) + 
               rgrid[i]*Qr1(i,j)*
                ((-4*Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/
                   ((-1 + rgrid[i])*pow(rgrid[i],3)*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3))) + 
                  (6*Q11.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                   ((-1 + rgrid[i])*pow(rgrid[i],3)*
                     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                       pow(mu,2)*pow(rgrid[i],3)))) + 
               (2*(1 + rgrid[i]*Q11(i,j))*
                  ((6*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
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
                         pow(mu,2)*pow(rgrid[i],3),2))))/
                pow(rgrid[i],2)))/pow(rgrid[i],2)))/
      (8.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
        (1 + rgrid[i]*Qrr(i,j)))) + 
  (h22.diff(d01,i,j)*(8*L*rgrid[i]*Qr1(i,j) + 
       (complex(0,48)*L*w*pow(rgrid[i],3)*Qr1(i,j))/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))) + 
       (pow(-1 + rgrid[i],2)*pow(rgrid[i],8)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          ((4*Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
               pow(1 + rgrid[i]*Qrr(i,j),2))/
             (pow(-1 + rgrid[i],2)*pow(rgrid[i],7)*
               pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3),2)) - 
            (2*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
               ((3*(1 + rgrid[i]*Q22(i,j))*Qr1(i,j)*
                    (-((L*(-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                      Q11.diff(d01,i,j)*Qr1(i,j)))/rgrid[i] + 
                 (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3)))))/
             ((-1 + rgrid[i])*pow(rgrid[i],4)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))) + 
            (pow(1 + rgrid[i]*Q11(i,j),2)*
               ((8*L*(1 + rgrid[i]*Q22(i,j))*
                    (Qr1.diff(d10,i,j)*rgrid[i] + Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],4)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 pow(rgrid[i],2)*pow(Qr1(i,j),2)*
                  ((4*Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],3)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) - 
                    (6*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],3)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3)))) + 
                 rgrid[i]*Qr1(i,j)*
                  ((6*L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],5)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) - 
                    (2*(1 + rgrid[i]*Q22(i,j))*
                       ((10*Qr1.diff(d01,i,j)*
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
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2))))/
                     pow(rgrid[i],2))))/pow(rgrid[i],4)))/
        (2.*pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j)))))/(4.*pow(L,2)) + 
  (2*A0.diff(d10,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) + 
  (2*A0(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*pow(rgrid[i],2)*
     (-12 + pow(mu,2) + (-12 + pow(mu,2))*rgrid[i] + 
       (-12 + pow(mu,2) + complex(0,6)*w)*pow(rgrid[i],2))*
     (-(L*(-(mu*a0(i,j)) - a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
       a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
   (L*(-12 + pow(mu,2))*(-1 + rgrid[i])*
     (1 + rgrid[i] + pow(rgrid[i],2))*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,4)*a0.diff(d01,i,j)*mu*muj*w*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) - 
  (complex(0,4)*a0.diff(d01,i,j)*mu*w*Ax(i,j)*exp(B1*mu*h(i,j)*rgrid[i])*
     pow(rgrid[i],2)*(1 + rgrid[i]*Qrr(i,j)))/
   (L*(-1 + rgrid[i])*pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3),2)*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  (2*A0.diff(d01,i,j)*exp(B1*mu*h(i,j)*rgrid[i])*(1 - rgrid[i])*
     pow(rgrid[i],4)*(((1 + rgrid[i]*Q11(i,j))*Qr1(i,j)*
          (-(L*(-(mu*a0(i,j)) - 
                 a0.diff(d10,i,j)*mu*(-1 + rgrid[i]))) - 
            a0.diff(d01,i,j)*mu*(-1 + rgrid[i])*rgrid[i]*Qr1(i,j)))/
        rgrid[i] - (2*a0.diff(d01,i,j)*mu*(1 + rgrid[i]*Qrr(i,j)))/
        (pow(rgrid[i],2)*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3)))))/
   (pow(L,2)*(-1 + rgrid[i])*
     (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
       pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
     (1 + rgrid[i]*Qtt(i,j))) + 
  Dh(i,j)*(mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
     (h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j))/L + 
     (complex(0,6)*w*pow(rgrid[i],3)*
        (L*mu*(h(i,j) + h.diff(d10,i,j)*rgrid[i]) - 
          h.diff(d01,i,j)*mu*pow(rgrid[i],2)*Qr1(i,j)))/
      (L*(-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) - 
     (exp((-2*mu*h(i,j)*rgrid[i])/(3.*A1))*pow(rgrid[i],5)*
        ((-2*pow(a0.diff(d01,i,j),2)*(2 + 3*pow(A1,2))*B1*
             pow(mu,2)*exp((2*mu*h(i,j)*rgrid[i])/(3.*A1) + 
               B1*mu*h(i,j)*rgrid[i])*(-1 + rgrid[i])*
             (1 + rgrid[i]*Qrr(i,j)))/
           (pow(rgrid[i],2)*
             (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))) - 
          ((1 + rgrid[i]*Q11(i,j))*
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
                pow(rgrid[i],4)))/pow(rgrid[i],2)))/
      ((2 + 3*pow(A1,2))*pow(L,2)*(-1 + rgrid[i])*
        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
          pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q11(i,j))*
        (1 + rgrid[i]*Qtt(i,j)))) + 
  (h11(i,j)*pow(rgrid[i],8)*((complex(0,-6)*L*w*
          (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
          (((1 + rgrid[i]*Q22(i,j))*
               (-((L*(-2 + Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j)))/pow(rgrid[i],3)) + 
                 (2*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Q11(i,j)))/
                  rgrid[i] + Q11.diff(d01,i,j)*Qr1(i,j)))/
             pow(rgrid[i],2) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j)))/pow(rgrid[i],3) - 
                 Q22.diff(d01,i,j)*Qr1(i,j)))/pow(rgrid[i],2)))/
        ((-12 + pow(mu,2))*pow(rgrid[i],2)*(-1 + pow(rgrid[i],3))) + 
       (pow(rgrid[i],4)*((Q11.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qrr(i,j))*
               ((2*Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],3)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))))*
               (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7) - 
            ((1 + rgrid[i]*Q11(i,j))*
               ((-2*pow(Q22.diff(d01,i,j),2)*
                    pow(1 + rgrid[i]*Qrr(i,j),2)*
                    (1 + rgrid[i]*Qtt(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],8)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 ((1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                    ((2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (4*Q22.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],3)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
                 (pow(1 + rgrid[i]*Q22(i,j),2)*
                    ((-2*pow(Qrr.diff(d01,i,j),2)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (4*Qrr.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j))*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (8*pow(1 + rgrid[i]*Qrr(i,j),2)*
                        (pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-1 + rgrid[i],2) + 
                        (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                       (pow(-1 + rgrid[i],2)*pow(rgrid[i],4)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2))))/
                  pow(rgrid[i],4)))/pow(rgrid[i],2)))/
        ((1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))))/
   (4.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
     pow(1 + rgrid[i]*Q22(i,j),2)) + 
  (h22(i,j)*((72*w*(complex(0,1)*pow(mu,2) + 2*(complex(0,-6) + w))*
          pow(rgrid[i],4))/
        (pow(-12 + pow(mu,2),2)*(-1 + rgrid[i])*
          pow(1 + rgrid[i] + pow(rgrid[i],2),2)) - 
       (complex(0,48)*w*rgrid[i])/
        ((-12 + pow(mu,2))*(1 + rgrid[i] + pow(rgrid[i],2))) + 
       (complex(0,24)*w*pow(rgrid[i],2)*
          (-2 + (pow(-1 + rgrid[i],2)*pow(rgrid[i],6)*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               ((6*(1 + rgrid[i]*Q11(i,j))*
                    (-((L*(-2 + 
                       Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q22(i,j)))/pow(rgrid[i],3)) + 
                      Q22.diff(d01,i,j)*Qr1(i,j))*
                    (1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],4)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 ((1 + rgrid[i]*Q22(i,j))*
                    ((-6*L*(-2 + 
                        Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                        rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],5)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      rgrid[i]*Qr1(i,j)*
                       ((-4*Qrr.diff(d01,i,j)*
                        (1 + rgrid[i]*Q11(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],3)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                        (6*Q11.diff(d01,i,j)*
                        (1 + rgrid[i]*Qrr(i,j)))/
                        ((-1 + rgrid[i])*pow(rgrid[i],3)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                      (2*(1 + rgrid[i]*Q11(i,j))*
                        ((6*Qr1.diff(d01,i,j)*
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
                        pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2))))/
                       pow(rgrid[i],2)))/pow(rgrid[i],2)))/
             (8.*L*(1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Q22(i,j))*
               (1 + rgrid[i]*Qrr(i,j)))))/
        ((-12 + pow(mu,2))*(-1 + pow(rgrid[i],3))) + 
       (pow(rgrid[i],12)*(-((Q11.diff(d01,i,j)*(-1 + rgrid[i])*
                 (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
                 ((2*Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j)))/
                    ((-1 + rgrid[i])*pow(rgrid[i],3)*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                   (2*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                    ((-1 + rgrid[i])*pow(rgrid[i],3)*
                      (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))))*
                 (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],7)) + 
            (L*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*
               pow(1 + rgrid[i]*Q11(i,j),2)*(1 + rgrid[i]*Q22(i,j))*
               ((-6*L*(-2 + Q22.diff(d10,i,j)*pow(rgrid[i],2) - 
                      rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j)))/
                  ((-1 + rgrid[i])*pow(rgrid[i],5)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 rgrid[i]*Qr1(i,j)*
                  ((-4*Qrr.diff(d01,i,j)*(1 + rgrid[i]*Q22(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],3)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                    (6*Q22.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                     ((-1 + rgrid[i])*pow(rgrid[i],3)*
                       (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3)))) + 
                 (2*(1 + rgrid[i]*Q22(i,j))*
                    ((6*Qr1.diff(d01,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (4*L*(((-1 + rgrid[i])*
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
                        pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3),2))))/
                  pow(rgrid[i],2))*(1 + rgrid[i]*Qtt(i,j)))/
             (2.*pow(rgrid[i],8)) + 
            ((1 + rgrid[i]*Q11(i,j))*
               ((-2*pow(Q22.diff(d01,i,j),2)*
                    pow(1 + rgrid[i]*Qrr(i,j),2)*
                    (1 + rgrid[i]*Qtt(i,j)))/
                  (pow(rgrid[i],8)*
                    (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                      pow(mu,2)*pow(rgrid[i],3))) + 
                 ((-1 + rgrid[i])*(1 + rgrid[i]*Q22(i,j))*
                    (1 + rgrid[i]*Qrr(i,j))*
                    ((2*Q22.diff(d01,i,j)*Qrr.diff(d01,i,j))/
                       ((-1 + rgrid[i])*pow(rgrid[i],2)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) + 
                      (4*Q22.diff(d02,i,j)*(1 + rgrid[i]*Qrr(i,j)))/
                       ((-1 + rgrid[i])*pow(rgrid[i],3)*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))))*
                    (1 + rgrid[i]*Qtt(i,j)))/pow(rgrid[i],6) + 
                 (pow(1 + rgrid[i]*Q22(i,j),2)*
                    ((-2*pow(Qrr.diff(d01,i,j),2)*
                        (1 + rgrid[i]*Qtt(i,j)))/
                       (pow(rgrid[i],4)*
                         (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3))) + 
                      (((4*Qrr.diff(d02,i,j))/
                        (rgrid[i]*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))) - 
                        (3*pow(L,2)*
                       (-2 + 
                       Q11.diff(d10,i,j)*pow(rgrid[i],2) - 
                       rgrid[i]*Q11(i,j)))/pow(rgrid[i],3) + 
                        3*Q11.diff(d01,i,j)*L*Qr1(i,j))*
                        (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j))\
)/pow(rgrid[i],4) + (8*pow(1 + rgrid[i]*Qrr(i,j),2)*
                         (pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                        exp(B1*mu*h(i,j)*rgrid[i])*
                        pow(-1 + rgrid[i],2) + 
                         (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                        (-1 + rgrid[i])*
                        (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                        pow(mu,2)*pow(rgrid[i],3))*
                        (1 + rgrid[i]*Qtt(i,j)))/2.))/
                       ((-1 + rgrid[i])*pow(rgrid[i],4)*
                         pow(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                         pow(mu,2)*pow(rgrid[i],3),2))))/
                  pow(rgrid[i],4)))/pow(rgrid[i],2)))/
        (pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),2)*
          pow(1 + rgrid[i]*Q22(i,j),2)*(1 + rgrid[i]*Qrr(i,j))*
          (1 + rgrid[i]*Qtt(i,j)))))/4.
;
