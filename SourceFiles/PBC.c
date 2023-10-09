
elem[0]=
(12*A0(i,j) - pow(mu,2)*A0(i,j) - complex(0,2)*w*A0(i,j))/
   (-12 + pow(mu,2)) - B1*mu*a0(i,j)*Dh(i,j) - (mu*a0(i,j)*htt(i,j))/2. + 
  (Ax.diff(d01,i,j)*(1 + Qtt(i,j)))/(L*(1 + Q11(i,j))) + 
  (muj*(-Q11.diff(d01,i,j) + Q22.diff(d01,i,j) + 
       2*h.diff(d01,i,j)*B1*mu + Q22.diff(d01,i,j)*Q11(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q11(i,j) - 
       Q11.diff(d01,i,j)*Q22(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q22(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q11(i,j)*Q22(i,j))*(1 + Qtt(i,j)))/
   (2.*L*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))) + 
  (Ax(i,j)*(-Q11.diff(d01,i,j) + Q22.diff(d01,i,j) + 
       2*h.diff(d01,i,j)*B1*mu + Q22.diff(d01,i,j)*Q11(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q11(i,j) - Q11.diff(d01,i,j)*Q22(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q22(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q11(i,j)*Q22(i,j))*(1 + Qtt(i,j)))/
   (2.*L*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[1]=
-A0.diff(d10,i,j) - (htt.diff(d10,i,j)*mu*a0(i,j))/2. - 
  Dh.diff(d10,i,j)*B1*mu*a0(i,j) + 
  (12*A0.diff(d10,i,j) - A0.diff(d10,i,j)*pow(mu,2) - 
     complex(0,2)*A0.diff(d10,i,j)*w - complex(0,2)*w*A0(i,j))/
   (-12 + pow(mu,2)) - 2*a0.diff(d10,i,j)*B1*mu*Dh(i,j) - 
  B1*mu*a0(i,j)*Dh(i,j) + mu*a0(i,j)*h22(i,j) - 
  a0.diff(d10,i,j)*mu*htt(i,j) + (A0.diff(d01,i,j)*Qr1(i,j))/L + 
  (a0.diff(d01,i,j)*B1*mu*Dh(i,j)*Qr1(i,j))/L + 
  (a0.diff(d01,i,j)*mu*htt(i,j)*Qr1(i,j))/(2.*L) - 
  (complex(0,0.5)*Ax.diff(d11,i,j)*(-12 + pow(mu,2))*(1 + Qtt(i,j)))/
   (L*w*(1 + Q11(i,j))) + (complex(0,0.5)*ht1.diff(d01,i,j)*mu*
     (-12 + pow(mu,2))*a0(i,j)*(1 + Qtt(i,j)))/(L*w*(1 + Q11(i,j))) - 
  (complex(0,0.25)*Ax.diff(d10,i,j)*(-12 + pow(mu,2))*
     (-Q11.diff(d01,i,j) + Q22.diff(d01,i,j) + 
       2*h.diff(d01,i,j)*B1*mu + Q22.diff(d01,i,j)*Q11(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q11(i,j) - 
       Q11.diff(d01,i,j)*Q22(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q22(i,j) + 
       2*h.diff(d01,i,j)*B1*mu*Q11(i,j)*Q22(i,j))*(1 + Qtt(i,j)))/
   (L*w*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))) + 
  (complex(0,0.5)*Ax.diff(d02,i,j)*(-12 + pow(mu,2))*Qr1(i,j)*
     (1 + Qtt(i,j)))/(pow(L,2)*w*(1 + Q11(i,j))) - 
  (complex(0,0.25)*(-12 + pow(mu,2))*ht1(i,j)*
     (Q11.diff(d01,i,j)*L*mu*a0(i,j)*(1 + Q22(i,j))*
        pow(1 + Qtt(i,j),2) + 
       (1 + Q11(i,j))*(-(Q22.diff(d01,i,j)*L*mu*a0(i,j)*
             pow(1 + Qtt(i,j),2)) - 
          2*(1 + Q22(i,j))*(a0.diff(d01,i,j)*L*mu + 
             h.diff(d01,i,j)*B1*L*pow(mu,2)*a0(i,j) + 
             2*a0.diff(d01,i,j)*L*mu*Qtt(i,j) + 
             2*h.diff(d01,i,j)*B1*L*pow(mu,2)*a0(i,j)*Qtt(i,j) + 
             a0.diff(d01,i,j)*L*mu*pow(Qtt(i,j),2) + 
             h.diff(d01,i,j)*B1*L*pow(mu,2)*a0(i,j)*pow(Qtt(i,j),2)))\
))/(pow(L,2)*w*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*(1 + Qtt(i,j))) + 
  (muj*((2*L*(((-12 + pow(mu,2))*(-Qrr.diff(d10,i,j) - Qtt(i,j)))/
             (pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*
               pow(1 + Qtt(i,j),2)) + 
            (((-12 + pow(mu,2))*(-Q22.diff(d10,i,j) - Q22(i,j)))/
                (pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),2)) + 
               ((-12 + 3*pow(mu,2) + 8*(-12 + pow(mu,2)))/
                   pow(1 + Q11(i,j),2) - 
                  (2*(-12 + pow(mu,2))*
                     (Q11.diff(d10,i,j) + Q11(i,j)))/
                   pow(1 + Q11(i,j),3))/(1 + Q22(i,j)))/(1 + Qtt(i,j))\
)*(Q11.diff(d01,i,j)*(-1 - Q22(i,j))*pow(1 + Qtt(i,j),2) + 
            (1 + Q11(i,j))*(Q22.diff(d01,i,j)*pow(1 + Qtt(i,j),2) + 
               2*(1 + Q22(i,j))*
                (h.diff(d01,i,j)*B1*mu + 
                  2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                  h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2)))))/
        (-12 + pow(mu,2)) + ((-12 + pow(mu,2))*
          ((complex(0,1)*(Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                  (-1 - Q22(i,j))*pow(1 + Qtt(i,j),2) + 
                 (1 + Q11(i,j))*
                  (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                     pow(1 + Qtt(i,j),2) + 
                    2*(1 + Q22(i,j))*
                     (Qr1.diff(d02,i,j) + 
                       h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu + 
                       2*Qr1.diff(d02,i,j)*Qtt(i,j) + 
                       2*h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                        Qtt(i,j) + 
                       Qr1.diff(d02,i,j)*pow(Qtt(i,j),2) + 
                       h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                        pow(Qtt(i,j),2)))))/w + 
            (6*L*((Q11.diff(d01,i,j)*(-1 - Q22(i,j))*
                     pow(1 + Qtt(i,j),2) + 
                    (1 + Q11(i,j))*
                     (Q22.diff(d01,i,j)*pow(1 + Qtt(i,j),2) + 
                       2*(1 + Q22(i,j))*
                        (h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2))))/3. \
+ (-(Q11.diff(d11,i,j)*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2)) + 
                    (Q11.diff(d10,i,j) + Q11(i,j) - 
                       2*(1 + Q11(i,j)))*
                     (Q22.diff(d01,i,j)*pow(1 + Qtt(i,j),2) + 
                       2*(1 + Q22(i,j))*
                        (h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2))) + 
                    Q11.diff(d01,i,j)*
                     (-((1 + Q22(i,j))*(1 + Qtt(i,j))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j))) - 
                       (1 + Qtt(i,j))*
                        ((Q22.diff(d10,i,j) + Q22(i,j) - 
                        7*(1 + Q22(i,j)))*(1 + Qtt(i,j)) + 
                        (1 + Q22(i,j))*(Qrr.diff(d10,i,j) + Qtt(i,j))\
)) + (1 + Q11(i,j))*(Q22.diff(d11,i,j)*pow(1 + Qtt(i,j),2) + 
                       2*(Q22.diff(d10,i,j) + Q22(i,j) - 
                        2*(1 + Q22(i,j)))*
                        (h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2)) + 
                       (1 + Q22(i,j))*
                        (-Qrr.diff(d11,i,j) + 
                        Qrr.diff(d10,i,j)*Qtt.diff(d01,i,j) - 
                        Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j) + Qtt\
.diff(d11,i,j) - 6*h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d11,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*Qrr.diff(d10,i,j)*B1*
                       mu + 2*h.diff(d01,i,j)*Qtt.diff(d10,i,j)*
                        B1*mu - Qrr.diff(d11,i,j)*Qtt(i,j) + 
                        Qtt.diff(d11,i,j)*Qtt(i,j) - 
                        8*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        4*h.diff(d11,i,j)*B1*mu*Qtt(i,j) + 
                        2*h.diff(d01,i,j)*Qrr.diff(d10,i,j)*B1*mu*
                        Qtt(i,j) + 
                        2*h.diff(d01,i,j)*Qtt.diff(d10,i,j)*B1*mu*
                        Qtt(i,j) - 
                        2*h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2) + 
                        2*h.diff(d11,i,j)*B1*mu*pow(Qtt(i,j),2)) + 
                       Q22.diff(d01,i,j)*
                        ((1 + Qtt(i,j))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j)) + 
                         (1 + Qtt(i,j))*
                         (Qrr.diff(d10,i,j) + Qtt(i,j) - 
                         5*(1 + Qtt(i,j))))))/3.))/(-12 + pow(mu,2))))/
        (pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*(1 + Qtt(i,j)))))/
   (4.*pow(L,2)) + ((2*L*(((-12 + pow(mu,2))*Ax(i,j)*
             (-Qrr.diff(d10,i,j) - Qtt(i,j)))/
           (pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2)) \
+ (((-12 + pow(mu,2))*Ax(i,j)*(-Q22.diff(d10,i,j) - Q22(i,j)))/
              (pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),2)) + 
             ((Ax.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                   (-12 + 3*pow(mu,2) + 8*(-12 + pow(mu,2)))*
                    Ax(i,j))/pow(1 + Q11(i,j),2) - 
                (2*(-12 + pow(mu,2))*Ax(i,j)*
                   (Q11.diff(d10,i,j) + Q11(i,j)))/
                 pow(1 + Q11(i,j),3))/(1 + Q22(i,j)))/(1 + Qtt(i,j)))*
        (Q11.diff(d01,i,j)*(-1 - Q22(i,j))*pow(1 + Qtt(i,j),2) + 
          (1 + Q11(i,j))*(Q22.diff(d01,i,j)*pow(1 + Qtt(i,j),2) + 
             2*(1 + Q22(i,j))*
              (h.diff(d01,i,j)*B1*mu + 
                2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2)))))/
      (-12 + pow(mu,2)) + ((-12 + pow(mu,2))*Ax(i,j)*
        ((complex(0,1)*(Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                (-1 - Q22(i,j))*pow(1 + Qtt(i,j),2) + 
               (1 + Q11(i,j))*
                (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                   pow(1 + Qtt(i,j),2) + 
                  2*(1 + Q22(i,j))*
                   (Qr1.diff(d02,i,j) + 
                     h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu + 
                     2*Qr1.diff(d02,i,j)*Qtt(i,j) + 
                     2*h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                      Qtt(i,j) + 
                     Qr1.diff(d02,i,j)*pow(Qtt(i,j),2) + 
                     h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                      pow(Qtt(i,j),2)))))/w + 
          (6*L*((Q11.diff(d01,i,j)*(-1 - Q22(i,j))*
                   pow(1 + Qtt(i,j),2) + 
                  (1 + Q11(i,j))*
                   (Q22.diff(d01,i,j)*pow(1 + Qtt(i,j),2) + 
                     2*(1 + Q22(i,j))*
                      (h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2))))/3. \
+ (-(Q11.diff(d11,i,j)*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2)) + 
                  (Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
                   (Q22.diff(d01,i,j)*pow(1 + Qtt(i,j),2) + 
                     2*(1 + Q22(i,j))*
                      (h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2))) + 
                  Q11.diff(d01,i,j)*
                   (-((1 + Q22(i,j))*(1 + Qtt(i,j))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j))) - 
                     (1 + Qtt(i,j))*
                      ((Q22.diff(d10,i,j) + Q22(i,j) - 
                        7*(1 + Q22(i,j)))*(1 + Qtt(i,j)) + 
                        (1 + Q22(i,j))*(Qrr.diff(d10,i,j) + Qtt(i,j))\
)) + (1 + Q11(i,j))*(Q22.diff(d11,i,j)*pow(1 + Qtt(i,j),2) + 
                     2*(Q22.diff(d10,i,j) + Q22(i,j) - 
                        2*(1 + Q22(i,j)))*
                      (h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2)) + 
                     (1 + Q22(i,j))*
                      (-Qrr.diff(d11,i,j) + 
                        Qrr.diff(d10,i,j)*Qtt.diff(d01,i,j) - 
                        Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j) + Qtt\
.diff(d11,i,j) - 6*h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d11,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*Qrr.diff(d10,i,j)*B1*
                       mu + 2*h.diff(d01,i,j)*Qtt.diff(d10,i,j)*
                        B1*mu - Qrr.diff(d11,i,j)*Qtt(i,j) + 
                        Qtt.diff(d11,i,j)*Qtt(i,j) - 
                        8*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        4*h.diff(d11,i,j)*B1*mu*Qtt(i,j) + 
                        2*h.diff(d01,i,j)*Qrr.diff(d10,i,j)*B1*mu*
                        Qtt(i,j) + 
                        2*h.diff(d01,i,j)*Qtt.diff(d10,i,j)*B1*mu*
                        Qtt(i,j) - 
                        2*h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2) + 
                        2*h.diff(d11,i,j)*B1*mu*pow(Qtt(i,j),2)) + 
                     Q22.diff(d01,i,j)*
                      ((1 + Qtt(i,j))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j)) + 
                        (1 + Qtt(i,j))*
                         (Qrr.diff(d10,i,j) + Qtt(i,j) - 
                         5*(1 + Qtt(i,j))))))/3.))/(-12 + pow(mu,2))))/
      (pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*(1 + Qtt(i,j))))/(4.*pow(L,2)) \
+ (2*L*(1 + Q11(i,j))*(Ax.diff(d11,i,j)/pow(1 + Q11(i,j),2) + 
        Ax.diff(d01,i,j)*(4/pow(1 + Q11(i,j),2) - 
           (2*(Q11.diff(d10,i,j) + Q11(i,j)))/pow(1 + Q11(i,j),3)))*
      (1 + Qtt(i,j)) + (Ax.diff(d01,i,j)*
        (2*L*(Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
           (1 + Qtt(i,j)) - (complex(0,0.5)*Q11.diff(d01,i,j)*
             (-12 + pow(mu,2))*Qr1(i,j)*(1 + Qtt(i,j)))/w + 
          (1 + Q11(i,j))*((complex(0,0.5)*Qtt.diff(d01,i,j)*
                (-12 + pow(mu,2))*Qr1(i,j))/w + 
             ((4*L*((-12 + 3*pow(mu,2) - 2*(-12 + pow(mu,2)))*
                      (1 + Qtt(i,j)) + 
                     (-12 + pow(mu,2))*
                      (Qtt.diff(d10,i,j) + Qtt(i,j))))/
                 (-12 + pow(mu,2)) + 
                (-12 + pow(mu,2))*(1 + Qtt(i,j))*
                 ((4*L)/(-12 + pow(mu,2)) + 
                   (complex(0,4)*Qr1.diff(d01,i,j))/w + 
                   (complex(0,0.5)*(-12 + pow(mu,2))*Qr1(i,j)*
                      ((2*Q22.diff(d01,i,j)*(1 + Qtt(i,j)))/
                         (-12 + pow(mu,2)) + 
                        (1 + Q22(i,j))*
                         ((-2*Qtt.diff(d01,i,j))/(-12 + pow(mu,2)) + 
                         (4*h.diff(d01,i,j)*B1*mu*(1 + Qtt(i,j)))/
                         (-12 + pow(mu,2)))))/
                    (w*(1 + Q22(i,j))*(1 + Qtt(i,j)))))/2.)))/
      pow(1 + Q11(i,j),2))/(2.*pow(L,2))
;
elem[2]=
Ax.diff(d10,i,j) + (complex(0,4)*A0.diff(d01,i,j)*w)/
   (L*pow(-12 + pow(mu,2),2)) + 
  (complex(0,4)*Ax.diff(d10,i,j)*w)/(-12 + pow(mu,2)) + 
  (complex(0,4)*a0.diff(d01,i,j)*B1*mu*w*Dh(i,j))/
   (L*pow(-12 + pow(mu,2),2)) - 
  (complex(0,4)*a0.diff(d01,i,j)*mu*w*h11(i,j))/
   (L*pow(-12 + pow(mu,2),2)) - mu*a0(i,j)*ht1(i,j) - 
  (complex(0,2)*mu*w*a0(i,j)*ht1(i,j))/(-12 + pow(mu,2)) + 
  (complex(0,2)*a0.diff(d01,i,j)*mu*w*htt(i,j))/
   (L*pow(-12 + pow(mu,2),2)) - 
  (Ax.diff(d01,i,j)*(-12 + pow(mu,2) + complex(0,4)*w)*Qr1(i,j))/
   (L*(-12 + pow(mu,2))) + (muj*
     ((complex(0,8)*w)/(-12 + pow(mu,2)) - 
       (complex(0,8)*(-12 + pow(mu,2) - complex(0,2)*w)*w)/
        pow(-12 + pow(mu,2),2) - 
       (2*(-1728*Qr1.diff(d01,i,j) + 
            432*Qr1.diff(d01,i,j)*pow(mu,2) - 
            36*Qr1.diff(d01,i,j)*pow(mu,4) + 
            Qr1.diff(d01,i,j)*pow(mu,6) - 96*L*pow(w,2) + 
            48*Qrr.diff(d10,i,j)*L*pow(w,2) - 
            48*Qtt.diff(d10,i,j)*L*pow(w,2) + 
            24*L*pow(mu,2)*pow(w,2) - 
            4*Qrr.diff(d10,i,j)*L*pow(mu,2)*pow(w,2) + 
            4*Qtt.diff(d10,i,j)*L*pow(mu,2)*pow(w,2) - 
            1728*Qr1.diff(d01,i,j)*Qtt(i,j) + 
            432*Qr1.diff(d01,i,j)*pow(mu,2)*Qtt(i,j) - 
            36*Qr1.diff(d01,i,j)*pow(mu,4)*Qtt(i,j) + 
            Qr1.diff(d01,i,j)*pow(mu,6)*Qtt(i,j) - 
            96*L*pow(w,2)*Qtt(i,j) + 
            24*L*pow(mu,2)*pow(w,2)*Qtt(i,j)))/
        (L*pow(-12 + pow(mu,2),3)*(1 + Qtt(i,j))) - 
       (complex(0,6)*w*(-2*L*(1 + Q11(i,j))*(1 + Q22(i,j))*
             pow(1 + Qtt(i,j),2)*
             ((-Qtt.diff(d10,i,j) - Qtt(i,j))/
                (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),3)) \
+ ((-Qrr.diff(d10,i,j) - Qtt(i,j))/
                   (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*
                     pow(1 + Qtt(i,j),2)) + 
                  ((-Q22.diff(d10,i,j) - Q22(i,j))/
                      (3.*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
                     ((-Q11.diff(d10,i,j) - Q11(i,j))/
                        (3.*pow(1 + Q11(i,j),2)) + 3/(1 + Q11(i,j)))/
                      (1 + Q22(i,j)))/(1 + Qtt(i,j)))/(1 + Qtt(i,j))) + 
            (-2*L*(Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
                (1 + Q22(i,j))*pow(1 + Qtt(i,j),2) - 
               (1 + Q22(i,j))*
                (2*L - Q11.diff(d10,i,j)*L + L*Q11(i,j) + 
                  Q11.diff(d01,i,j)*Qr1(i,j))*pow(1 + Qtt(i,j),2) + 
               (1 + Q11(i,j))*
                (-2*L*(Q22.diff(d10,i,j) + Q22(i,j) - 
                     2*(1 + Q22(i,j)))*pow(1 + Qtt(i,j),2) + 
                  (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
                     Q22.diff(d01,i,j)*Qr1(i,j))*pow(1 + Qtt(i,j),2) \
+ (1 + Q22(i,j))*(-((L*(1 + Qtt(i,j))*
                        (12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                        Qtt.diff(d10,i,j)*pow(mu,2) + 
                        2*pow(mu,2)*Qtt(i,j)))/(-12 + pow(mu,2))) + 
                     2*Qr1(i,j)*
                      (h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2)) + 
                     2*(-(L*(-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j))*
                        ((-2/(-12 + pow(mu,2)) - 
                        (3*(-4 + pow(mu,2)))/
                        pow(-12 + pow(mu,2),2))*(1 + Qtt(i,j)) + 
                        (Qrr.diff(d10,i,j) + Qtt(i,j))/
                        (-12 + pow(mu,2))))/2. + 
                        ((1 + Qtt(i,j))*
                         ((-12 + pow(mu,2))*
                        (Qr1.diff(d01,i,j) - 
                        h.diff(d10,i,j)*B1*L*mu - B1*L*mu*h(i,j))*
                        (1 + Qtt(i,j)) - 
                         L*(12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                         Qtt.diff(d10,i,j)*pow(mu,2) + 
                         2*pow(mu,2)*Qtt(i,j))))/(-12 + pow(mu,2)))))\
)/(3.*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2))))/
        (L*(-12 + pow(mu,2)))))/2. + 
  (Ax(i,j)*((complex(0,8)*w)/(-12 + pow(mu,2)) - 
       (complex(0,8)*(-12 + pow(mu,2) - complex(0,2)*w)*w)/
        pow(-12 + pow(mu,2),2) - 
       (2*(-1728*Qr1.diff(d01,i,j) + 
            432*Qr1.diff(d01,i,j)*pow(mu,2) - 
            36*Qr1.diff(d01,i,j)*pow(mu,4) + 
            Qr1.diff(d01,i,j)*pow(mu,6) - 96*L*pow(w,2) + 
            48*Qrr.diff(d10,i,j)*L*pow(w,2) - 
            48*Qtt.diff(d10,i,j)*L*pow(w,2) + 
            24*L*pow(mu,2)*pow(w,2) - 
            4*Qrr.diff(d10,i,j)*L*pow(mu,2)*pow(w,2) + 
            4*Qtt.diff(d10,i,j)*L*pow(mu,2)*pow(w,2) - 
            1728*Qr1.diff(d01,i,j)*Qtt(i,j) + 
            432*Qr1.diff(d01,i,j)*pow(mu,2)*Qtt(i,j) - 
            36*Qr1.diff(d01,i,j)*pow(mu,4)*Qtt(i,j) + 
            Qr1.diff(d01,i,j)*pow(mu,6)*Qtt(i,j) - 
            96*L*pow(w,2)*Qtt(i,j) + 
            24*L*pow(mu,2)*pow(w,2)*Qtt(i,j)))/
        (L*pow(-12 + pow(mu,2),3)*(1 + Qtt(i,j))) - 
       (complex(0,6)*w*(-2*L*(1 + Q11(i,j))*(1 + Q22(i,j))*
             pow(1 + Qtt(i,j),2)*
             ((-Qtt.diff(d10,i,j) - Qtt(i,j))/
                (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),3)) \
+ ((-Qrr.diff(d10,i,j) - Qtt(i,j))/
                   (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*
                     pow(1 + Qtt(i,j),2)) + 
                  ((-Q22.diff(d10,i,j) - Q22(i,j))/
                      (3.*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
                     ((-Q11.diff(d10,i,j) - Q11(i,j))/
                        (3.*pow(1 + Q11(i,j),2)) + 3/(1 + Q11(i,j)))/
                      (1 + Q22(i,j)))/(1 + Qtt(i,j)))/(1 + Qtt(i,j))) + 
            (-2*L*(Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
                (1 + Q22(i,j))*pow(1 + Qtt(i,j),2) - 
               (1 + Q22(i,j))*
                (2*L - Q11.diff(d10,i,j)*L + L*Q11(i,j) + 
                  Q11.diff(d01,i,j)*Qr1(i,j))*pow(1 + Qtt(i,j),2) + 
               (1 + Q11(i,j))*
                (-2*L*(Q22.diff(d10,i,j) + Q22(i,j) - 
                     2*(1 + Q22(i,j)))*pow(1 + Qtt(i,j),2) + 
                  (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
                     Q22.diff(d01,i,j)*Qr1(i,j))*pow(1 + Qtt(i,j),2) \
+ (1 + Q22(i,j))*(-((L*(1 + Qtt(i,j))*
                         (12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                         Qtt.diff(d10,i,j)*pow(mu,2) + 
                         2*pow(mu,2)*Qtt(i,j)))/(-12 + pow(mu,2))) + 
                     2*Qr1(i,j)*
                      (h.diff(d01,i,j)*B1*mu + 
                        2*h.diff(d01,i,j)*B1*mu*Qtt(i,j) + 
                        h.diff(d01,i,j)*B1*mu*pow(Qtt(i,j),2)) + 
                     2*(-(L*(-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j))*
                         ((-2/(-12 + pow(mu,2)) - 
                        (3*(-4 + pow(mu,2)))/
                        pow(-12 + pow(mu,2),2))*(1 + Qtt(i,j)) + 
                         (Qrr.diff(d10,i,j) + Qtt(i,j))/
                         (-12 + pow(mu,2))))/2. + 
                        ((1 + Qtt(i,j))*
                         ((-12 + pow(mu,2))*
                         (Qr1.diff(d01,i,j) - 
                        h.diff(d10,i,j)*B1*L*mu - B1*L*mu*h(i,j))*
                         (1 + Qtt(i,j)) - 
                         L*(12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                         Qtt.diff(d10,i,j)*pow(mu,2) + 
                         2*pow(mu,2)*Qtt(i,j))))/(-12 + pow(mu,2))))))/
             (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2))))/
        (L*(-12 + pow(mu,2)))))/2.
;
elem[3]=
Dh.diff(d02,i,j)/(pow(L,2)*(1 + Q11(i,j))) - 
  (h11.diff(d01,i,j)*h.diff(d01,i,j)*mu)/(pow(L,2)*(1 + Q11(i,j))) - 
  (h.diff(d01,i,j)*htt.diff(d01,i,j)*mu)/
   (2.*pow(L,2)*(1 + Q11(i,j))) - 
  (complex(0,1)*h.diff(d01,i,j)*mu*w*ht1(i,j))/(L*(1 + Q11(i,j))) + 
  (B1*mu*(-12 + pow(mu,2) + complex(0,2)*w)*a0(i,j)*A0(i,j)*
     exp(B1*mu*h(i,j)))/((-12 + pow(mu,2))*pow(1 + Qtt(i,j),2)) + 
  (Dh.diff(d10,i,j)*(-12 + pow(mu,2) + complex(0,4)*w))/
   (2.*(1 + Qtt(i,j))) - (complex(0,2)*a0.diff(d01,i,j)*B1*mu*muj*w*
     exp(B1*mu*h(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) - 
  (complex(0,2)*a0.diff(d01,i,j)*B1*mu*w*Ax(i,j)*exp(B1*mu*h(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) + 
  ((-12 + pow(mu,2))*htt(i,j)*
     ((complex(0,-2)*L*w*(h.diff(d10,i,j)*L*mu + L*mu*h(i,j) - 
            h.diff(d01,i,j)*mu*Qr1(i,j)))/(-12 + pow(mu,2)) + 
       (2*B1*pow(L,2)*pow(mu,2)*pow(a0(i,j),2)*exp(B1*mu*h(i,j)))/
        ((-12 + pow(mu,2))*(1 + Qtt(i,j)))))/
   (4.*pow(L,2)*(1 + Qtt(i,j))) + 
  (Dh.diff(d01,i,j)*(144 - 24*pow(mu,2) + pow(mu,4))*
     ((-4*L*(-12 + pow(mu,2) + complex(0,4)*w)*Qr1(i,j)*(1 + Qtt(i,j)))/
        pow(-12 + pow(mu,2),2) + 
       (4*(Q22.diff(d01,i,j)*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2) + 
            (1 + Q22(i,j))*(Qtt.diff(d01,i,j)*(1 + Q11(i,j))*
                (1 + Qtt(i,j)) - 
               (1 + Qtt(i,j))*
                (Q11.diff(d01,i,j) - Qtt.diff(d01,i,j) - 
                  Qtt.diff(d01,i,j)*Q11(i,j) + 
                  Q11.diff(d01,i,j)*Qtt(i,j)))))/
        (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))))/
   (8.*pow(L,2)*pow(1 + Qtt(i,j),2)) + 
  (h11(i,j)*(h.diff(d01,i,j)*Q11.diff(d01,i,j)*mu*(1 + Q22(i,j))*
        pow(1 + Qtt(i,j),2) - 
       (1 + Q11(i,j))*(h.diff(d01,i,j)*Q22.diff(d01,i,j)*mu*
           pow(1 + Qtt(i,j),2) + 
          (1 + Q22(i,j))*(h.diff(d01,i,j)*Qtt.diff(d01,i,j)*mu*
              (1 + Qtt(i,j)) + 
             (1 + Qtt(i,j))*(2*h.diff(d02,i,j)*mu + 
                h.diff(d01,i,j)*Qtt.diff(d01,i,j)*mu + 
                2*h.diff(d02,i,j)*mu*Qtt(i,j))))))/
   (2.*pow(L,2)*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*
     pow(1 + Qtt(i,j),2)) + (Dh(i,j)*
     ((16*exp((-2*mu*h(i,j))/(3.*A1)) + 24*pow(A1,2)*exp(A1*mu*h(i,j)))/
        (2 + 3*pow(A1,2)) + (-12 + pow(mu,2))/(1 + Qtt(i,j)) + 
       (complex(0,4)*w)/(1 + Qtt(i,j)) + 
       (4*pow(w,2))/((-12 + pow(mu,2))*(1 + Qtt(i,j))) - 
       (complex(0,2)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
          (-72 + 12*Qrr.diff(d10,i,j) + 8*pow(mu,2) - 
            Qrr.diff(d10,i,j)*pow(mu,2) - 60*Qtt(i,j) + 
            7*pow(mu,2)*Qtt(i,j)))/
        (pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2)) - 
       (48*pow(w,2) - 48*Qtt.diff(d10,i,j)*pow(w,2) + 
          4*pow(mu,2)*pow(w,2) + 
          4*Qtt.diff(d10,i,j)*pow(mu,2)*pow(w,2) - 
          144*pow(B1,2)*pow(mu,2)*pow(a0(i,j),2)*
           exp(B1*mu*h(i,j)) + 
          24*pow(B1,2)*pow(mu,4)*pow(a0(i,j),2)*exp(B1*mu*h(i,j)) - 
          pow(B1,2)*pow(mu,6)*pow(a0(i,j),2)*exp(B1*mu*h(i,j)) + 
          8*pow(mu,2)*pow(w,2)*Qtt(i,j))/
        (pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2)) + 
       (complex(0,3)*w*(2*L*(1 + Q11(i,j))*(1 + Q22(i,j))*
             pow(1 + Qtt(i,j),2)*
             (((-12 + pow(mu,2))*(-Qtt.diff(d10,i,j) - Qtt(i,j)))/
                (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),4)) \
+ ((((-12 + pow(mu,2))*(-Q22.diff(d10,i,j) - Q22(i,j)))/
                      (3.*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
                     (((-12 + pow(mu,2))*
                       (-Q11.diff(d10,i,j) - Q11(i,j)))/
                        (3.*pow(1 + Q11(i,j),2)) + 
                        ((11*(-12 + pow(mu,2)))/3. + 
                        (-12 + 3*pow(mu,2))/3.)/(1 + Q11(i,j)))/
                      (1 + Q22(i,j)))/pow(1 + Qtt(i,j),2) - 
                  (2*(-12 + pow(mu,2))*
                     (Qrr.diff(d10,i,j) + Qtt(i,j)))/
                   (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*
                     pow(1 + Qtt(i,j),3)))/(1 + Qtt(i,j))) + 
            ((-12 + pow(mu,2))*
               (2*L*(Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
                  (1 + Q22(i,j))*pow(1 + Qtt(i,j),2) - 
                 (1 + Q22(i,j))*
                  (2*L - Q11.diff(d10,i,j)*L + L*Q11(i,j) + 
                    Q11.diff(d01,i,j)*Qr1(i,j))*pow(1 + Qtt(i,j),2) + 
                 (1 + Q11(i,j))*
                  (2*L*(Q22.diff(d10,i,j) + Q22(i,j) - 
                       2*(1 + Q22(i,j)))*pow(1 + Qtt(i,j),2) - 
                    (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
                       Q22.diff(d01,i,j)*Qr1(i,j))*
                     pow(1 + Qtt(i,j),2) + 
                    ((1 + Q22(i,j))*
                       (24*Qr1.diff(d01,i,j) - 
                         12*Qrr.diff(d10,i,j)*L - 
                         36*Qtt.diff(d10,i,j)*L - 
                         2*Qr1.diff(d01,i,j)*pow(mu,2) + 
                         4*L*pow(mu,2) + 
                         Qrr.diff(d10,i,j)*L*pow(mu,2) + 
                         3*Qtt.diff(d10,i,j)*L*pow(mu,2) + 
                         48*Qr1.diff(d01,i,j)*Qtt(i,j) - 
                         48*L*Qtt(i,j) - 
                         12*Qrr.diff(d10,i,j)*L*Qtt(i,j) - 
                         36*Qtt.diff(d10,i,j)*L*Qtt(i,j) - 
                         4*Qr1.diff(d01,i,j)*pow(mu,2)*Qtt(i,j) + 
                         12*L*pow(mu,2)*Qtt(i,j) + 
                         Qrr.diff(d10,i,j)*L*pow(mu,2)*Qtt(i,j) + 
                         3*Qtt.diff(d10,i,j)*L*pow(mu,2)*Qtt(i,j) + 
                         24*Qr1.diff(d01,i,j)*pow(Qtt(i,j),2) - 
                         48*L*pow(Qtt(i,j),2) - 
                         2*Qr1.diff(d01,i,j)*pow(mu,2)*
                         pow(Qtt(i,j),2) + 
                         8*L*pow(mu,2)*pow(Qtt(i,j),2)))/
                     (-12 + pow(mu,2)))))/
             (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),3))))/
        (L*(-12 + pow(mu,2)))))/2.
;
elem[4]=
h11.diff(d10,i,j) + (complex(0,4)*h11.diff(d10,i,j)*w)/
   (-12 + pow(mu,2)) + (complex(0,8)*a0.diff(d01,i,j)*mu*muj*w*
     exp(B1*mu*h(i,j)))/(L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))) + 
  (complex(0,8)*a0.diff(d01,i,j)*mu*w*Ax(i,j)*exp(B1*mu*h(i,j)))/
   (L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))) - 
  (h11.diff(d01,i,j)*(-12 + pow(mu,2) + complex(0,4)*w)*Qr1(i,j))/
   (L*(-12 + pow(mu,2))) - (htt.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (4*Dh.diff(d01,i,j)*h.diff(d01,i,j)*mu*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (complex(0,2)*ht1.diff(d01,i,j)*w*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (htt(i,j)*(complex(0,1.5)*L*(-12 + pow(mu,2))*w*(1 + Q22(i,j))*
        (2*Qr1.diff(d01,i,j) + 2*L - Q11.diff(d10,i,j)*L + 
          2*Qr1.diff(d01,i,j)*Q11(i,j) + L*Q11(i,j) + 
          Q11.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)) - 
       complex(0,1.5)*L*(-12 + pow(mu,2))*w*(1 + Q11(i,j))*
        (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
          Q22.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j))))/
   (3.*pow(L,2)*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*(1 + Q22(i,j))*
     (1 + Qtt(i,j))) + (complex(0,2)*w*ht1(i,j)*
     ((Q22.diff(d01,i,j)*(-12 + pow(mu,2))*(1 + Q11(i,j))*
          (1 + Qtt(i,j)))/2. + 
       (1 + Q22(i,j))*(-(Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
             (1 + Q11(i,j))) + 
          (Q11.diff(d01,i,j)*(-12 + pow(mu,2))*(1 + Qtt(i,j)))/2.)))/
   (L*pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))) + 
  (htt.diff(d01,i,j)*(Q11.diff(d01,i,j)*(1 + Q22(i,j))*
        pow(1 + Qtt(i,j),2) + 
       (1 + Q11(i,j))*(Q22.diff(d01,i,j)*pow(1 + Qtt(i,j),2) - 
          2*(1 + Q22(i,j))*(Qtt.diff(d01,i,j) + 
             Qtt.diff(d01,i,j)*Qtt(i,j)))))/
   (2.*pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*
     (1 + Qtt(i,j))) + (h11(i,j)*
     ((complex(0,16)*w)/(-12 + pow(mu,2)) - 
       (complex(0,16)*(-12 + pow(mu,2) - complex(0,2)*w)*w)/
        pow(-12 + pow(mu,2),2) + 
       4*((4*pow(w,2)*((-2/(-12 + pow(mu,2)) - 
                  (3*(-4 + pow(mu,2)))/pow(-12 + pow(mu,2),2))*
                (1 + Qtt(i,j)) + 
               (Qrr.diff(d10,i,j) + Qtt(i,j))/(-12 + pow(mu,2))))/
           ((-12 + pow(mu,2))*(1 + Qtt(i,j))) + 
          ((1 + Qtt(i,j))*((72*pow(A1,2)*
                   exp((-2*mu*h(i,j))/(3.*A1)) + 48*exp(A1*mu*h(i,j)))/
                (2 + 3*pow(A1,2)) - 
               (2*pow(h.diff(d01,i,j),2)*pow(mu,2))/
                (pow(L,2)*(1 + Q11(i,j))) + 
               pow(Qtt.diff(d01,i,j),2)/
                (pow(L,2)*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2)) + 
               (2*(((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-12 + pow(mu,2)))/2. - 
                       (-12*Qtt.diff(d02,i,j) + 
                        Qtt.diff(d02,i,j)*pow(mu,2))*
                        (1 + Q11(i,j)) - 
                       4*pow(L,2)*pow(w,2)*
                        (2 - Q11.diff(d10,i,j) + 3*Q11(i,j) - 
                        Q11.diff(d10,i,j)*Q11(i,j) + 
                        pow(Q11(i,j),2)))/
                     ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*
                       (1 + Qtt(i,j))) + 
                    2*pow(L,2)*pow(w,2)*
                     (1 + 2*Q11(i,j) + pow(Q11(i,j),2))*
                     ((-Qtt.diff(d10,i,j) - Qtt(i,j))/
                        ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*
                         pow(1 + Qtt(i,j),2)) + 
                       ((6/(-12 + pow(mu,2)) - 
                        (3*(-4 + pow(mu,2)))/
                        pow(-12 + pow(mu,2),2))/
                        pow(1 + Q11(i,j),2) - 
                         (2*(Q11.diff(d10,i,j) + Q11(i,j)))/
                         ((-12 + pow(mu,2))*pow(1 + Q11(i,j),3)))/
                        (1 + Qtt(i,j)))))/pow(L,2)))/(-12 + pow(mu,2))) \
+ (-2*pow(L,2)*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2)*
           ((-Qtt.diff(d10,i,j) - Qtt(i,j))/
              ((1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),3)) + 
             ((-Qrr.diff(d10,i,j) - Qtt(i,j))/
                 ((1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2)) \
+ ((-Q22.diff(d10,i,j) - Q22(i,j))/
                    ((1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
                   ((-Q11.diff(d10,i,j) - Q11(i,j))/
                       pow(1 + Q11(i,j),2) + 8/(1 + Q11(i,j)))/
                    (1 + Q22(i,j)))/(1 + Qtt(i,j)))/(1 + Qtt(i,j))) + 
          (-2*pow(L,2)*(1 + Q11(i,j))*
              (Q22.diff(d10,i,j) + Q22(i,j) - 2*(1 + Q22(i,j)))*
              pow(1 + Qtt(i,j),2) + 
             L*(1 + Q11(i,j))*
              (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
                Q22.diff(d01,i,j)*Qr1(i,j))*pow(1 + Qtt(i,j),2) + 
             (1 + Q22(i,j))*(Q11.diff(d01,i,j)*L*Qr1(i,j)*
                 pow(1 + Qtt(i,j),2) - 
                (2*pow(L,2)*(1 + Qtt(i,j))*
                   (48 - 12*Q11.diff(d10,i,j) - 
                     24*Qtt.diff(d10,i,j) + 
                     Q11.diff(d10,i,j)*pow(mu,2) + 
                     2*Qtt.diff(d10,i,j)*pow(mu,2) + 36*Q11(i,j) - 
                     24*Qtt.diff(d10,i,j)*Q11(i,j) + 
                     pow(mu,2)*Q11(i,j) + 
                     2*Qtt.diff(d10,i,j)*pow(mu,2)*Q11(i,j) + 
                     24*Qtt(i,j) - 12*Q11.diff(d10,i,j)*Qtt(i,j) + 
                     2*pow(mu,2)*Qtt(i,j) + 
                     Q11.diff(d10,i,j)*pow(mu,2)*Qtt(i,j) + 
                     12*Q11(i,j)*Qtt(i,j) + 
                     3*pow(mu,2)*Q11(i,j)*Qtt(i,j)))/
                 (-12 + pow(mu,2)) + 
                ((-12 + pow(mu,2))*(1 + Qtt(i,j))*
                   ((8*pow(Qtt.diff(d01,i,j),2))/
                      pow(-12 + pow(mu,2),2) + 
                     (2*L*(2*Qr1.diff(d01,i,j) + 2*L - 
                        Q11.diff(d10,i,j)*L + 
                        2*Qr1.diff(d01,i,j)*Q11(i,j) + L*Q11(i,j))*
                        (1 + Qtt(i,j)))/(-12 + pow(mu,2))))/2.))/
           ((1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2)))/
        pow(L,2) - (complex(0,12)*w*
          (-2*L*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2)*
             ((-Qtt.diff(d10,i,j) - Qtt(i,j))/
                (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),3)) \
+ ((-Qrr.diff(d10,i,j) - Qtt(i,j))/
                   (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*
                     pow(1 + Qtt(i,j),2)) + 
                  ((-Q22.diff(d10,i,j) - Q22(i,j))/
                      (3.*(1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
                     ((-Q11.diff(d10,i,j) - Q11(i,j))/
                        (3.*pow(1 + Q11(i,j),2)) + 3/(1 + Q11(i,j)))/
                      (1 + Q22(i,j)))/(1 + Qtt(i,j)))/(1 + Qtt(i,j))) + 
            (-2*L*(1 + Q11(i,j))*
                (Q22.diff(d10,i,j) + Q22(i,j) - 2*(1 + Q22(i,j)))*
                pow(1 + Qtt(i,j),2) + 
               (1 + Q11(i,j))*
                (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
                  Q22.diff(d01,i,j)*Qr1(i,j))*pow(1 + Qtt(i,j),2) + 
               (1 + Q22(i,j))*
                (L*(2 - Q11.diff(d10,i,j) + Q11(i,j))*
                   pow(1 + Qtt(i,j),2) - 
                  2*L*(Q11.diff(d10,i,j) + Q11(i,j) - 
                     2*(1 + Q11(i,j)))*pow(1 + Qtt(i,j),2) + 
                  Qr1(i,j)*(Q11.diff(d01,i,j) + 
                     2*Q11.diff(d01,i,j)*Qtt(i,j) + 
                     Q11.diff(d01,i,j)*pow(Qtt(i,j),2)) + 
                  (1 + Q11(i,j))*
                   (-((L*(1 + Qtt(i,j))*
                        (12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                        Qtt.diff(d10,i,j)*pow(mu,2) + 
                        2*pow(mu,2)*Qtt(i,j)))/(-12 + pow(mu,2))) + 
                     ((1 + Qtt(i,j))*
                        (-24*Qr1.diff(d01,i,j) - 60*L + 
                         12*Qrr.diff(d10,i,j)*L + 
                         24*Qtt.diff(d10,i,j)*L + 
                         2*Qr1.diff(d01,i,j)*pow(mu,2) + 
                         3*L*pow(mu,2) - 
                         Qrr.diff(d10,i,j)*L*pow(mu,2) - 
                         2*Qtt.diff(d10,i,j)*L*pow(mu,2) - 
                         24*Qr1.diff(d01,i,j)*Qtt(i,j) - 
                         24*L*Qtt(i,j) + 
                         2*Qr1.diff(d01,i,j)*pow(mu,2)*Qtt(i,j)))/
                      (-12 + pow(mu,2)))))/
             (3.*(1 + Q11(i,j))*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2))))/
        (L*(-12 + pow(mu,2))) + 
       (4*((pow(L,2)*(1 + 2*Q11(i,j) + pow(Q11(i,j),2))*
               (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
               pow(-12 + pow(mu,2) - 12*Qtt(i,j) + 
                 pow(mu,2)*Qtt(i,j),2)*
               ((((12/pow(-12 + pow(mu,2),2) - 
                       (6*(-4 + pow(mu,2)))/
                       pow(-12 + pow(mu,2),3))/
                       pow(1 + Q11(i,j),2) - 
                       (2*(Q11.diff(d10,i,j) + Q11(i,j)))/
                        (pow(-12 + pow(mu,2),2)*
                        pow(1 + Q11(i,j),3)))/pow(1 + Q22(i,j),2) \
- (2*(Q22.diff(d10,i,j) + Q22(i,j)))/
                     (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
                       pow(1 + Q22(i,j),3)))/pow(1 + Qtt(i,j),2) - 
                 (2*(Qtt.diff(d10,i,j) + Qtt(i,j)))/
                  (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
                    pow(1 + Q22(i,j),2)*pow(1 + Qtt(i,j),3))))/2. + 
            ((pow(L,2)*(-4*(1 + 2*Q11(i,j) + pow(Q11(i,j),2)) + 
                    2*(Q11.diff(d10,i,j) + Q11(i,j) + 
                       Q11.diff(d10,i,j)*Q11(i,j) + pow(Q11(i,j),2)\
))*(1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                  pow(-12 + pow(mu,2) - 12*Qtt(i,j) + 
                    pow(mu,2)*Qtt(i,j),2))/2. + 
               Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                (-12 + pow(mu,2))*(1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)) - 
               ((-12 + pow(mu,2))*(1 + Q11(i,j))*
                  (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*(1 + Qtt(i,j))*
                  (2*pow(Qtt.diff(d01,i,j),2) + 
                    4*(Qtt.diff(d02,i,j) + 
                       Qtt.diff(d02,i,j)*Qtt(i,j)) - 
                    (pow(L,2)*(2 - Q11.diff(d10,i,j) + Q11(i,j))*
                       (-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j)))/2. - 
                    (Q11.diff(d01,i,j)*L*Qr1(i,j)*
                       (-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j)))/2.))/2. + 
               (1 + 2*Q11(i,j) + pow(Q11(i,j),2))*
                (((-12 + pow(mu,2))*(1 + Q22(i,j))*(1 + Qtt(i,j))*
                     ((pow(L,2)*
                        (2 - Q22.diff(d10,i,j) + Q22(i,j))*
                        (-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j)))/2. + 
                       (Q22.diff(d01,i,j)*L*Qr1(i,j)*
                        (-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j)))/2.))/2. + 
                  2*((pow(L,2)*
                        (-4*(1 + 2*Q22(i,j) + pow(Q22(i,j),2)) + 
                        2*(Q22.diff(d10,i,j) + Q22(i,j) + 
                        Q22.diff(d10,i,j)*Q22(i,j) + 
                        pow(Q22(i,j),2)))*
                        pow(-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j),2))/4. + 
                     (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                      (-(Qr1(i,j)*
                         (-(L*(-12 + pow(mu,2))*
                         (-12*Qtt.diff(d01,i,j) + 
                        Qtt.diff(d01,i,j)*pow(mu,2))*(1 + Qtt(i,j))) \
+ (Qtt.diff(d01,i,j)*L*(-12 + pow(mu,2))*
                         (-12 + pow(mu,2) - 12*Qtt(i,j) + 
                         pow(mu,2)*Qtt(i,j)))/2.)) + 
                        (L*(-12 + pow(mu,2))*(1 + Qtt(i,j))*
                         (2*L*pow(mu,2)*pow(a0(i,j),2)*
                        exp(B1*mu*h(i,j)) + 
                         (Qr1.diff(d01,i,j)*
                        (-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j)))/2. - 
                         2*L*
                         (12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                         Qtt.diff(d10,i,j)*pow(mu,2) + 
                         2*pow(mu,2)*Qtt(i,j))))/2. + 
                        pow(L,2)*
                         ((-3*
                         pow(-12 + pow(mu,2) - 12*Qtt(i,j) + 
                        pow(mu,2)*Qtt(i,j),2))/2. + 
                         ((-12 + pow(mu,2) - 12*Qtt(i,j) + 
                         pow(mu,2)*Qtt(i,j))*
                         (-12 - 24*Qtt.diff(d10,i,j) + 5*pow(mu,2) + 
                         2*Qtt.diff(d10,i,j)*pow(mu,2) - 
                         36*Qtt(i,j) + 7*pow(mu,2)*Qtt(i,j)))/2.)))))/
             (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
               pow(1 + Q22(i,j),2)*pow(1 + Qtt(i,j),2))))/pow(L,2)))/4.
;
elem[5]=
(w*(complex(0,-2) + (8*w)/(-12 + pow(mu,2)))*h22(i,j))/4. - 
  (complex(0,2)*a0.diff(d01,i,j)*mu*muj*w*exp(B1*mu*h(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (complex(0,2)*a0.diff(d01,i,j)*mu*w*Ax(i,j)*exp(B1*mu*h(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (complex(0,1)*w*Dh(i,j)*(h.diff(d10,i,j)*L*mu + L*mu*h(i,j) - 
       h.diff(d01,i,j)*mu*Qr1(i,j)))/L - 
  (complex(0,0.25)*w*h11(i,j)*((1 + Q22(i,j))*
        (2*Qr1.diff(d01,i,j) + 2*L - Q11.diff(d10,i,j)*L + 
          2*Qr1.diff(d01,i,j)*Q11(i,j) + L*Q11(i,j) + 
          Q11.diff(d01,i,j)*Qr1(i,j)) + 
       (-1 - Q11(i,j))*(2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
          Q22.diff(d01,i,j)*Qr1(i,j))))/(L*(1 + Q11(i,j))*(1 + Q22(i,j))) \
+ (complex(0,0.5)*ht1.diff(d01,i,j)*w*(1 + Qtt(i,j)))/(L*(1 + Q11(i,j))) - 
  (complex(0,0.25)*w*ht1(i,j)*(Q22.diff(d01,i,j)*(-1 - Q11(i,j))*
        pow(1 + Qtt(i,j),2) + 
       (1 + Q22(i,j))*(Qtt.diff(d01,i,j)*(1 + Q11(i,j))*(1 + Qtt(i,j)) + 
          (1 + Qtt(i,j))*(Q11.diff(d01,i,j) - 3*Qtt.diff(d01,i,j) - 
             3*Qtt.diff(d01,i,j)*Q11(i,j) + Q11.diff(d01,i,j)*Qtt(i,j)))\
))/(L*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*(1 + Qtt(i,j)))
;
elem[6]=
complex(0,-1)*h22.diff(d10,i,j)*w - 
  (a0.diff(d01,i,j)*Ax.diff(d10,i,j)*mu*exp(B1*mu*h(i,j)))/
   (L*(1 + Q11(i,j))) + (mu*muj*
     (-(a0.diff(d01,i,j)*(-3*Qr1.diff(d01,i,j)*
              (-12 + pow(mu,2)) + complex(0,12)*L*w)*exp(B1*mu*h(i,j)))/
        (3.*(1 + Q11(i,j))) - complex(0,6)*L*w*
        ((a0.diff(d11,i,j)*exp(B1*mu*h(i,j)))/(3.*(1 + Q11(i,j))) + 
          a0.diff(d01,i,j)*((exp(B1*mu*h(i,j))*
                (-Q11.diff(d10,i,j) - Q11(i,j)))/
              (3.*pow(1 + Q11(i,j),2)) + 
             (-exp(B1*mu*h(i,j))/3. + 
                (2*exp(B1*mu*h(i,j)) + 
                   B1*mu*exp(B1*mu*h(i,j))*(h.diff(d10,i,j) + h(i,j)))/
                 3.)/(1 + Q11(i,j))))))/(pow(L,2)*(-12 + pow(mu,2))) + 
  (mu*(-(a0.diff(d01,i,j)*(-3*Qr1.diff(d01,i,j)*
              (-12 + pow(mu,2)) + complex(0,12)*L*w)*Ax(i,j)*
           exp(B1*mu*h(i,j)))/(3.*(1 + Q11(i,j))) - 
       complex(0,6)*L*w*((a0.diff(d11,i,j)*Ax(i,j)*exp(B1*mu*h(i,j)))/
           (3.*(1 + Q11(i,j))) + 
          a0.diff(d01,i,j)*((Ax(i,j)*exp(B1*mu*h(i,j))*
                (-Q11.diff(d10,i,j) - Q11(i,j)))/
              (3.*pow(1 + Q11(i,j),2)) + 
             ((Ax.diff(d10,i,j)*exp(B1*mu*h(i,j)))/3. + 
                Ax(i,j)*(-exp(B1*mu*h(i,j))/3. + 
                   (2*exp(B1*mu*h(i,j)) + 
                      B1*mu*exp(B1*mu*h(i,j))*
                       (h.diff(d10,i,j) + h(i,j)))/3.))/(1 + Q11(i,j)))))\
)/(pow(L,2)*(-12 + pow(mu,2))) + 
  (complex(0,1)*h22.diff(d01,i,j)*w*Qr1(i,j))/L + 
  (a0.diff(d01,i,j)*Ax.diff(d01,i,j)*mu*exp(B1*mu*h(i,j))*Qr1(i,j))/
   (pow(L,2)*(1 + Q11(i,j))) + 
  (complex(0,1)*w*((Dh.diff(d10,i,j) + Dh(i,j))*
        (h.diff(d10,i,j)*L*mu + L*mu*h(i,j) - 
          h.diff(d01,i,j)*mu*Qr1(i,j)) + 
       Dh(i,j)*(-(h.diff(d01,i,j)*Qr1.diff(d10,i,j)*mu) + 
          2*h.diff(d10,i,j)*L*mu + h.diff(d20,i,j)*L*mu - 
          2*h.diff(d01,i,j)*mu*Qr1(i,j) - h.diff(d11,i,j)*mu*Qr1(i,j)))\
)/L - (complex(0,0.25)*w*(((h11(i,j)*(-Q22.diff(d10,i,j) - Q22(i,j)))/
           ((1 + Q11(i,j))*pow(1 + Q22(i,j),2)) + 
          ((h11(i,j)*(-Q11.diff(d10,i,j) - Q11(i,j)))/
              pow(1 + Q11(i,j),2) + 
             (h11.diff(d10,i,j) + 4*h11(i,j))/(1 + Q11(i,j)))/
           (1 + Q22(i,j)))*((1 + Q22(i,j))*
           (2*Qr1.diff(d01,i,j) + 2*L - Q11.diff(d10,i,j)*L + 
             2*Qr1.diff(d01,i,j)*Q11(i,j) + L*Q11(i,j) + 
             Q11.diff(d01,i,j)*Qr1(i,j)) + 
          (-1 - Q11(i,j))*(2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
             Q22.diff(d01,i,j)*Qr1(i,j))) + 
       (h11(i,j)*((Q22.diff(d10,i,j) + Q22(i,j) - 2*(1 + Q22(i,j)))*
             (2*Qr1.diff(d01,i,j) + 2*L - Q11.diff(d10,i,j)*L + 
               2*Qr1.diff(d01,i,j)*Q11(i,j) + L*Q11(i,j) + 
               Q11.diff(d01,i,j)*Qr1(i,j)) + 
            (1 + Q22(i,j))*(-2*Qr1.diff(d01,i,j) + 
               2*Q11.diff(d10,i,j)*Qr1.diff(d01,i,j) + 
               Q11.diff(d01,i,j)*Qr1.diff(d10,i,j) + 
               2*Qr1.diff(d11,i,j) - 6*L + 2*Q11.diff(d10,i,j)*L - 
               Q11.diff(d20,i,j)*L + 2*Qr1.diff(d11,i,j)*Q11(i,j) - 
               2*L*Q11(i,j) + Q11.diff(d11,i,j)*Qr1(i,j)) + 
            (Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
             (-2*L + Q22.diff(d10,i,j)*L - L*Q22(i,j) - 
               Q22.diff(d01,i,j)*Qr1(i,j)) + 
            (1 + Q11(i,j))*(-(Q22.diff(d01,i,j)*Qr1.diff(d10,i,j)) + 
               6*L - 2*Q22.diff(d10,i,j)*L + Q22.diff(d20,i,j)*L + 
               2*L*Q22(i,j) - Q22.diff(d11,i,j)*Qr1(i,j))))/
        ((1 + Q11(i,j))*(1 + Q22(i,j)))))/L + 
  (ht1.diff(d11,i,j)*(-12 + pow(mu,2))*(1 + Qtt(i,j)))/
   (4.*L*(1 + Q11(i,j))) - (ht1.diff(d02,i,j)*(-12 + pow(mu,2))*
     Qr1(i,j)*(1 + Qtt(i,j)))/(4.*pow(L,2)*(1 + Q11(i,j))) + 
  (ht1.diff(d10,i,j)*(-12 + pow(mu,2))*
     (Q22.diff(d01,i,j)*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2) - 
       (1 + Q22(i,j))*(Q11.diff(d01,i,j) - 2*Qtt.diff(d01,i,j) - 
          2*Qtt.diff(d01,i,j)*Q11(i,j) + 
          2*Q11.diff(d01,i,j)*Qtt(i,j) - 
          2*Qtt.diff(d01,i,j)*Qtt(i,j) - 
          2*Qtt.diff(d01,i,j)*Q11(i,j)*Qtt(i,j) + 
          Q11.diff(d01,i,j)*pow(Qtt(i,j),2))))/
   (8.*L*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*(1 + Qtt(i,j))) + 
  (complex(0,2)*L*w*(1 + Q11(i,j))*
      (ht1.diff(d11,i,j)/pow(1 + Q11(i,j),2) + 
        ht1.diff(d01,i,j)*(4/pow(1 + Q11(i,j),2) - 
           (2*(Q11.diff(d10,i,j) + Q11(i,j)))/pow(1 + Q11(i,j),3)))*
      (1 + Qtt(i,j)) + (ht1.diff(d01,i,j)*
        (complex(0,2)*L*w*(Q11.diff(d10,i,j) + Q11(i,j) - 
             2*(1 + Q11(i,j)))*(1 + Qtt(i,j)) + 
          (Q11.diff(d01,i,j)*(-12 + pow(mu,2))*Qr1(i,j)*
             (1 + Qtt(i,j)))/2. + 
          (1 + Q11(i,j))*(18*Qtt.diff(d01,i,j)*Qr1(i,j) - 
             (3*Qtt.diff(d01,i,j)*pow(mu,2)*Qr1(i,j))/2. + 
             ((-12 + pow(mu,2))*(1 + Qtt(i,j))*
                 (-4*Qr1.diff(d01,i,j) + 
                   (complex(0,4)*L*w)/(-12 + pow(mu,2)) - 
                   (Q22.diff(d01,i,j)*Qr1(i,j))/(1 + Q22(i,j)) + 
                   (Qtt.diff(d01,i,j)*Qr1(i,j))/(1 + Qtt(i,j))) + 
                (complex(0,4)*L*w*
                   ((-12 + 3*pow(mu,2) - 2*(-12 + pow(mu,2)))*
                      (1 + Qtt(i,j)) + 
                     (-12 + pow(mu,2))*
                      (Qtt.diff(d10,i,j) + Qtt(i,j))))/
                 (-12 + pow(mu,2)))/2.)))/pow(1 + Q11(i,j),2))/
   (4.*pow(L,2)) + ((complex(0,-2)*L*w*
        (((-12 + pow(mu,2))*ht1(i,j)*
             (-Qrr.diff(d10,i,j) - Qtt(i,j)))/
           (pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2)) \
+ (((-12 + pow(mu,2))*ht1(i,j)*(-Q22.diff(d10,i,j) - Q22(i,j)))/
              (pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),2)) + 
             ((ht1.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                   (-12 + 3*pow(mu,2) + 8*(-12 + pow(mu,2)))*
                    ht1(i,j))/pow(1 + Q11(i,j),2) - 
                (2*(-12 + pow(mu,2))*ht1(i,j)*
                   (Q11.diff(d10,i,j) + Q11(i,j)))/
                 pow(1 + Q11(i,j),3))/(1 + Q22(i,j)))/(1 + Qtt(i,j)))*
        (Q22.diff(d01,i,j)*(-1 - Q11(i,j))*pow(1 + Qtt(i,j),2) + 
          (1 + Q22(i,j))*(Qtt.diff(d01,i,j)*(1 + Q11(i,j))*
              (1 + Qtt(i,j)) + 
             (1 + Qtt(i,j))*(Q11.diff(d01,i,j) - 
                3*Qtt.diff(d01,i,j) - 3*Qtt.diff(d01,i,j)*Q11(i,j) + 
                Q11.diff(d01,i,j)*Qtt(i,j)))))/(-12 + pow(mu,2)) + 
     ((-12 + pow(mu,2))*ht1(i,j)*
        (Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*(1 + Q22(i,j))*
           pow(1 + Qtt(i,j),2) - 
          (1 + Q11(i,j))*(Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*
              pow(1 + Qtt(i,j),2) + 
             (1 + Q22(i,j))*((2*
                   ((3*Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-12 + pow(mu,2)))/2. - 
                     4*a0.diff(d01,i,j)*L*pow(mu,2)*a0(i,j)*
                      exp(B1*mu*h(i,j)))*(1 + Qtt(i,j)))/
                 (-12 + pow(mu,2)) - 
                (1 + Qtt(i,j))*
                 (-2*Qr1.diff(d02,i,j) + 
                   Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j) - 
                   2*Qr1.diff(d02,i,j)*Qtt(i,j)))) - 
          (complex(0,6)*L*w*((Q22.diff(d01,i,j)*(-1 - Q11(i,j))*
                   pow(1 + Qtt(i,j),2) + 
                  (1 + Q22(i,j))*
                   (Qtt.diff(d01,i,j)*(1 + Q11(i,j))*
                      (1 + Qtt(i,j)) + 
                     (1 + Qtt(i,j))*
                      (Q11.diff(d01,i,j) - 3*Qtt.diff(d01,i,j) - 
                        3*Qtt.diff(d01,i,j)*Q11(i,j) + 
                        Q11.diff(d01,i,j)*Qtt(i,j))))/3. + 
               (-(Q22.diff(d11,i,j)*(1 + Q11(i,j))*
                     pow(1 + Qtt(i,j),2)) + 
                  (Q22.diff(d10,i,j) + Q22(i,j) - 2*(1 + Q22(i,j)))*
                   (Qtt.diff(d01,i,j)*(1 + Q11(i,j))*
                      (1 + Qtt(i,j)) + 
                     (1 + Qtt(i,j))*
                      (Q11.diff(d01,i,j) - 3*Qtt.diff(d01,i,j) - 
                        3*Qtt.diff(d01,i,j)*Q11(i,j) + 
                        Q11.diff(d01,i,j)*Qtt(i,j))) + 
                  Q22.diff(d01,i,j)*
                   (-((1 + Q11(i,j))*(1 + Qtt(i,j))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j))) - 
                     (1 + Qtt(i,j))*
                      ((Q11.diff(d10,i,j) + Q11(i,j) - 
                        7*(1 + Q11(i,j)))*(1 + Qtt(i,j)) + 
                        (1 + Q11(i,j))*(Qrr.diff(d10,i,j) + Qtt(i,j))\
)) + (1 + Q22(i,j))*(-5*Q11.diff(d01,i,j) + Q11.diff(d11,i,j) + 
                     Q11.diff(d01,i,j)*Qrr.diff(d10,i,j) + 
                     15*Qtt.diff(d01,i,j) - 
                     3*Q11.diff(d10,i,j)*Qtt.diff(d01,i,j) - 
                     3*Qrr.diff(d10,i,j)*Qtt.diff(d01,i,j) + 
                     Q11.diff(d01,i,j)*Qtt.diff(d10,i,j) - 
                     3*Qtt.diff(d11,i,j) + 
                     12*Qtt.diff(d01,i,j)*Q11(i,j) - 
                     3*Qrr.diff(d10,i,j)*Qtt.diff(d01,i,j)*
                      Q11(i,j) - 3*Qtt.diff(d11,i,j)*Q11(i,j) - 
                     8*Q11.diff(d01,i,j)*Qtt(i,j) + 
                     2*Q11.diff(d11,i,j)*Qtt(i,j) + 
                     Q11.diff(d01,i,j)*Qrr.diff(d10,i,j)*Qtt(i,j) + 
                     12*Qtt.diff(d01,i,j)*Qtt(i,j) - 
                     3*Q11.diff(d10,i,j)*Qtt.diff(d01,i,j)*
                      Qtt(i,j) + 
                     Q11.diff(d01,i,j)*Qtt.diff(d10,i,j)*Qtt(i,j) - 
                     3*Qtt.diff(d11,i,j)*Qtt(i,j) + 
                     9*Qtt.diff(d01,i,j)*Q11(i,j)*Qtt(i,j) - 
                     3*Qtt.diff(d11,i,j)*Q11(i,j)*Qtt(i,j) - 
                     3*Q11.diff(d01,i,j)*pow(Qtt(i,j),2) + 
                     Q11.diff(d11,i,j)*pow(Qtt(i,j),2) + 
                     Qrr.diff(d11,i,j)*(1 + Q11(i,j))*
                      (1 + Qtt(i,j)) + 
                     Qtt.diff(d01,i,j)*
                      ((Q11.diff(d10,i,j) + Q11(i,j) - 
                        5*(1 + Q11(i,j)))*(1 + Qtt(i,j)) + 
                        (1 + Q11(i,j))*(Qtt.diff(d10,i,j) + Qtt(i,j)))\
))/3.))/(-12 + pow(mu,2))))/
      (pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*(1 + Qtt(i,j))))/(8.*pow(L,2)) \
+ (w*(-(h22.diff(d10,i,j)*(complex(0,2) - (8*w)/(-12 + pow(mu,2)))) + 
       h22(i,j)*((8*w)/(-12 + pow(mu,2)) - 
          (complex(0,2)*(L*(1 + Q11(i,j))*(1 + Q22(i,j))*
                (-12 + pow(mu,2) - 12*Qtt(i,j) + pow(mu,2)*Qtt(i,j))*
                (-((-Qtt.diff(d10,i,j) - Qtt(i,j))/
                     ((-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
                       pow(1 + Qtt(i,j),2))) + 
                  (-((-Q22.diff(d10,i,j) - Q22(i,j))/
                        ((-12 + pow(mu,2))*(1 + Q11(i,j))*
                         pow(1 + Q22(i,j),2))) + 
                     (-((-Q11.diff(d10,i,j) - Q11(i,j))/
                        ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2))) + 
                        (-6/(-12 + pow(mu,2)) + 
                        (3*(-4 + pow(mu,2)))/
                        pow(-12 + pow(mu,2),2))/(1 + Q11(i,j)))/
                      (1 + Q22(i,j)))/(1 + Qtt(i,j))) - 
               (((-12 + pow(mu,2))*(1 + Q11(i,j))*
                     (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
                       Q22.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)))/2. + 
                  L*(1 + Q11(i,j))*
                   (Q22.diff(d10,i,j) + Q22(i,j) - 2*(1 + Q22(i,j)))*
                   (-12 + pow(mu,2) - 12*Qtt(i,j) + 
                     pow(mu,2)*Qtt(i,j)) + 
                  (1 + Q22(i,j))*
                   ((L*(-12 + pow(mu,2))*
                        (2 - Q11.diff(d10,i,j) + Q11(i,j))*
                        (1 + Qtt(i,j)))/2. + 
                     Qr1(i,j)*
                      (12*Qtt.diff(d01,i,j) - 
                        Qtt.diff(d01,i,j)*pow(mu,2) + 
                        12*Qtt.diff(d01,i,j)*Q11(i,j) - 
                        Qtt.diff(d01,i,j)*pow(mu,2)*Q11(i,j) + 
                        (Q11.diff(d01,i,j)*(-12 + pow(mu,2))*
                         (1 + Qtt(i,j)))/2.) + 
                     2*((L*(Q11.diff(d10,i,j) + Q11(i,j) - 
                        2*(1 + Q11(i,j)))*
                         (-12 + pow(mu,2) - 12*Qtt(i,j) + 
                         pow(mu,2)*Qtt(i,j)))/2. + 
                        (1 + Q11(i,j))*
                         ((Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*
                         (1 + Qtt(i,j)))/2. + 
                         L*(12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                         Qtt.diff(d10,i,j)*pow(mu,2) + 
                         2*pow(mu,2)*Qtt(i,j))))))/
                ((-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
                  (1 + Qtt(i,j)))))/L)))/4.
;
elem[7]=
((-2 - (complex(0,4)*w)/(-12 + pow(mu,2)) - 
       (16*pow(w,2))/pow(-12 + pow(mu,2),2))*h22(i,j))/4. + 
  (2*mu*(-12 + pow(mu,2) + complex(0,2)*w)*a0(i,j)*A0(i,j)*
     exp(B1*mu*h(i,j)))/(pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j))) - 
  (h11.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (htt.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (2*Dh.diff(d01,i,j)*h.diff(d01,i,j)*mu*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) - 
  (complex(0,2)*ht1.diff(d01,i,j)*w*(1 + Qtt(i,j)))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (htt(i,j)*((complex(0,1)*w*((2*Qr1.diff(d01,i,j))/L + 
            (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
               (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
            (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
               (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))))/
        (-12 + pow(mu,2)) + (2*pow(mu,2)*pow(a0(i,j),2)*
          exp(B1*mu*h(i,j)))/((-12 + pow(mu,2))*(1 + Qtt(i,j)))))/2. - 
  (complex(0,2)*muj*w*exp(B1*mu*h(i,j))*
     (2*L*mu*a0(i,j)*(1 + Q11(i,j))*Qr1(i,j) - 
       (2*a0.diff(d01,i,j)*mu*(1 + Qtt(i,j)))/(-12 + pow(mu,2))))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) - 
  (complex(0,2)*w*Ax(i,j)*exp(B1*mu*h(i,j))*
     (2*L*mu*a0(i,j)*(1 + Q11(i,j))*Qr1(i,j) - 
       (2*a0.diff(d01,i,j)*mu*(1 + Qtt(i,j)))/(-12 + pow(mu,2))))/
   (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) + 
  (Dh(i,j)*((complex(0,-4)*L*w*
          (h.diff(d10,i,j)*L*mu + L*mu*h(i,j) + 
            h.diff(d01,i,j)*mu*Qr1(i,j)))/(-12 + pow(mu,2)) + 
       (2*exp((-2*mu*h(i,j))/(3.*A1))*
          ((2 + 3*pow(A1,2))*B1*pow(L,2)*pow(mu,2)*
             pow(a0(i,j),2)*exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j)) - 
            24*A1*pow(L,2)*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
             pow(1 + Qtt(i,j),2)))/
        ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))))/
   (2.*pow(L,2)) + (h11.diff(d01,i,j)*
     ((complex(0,8)*L*w*Qr1(i,j))/(-12 + pow(mu,2)) + 
       (2*(1 + Qtt(i,j))*(Q11.diff(d01,i,j) + 
            (1 + Q11(i,j))*((-3*Q22.diff(d01,i,j))/(1 + Q22(i,j)) - 
               (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))))/
        ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2))))/(4.*pow(L,2)) + 
  (htt.diff(d01,i,j)*(L*(2 + (complex(0,8)*w)/(-12 + pow(mu,2)))*
        Qr1(i,j) + (2*(1 + Qtt(i,j))*
          (Q11.diff(d01,i,j) + 
            (1 + Q11(i,j))*(-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) - 
               (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))))/
        ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2))))/(4.*pow(L,2)) + 
  (w*ht1(i,j)*((complex(0,2) - (4*w)/(-12 + pow(mu,2)))*Qr1(i,j) + 
       (complex(0,4)*(-(Q22.diff(d01,i,j)*(-12 + pow(mu,2))*
                (1 + Q11(i,j))*(1 + Qtt(i,j)))/2. + 
            (1 + Q22(i,j))*(-(Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                  (1 + Q11(i,j))) + 
               (Q11.diff(d01,i,j)*(-12 + pow(mu,2))*(1 + Qtt(i,j)))/2.\
)))/(L*pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))))/2. \
+ (h11(i,j)*((complex(0,2)*w*((2*Qr1.diff(d01,i,j))/L + 
            (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
               (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
            (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
               (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))))/
        (-12 + pow(mu,2)) + (8*
          ((Q11.diff(d01,i,j)*(-12 + pow(mu,2))*(1 + Q22(i,j))*
               (1 + Qtt(i,j))*
               ((Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                    (1 + Q22(i,j)))/2. + 
                 (Q22.diff(d01,i,j)*(-12 + pow(mu,2))*
                    (1 + Qtt(i,j)))/2.))/2. + 
            (1 + Q11(i,j))*((pow(Q22.diff(d01,i,j),2)*
                  (144 - 24*pow(mu,2) + pow(mu,4))*
                  (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/4. - 
               ((-12 + pow(mu,2))*(1 + Q22(i,j))*(1 + Qtt(i,j))*
                  ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-12 + pow(mu,2)))/2. + 
                    Q22.diff(d02,i,j)*(-12 + pow(mu,2))*
                     (1 + Qtt(i,j))))/2. + 
               (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                ((pow(Qtt.diff(d01,i,j),2)*
                     (144 - 24*pow(mu,2) + pow(mu,4)))/4. + 
                  (-12 + pow(mu,2))*
                   (6*Qtt.diff(d02,i,j) - 
                     (Qtt.diff(d02,i,j)*pow(mu,2))/2.)*(1 + Qtt(i,j)) \
- (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                     (144 - 24*pow(mu,2) + pow(mu,4))*
                     (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/2.))))/
        (pow(L,2)*pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
          pow(1 + Q22(i,j),2)*(1 + Qtt(i,j)))))/4.
;
elem[8]=
-h22.diff(d10,i,j)/2. - Dh.diff(d10,i,j)*h.diff(d10,i,j)*mu - 
  Dh.diff(d10,i,j)*mu*h(i,j) + (h11.diff(d11,i,j)*Qr1(i,j))/L + 
  (htt.diff(d11,i,j)*Qr1(i,j))/L - 
  (Dh.diff(d10,i,j)*h.diff(d01,i,j)*mu*Qr1(i,j))/L + 
  complex(0,1)*ht1.diff(d10,i,j)*w*Qr1(i,j) + 
  (htt.diff(d10,i,j)*((2*Qr1.diff(d01,i,j))/L + 
       (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
          (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
       (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
          (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))))/4. + 
  (h11.diff(d10,i,j)*((2*Qr1.diff(d01,i,j))/L + 
       (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
          (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
       (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
          (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))))/4. + 
  (2*A0.diff(d10,i,j)*mu*a0(i,j)*exp(B1*mu*h(i,j)))/
   ((-12 + pow(mu,2))*(1 + Qtt(i,j))) - 
  (h22.diff(d02,i,j)*(1 + Qtt(i,j)))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))) + 
  (2*A0.diff(d01,i,j)*exp(B1*mu*h(i,j))*
     (L*mu*a0(i,j)*(1 + Q11(i,j))*Qr1(i,j) - 
       (2*a0.diff(d01,i,j)*mu*(1 + Qtt(i,j)))/(-12 + pow(mu,2))))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))) - 
  (2*(((-36 + 3*pow(mu,2) + complex(0,6)*w)*A0(i,j)*exp(B1*mu*h(i,j))*
          (-2*a0.diff(d10,i,j)*L*mu - a0.diff(d01,i,j)*mu*Qr1(i,j)))/
        (3.*(-12 + pow(mu,2))*(1 + Qtt(i,j))) - 
       L*mu*a0(i,j)*(((-36 + 3*pow(mu,2) + complex(0,6)*w)*A0(i,j)*
             exp(B1*mu*h(i,j))*(-Qtt.diff(d10,i,j) - Qtt(i,j)))/
           (3.*(-12 + pow(mu,2))*pow(1 + Qtt(i,j),2)) + 
          ((A0.diff(d10,i,j)*(-36 + 3*pow(mu,2) + complex(0,6)*w)*
                exp(B1*mu*h(i,j)))/(3.*(-12 + pow(mu,2))) + 
             A0(i,j)*(((-36 + 3*pow(mu,2) + complex(0,12)*w)*
                   exp(B1*mu*h(i,j)))/(3.*(-12 + pow(mu,2))) + 
                (-36 + 3*pow(mu,2) + complex(0,6)*w)*
                 (-(((-4 + pow(mu,2))*exp(B1*mu*h(i,j)))/
                      pow(-12 + pow(mu,2),2)) + 
                   (-exp(B1*mu*h(i,j))/3. + 
                      (2*exp(B1*mu*h(i,j)) + 
                         B1*mu*exp(B1*mu*h(i,j))*
                         (h.diff(d10,i,j) + h(i,j)))/3.)/
                    (-12 + pow(mu,2)))))/(1 + Qtt(i,j)))))/
   (L*(-12 + pow(mu,2))) + (h22.diff(d01,i,j)*
     (4*L*Qr1(i,j) + (complex(0,8)*L*w*Qr1(i,j))/(-12 + pow(mu,2)) + 
       (2*L*(1 + 2*Q11(i,j) + pow(Q11(i,j),2))*(1 + Q22(i,j))*Qr1(i,j)*
           pow(1 + Qtt(i,j),2) - 
          (2*Q22.diff(d01,i,j)*(1 + Q11(i,j))*pow(1 + Qtt(i,j),3))/
           (-12 + pow(mu,2)) + 
          (2*Q11.diff(d01,i,j)*(1 + Q22(i,j))*pow(1 + Qtt(i,j),3))/
           (-12 + pow(mu,2)))/
        (pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2))))/
   (4.*pow(L,2)) + ((-2*h11.diff(d12,i,j)*(1 + Qtt(i,j)))/
      ((-12 + pow(mu,2))*(1 + Q11(i,j))) + 
     h11.diff(d02,i,j)*(-2*pow(Qr1(i,j),2) - 
        2*(((-Q11.diff(d10,i,j) - Q11(i,j))/
               ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
              (3*(-4 + pow(mu,2)))/
               (pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))))*(1 + Qtt(i,j)) \
+ (Qrr.diff(d10,i,j) + Qtt(i,j))/((-12 + pow(mu,2))*(1 + Q11(i,j))))))/
   (2.*pow(L,2)) + ((-2*htt.diff(d12,i,j)*(1 + Qtt(i,j)))/
      ((-12 + pow(mu,2))*(1 + Q11(i,j))) + 
     htt.diff(d02,i,j)*(-2*pow(Qr1(i,j),2) - 
        2*(((-Q11.diff(d10,i,j) - Q11(i,j))/
               ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) - 
              (3*(-4 + pow(mu,2)))/
               (pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))))*(1 + Qtt(i,j)) \
+ (Qrr.diff(d10,i,j) + Qtt(i,j))/((-12 + pow(mu,2))*(1 + Q11(i,j))))))/
   (2.*pow(L,2)) - (complex(0,1)*w*
     ((2*ht1.diff(d11,i,j)*(1 + Qtt(i,j)))/
        ((-12 + pow(mu,2))*(1 + Q11(i,j))) + 
       ht1.diff(d01,i,j)*((2*
             ((-Q11.diff(d10,i,j) - Q11(i,j))/pow(1 + Q11(i,j),2) + 
               2/(1 + Q11(i,j)))*(1 + Qtt(i,j)))/(-12 + pow(mu,2)) + 
          (pow(Qr1(i,j),2) + Q11(i,j)*pow(Qr1(i,j),2) - 
             (2*(-36 + 12*Qrr.diff(d10,i,j) + 5*pow(mu,2) - 
                  Qrr.diff(d10,i,j)*pow(mu,2) - 24*Qtt(i,j) + 
                  4*pow(mu,2)*Qtt(i,j)))/pow(-12 + pow(mu,2),2))/
           (1 + Q11(i,j)))))/L - 
  (complex(0,2)*muj*w*(((exp(B1*mu*h(i,j))*
             (-Qtt.diff(d10,i,j) - Qtt(i,j)))/
           ((-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2)) + 
          ((exp(B1*mu*h(i,j))*(-Q11.diff(d10,i,j) - Q11(i,j)))/
              ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
             ((-3*(-4 + pow(mu,2))*exp(B1*mu*h(i,j)))/
                 pow(-12 + pow(mu,2),2) + 
                (4*exp(B1*mu*h(i,j)) + 
                   B1*mu*exp(B1*mu*h(i,j))*
                    (h.diff(d10,i,j) + h(i,j)))/(-12 + pow(mu,2)))/
              (1 + Q11(i,j)))/(1 + Qtt(i,j)))*
        (2*L*mu*a0(i,j)*(1 + Q11(i,j))*Qr1(i,j) - 
          (2*a0.diff(d01,i,j)*mu*(1 + Qtt(i,j)))/(-12 + pow(mu,2))) + 
       (exp(B1*mu*h(i,j))*(2*(L*mu*a0(i,j)*
                (Qr1.diff(d10,i,j)*(1 + Q11(i,j)) + 
                  (-1 + Q11.diff(d10,i,j))*Qr1(i,j)) + 
               (1 + Q11(i,j))*Qr1(i,j)*
                (2*a0.diff(d10,i,j)*L*mu - 
                  a0.diff(d01,i,j)*mu*Qr1(i,j))) - 
            2*mu*((a0.diff(d11,i,j)*(1 + Qtt(i,j)))/
                (-12 + pow(mu,2)) + 
               a0.diff(d01,i,j)*
                ((-2/(-12 + pow(mu,2)) - 
                     (3*(-4 + pow(mu,2)))/pow(-12 + pow(mu,2),2))*
                   (1 + Qtt(i,j)) + 
                  (Qrr.diff(d10,i,j) + Qtt(i,j))/(-12 + pow(mu,2))))))/
        ((-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j)))))/L - 
  (complex(0,2)*w*(((Ax(i,j)*exp(B1*mu*h(i,j))*
             (-Qtt.diff(d10,i,j) - Qtt(i,j)))/
           ((-12 + pow(mu,2))*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2)) + 
          ((Ax(i,j)*exp(B1*mu*h(i,j))*
                (-Q11.diff(d10,i,j) - Q11(i,j)))/
              ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
             ((Ax.diff(d10,i,j)*exp(B1*mu*h(i,j)))/
                 (-12 + pow(mu,2)) + 
                Ax(i,j)*((-3*(-4 + pow(mu,2))*exp(B1*mu*h(i,j)))/
                    pow(-12 + pow(mu,2),2) + 
                   (4*exp(B1*mu*h(i,j)) + 
                      B1*mu*exp(B1*mu*h(i,j))*
                       (h.diff(d10,i,j) + h(i,j)))/
                    (-12 + pow(mu,2))))/(1 + Q11(i,j)))/(1 + Qtt(i,j)))*
        (2*L*mu*a0(i,j)*(1 + Q11(i,j))*Qr1(i,j) - 
          (2*a0.diff(d01,i,j)*mu*(1 + Qtt(i,j)))/(-12 + pow(mu,2))) + 
       (Ax(i,j)*exp(B1*mu*h(i,j))*
          (2*(L*mu*a0(i,j)*(Qr1.diff(d10,i,j)*(1 + Q11(i,j)) + 
                  (-1 + Q11.diff(d10,i,j))*Qr1(i,j)) + 
               (1 + Q11(i,j))*Qr1(i,j)*
                (2*a0.diff(d10,i,j)*L*mu - 
                  a0.diff(d01,i,j)*mu*Qr1(i,j))) - 
            2*mu*((a0.diff(d11,i,j)*(1 + Qtt(i,j)))/
                (-12 + pow(mu,2)) + 
               a0.diff(d01,i,j)*
                ((-2/(-12 + pow(mu,2)) - 
                     (3*(-4 + pow(mu,2)))/pow(-12 + pow(mu,2),2))*
                   (1 + Qtt(i,j)) + 
                  (Qrr.diff(d10,i,j) + Qtt(i,j))/(-12 + pow(mu,2))))))/
        ((-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j)))))/L + 
  ((2*h.diff(d01,i,j)*mu*(Dh.diff(d11,i,j)/(1 + Q11(i,j)) + 
          Dh.diff(d01,i,j)*
           ((-Q11.diff(d10,i,j) - Q11(i,j))/pow(1 + Q11(i,j),2) + 
             3/(1 + Q11(i,j))))*(1 + Qtt(i,j)))/(-12 + pow(mu,2)) + 
     (Dh.diff(d01,i,j)*((1 + Q11(i,j))*Qr1(i,j)*
           (-(h.diff(d10,i,j)*L*mu) - L*mu*h(i,j) + 
             3*h.diff(d01,i,j)*mu*Qr1(i,j)) + 
          2*mu*((h.diff(d11,i,j)*(1 + Qtt(i,j)))/(-12 + pow(mu,2)) + 
             h.diff(d01,i,j)*
              ((-(1/(-12 + pow(mu,2))) - 
                   (3*(-4 + pow(mu,2)))/pow(-12 + pow(mu,2),2))*
                 (1 + Qtt(i,j)) + 
                (Qrr.diff(d10,i,j) + Qtt(i,j))/(-12 + pow(mu,2))))))/
      (1 + Q11(i,j)))/pow(L,2) + 
  (htt.diff(d11,i,j)*(L*(2 + (complex(0,8)*w)/(-12 + pow(mu,2)))*
         Qr1(i,j) + (2*(1 + Qtt(i,j))*
           (Q11.diff(d01,i,j) + 
             (1 + Q11(i,j))*(-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) - 
                (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))))/
         ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2))) + 
     htt.diff(d01,i,j)*(L*(2 + (complex(0,8)*w)/(-12 + pow(mu,2)))*
         (Qr1.diff(d10,i,j) + Qr1(i,j)) + 
        pow(Qr1(i,j),2)*(Q11.diff(d01,i,j)/(1 + Q11(i,j)) - 
           Q22.diff(d01,i,j)/(1 + Q22(i,j)) - 
           (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j))) + 
        Qr1(i,j)*(-2*Qr1.diff(d01,i,j) + 
           L*((complex(0,8)*w)/(-12 + pow(mu,2)) + 
              (2 - Q11.diff(d10,i,j) + Q11(i,j))/(1 + Q11(i,j)) + 
              (-2 + Q22.diff(d10,i,j) - Q22(i,j))/(1 + Q22(i,j)) + 
              (2*(12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                   Qtt.diff(d10,i,j)*pow(mu,2) + 
                   2*pow(mu,2)*Qtt(i,j)))/
               ((-12 + pow(mu,2))*(1 + Qtt(i,j))))) + 
        2*((((2/(-12 + pow(mu,2)) - 
                    (3*(-4 + pow(mu,2)))/pow(-12 + pow(mu,2),2)\
)/pow(1 + Q11(i,j),2) - (2*(Q11.diff(d10,i,j) + Q11(i,j)))/
                  ((-12 + pow(mu,2))*pow(1 + Q11(i,j),3)))*
               (1 + Qtt(i,j)) + 
              (Qrr.diff(d10,i,j) + Qtt(i,j))/
               ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))*
            (Q11.diff(d01,i,j) + 
              (1 + Q11(i,j))*(-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) - 
                 (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))) + 
           ((1 + Qtt(i,j))*(-Q11.diff(d01,i,j) + Q11.diff(d11,i,j) + 
                (Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
                 (-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) - 
                   (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j))) + 
                (1 + Q11(i,j))*
                 ((-Q22.diff(d01,i,j) + 
                      Q22.diff(d01,i,j)*Q22.diff(d10,i,j) - Q22.di\
ff(d11,i,j) - Q22.diff(d11,i,j)*Q22(i,j))/pow(1 + Q22(i,j),2) + 
                   (2*(-Qtt.diff(d01,i,j) + 
                        Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j) - Qtt.d\
iff(d11,i,j) - Qtt.diff(d11,i,j)*Qtt(i,j)))/pow(1 + Qtt(i,j),2))))/
            ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))))/(4.*pow(L,2)) + 
  (h11.diff(d11,i,j)*((complex(0,8)*L*w*Qr1(i,j))/(-12 + pow(mu,2)) + 
        (2*(1 + Qtt(i,j))*(Q11.diff(d01,i,j) + 
             (1 + Q11(i,j))*((-3*Q22.diff(d01,i,j))/(1 + Q22(i,j)) - 
                (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))))/
         ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2))) + 
     h11.diff(d01,i,j)*((-2*Qr1.diff(d01,i,j) + 
           L*((complex(0,8)*w)/(-12 + pow(mu,2)) + 
              (2 - Q11.diff(d10,i,j) + Q11(i,j))/(1 + Q11(i,j)) + 
              (-2 + Q22.diff(d10,i,j) - Q22(i,j))/(1 + Q22(i,j))))*
         Qr1(i,j) + (Q11.diff(d01,i,j)*pow(Qr1(i,j),2))/
         (1 + Q11(i,j)) - (5*Q22.diff(d01,i,j)*pow(Qr1(i,j),2))/
         (1 + Q22(i,j)) + (complex(0,8)*L*w*
           (Qr1.diff(d10,i,j) + Qr1(i,j)))/(-12 + pow(mu,2)) + 
        2*((((2/(-12 + pow(mu,2)) - 
                    (3*(-4 + pow(mu,2)))/pow(-12 + pow(mu,2),2)\
)/pow(1 + Q11(i,j),2) - (2*(Q11.diff(d10,i,j) + Q11(i,j)))/
                  ((-12 + pow(mu,2))*pow(1 + Q11(i,j),3)))*
               (1 + Qtt(i,j)) + 
              (Qrr.diff(d10,i,j) + Qtt(i,j))/
               ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))*
            (Q11.diff(d01,i,j) + 
              (1 + Q11(i,j))*((-3*Q22.diff(d01,i,j))/(1 + Q22(i,j)) - 
                 (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))) + 
           ((1 + Qtt(i,j))*(-Q11.diff(d01,i,j) + Q11.diff(d11,i,j) + 
                (Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
                 ((-3*Q22.diff(d01,i,j))/(1 + Q22(i,j)) - 
                   (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j))) + 
                (1 + Q11(i,j))*
                 ((3*(-Q22.diff(d01,i,j) + 
                        Q22.diff(d01,i,j)*Q22.diff(d10,i,j) - Q22.\
diff(d11,i,j) - Q22.diff(d11,i,j)*Q22(i,j)))/pow(1 + Q22(i,j),2) + 
                   (2*(-Qtt.diff(d01,i,j) + 
                        Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j) - Qtt.d\
iff(d11,i,j) - Qtt.diff(d11,i,j)*Qtt(i,j)))/pow(1 + Qtt(i,j),2))))/
            ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))))/(4.*pow(L,2)) + 
  (htt.diff(d10,i,j)*((complex(0,1)*w*
           ((2*Qr1.diff(d01,i,j))/L + 
             (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
                (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
             (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
                (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))))/
         (-12 + pow(mu,2)) + 
        (2*pow(mu,2)*pow(a0(i,j),2)*exp(B1*mu*h(i,j)))/
         ((-12 + pow(mu,2))*(1 + Qtt(i,j)))) + 
     htt(i,j)*((2*(pow(L,2)*pow(mu,2)*pow(a0(i,j),2)*
              (1 + Q11(i,j))*
              ((exp(B1*mu*h(i,j))*(-Qtt.diff(d10,i,j) - Qtt(i,j)))/
                 ((-12 + pow(mu,2))*(1 + Q11(i,j))*
                   pow(1 + Qtt(i,j),2)) + 
                ((exp(B1*mu*h(i,j))*
                      (-Q11.diff(d10,i,j) - Q11(i,j)))/
                    ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
                   ((-3*(-4 + pow(mu,2))*exp(B1*mu*h(i,j)))/
                       pow(-12 + pow(mu,2),2) + 
                      (4*exp(B1*mu*h(i,j)) + 
                        B1*mu*exp(B1*mu*h(i,j))*
                        (h.diff(d10,i,j) + h(i,j)))/
                       (-12 + pow(mu,2)))/(1 + Q11(i,j)))/
                 (1 + Qtt(i,j))) + 
             (exp(B1*mu*h(i,j))*
                (pow(L,2)*pow(mu,2)*pow(a0(i,j),2)*
                   (Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j))) \
+ (1 + Q11(i,j))*(4*a0.diff(d10,i,j)*pow(L,2)*pow(mu,2)*
                      a0(i,j) + 
                     2*a0.diff(d01,i,j)*L*pow(mu,2)*a0(i,j)*
                      Qr1(i,j)) - 
                  (2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                     (1 + Qtt(i,j)))/(-12 + pow(mu,2))))/
              ((-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j)))))/
         pow(L,2) + (complex(0,3)*w*
           (((2*Qr1.diff(d01,i,j))/L + 
                (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
                   (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
                (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
                   (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)))/3. + 
             (((-Q11.diff(d10,i,j) - Q11(i,j))/
                    pow(1 + Q11(i,j),2) + 2/(1 + Q11(i,j)))*
                 (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
                   (Q11.diff(d01,i,j)*Qr1(i,j))/L) + 
                (-6 + 2*Q11.diff(d10,i,j) - Q11.diff(d20,i,j) + 
                   (Q11.diff(d01,i,j)*Qr1.diff(d10,i,j))/L - 
                   2*Q11(i,j) + (Q11.diff(d11,i,j)*Qr1(i,j))/L)/
                 (1 + Q11(i,j)) + 
                ((-Q22.diff(d10,i,j) - Q22(i,j))/
                    pow(1 + Q22(i,j),2) + 2/(1 + Q22(i,j)))*
                 (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
                   (Q22.diff(d01,i,j)*Qr1(i,j))/L) + 
                (-6 + 2*Q22.diff(d10,i,j) - Q22.diff(d20,i,j) + 
                   (Q22.diff(d01,i,j)*Qr1.diff(d10,i,j))/L - 
                   2*Q22(i,j) + (Q22.diff(d11,i,j)*Qr1(i,j))/L)/
                 (1 + Q22(i,j)) + 
                (2*(Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j) + 
                     (Qr1(i,j)*
                        (-Qrr.diff(d11,i,j) + 
                        Qrr.diff(d10,i,j)*Qtt.diff(d01,i,j) - 
                        Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j) + Qtt\
.diff(d11,i,j) - Qrr.diff(d11,i,j)*Qtt(i,j) + 
                        Qtt.diff(d11,i,j)*Qtt(i,j)))/
                      pow(1 + Qtt(i,j),2)))/L)/3.))/(-12 + pow(mu,2))))/
   2. + (-(h22.diff(d10,i,j)*
        (2 + (complex(0,4)*w)/(-12 + pow(mu,2)) + 
          (16*pow(w,2))/pow(-12 + pow(mu,2),2))) + 
     h22(i,j)*((2*Qr1.diff(d01,i,j))/L + 
        (2*L - Q11.diff(d10,i,j)*L + L*Q11(i,j) + 
           Q11.diff(d01,i,j)*Qr1(i,j))/(L*(1 + Q11(i,j))) + 
        (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
           Q22.diff(d01,i,j)*Qr1(i,j))/(L*(1 + Q22(i,j))) - 
        (2*Qtt.diff(d01,i,j)*Qr1(i,j))/(L*(1 + Qtt(i,j))) - 
        (2*(12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
             Qtt.diff(d10,i,j)*pow(mu,2) + 2*pow(mu,2)*Qtt(i,j)))/
         ((-12 + pow(mu,2))*(1 + Qtt(i,j))) + 
        (complex(0,6)*w*(-(lType(2)/lType(3)) +
             ((2*Qr1.diff(d01,i,j))/L - 
                (-2 + Q11.diff(d10,i,j) - Q11(i,j) - 
                   (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) - 
                (-2 + Q22.diff(d10,i,j) - Q22(i,j) - 
                   (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)) - 
                (2*Qtt.diff(d01,i,j)*Qr1(i,j))/(L*(1 + Qtt(i,j))) - 
                (2*(12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                     Qtt.diff(d10,i,j)*pow(mu,2) + 
                     2*pow(mu,2)*Qtt(i,j)))/
                 ((-12 + pow(mu,2))*(1 + Qtt(i,j))))/3.))/
         (-12 + pow(mu,2)) + 
        (8*(2*pow(L,2)*(-12 + pow(mu,2))*pow(w,2)*
              (1 + 2*Q11(i,j) + pow(Q11(i,j),2))*
              (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*(1 + Qtt(i,j))*
              ((2*(Qtt.diff(d10,i,j) + Qtt(i,j)))/
                 (pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
                   pow(1 + Q22(i,j),2)*pow(1 + Qtt(i,j),2)) + 
                ((((-10/pow(-12 + pow(mu,2),3) + 
                       (9*(-4 + pow(mu,2)))/
                       pow(-12 + pow(mu,2),4))/
                       pow(1 + Q11(i,j),2) + 
                       (2*(Q11.diff(d10,i,j) + Q11(i,j)))/
                       (pow(-12 + pow(mu,2),3)*
                       pow(1 + Q11(i,j),3)))/pow(1 + Q22(i,j),2) \
+ (2*(Q22.diff(d10,i,j) + Q22(i,j)))/
                       (pow(-12 + pow(mu,2),3)*
                        pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),3)))*
                    (1 + Qtt(i,j)) - 
                   (Qrr.diff(d10,i,j) + Qtt(i,j))/
                    (pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
                      pow(1 + Q22(i,j),2)))/pow(1 + Qtt(i,j),2)) - 
             ((Q11.diff(d01,i,j)*(-12 + pow(mu,2))*(1 + Q22(i,j))*
                   (1 + Qtt(i,j))*
                   ((Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Q22(i,j)))/2. + 
                     (Q22.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j)))/2.))/2. + 
                2*pow(L,2)*pow(w,2)*
                 ((((-12 + 3*pow(mu,2) - 10*(-12 + pow(mu,2)))*
                       (1 + 2*Q11(i,j) + pow(Q11(i,j),2)) + 
                       2*(-12 + pow(mu,2))*
                       (Q11.diff(d10,i,j) + Q11(i,j) + 
                       Q11.diff(d10,i,j)*Q11(i,j) + 
                       pow(Q11(i,j),2)))*
                       (1 + 2*Q22(i,j) + pow(Q22(i,j),2)) + 
                      2*(-12 + pow(mu,2))*
                       (1 + 2*Q11(i,j) + pow(Q11(i,j),2))*
                       (Q22.diff(d10,i,j) + Q22(i,j) + 
                        Q22.diff(d10,i,j)*Q22(i,j) + 
                        pow(Q22(i,j),2)))*(1 + Qtt(i,j)) + 
                   (-12 + pow(mu,2))*
                    (1 + 2*Q11(i,j) + pow(Q11(i,j),2))*
                    (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                    (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
                (1 + Q11(i,j))*
                 ((pow(Q22.diff(d01,i,j),2)*
                      (144 - 24*pow(mu,2) + pow(mu,4))*
                      (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/4. - 
                   ((-12 + pow(mu,2))*(1 + Q22(i,j))*(1 + Qtt(i,j))*
                      ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-12 + pow(mu,2)))/2. + 
                        Q22.diff(d02,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j))))/2. + 
                   (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                    (36*pow(Qtt.diff(d01,i,j),2) - 
                      6*pow(Qtt.diff(d01,i,j),2)*pow(mu,2) + 
                      (pow(Qtt.diff(d01,i,j),2)*pow(mu,4))/4. + 
                      (-12 + pow(mu,2))*
                       (6*Qtt.diff(d02,i,j) - 
                        (Qtt.diff(d02,i,j)*pow(mu,2))/2.)*
                       (1 + Qtt(i,j)) - 
                      (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                         (144 - 24*pow(mu,2) + pow(mu,4))*
                         (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/2.)))/
              (pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
                pow(1 + Q22(i,j),2)*(1 + Qtt(i,j)))))/pow(L,2)))/4. + 
  (Dh.diff(d10,i,j)*((complex(0,-4)*L*w*
           (h.diff(d10,i,j)*L*mu + L*mu*h(i,j) + 
             h.diff(d01,i,j)*mu*Qr1(i,j)))/(-12 + pow(mu,2)) + 
        (2*exp((-2*mu*h(i,j))/(3.*A1))*
           ((2 + 3*pow(A1,2))*B1*pow(L,2)*pow(mu,2)*
              pow(a0(i,j),2)*exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j)) \
- 24*A1*pow(L,2)*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
              pow(1 + Qtt(i,j),2)))/
         ((2 + 3*pow(A1,2))*(-12 + pow(mu,2))*(1 + Qtt(i,j)))) + 
     Dh(i,j)*(-2*h.diff(d10,i,j)*pow(L,2)*mu - 
        2*pow(L,2)*mu*h(i,j) - 2*h.diff(d01,i,j)*L*mu*Qr1(i,j) - 
        (complex(0,4)*L*w*(h.diff(d01,i,j)*Qr1.diff(d10,i,j)*mu + 
             4*h.diff(d10,i,j)*L*mu + h.diff(d20,i,j)*L*mu + 
             2*L*mu*h(i,j) + 4*h.diff(d01,i,j)*mu*Qr1(i,j) + 
             h.diff(d11,i,j)*mu*Qr1(i,j)))/(-12 + pow(mu,2)) + 
        (2*((1 + Q11(i,j))*((exp((-2*mu*h(i,j))/(3.*A1))*
                   (-Qtt.diff(d10,i,j) - Qtt(i,j)))/
                 ((-12 + pow(mu,2))*(1 + Q11(i,j))*
                   pow(1 + Qtt(i,j),2)) + 
                ((exp((-2*mu*h(i,j))/(3.*A1))*
                      (-Q11.diff(d10,i,j) - Q11(i,j)))/
                    ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)) + 
                   ((-3*(-4 + pow(mu,2))*
                       exp((-2*mu*h(i,j))/(3.*A1)))/
                       pow(-12 + pow(mu,2),2) + 
                      (5*exp((-2*mu*h(i,j))/(3.*A1)) - 
                        (2*mu*exp((-2*mu*h(i,j))/(3.*A1))*
                       (h.diff(d10,i,j) + h(i,j)))/(3.*A1))/
                       (-12 + pow(mu,2)))/(1 + Q11(i,j)))/
                 (1 + Qtt(i,j)))*
              ((2 + 3*pow(A1,2))*B1*pow(L,2)*pow(mu,2)*
                 pow(a0(i,j),2)*
                 exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j)) - 
                24*A1*pow(L,2)*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
                 pow(1 + Qtt(i,j),2)) + 
             (exp((-2*mu*h(i,j))/(3.*A1))*
                ((-2*pow(a0.diff(d01,i,j),2)*(2 + 3*pow(A1,2))*
                     B1*pow(mu,2)*
                     exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j))*
                     (1 + Qtt(i,j)))/(-12 + pow(mu,2)) + 
                  (Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
                   ((2 + 3*pow(A1,2))*B1*pow(L,2)*pow(mu,2)*
                      pow(a0(i,j),2)*
                      exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j)) - 
                     24*A1*pow(L,2)*
                      (-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
                      pow(1 + Qtt(i,j),2)) + 
                  (1 + Q11(i,j))*
                   (((2 + 3*pow(A1,2))*B1*L*pow(mu,2)*a0(i,j)*
                        exp((2*mu*h(i,j))/(3.*A1) + B1*mu*h(i,j))*
                        (12*a0.diff(d10,i,j)*A1*L + 
                        2*h.diff(d10,i,j)*L*mu*a0(i,j) + 
                        3*h.diff(d10,i,j)*A1*B1*L*mu*a0(i,j) + 
                        2*L*mu*a0(i,j)*h(i,j) + 
                        3*A1*B1*L*mu*a0(i,j)*h(i,j) + 
                        6*a0.diff(d01,i,j)*A1*Qr1(i,j)))/(3.*A1) - 
                     24*A1*pow(L,2)*
                      ((-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
                         (1 + Qtt(i,j))*(Qtt.diff(d10,i,j) + Qtt(i,j)) \
+ (1 + Qtt(i,j))*((-4*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j))) + 
                        (2/(3.*A1) + A1)*mu*
                        exp((2/(3.*A1) + A1)*mu*h(i,j))*
                        (h.diff(d10,i,j) + h(i,j)))*(1 + Qtt(i,j)) + 
                         (-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
                         (Qrr.diff(d10,i,j) + Qtt(i,j)))))))/
              ((-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j)))))/
         (2 + 3*pow(A1,2))))/(2.*pow(L,2)) + 
  (w*(ht1.diff(d10,i,j)*((complex(0,2) - (4*w)/(-12 + pow(mu,2)))*
           Qr1(i,j) + (complex(0,4)*
             (-(Q22.diff(d01,i,j)*(-12 + pow(mu,2))*(1 + Q11(i,j))*
                   (1 + Qtt(i,j)))/2. + 
               (1 + Q22(i,j))*
                (-(Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                     (1 + Q11(i,j))) + 
                  (Q11.diff(d01,i,j)*(-12 + pow(mu,2))*
                     (1 + Qtt(i,j)))/2.)))/
           (L*pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
             (1 + Q22(i,j)))) + 
       ht1(i,j)*((complex(0,2) - (4*w)/(-12 + pow(mu,2)))*
           (Qr1.diff(d10,i,j) + Qr1(i,j)) + 
          (complex(0,2)*pow(Qr1(i,j),2)*
             (Q11.diff(d01,i,j)/(1 + Q11(i,j)) - 
               Qtt.diff(d01,i,j)/(1 + Qtt(i,j))))/L + 
          Qr1(i,j)*((complex(0,2)*Qr1.diff(d01,i,j))/L - 
             (4*w)/(-12 + pow(mu,2)) + 
             (complex(0,2)*(2 - Q11.diff(d10,i,j) + Q11(i,j)))/
              (1 + Q11(i,j)) + 
             (complex(0,2)*(12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
                  Qtt.diff(d10,i,j)*pow(mu,2) + 
                  2*pow(mu,2)*Qtt(i,j)))/
              ((-12 + pow(mu,2))*(1 + Qtt(i,j)))) + 
          (complex(0,4)*((-(Q22.diff(d01,i,j)*(-12 + pow(mu,2))*
                      (1 + Q11(i,j))*(1 + Qtt(i,j)))/2. + 
                  (1 + Q22(i,j))*
                   (-(Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Q11(i,j))) + 
                     (Q11.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j)))/2.))*
                ((-Qtt.diff(d10,i,j) - Qtt(i,j))/
                   (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
                     (1 + Q22(i,j))*(1 + Qtt(i,j))) + 
                  (((-Q22.diff(d10,i,j) - Q22(i,j))/
                        (pow(-12 + pow(mu,2),2)*
                        pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),2)) \
+ ((6/pow(-12 + pow(mu,2),2) - 
                       (6*(-4 + pow(mu,2)))/
                       pow(-12 + pow(mu,2),3))/
                       pow(1 + Q11(i,j),2) - 
                        (2*(Q11.diff(d10,i,j) + Q11(i,j)))/
                        (pow(-12 + pow(mu,2),2)*
                        pow(1 + Q11(i,j),3)))/(1 + Q22(i,j)))*
                      (1 + Qtt(i,j)) + 
                     (Qrr.diff(d10,i,j) + Qtt(i,j))/
                      (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
                        (1 + Q22(i,j))))/(1 + Qtt(i,j))) + 
               ((Q22.diff(d10,i,j) + Q22(i,j) - 2*(1 + Q22(i,j)))*
                   (-(Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Q11(i,j))) + 
                     (Q11.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j)))/2.) + 
                  (-(Q22.diff(d11,i,j)*(-12 + pow(mu,2))*
                        (1 + Q11(i,j))*(1 + Qtt(i,j))) + 
                     Q22.diff(d01,i,j)*
                      (-(((-12 + 3*pow(mu,2) - 
                       5*(-12 + pow(mu,2)))*(1 + Q11(i,j)) + 
                        (-12 + pow(mu,2))*
                        (Q11.diff(d10,i,j) + Q11(i,j)))*
                        (1 + Qtt(i,j))) - 
                        (-12 + pow(mu,2))*(1 + Q11(i,j))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j))))/2. + 
                  (1 + Q22(i,j))*
                   (-(Qtt.diff(d11,i,j)*(-12 + pow(mu,2))*
                        (1 + Q11(i,j))) - 
                     Qtt.diff(d01,i,j)*
                      ((-12 + 3*pow(mu,2) - 3*(-12 + pow(mu,2)))*
                        (1 + Q11(i,j)) + 
                        (-12 + pow(mu,2))*
                         (Q11.diff(d10,i,j) + Q11(i,j))) + 
                     (Q11.diff(d11,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j)) + 
                        Q11.diff(d01,i,j)*
                         ((-12 + 3*pow(mu,2) - 
                        3*(-12 + pow(mu,2)))*(1 + Qtt(i,j)) + 
                         (-12 + pow(mu,2))*
                         (Qtt.diff(d10,i,j) + Qtt(i,j))))/2.))/
                (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
                  (1 + Q22(i,j)))))/L)))/2. + 
  (h11.diff(d10,i,j)*((complex(0,2)*w*
           ((2*Qr1.diff(d01,i,j))/L + 
             (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
                (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
             (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
                (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))))/
         (-12 + pow(mu,2)) + 
        (8*((Q11.diff(d01,i,j)*(-12 + pow(mu,2))*(1 + Q22(i,j))*
                (1 + Qtt(i,j))*
                ((Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                     (1 + Q22(i,j)))/2. + 
                  (Q22.diff(d01,i,j)*(-12 + pow(mu,2))*
                     (1 + Qtt(i,j)))/2.))/2. + 
             (1 + Q11(i,j))*((pow(Q22.diff(d01,i,j),2)*
                   (144 - 24*pow(mu,2) + pow(mu,4))*
                   (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/4. - 
                ((-12 + pow(mu,2))*(1 + Q22(i,j))*(1 + Qtt(i,j))*
                   ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-12 + pow(mu,2)))/2. + 
                     Q22.diff(d02,i,j)*(-12 + pow(mu,2))*
                      (1 + Qtt(i,j))))/2. + 
                (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                 ((pow(Qtt.diff(d01,i,j),2)*
                      (144 - 24*pow(mu,2) + pow(mu,4)))/4. + 
                   (-12 + pow(mu,2))*
                    (6*Qtt.diff(d02,i,j) - 
                      (Qtt.diff(d02,i,j)*pow(mu,2))/2.)*
                    (1 + Qtt(i,j)) - 
                   (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                      (144 - 24*pow(mu,2) + pow(mu,4))*
                      (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/2.))))/
         (pow(L,2)*pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
           pow(1 + Q22(i,j),2)*(1 + Qtt(i,j)))) + 
     h11(i,j)*((complex(0,6)*w*
           (((2*Qr1.diff(d01,i,j))/L + 
                (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
                   (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
                (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
                   (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)))/3. \
+ (((-Q11.diff(d10,i,j) - Q11(i,j))/pow(1 + Q11(i,j),2) + 
                   2/(1 + Q11(i,j)))*
                 (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
                   (Q11.diff(d01,i,j)*Qr1(i,j))/L) + 
                (-6 + 2*Q11.diff(d10,i,j) - Q11.diff(d20,i,j) + 
                   (Q11.diff(d01,i,j)*Qr1.diff(d10,i,j))/L - 
                   2*Q11(i,j) + (Q11.diff(d11,i,j)*Qr1(i,j))/L)/
                 (1 + Q11(i,j)) + 
                ((-Q22.diff(d10,i,j) - Q22(i,j))/
                    pow(1 + Q22(i,j),2) + 2/(1 + Q22(i,j)))*
                 (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
                   (3*Q22.diff(d01,i,j)*Qr1(i,j))/L) + 
                (6 - 2*Q22.diff(d10,i,j) + Q22.diff(d20,i,j) + 
                   (3*Q22.diff(d01,i,j)*Qr1.diff(d10,i,j))/L + 
                   2*Q22(i,j) + (3*Q22.diff(d11,i,j)*Qr1(i,j))/L)/
                 (1 + Q22(i,j)) + 
                (2*(Qr1.diff(d01,i,j) + Qr1.diff(d11,i,j) + 
                     (Qr1(i,j)*
                        (-Qrr.diff(d11,i,j) + 
                        Qrr.diff(d10,i,j)*Qtt.diff(d01,i,j) - 
                        Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j) + Qtt\
.diff(d11,i,j) - Qrr.diff(d11,i,j)*Qtt(i,j) + 
                        Qtt.diff(d11,i,j)*Qtt(i,j)))/
                      pow(1 + Qtt(i,j),2)))/L)/3.))/(-12 + pow(mu,2)) + 
        (8*(((-2*(Qtt.diff(d10,i,j) + Qtt(i,j)))/
                 (pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
                   pow(1 + Q22(i,j),2)*pow(1 + Qtt(i,j),2)) + 
                ((((10/pow(-12 + pow(mu,2),3) - 
                       (9*(-4 + pow(mu,2)))/
                       pow(-12 + pow(mu,2),4))/
                       pow(1 + Q11(i,j),2) - 
                       (2*(Q11.diff(d10,i,j) + Q11(i,j)))/
                       (pow(-12 + pow(mu,2),3)*
                       pow(1 + Q11(i,j),3)))/pow(1 + Q22(i,j),2) \
- (2*(Q22.diff(d10,i,j) + Q22(i,j)))/
                       (pow(-12 + pow(mu,2),3)*
                        pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),3)))*
                    (1 + Qtt(i,j)) + 
                   (Qrr.diff(d10,i,j) + Qtt(i,j))/
                    (pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
                      pow(1 + Q22(i,j),2)))/pow(1 + Qtt(i,j),2))*
              ((Q11.diff(d01,i,j)*(-12 + pow(mu,2))*(1 + Q22(i,j))*
                   (1 + Qtt(i,j))*
                   ((Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Q22(i,j)))/2. + 
                     (Q22.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j)))/2.))/2. + 
                (1 + Q11(i,j))*
                 ((pow(Q22.diff(d01,i,j),2)*
                      (144 - 24*pow(mu,2) + pow(mu,4))*
                      (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/4. - 
                   ((-12 + pow(mu,2))*(1 + Q22(i,j))*(1 + Qtt(i,j))*
                      ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-12 + pow(mu,2)))/2. + 
                        Q22.diff(d02,i,j)*(-12 + pow(mu,2))*
                         (1 + Qtt(i,j))))/2. + 
                   (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                    ((pow(Qtt.diff(d01,i,j),2)*
                         (144 - 24*pow(mu,2) + pow(mu,4)))/4. + 
                      (-12 + pow(mu,2))*
                       (6*Qtt.diff(d02,i,j) - 
                         (Qtt.diff(d02,i,j)*pow(mu,2))/2.)*
                       (1 + Qtt(i,j)) - 
                      (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                         (144 - 24*pow(mu,2) + pow(mu,4))*
                         (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/2.))) + 
             ((Q11.diff(d10,i,j) + Q11(i,j) - 2*(1 + Q11(i,j)))*
                 ((pow(Q22.diff(d01,i,j),2)*
                      (144 - 24*pow(mu,2) + pow(mu,4))*
                      (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/4. - 
                   ((-12 + pow(mu,2))*(1 + Q22(i,j))*(1 + Qtt(i,j))*
                      ((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                        (-12 + pow(mu,2)))/2. + 
                        Q22.diff(d02,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j))))/2. + 
                   (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                    ((pow(Qtt.diff(d01,i,j),2)*
                        (144 - 24*pow(mu,2) + pow(mu,4)))/4. + 
                      (-12 + pow(mu,2))*
                       (6*Qtt.diff(d02,i,j) - 
                        (Qtt.diff(d02,i,j)*pow(mu,2))/2.)*
                       (1 + Qtt(i,j)) - 
                      (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                         (144 - 24*pow(mu,2) + pow(mu,4))*
                         (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/2.)) + 
                (((Qtt.diff(d01,i,j)*(-12 + pow(mu,2))*
                       (1 + Q22(i,j)))/2. + 
                      (Q22.diff(d01,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j)))/2.)*
                    (Q11.diff(d11,i,j)*(-12 + pow(mu,2))*
                       (1 + Q22(i,j))*(1 + Qtt(i,j)) + 
                      Q11.diff(d01,i,j)*
                       (((-12 + 3*pow(mu,2) - 
                       5*(-12 + pow(mu,2)))*(1 + Q22(i,j)) + 
                        (-12 + pow(mu,2))*
                        (Q22.diff(d10,i,j) + Q22(i,j)))*
                        (1 + Qtt(i,j)) + 
                        (-12 + pow(mu,2))*(1 + Q22(i,j))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j)))) + 
                   Q11.diff(d01,i,j)*(-12 + pow(mu,2))*
                    (1 + Q22(i,j))*(1 + Qtt(i,j))*
                    ((Qtt.diff(d11,i,j)*(-12 + pow(mu,2))*
                        (1 + Q22(i,j)) + 
                        Qtt.diff(d01,i,j)*
                        ((-12 + 3*pow(mu,2) - 
                       3*(-12 + pow(mu,2)))*(1 + Q22(i,j)) + 
                        (-12 + pow(mu,2))*
                        (Q22.diff(d10,i,j) + Q22(i,j))))/2. + 
                      (Q22.diff(d11,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j)) + 
                        Q22.diff(d01,i,j)*
                        ((-12 + 3*pow(mu,2) - 
                       3*(-12 + pow(mu,2)))*(1 + Qtt(i,j)) + 
                        (-12 + pow(mu,2))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j))))/2.))/2. + 
                (1 + Q11(i,j))*
                 ((-4*(1 + 2*Q22(i,j) + pow(Q22(i,j),2)) + 
                      2*(Q22.diff(d10,i,j) + Q22(i,j) + 
                         Q22.diff(d10,i,j)*Q22(i,j) + 
                         pow(Q22(i,j),2)))*
                    ((pow(Qtt.diff(d01,i,j),2)*
                        (144 - 24*pow(mu,2) + pow(mu,4)))/4. + 
                      (-12 + pow(mu,2))*
                       (6*Qtt.diff(d02,i,j) - 
                        (Qtt.diff(d02,i,j)*pow(mu,2))/2.)*
                       (1 + Qtt(i,j)) - 
                      (pow(h.diff(d01,i,j),2)*pow(mu,2)*
                         (144 - 24*pow(mu,2) + pow(mu,4))*
                         (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)))/2.) + 
                   (2*Q22.diff(d01,i,j)*Q22.diff(d11,i,j)*
                       (144 - 24*pow(mu,2) + pow(mu,4))*
                       (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)) + 
                      pow(Q22.diff(d01,i,j),2)*
                       ((288 - 96*pow(mu,2) + 6*pow(mu,4) - 
                        6*(144 - 24*pow(mu,2) + pow(mu,4)))*
                        (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)) + 
                         2*(144 - 24*pow(mu,2) + pow(mu,4))*
                         (Qtt.diff(d10,i,j) + Qtt(i,j) + 
                         Qtt.diff(d10,i,j)*Qtt(i,j) + 
                         pow(Qtt(i,j),2))))/4. + 
                   (((Q22.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-12 + pow(mu,2)))/2. + 
                        Q22.diff(d02,i,j)*(-12 + pow(mu,2))*
                        (1 + Qtt(i,j)))*
                       (-(((-12 + 3*pow(mu,2) - 
                       4*(-12 + pow(mu,2)))*(1 + Q22(i,j)) + 
                        (-12 + pow(mu,2))*
                        (Q22.diff(d10,i,j) + Q22(i,j)))*
                        (1 + Qtt(i,j))) - 
                        (-12 + pow(mu,2))*(1 + Q22(i,j))*
                        (Qtt.diff(d10,i,j) + Qtt(i,j))) - 
                      (-12 + pow(mu,2))*(1 + Q22(i,j))*(1 + Qtt(i,j))*
                       ((Q22.diff(d01,i,j)*Qtt.diff(d11,i,j)*
                        (-12 + pow(mu,2)) + 
                        Qtt.diff(d01,i,j)*
                        (Q22.diff(d11,i,j)*(-12 + pow(mu,2)) + 
                        Q22.diff(d01,i,j)*
                        (-12 + 3*pow(mu,2) - 2*(-12 + pow(mu,2))))\
)/2. + Q22.diff(d12,i,j)*(-12 + pow(mu,2))*(1 + Qtt(i,j)) + 
                         Q22.diff(d02,i,j)*
                         ((-12 + 3*pow(mu,2) - 
                        3*(-12 + pow(mu,2)))*(1 + Qtt(i,j)) + 
                         (-12 + pow(mu,2))*
                         (Qtt.diff(d10,i,j) + Qtt(i,j)))))/2. + 
                   (1 + 2*Q22(i,j) + pow(Q22(i,j),2))*
                    ((2*Qtt.diff(d01,i,j)*Qtt.diff(d11,i,j)*
                        (144 - 24*pow(mu,2) + pow(mu,4)) + 
                         pow(Qtt.diff(d01,i,j),2)*
                         (288 - 96*pow(mu,2) + 6*pow(mu,4) - 
                         2*(144 - 24*pow(mu,2) + pow(mu,4))))/4. + 
                      (-12 + pow(mu,2))*
                       (6*Qtt.diff(d12,i,j) - 
                         Qtt.diff(d02,i,j)*pow(mu,2) - 
                         (Qtt.diff(d12,i,j)*pow(mu,2))/2. + 
                         pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                         exp(B1*mu*h(i,j)))*(1 + Qtt(i,j)) + 
                      (6*Qtt.diff(d02,i,j) - 
                         (Qtt.diff(d02,i,j)*pow(mu,2))/2.)*
                       ((-12 + 3*pow(mu,2) - 2*(-12 + pow(mu,2)))*
                         (1 + Qtt(i,j)) + 
                         (-12 + pow(mu,2))*
                         (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
                      (pow(mu,2)*
                         (-2*h.diff(d01,i,j)*h.diff(d11,i,j)*
                         (144 - 24*pow(mu,2) + pow(mu,4))*
                         (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2)) + 
                         pow(h.diff(d01,i,j),2)*
                         (-((288 - 96*pow(mu,2) + 6*pow(mu,4) - 
                        2*(144 - 24*pow(mu,2) + pow(mu,4)))*
                         (1 + 2*Qtt(i,j) + pow(Qtt(i,j),2))) - 
                         2*(144 - 24*pow(mu,2) + pow(mu,4))*
                         (Qtt.diff(d10,i,j) + Qtt(i,j) + 
                         Qtt.diff(d10,i,j)*Qtt(i,j) + pow(Qtt(i,j),2)\
))))/2.)))/(pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
                pow(1 + Q22(i,j),2)*(1 + Qtt(i,j)))))/pow(L,2)))/4.
;
elem[9]=
(htt.diff(d01,i,j)*(-12 + pow(mu,2) + complex(0,4)*w))/
   (2.*pow(-12 + pow(mu,2),2)) + 
  (complex(0,2)*h11.diff(d01,i,j)*w)/pow(-12 + pow(mu,2),2) - 
  (complex(0,4)*h.diff(d01,i,j)*mu*w*Dh(i,j))/pow(-12 + pow(mu,2),2) + 
  (complex(0,1)*L*(-12 + pow(mu,2) + complex(0,2)*w)*w*ht1(i,j))/
   pow(-12 + pow(mu,2),2) + 
  (complex(0,2)*Q22.diff(d01,i,j)*w*h11(i,j))/
   (pow(-12 + pow(mu,2),2)*(1 + Q22(i,j))) - 
  (complex(0,4)*L*mu*muj*w*a0(i,j)*exp(B1*mu*h(i,j)))/
   (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j))) - 
  (complex(0,4)*L*mu*w*a0(i,j)*Ax(i,j)*exp(B1*mu*h(i,j)))/
   (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j)))
;
elem[10]=
h11.diff(d11,i,j)/(-12 + pow(mu,2)) + 
  htt.diff(d11,i,j)/(-12 + pow(mu,2)) - 
  (2*Dh.diff(d10,i,j)*h.diff(d01,i,j)*mu)/(-12 + pow(mu,2)) + 
  (h22.diff(d01,i,j)*(-12 + pow(mu,2) + complex(0,2)*w))/
   pow(-12 + pow(mu,2),2) + 
  (complex(0,1)*ht1.diff(d10,i,j)*L*w)/(-12 + pow(mu,2)) + 
  (h11.diff(d10,i,j)*Q22.diff(d01,i,j))/
   ((-12 + pow(mu,2))*(1 + Q22(i,j))) - 
  (h11.diff(d02,i,j)*Qr1(i,j))/(L*(-12 + pow(mu,2))) - 
  (htt.diff(d02,i,j)*Qr1(i,j))/(L*(-12 + pow(mu,2))) - 
  (complex(0,1)*ht1.diff(d01,i,j)*w*Qr1(i,j))/(-12 + pow(mu,2)) + 
  (2*Dh.diff(d01,i,j)*(-(h.diff(d10,i,j)*L*mu) - L*mu*h(i,j) + 
       2*h.diff(d01,i,j)*mu*Qr1(i,j)))/(L*(-12 + pow(mu,2))) + 
  ((complex(0,4)*(h11.diff(d11,i,j)/(-12 + pow(mu,2)) - 
          (3*h11.diff(d01,i,j)*(-4 + pow(mu,2)))/
           pow(-12 + pow(mu,2),2))*w)/(-12 + pow(mu,2)) + 
     (h11.diff(d01,i,j)*((complex(0,4)*w)/(-12 + pow(mu,2)) + 
          (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
             (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
          (-2 + Q22.diff(d10,i,j) - Q22(i,j) - 
             (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))))/
      (-12 + pow(mu,2)))/2. + 
  (4*A0.diff(d01,i,j)*mu*a0(i,j)*exp(B1*mu*h(i,j)))/
   (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j))) + 
  (4*a0.diff(d01,i,j)*mu*(-12 + pow(mu,2) + complex(0,2)*w)*A0(i,j)*
     exp(B1*mu*h(i,j)))/(pow(-12 + pow(mu,2),3)*(1 + Qtt(i,j))) - 
  (Qtt.diff(d01,i,j)*(-12 + pow(mu,2) + complex(0,2)*w)*h22(i,j))/
   (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j))) + 
  (htt(i,j)*(complex(0,-1)*Qrr.diff(d11,i,j)*w + 
       complex(0,1)*Qrr.diff(d10,i,j)*Qtt.diff(d01,i,j)*w - 
       complex(0,1)*Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j)*w + 
       complex(0,1)*Qtt.diff(d11,i,j)*w + 
       4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j)) - 
       complex(0,1)*Qrr.diff(d11,i,j)*w*Qtt(i,j) + 
       complex(0,1)*Qtt.diff(d11,i,j)*w*Qtt(i,j) + 
       4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j))*Qtt(i,j))\
)/(pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2)) - 
  complex(0,4)*muj*w*((exp(B1*mu*h(i,j))*
        (2*a0.diff(d10,i,j)*L*mu - a0.diff(d01,i,j)*mu*Qr1(i,j)))/
      (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j))) + 
     L*mu*a0(i,j)*((exp(B1*mu*h(i,j))*(-Qtt.diff(d10,i,j) - Qtt(i,j)))/
         (pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2)) + 
        ((-6*(-4 + pow(mu,2))*exp(B1*mu*h(i,j)))/
            pow(-12 + pow(mu,2),3) + 
           (2*exp(B1*mu*h(i,j)) + 
              B1*mu*exp(B1*mu*h(i,j))*(h.diff(d10,i,j) + h(i,j)))/
            pow(-12 + pow(mu,2),2))/(1 + Qtt(i,j)))) - 
  complex(0,4)*w*((Ax(i,j)*exp(B1*mu*h(i,j))*
        (2*a0.diff(d10,i,j)*L*mu - a0.diff(d01,i,j)*mu*Qr1(i,j)))/
      (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j))) + 
     L*mu*a0(i,j)*((Ax(i,j)*exp(B1*mu*h(i,j))*
           (-Qtt.diff(d10,i,j) - Qtt(i,j)))/
         (pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2)) + 
        ((Ax.diff(d10,i,j)*exp(B1*mu*h(i,j)))/
            pow(-12 + pow(mu,2),2) + 
           Ax(i,j)*((-6*(-4 + pow(mu,2))*exp(B1*mu*h(i,j)))/
               pow(-12 + pow(mu,2),3) + 
              (2*exp(B1*mu*h(i,j)) + 
                 B1*mu*exp(B1*mu*h(i,j))*(h.diff(d10,i,j) + h(i,j)))/
               pow(-12 + pow(mu,2),2)))/(1 + Qtt(i,j)))) + 
  (complex(0,3)*w*((2*Q22.diff(d01,i,j)*
          (h11.diff(d10,i,j)/(3.*(-12 + pow(mu,2))) + 
            (1/(3.*(-12 + pow(mu,2))) - 
               (-4 + pow(mu,2))/pow(-12 + pow(mu,2),2))*h11(i,j)))/
        (1 + Q22(i,j)) + (h11(i,j)*
          ((-2*(-Q22.diff(d01,i,j) + 
                 Q22.diff(d01,i,j)*Q22.diff(d10,i,j) - Q22.diff(d1\
1,i,j) - Q22.diff(d11,i,j)*Q22(i,j)))/pow(1 + Q22(i,j),2) + 
            (-Qrr.diff(d11,i,j) - Qtt.diff(d01,i,j) + 
               Qrr.diff(d10,i,j)*Qtt.diff(d01,i,j) - 
               Qrr.diff(d11,i,j)*Qtt(i,j))/pow(1 + Qtt(i,j),2) + 
            (Qtt.diff(d01,i,j) - 
               Qtt.diff(d01,i,j)*Qtt.diff(d10,i,j) + Qtt.diff(d11,i\
,j) + Qtt.diff(d11,i,j)*Qtt(i,j))/pow(1 + Qtt(i,j),2)))/
        (3.*(-12 + pow(mu,2)))))/(-12 + pow(mu,2)) + 
  ((htt.diff(d11,i,j)/(-12 + pow(mu,2)) - 
        (3*htt.diff(d01,i,j)*(-4 + pow(mu,2)))/
         pow(-12 + pow(mu,2),2))*
      (1 + (complex(0,4)*w)/(-12 + pow(mu,2))) + 
     (htt.diff(d01,i,j)*((complex(0,4)*w)/(-12 + pow(mu,2)) + 
          (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
             (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) - 
          (Qtt.diff(d01,i,j)*Qr1(i,j))/(L*(1 + Qtt(i,j))) + 
          (12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
             Qtt.diff(d10,i,j)*pow(mu,2) + 2*pow(mu,2)*Qtt(i,j))/
           ((-12 + pow(mu,2))*(1 + Qtt(i,j)))))/(-12 + pow(mu,2)))/2. + 
  w*((complex(0,1)*L - (2*L*w)/(-12 + pow(mu,2)))*
      (ht1.diff(d10,i,j)/(-12 + pow(mu,2)) - 
        (3*(-4 + pow(mu,2))*ht1(i,j))/pow(-12 + pow(mu,2),2)) + 
     (ht1(i,j)*(complex(0,1)*Qr1.diff(d01,i,j) - 
          (2*L*w)/(-12 + pow(mu,2)) + 
          (complex(0,1)*(2*L - Q11.diff(d10,i,j)*L + L*Q11(i,j) + 
               Q11.diff(d01,i,j)*Qr1(i,j)))/(1 + Q11(i,j)) - 
          (complex(0,1)*Qtt.diff(d01,i,j)*Qr1(i,j))/(1 + Qtt(i,j)) + 
          (complex(0,1)*L*(12 - 12*Qtt.diff(d10,i,j) + pow(mu,2) + 
               Qtt.diff(d10,i,j)*pow(mu,2) + 2*pow(mu,2)*Qtt(i,j)))/
           ((-12 + pow(mu,2))*(1 + Qtt(i,j)))))/(-12 + pow(mu,2))) - 
  (4*(complex(0,3)*h.diff(d01,i,j)*L*mu*(-12 + pow(mu,2))*w*
        (1 + Qtt(i,j))*((Dh(i,j)*(-Qtt.diff(d10,i,j) - Qtt(i,j)))/
           (3.*pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2)) + 
          (Dh.diff(d10,i,j)/(3.*pow(-12 + pow(mu,2),2)) + 
             (1/(3.*pow(-12 + pow(mu,2),2)) - 
                (2*(-4 + pow(mu,2)))/pow(-12 + pow(mu,2),3))*Dh(i,j)\
)/(1 + Qtt(i,j))) + (L*Dh(i,j)*
          (-3*a0.diff(d01,i,j)*B1*pow(mu,2)*(-12 + pow(mu,2))*
             a0(i,j)*exp(B1*mu*h(i,j)) + 
            (mu*(complex(0,6)*h.diff(d11,i,j)*(-12 + pow(mu,2))*w*
                  (1 + Qtt(i,j)) + 
                 h.diff(d01,i,j)*
                  (((-12 + pow(mu,2))*
                        (-36 + 3*pow(mu,2) + complex(0,18)*w) + 
                       complex(0,12)*pow(mu,2)*w)*(1 + Qtt(i,j)) + 
                    complex(0,6)*(-12 + pow(mu,2))*w*
                     (Qtt.diff(d10,i,j) + Qtt(i,j)))))/2.))/
        (3.*pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j)))))/
   (L*(-12 + pow(mu,2)))
;
