
elem[0]=
(((1 + rgrid[i]*Q11(i,j))*(2*(-1 + rgrid[i])*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*
          (rgrid[i]*(L*(Q11.diff(d10,i,j)*rgrid[i] + Q11(i,j)) - 
               2*Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)) + 
            2*L*(1 + rgrid[i]*Qrr(i,j)))*(1 + rgrid[i]*Qtt(i,j)) + 
         (-1 + rgrid[i])*pow(rgrid[i],2)*
          (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*
          pow(1 + rgrid[i]*Q11(i,j),2)*Qr1(i,j)*
          (2*L*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
               pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j)*
             (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
            (1 + rgrid[i]*Q22(i,j))*
             (2*Qr1.diff(d10,i,j)*L*rgrid[i]*
                (-4 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                  pow(mu,2)*pow(rgrid[i],4))*
                (1 + rgrid[i]*Qtt(i,j)) + 
               L*pow(rgrid[i],2)*
                pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                  pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),3)*
                (1 + rgrid[i]*Qtt(i,j)) + 
               Qr1(i,j)*(L*(8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                     2*pow(mu,2)*pow(rgrid[i],4))*
                   (1 + rgrid[i]*Qrr(i,j)) + 
                  pow(rgrid[i],2)*
                   (L*(12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
                     2*Qr1.diff(d01,i,j)*
                      (4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                        pow(mu,2)*pow(rgrid[i],4)))*
                   (1 + rgrid[i]*Qtt(i,j))))) - 
         2*(1 + rgrid[i]*Q11(i,j))*
          (-(L*(-1 + rgrid[i])*
               (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                 pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))*
               (rgrid[i]*(Q22.diff(d10,i,j)*rgrid[i] + Q22(i,j)) + 
                 2*(1 + rgrid[i]*Qrr(i,j)))*(1 + rgrid[i]*Qtt(i,j))) + 
            (1 + rgrid[i]*Q22(i,j))*
             (L*(-8 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                  2*pow(mu,2)*pow(rgrid[i],4))*
                pow(1 + rgrid[i]*Qrr(i,j),2) - 
               (-1 + rgrid[i])*rgrid[i]*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*
                (2*Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j) - 
                  L*(Qrr.diff(d10,i,j)*rgrid[i] + Qrr(i,j)))*
                (1 + rgrid[i]*Qtt(i,j)) - 
               (1 + rgrid[i]*Qrr(i,j))*
                (L*rgrid[i]*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                     pow(mu,2)*pow(rgrid[i],4))*
                   (Qtt.diff(d10,i,j)*rgrid[i] + Qtt(i,j)) + 
                  (L*(-24 + 3*(4 + pow(mu,2))*pow(rgrid[i],3) - 
                        2*pow(mu,2)*pow(rgrid[i],4)) + 
                     Qr1.diff(d01,i,j)*
                      (-8*pow(rgrid[i],2) + 
                        2*(4 + pow(mu,2))*pow(rgrid[i],5) - 
                        2*pow(mu,2)*pow(rgrid[i],6)) + 
                     2*L*pow(rgrid[i],2)*
                      pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                        pow(mu,2)*pow(rgrid[i],4),2)*
                      pow(Qr1(i,j),2))*(1 + rgrid[i]*Qtt(i,j))))))*
       ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
            pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Q22(i,j))*
          (1 + rgrid[i]*Qrr(i,j))*
          (rgrid[i]*(L*(Q11.diff(d10,i,j)*rgrid[i] + Q11(i,j)) - 
               Q11.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)) + 
            2*L*(1 + rgrid[i]*Qrr(i,j)))*(1 + rgrid[i]*Qtt(i,j)) + 
         (1 + rgrid[i]*Q11(i,j))*
          ((-1 + rgrid[i])*(-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
               pow(mu,2)*pow(rgrid[i],3))*(1 + rgrid[i]*Qrr(i,j))*
             (rgrid[i]*(L*(Q22.diff(d10,i,j)*rgrid[i] + Q22(i,j)) - 
                  Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j)) + 
               2*L*(1 + rgrid[i]*Qrr(i,j)))*(1 + rgrid[i]*Qtt(i,j)) + 
            (1 + rgrid[i]*Q22(i,j))*
             (L*(8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                  2*pow(mu,2)*pow(rgrid[i],4))*
                pow(1 + rgrid[i]*Qrr(i,j),2) + 
               (-1 + rgrid[i])*rgrid[i]*
                (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                  pow(mu,2)*pow(rgrid[i],3))*
                (Qrr.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j) - 
                  L*(Qrr.diff(d10,i,j)*rgrid[i] + Qrr(i,j)))*
                (1 + rgrid[i]*Qtt(i,j)) + 
               (1 + rgrid[i]*Qrr(i,j))*
                ((L*(-24 + 3*(4 + pow(mu,2))*pow(rgrid[i],3) - 
                        2*pow(mu,2)*pow(rgrid[i],4)) + 
                     Qr1.diff(d01,i,j)*
                      (-8*pow(rgrid[i],2) + 
                        2*(4 + pow(mu,2))*pow(rgrid[i],5) - 
                        2*pow(mu,2)*pow(rgrid[i],6)) + 
                     L*pow(rgrid[i],2)*
                      pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                        pow(mu,2)*pow(rgrid[i],4),2)*
                      pow(Qr1(i,j),2))*(1 + rgrid[i]*Qtt(i,j)) - 
                  (-1 + rgrid[i])*rgrid[i]*
                   (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                     pow(mu,2)*pow(rgrid[i],3))*
                   (Qtt.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j) - 
                     L*(Qtt.diff(d10,i,j)*rgrid[i] + Qtt(i,j))))))))/
     (4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
       pow(mu,2)*pow(rgrid[i],4)) + 
    pow(rgrid[i],2)*(-2*Q11.diff(d01,i,j)*rgrid[i]*
        (1 + rgrid[i]*Q22(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (1 + rgrid[i]*Qtt(i,j)) + 
       2*(1 + rgrid[i]*Q11(i,j))*
        (Q22.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j))*
           (1 + rgrid[i]*Qtt(i,j)) + 
          (1 + rgrid[i]*Q22(i,j))*
           (Qtt.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)) + 
             Qrr.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qtt(i,j)) + 
             L*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j)*
              (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)))) + 
       pow(1 + rgrid[i]*Q11(i,j),2)*
        (2*L*(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
             pow(mu,2)*pow(rgrid[i],4))*Qr1(i,j)*
           (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
          (1 + rgrid[i]*Q22(i,j))*
           (2*Qr1.diff(d10,i,j)*L*rgrid[i]*
              (-4 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Qtt(i,j)) + 
             L*pow(rgrid[i],2)*
              pow(4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                pow(mu,2)*pow(rgrid[i],4),2)*pow(Qr1(i,j),3)*
              (1 + rgrid[i]*Qtt(i,j)) + 
             Qr1(i,j)*(L*(8 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                   2*pow(mu,2)*pow(rgrid[i],4))*
                 (1 + rgrid[i]*Qrr(i,j)) + 
                pow(rgrid[i],2)*
                 (L*(12 + pow(mu,2)*(3 - 4*rgrid[i]))*rgrid[i] + 
                   2*Qr1.diff(d01,i,j)*
                    (4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                      pow(mu,2)*pow(rgrid[i],4)))*
                 (1 + rgrid[i]*Qtt(i,j))))))*
     (-2*Q11.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Q22(i,j))*
        pow(1 + rgrid[i]*Qrr(i,j),2)*(1 + rgrid[i]*Qtt(i,j)) + 
       (1 + rgrid[i]*Q11(i,j))*(1 + rgrid[i]*Qrr(i,j))*
        (2*Q22.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j))*
           (1 + rgrid[i]*Qtt(i,j)) + 
          (1 + rgrid[i]*Q22(i,j))*
           (L*rgrid[i]*(-4 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                pow(mu,2)*pow(rgrid[i],4))*
              (Q11.diff(d10,i,j)*rgrid[i] + Q11(i,j))*Qr1(i,j)*
              (1 + rgrid[i]*Qtt(i,j)) + 
             Q11.diff(d01,i,j)*pow(rgrid[i],3)*
              (4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                pow(mu,2)*pow(rgrid[i],4))*pow(Qr1(i,j),2)*
              (1 + rgrid[i]*Qtt(i,j)) + 
             2*(Qtt.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j)) + 
                Qrr.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qtt(i,j))))) + 
       pow(1 + rgrid[i]*Q11(i,j),2)*
        ((-1 + rgrid[i])*rgrid[i]*
           (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
             pow(mu,2)*pow(rgrid[i],3))*Qr1(i,j)*
           (-(L*(Q22.diff(d10,i,j)*rgrid[i] + Q22(i,j))) + 
             Q22.diff(d01,i,j)*pow(rgrid[i],2)*Qr1(i,j))*
           (1 + rgrid[i]*Qrr(i,j))*(1 + rgrid[i]*Qtt(i,j)) + 
          (1 + rgrid[i]*Q22(i,j))*
           (2*Qr1.diff(d10,i,j)*L*rgrid[i]*
              (-4 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                pow(mu,2)*pow(rgrid[i],4))*(1 + rgrid[i]*Qrr(i,j))*
              (1 + rgrid[i]*Qtt(i,j)) - 
             (-1 + rgrid[i])*pow(rgrid[i],2)*
              (-4 - 4*rgrid[i] - 4*pow(rgrid[i],2) + 
                pow(mu,2)*pow(rgrid[i],3))*pow(Qr1(i,j),2)*
              (-(Qtt.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qrr(i,j))) + 
                Qrr.diff(d01,i,j)*rgrid[i]*(1 + rgrid[i]*Qtt(i,j))) + 
             Qr1(i,j)*(L*rgrid[i]*
                 (4 - (4 + pow(mu,2))*pow(rgrid[i],3) + 
                   pow(mu,2)*pow(rgrid[i],4))*
                 (Qrr.diff(d10,i,j)*rgrid[i] + Qrr(i,j))*
                 (1 + rgrid[i]*Qtt(i,j)) + 
                (1 + rgrid[i]*Qrr(i,j))*
                 (L*rgrid[i]*(-4 + (4 + pow(mu,2))*pow(rgrid[i],3) - 
                      pow(mu,2)*pow(rgrid[i],4))*
                    (Qtt.diff(d10,i,j)*rgrid[i] + Qtt(i,j)) - 
                   2*(L*(-12 + pow(mu,2)*pow(rgrid[i],4)) + 
                      Qr1.diff(d01,i,j)*
                       (-8*pow(rgrid[i],2) + 
                         2*(4 + pow(mu,2))*pow(rgrid[i],5) - 
                         2*pow(mu,2)*pow(rgrid[i],6)))*
                    (1 + rgrid[i]*Qtt(i,j))))))))/
  (16.*pow(L,2)*pow(1 + rgrid[i]*Q11(i,j),3)*
    pow(1 + rgrid[i]*Q22(i,j),2)*pow(1 + rgrid[i]*Qrr(i,j),3)*
    pow(1 + rgrid[i]*Qtt(i,j),2))
;
