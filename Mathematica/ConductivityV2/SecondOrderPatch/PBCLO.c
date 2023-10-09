
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
((complex(0,-1)*pow(mu,2) + 2*(complex(0,6) + w))*(1 + Qtt(i,j))*
    D[-1 + d11](i,j,i1,j1))/(2.*L*w*(1 + Q11(i,j)))
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
((-12 + pow(mu,2) + complex(0,2)*w)*(1 + Qtt(i,j))*D[-1 + d11](i,j,i1,j1))/
  (4.*L*(1 + Q11(i,j)))
;
elem[6][6]=
0
;
elem[7][0]=
0
;
elem[7][1]=
0
;
elem[7][2]=
0
;
elem[7][3]=
0
;
elem[7][4]=
0
;
elem[7][5]=
0
;
elem[7][6]=
0
;
elem[8][0]=
0
;
elem[8][1]=
((4*L*Qr1(i,j) + (complex(0,8)*L*w*Qr1(i,j))/(-12 + pow(mu,2)) + 
       (2*(1 + Qtt(i,j))*(Q11.diff(d01,i,j) + 
            (1 + Q11(i,j))*((-3*Q22.diff(d01,i,j))/(1 + Q22(i,j)) - 
               (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))))/
        ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))*D[-1 + d11](i,j,i1,j1)\
)/(4.*pow(L,2)) - ((1 + Qtt(i,j))*D[-1 + d12](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[8][2]=
(2*h.diff(d01,i,j)*mu*(1 + Qtt(i,j))*D[-1 + d11](i,j,i1,j1))/
  (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[8][3]=
((4*L*Qr1(i,j) + L*(2 + (complex(0,8)*w)/(-12 + pow(mu,2)))*Qr1(i,j) + 
       (2*(1 + Qtt(i,j))*(Q11.diff(d01,i,j) + 
            (1 + Q11(i,j))*(-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) - 
               (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))))/
        ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))*D[-1 + d11](i,j,i1,j1)\
)/(4.*pow(L,2)) - ((1 + Qtt(i,j))*D[-1 + d12](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[8][4]=
0
;
elem[8][5]=
(complex(0,-2)*w*(1 + Qtt(i,j))*D[-1 + d11](i,j,i1,j1))/
  (L*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[8][6]=
0
;
elem[9][0]=
0
;
elem[9][1]=
0
;
elem[9][2]=
0
;
elem[9][3]=
0
;
elem[9][4]=
0
;
elem[9][5]=
0
;
elem[9][6]=
0
;
elem[10][0]=
0
;
elem[10][1]=
((-12 + pow(mu,2) + complex(0,2)*w)*D[-1 + d11](i,j,i1,j1))/
  pow(-12 + pow(mu,2),2)
;
elem[10][2]=
0
;
elem[10][3]=
((-36 + 3*pow(mu,2) + complex(0,4)*w)*D[-1 + d11](i,j,i1,j1))/
  (2.*pow(-12 + pow(mu,2),2))
;
elem[10][4]=
0
;
elem[10][5]=
0
;
elem[10][6]=
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
((complex(0,-1)*pow(mu,2) + 2*(complex(0,6) + w))*
    (Q22.diff(d01,i,j)*(1 + Q11(i,j)) - 
      Q11.diff(d01,i,j)*(1 + Q22(i,j)) + 
      2*h.diff(d01,i,j)*B1*mu*(1 + Q11(i,j))*(1 + Q22(i,j)))*
    (1 + Qtt(i,j))*D[-1 + d10](i,j,i1,j1))/
  (4.*L*w*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[1][1]+=
0
;
elem[1][2]+=
-(B1*mu*a0(i,j)*D[-1 + d10](i,j,i1,j1))
;
elem[1][3]+=
-(mu*a0(i,j)*D[-1 + d10](i,j,i1,j1))/2.
;
elem[1][4]+=
0
;
elem[1][5]+=
0
;
elem[1][6]+=
(-2*(-12 + pow(mu,2) + complex(0,1)*w)*D[-1 + d10](i,j,i1,j1))/
  (-12 + pow(mu,2))
;
elem[2][0]+=
(1 + (complex(0,4)*w)/(-12 + pow(mu,2)))*D[-1 + d10](i,j,i1,j1)
;
elem[2][1]+=
0
;
elem[2][2]+=
0
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
((-12 + pow(mu,2) + complex(0,4)*w)*D[-1 + d10](i,j,i1,j1))/
  (2.*(1 + Qtt(i,j)))
;
elem[3][3]+=
0
;
elem[3][4]+=
0
;
elem[3][5]+=
0
;
elem[3][6]+=
0
;
elem[4][0]+=
0
;
elem[4][1]+=
(1 + (complex(0,4)*w)/(-12 + pow(mu,2)))*D[-1 + d10](i,j,i1,j1)
;
elem[4][2]+=
0
;
elem[4][3]+=
0
;
elem[4][4]+=
0
;
elem[4][5]+=
0
;
elem[4][6]+=
0
;
elem[5][0]+=
0
;
elem[5][1]+=
0
;
elem[5][2]+=
0
;
elem[5][3]+=
0
;
elem[5][4]+=
0
;
elem[5][5]+=
0
;
elem[5][6]+=
0
;
elem[6][0]+=
-((a0.diff(d01,i,j)*mu*(-12 + pow(mu,2) + complex(0,2)*w)*
      exp(B1*mu*h(i,j))*D[-1 + d10](i,j,i1,j1))/
    (L*(-12 + pow(mu,2))*(1 + Q11(i,j))))
;
elem[6][1]+=
(complex(0,-0.25)*w*((1 + Q22(i,j))*
       (2*Qr1.diff(d01,i,j) + 2*L - Q11.diff(d10,i,j)*L + 
         (2*Qr1.diff(d01,i,j) + L)*Q11(i,j) + 
         Q11.diff(d01,i,j)*Qr1(i,j)) - 
      (1 + Q11(i,j))*(-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) + 
         Q22.diff(d01,i,j)*Qr1(i,j)))*D[-1 + d10](i,j,i1,j1))/
  (L*(1 + Q11(i,j))*(1 + Q22(i,j)))
;
elem[6][2]+=
(complex(0,1)*mu*w*(h.diff(d10,i,j)*L + L*h(i,j) - 
      h.diff(d01,i,j)*Qr1(i,j))*D[-1 + d10](i,j,i1,j1))/L
;
elem[6][3]+=
0
;
elem[6][4]+=
(w*(complex(0,-3) + (4*w)/(-12 + pow(mu,2)))*D[-1 + d10](i,j,i1,j1))/2.
;
elem[6][5]+=
((-12 + pow(mu,2) + complex(0,2)*w)*
    (-(Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j))) + 
      (1 + Q11(i,j))*(2*Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
         Q22.diff(d01,i,j)*(1 + Qtt(i,j))))*D[-1 + d10](i,j,i1,j1))/
  (8.*L*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[6][6]+=
0
;
elem[7][0]+=
0
;
elem[7][1]+=
0
;
elem[7][2]+=
0
;
elem[7][3]+=
0
;
elem[7][4]+=
0
;
elem[7][5]+=
0
;
elem[7][6]+=
0
;
elem[8][0]+=
(complex(0,-4)*mu*w*exp(B1*mu*h(i,j))*
    (L*(-12 + pow(mu,2))*a0(i,j)*(1 + Q11(i,j))*Qr1(i,j) - 
      a0.diff(d01,i,j)*(1 + Qtt(i,j)))*D[-1 + d10](i,j,i1,j1))/
  (L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*(1 + Qtt(i,j)))
;
elem[8][1]+=
(((2*Qr1.diff(d01,i,j))/L + (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
         (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
      (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
         (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)) + 
      (complex(0,2)*w*((2*Qr1.diff(d01,i,j))/L + 
           (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
              (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
           (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
              (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))))/
       (-12 + pow(mu,2)) + (2*
         (Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j))*
            (Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
              Q22.diff(d01,i,j)*(1 + Qtt(i,j))) + 
           (1 + Q11(i,j))*(pow(Q22.diff(d01,i,j),2)*
               pow(1 + Qtt(i,j),2) - 
              (1 + Q22(i,j))*(1 + Qtt(i,j))*
               (Q22.diff(d01,i,j)*Qtt.diff(d01,i,j) + 
                 2*Q22.diff(d02,i,j)*(1 + Qtt(i,j))) - 
              pow(1 + Q22(i,j),2)*
               (-pow(Qtt.diff(d01,i,j),2) + 
                 2*Qtt.diff(d02,i,j)*(1 + Qtt(i,j)) + 
                 2*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                  pow(1 + Qtt(i,j),2)))))/
       (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*
         pow(1 + Q22(i,j),2)*(1 + Qtt(i,j))))*D[-1 + d10](i,j,i1,j1))/4.
;
elem[8][2]+=
(-(h.diff(d10,i,j)*mu) - mu*h(i,j) - (h.diff(d01,i,j)*mu*Qr1(i,j))/L + 
    (complex(0,-2)*mu*w*(h.diff(d10,i,j)*L + L*h(i,j) + 
          h.diff(d01,i,j)*Qr1(i,j)) + 
       (L*exp((-2*mu*h(i,j))/(3.*A1))*
          ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
             exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) - 
            24*A1*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
             pow(1 + Qtt(i,j),2)))/((2 + 3*pow(A1,2))*(1 + Qtt(i,j))))/
     (L*(-12 + pow(mu,2))))*D[-1 + d10](i,j,i1,j1)
;
elem[8][3]+=
(((2*Qr1.diff(d01,i,j))/L + (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
         (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
      (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
         (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)) + 
      (2*(complex(0,1)*w*((2*Qr1.diff(d01,i,j))/L + 
              (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
                 (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
              (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
                 (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))) + 
           (2*pow(mu,2)*pow(a0(i,j),2)*exp(B1*mu*h(i,j)))/
            (1 + Qtt(i,j))))/(-12 + pow(mu,2)))*D[-1 + d10](i,j,i1,j1))/4.
;
elem[8][4]+=
(-1 - (complex(0,1)*w)/(-12 + pow(mu,2)) - 
    (4*pow(w,2))/pow(-12 + pow(mu,2),2))*D[-1 + d10](i,j,i1,j1)
;
elem[8][5]+=
(w*(complex(0,2)*Qr1(i,j) + (complex(0,2) - (4*w)/(-12 + pow(mu,2)))*
       Qr1(i,j) - (complex(0,2)*
         (-(Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j))) + 
           (1 + Q11(i,j))*(2*Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
              Q22.diff(d01,i,j)*(1 + Qtt(i,j)))))/
       (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))))*
    D[-1 + d10](i,j,i1,j1))/2.
;
elem[8][6]+=
(4*mu*(-12 + pow(mu,2) + complex(0,1)*w)*a0(i,j)*exp(B1*mu*h(i,j))*
    D[-1 + d10](i,j,i1,j1))/(pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j)))
;
elem[9][0]+=
0
;
elem[9][1]+=
0
;
elem[9][2]+=
0
;
elem[9][3]+=
0
;
elem[9][4]+=
0
;
elem[9][5]+=
0
;
elem[9][6]+=
0
;
elem[10][0]+=
(complex(0,-4)*L*mu*w*a0(i,j)*exp(B1*mu*h(i,j))*D[-1 + d10](i,j,i1,j1))/
  (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j)))
;
elem[10][1]+=
(Q22.diff(d01,i,j)*(-12 + pow(mu,2) + complex(0,2)*w)*
    D[-1 + d10](i,j,i1,j1))/(pow(-12 + pow(mu,2),2)*(1 + Q22(i,j)))
;
elem[10][2]+=
(-2*h.diff(d01,i,j)*mu*(-12 + pow(mu,2) + complex(0,2)*w)*
    D[-1 + d10](i,j,i1,j1))/pow(-12 + pow(mu,2),2)
;
elem[10][3]+=
0
;
elem[10][4]+=
0
;
elem[10][5]+=
(complex(0,2)*L*(-12 + pow(mu,2) + complex(0,1)*w)*w*
    D[-1 + d10](i,j,i1,j1))/pow(-12 + pow(mu,2),2)
;
elem[10][6]+=
0
;
}
if(i==i1){
elem[0][0]+=
((1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/(L + L*Q11(i,j))
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
((4*L*(2 - Q11.diff(d10,i,j) + Q11(i,j))*(1 + Qtt(i,j)) - 
       (complex(0,1)*Q11.diff(d01,i,j)*(-12 + pow(mu,2))*Qr1(i,j)*
          (1 + Qtt(i,j)))/w + 
       (1 + Q11(i,j))*((complex(0,1)*Qtt.diff(d01,i,j)*
             (-12 + pow(mu,2))*Qr1(i,j))/w + 
          (4*L*(12 + pow(mu,2) + 
               Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
               2*pow(mu,2)*Qtt(i,j)))/(-12 + pow(mu,2)) + 
          (-12 + pow(mu,2))*(1 + Qtt(i,j))*
           ((4*L)/(-12 + pow(mu,2)) + 
             (complex(0,4)*Qr1.diff(d01,i,j))/w + 
             (complex(0,1)*Qr1(i,j)*
                (-(Qtt.diff(d01,i,j)*(1 + Q22(i,j))) + 
                  Q22.diff(d01,i,j)*(1 + Qtt(i,j)) + 
                  2*h.diff(d01,i,j)*B1*mu*(1 + Q22(i,j))*(1 + Qtt(i,j))\
))/(w*(1 + Q22(i,j))*(1 + Qtt(i,j))))))*D[-1 + d01](i,j,i1,j1))/
   (4.*pow(L,2)*pow(1 + Q11(i,j),2)) + 
  (complex(0,0.5)*(-12 + pow(mu,2))*Qr1(i,j)*(1 + Qtt(i,j))*
     D[-1 + d02](i,j,i1,j1))/(pow(L,2)*w*(1 + Q11(i,j)))
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
(complex(0,0.5)*mu*(-12 + pow(mu,2))*a0(i,j)*(1 + Qtt(i,j))*
    D[-1 + d01](i,j,i1,j1))/(L*w*(1 + Q11(i,j)))
;
elem[1][6]+=
(Qr1(i,j)*D[-1 + d01](i,j,i1,j1))/L
;
elem[2][0]+=
-(((-12 + pow(mu,2) + complex(0,4)*w)*Qr1(i,j)*D[-1 + d01](i,j,i1,j1))/
    (L*(-12 + pow(mu,2))))
;
elem[2][1]+=
0
;
elem[2][2]+=
0
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
(complex(0,4)*w*D[-1 + d01](i,j,i1,j1))/(L*pow(-12 + pow(mu,2),2))
;
elem[3][0]+=
0
;
elem[3][1]+=
-((h.diff(d01,i,j)*mu*D[-1 + d01](i,j,i1,j1))/(pow(L,2)*(1 + Q11(i,j))))
;
elem[3][2]+=
((-(L*(-12 + pow(mu,2) + complex(0,4)*w)*Qr1(i,j)) + 
       (Q22.diff(d01,i,j)*(1 + Q11(i,j))*(1 + Qtt(i,j)) - 
          (1 + Q22(i,j))*(-2*Qtt.diff(d01,i,j)*(1 + Q11(i,j)) + 
             Q11.diff(d01,i,j)*(1 + Qtt(i,j))))/
        (pow(1 + Q11(i,j),2)*(1 + Q22(i,j))))*D[-1 + d01](i,j,i1,j1))/
   (2.*pow(L,2)*(1 + Qtt(i,j))) + 
  D[-1 + d02](i,j,i1,j1)/(pow(L,2)*(1 + Q11(i,j)))
;
elem[3][3]+=
-(h.diff(d01,i,j)*mu*D[-1 + d01](i,j,i1,j1))/(2.*pow(L,2)*(1 + Q11(i,j)))
;
elem[3][4]+=
0
;
elem[3][5]+=
0
;
elem[3][6]+=
0
;
elem[4][0]+=
0
;
elem[4][1]+=
-(((-12 + pow(mu,2) + complex(0,4)*w)*Qr1(i,j)*D[-1 + d01](i,j,i1,j1))/
    (L*(-12 + pow(mu,2))))
;
elem[4][2]+=
(4*h.diff(d01,i,j)*mu*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[4][3]+=
((Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j)) + 
       (1 + Q11(i,j))*(-2*Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
          Q22.diff(d01,i,j)*(1 + Qtt(i,j))))*D[-1 + d01](i,j,i1,j1))/
   (2.*pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))) \
- ((1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[4][4]+=
0
;
elem[4][5]+=
(complex(0,-2)*w*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
  (L*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[4][6]+=
0
;
elem[5][0]+=
0
;
elem[5][1]+=
0
;
elem[5][2]+=
0
;
elem[5][3]+=
0
;
elem[5][4]+=
0
;
elem[5][5]+=
(complex(0,0.5)*w*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/(L*(1 + Q11(i,j)))
;
elem[5][6]+=
0
;
elem[6][0]+=
(a0.diff(d01,i,j)*mu*exp(B1*mu*h(i,j))*Qr1(i,j)*D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*(1 + Q11(i,j)))
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
(complex(0,1)*w*Qr1(i,j)*D[-1 + d01](i,j,i1,j1))/L
;
elem[6][5]+=
((complex(0,4)*L*w*(2 - Q11.diff(d10,i,j) + Q11(i,j))*(1 + Qtt(i,j)) + 
       Q11.diff(d01,i,j)*(-12 + pow(mu,2))*Qr1(i,j)*(1 + Qtt(i,j)) + 
       (1 + Q11(i,j))*(36*Qtt.diff(d01,i,j)*Qr1(i,j) - 
          3*Qtt.diff(d01,i,j)*pow(mu,2)*Qr1(i,j) + 
          (complex(0,4)*L*w*(12 + pow(mu,2) + 
               Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
               2*pow(mu,2)*Qtt(i,j)))/(-12 + pow(mu,2)) + 
          (-12 + pow(mu,2))*(1 + Qtt(i,j))*
           (-4*Qr1.diff(d01,i,j) + 
             (complex(0,4)*L*w)/(-12 + pow(mu,2)) + 
             Qr1(i,j)*(-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) + 
                Qtt.diff(d01,i,j)/(1 + Qtt(i,j))))))*
     D[-1 + d01](i,j,i1,j1))/(8.*pow(L,2)*pow(1 + Q11(i,j),2)) - 
  ((-12 + pow(mu,2))*Qr1(i,j)*(1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (4.*pow(L,2)*(1 + Q11(i,j)))
;
elem[6][6]+=
0
;
elem[7][0]+=
0
;
elem[7][1]+=
((complex(0,4)*L*w*Qr1(i,j) + 
       ((1 + Qtt(i,j))*(Q11.diff(d01,i,j) + 
            (1 + Q11(i,j))*((-3*Q22.diff(d01,i,j))/(1 + Q22(i,j)) - 
               (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))))/
        pow(1 + Q11(i,j),2))*D[-1 + d01](i,j,i1,j1))/
   (2.*pow(L,2)*(-12 + pow(mu,2))) - 
  ((1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[7][2]+=
(2*h.diff(d01,i,j)*mu*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[7][3]+=
((L*(2 + (complex(0,8)*w)/(-12 + pow(mu,2)))*Qr1(i,j) + 
       (2*(1 + Qtt(i,j))*(Q11.diff(d01,i,j) + 
            (1 + Q11(i,j))*(-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) - 
               (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))))/
        ((-12 + pow(mu,2))*pow(1 + Q11(i,j),2)))*D[-1 + d01](i,j,i1,j1)\
)/(4.*pow(L,2)) - ((1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[7][4]+=
0
;
elem[7][5]+=
(complex(0,-2)*w*(1 + Qtt(i,j))*D[-1 + d01](i,j,i1,j1))/
  (L*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[7][6]+=
0
;
elem[8][0]+=
0
;
elem[8][1]+=
(((-2*Qr1.diff(d01,i,j) + L*
           ((complex(0,8)*w)/(-12 + pow(mu,2)) + 
             (2 - Q11.diff(d10,i,j) + Q11(i,j))/(1 + Q11(i,j)) + 
             (-2 + Q22.diff(d10,i,j) - Q22(i,j))/(1 + Q22(i,j))))*
        Qr1(i,j) + (Q11.diff(d01,i,j)*pow(Qr1(i,j),2))/
        (1 + Q11(i,j)) - (5*Q22.diff(d01,i,j)*pow(Qr1(i,j),2))/
        (1 + Q22(i,j)) + (complex(0,8)*L*w*
          (Qr1.diff(d10,i,j) + Qr1(i,j)))/(-12 + pow(mu,2)) + 
       (2*((-((12 + pow(mu,2) + 
                    2*Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                    3*(-4 + pow(mu,2))*Q11(i,j))*(1 + Qtt(i,j))) + 
               (-12 + pow(mu,2))*(1 + Q11(i,j))*
                (Qrr.diff(d10,i,j) + Qtt(i,j)))*
             (Q11.diff(d01,i,j) + 
               (1 + Q11(i,j))*
                ((-3*Q22.diff(d01,i,j))/(1 + Q22(i,j)) - 
                  (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))) + 
            (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))*
             (-Q11.diff(d01,i,j) + Q11.diff(d11,i,j) + 
               (-2 + Q11.diff(d10,i,j) - Q11(i,j))*
                ((-3*Q22.diff(d01,i,j))/(1 + Q22(i,j)) - 
                  (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j))) + 
               (1 + Q11(i,j))*
                ((-3*(-(Q22.diff(d01,i,j)*
                        (-1 + Q22.diff(d10,i,j))) + 
                       Q22.diff(d11,i,j)*(1 + Q22(i,j))))/
                   pow(1 + Q22(i,j),2) - 
                  (2*(-(Qtt.diff(d01,i,j)*
                        (-1 + Qtt.diff(d10,i,j))) + 
                       Qtt.diff(d11,i,j)*(1 + Qtt(i,j))))/
                   pow(1 + Qtt(i,j),2)))))/
        (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),3)))*
     D[-1 + d01](i,j,i1,j1))/(4.*pow(L,2)) - 
  ((pow(Qr1(i,j),2) + (-((Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
               3*(-4 + pow(mu,2)) + 4*(-6 + pow(mu,2))*Q11(i,j))*
             (1 + Qtt(i,j))) + 
          (-12 + pow(mu,2))*(1 + Q11(i,j))*
           (Qrr.diff(d10,i,j) + Qtt(i,j)))/
        (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)))*
     D[-1 + d02](i,j,i1,j1))/pow(L,2)
;
elem[8][2]+=
(((2*h.diff(d01,i,j)*mu*(3 - Q11.diff(d10,i,j) + 2*Q11(i,j))*
         (1 + Qtt(i,j)))/(-12 + pow(mu,2)) - 
      (1 + Q11(i,j))*(mu*(1 + Q11(i,j))*Qr1(i,j)*
          (h.diff(d10,i,j)*L + L*h(i,j) - 3*h.diff(d01,i,j)*Qr1(i,j)) \
+ (2*mu*(-(h.diff(d11,i,j)*(-12 + pow(mu,2))*(1 + Qtt(i,j))) + 
              h.diff(d01,i,j)*
               (-(Qrr.diff(d10,i,j)*(-12 + pow(mu,2))) + 
                 4*(-6 + pow(mu,2)) + 3*(-4 + pow(mu,2))*Qtt(i,j))))/
          pow(-12 + pow(mu,2),2)))*D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*pow(1 + Q11(i,j),2))
;
elem[8][3]+=
((L*(2 + (complex(0,8)*w)/(-12 + pow(mu,2)))*
        (Qr1.diff(d10,i,j) + Qr1(i,j)) + 
       pow(Qr1(i,j),2)*(Q11.diff(d01,i,j)/(1 + Q11(i,j)) - 
          Q22.diff(d01,i,j)/(1 + Q22(i,j)) - 
          (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j))) + 
       Qr1(i,j)*(-2*Qr1.diff(d01,i,j) + 
          L*((complex(0,8)*w)/(-12 + pow(mu,2)) + 
             (2 - Q11.diff(d10,i,j) + Q11(i,j))/(1 + Q11(i,j)) + 
             (-2 + Q22.diff(d10,i,j) - Q22(i,j))/(1 + Q22(i,j)) + 
             (2*(12 + pow(mu,2) + 
                  Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                  2*pow(mu,2)*Qtt(i,j)))/
              ((-12 + pow(mu,2))*(1 + Qtt(i,j))))) + 
       (2*((-((12 + pow(mu,2) + 
                    2*Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                    3*(-4 + pow(mu,2))*Q11(i,j))*(1 + Qtt(i,j))) + 
               (-12 + pow(mu,2))*(1 + Q11(i,j))*
                (Qrr.diff(d10,i,j) + Qtt(i,j)))*
             (Q11.diff(d01,i,j) + 
               (1 + Q11(i,j))*
                (-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) - 
                  (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j)))) + 
            (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))*
             (-Q11.diff(d01,i,j) + Q11.diff(d11,i,j) + 
               (-2 + Q11.diff(d10,i,j) - Q11(i,j))*
                (-(Q22.diff(d01,i,j)/(1 + Q22(i,j))) - 
                  (2*Qtt.diff(d01,i,j))/(1 + Qtt(i,j))) + 
               (1 + Q11(i,j))*
                ((Q22.diff(d01,i,j)*(-1 + Q22.diff(d10,i,j)) - 
                     Q22.diff(d11,i,j)*(1 + Q22(i,j)))/
                   pow(1 + Q22(i,j),2) - 
                  (2*(-(Qtt.diff(d01,i,j)*
                        (-1 + Qtt.diff(d10,i,j))) + 
                       Qtt.diff(d11,i,j)*(1 + Qtt(i,j))))/
                   pow(1 + Qtt(i,j),2)))))/
        (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),3)))*
     D[-1 + d01](i,j,i1,j1))/(4.*pow(L,2)) - 
  ((pow(Qr1(i,j),2) + (-((Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
               3*(-4 + pow(mu,2)) + 4*(-6 + pow(mu,2))*Q11(i,j))*
             (1 + Qtt(i,j))) + 
          (-12 + pow(mu,2))*(1 + Q11(i,j))*
           (Qrr.diff(d10,i,j) + Qtt(i,j)))/
        (pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)))*
     D[-1 + d02](i,j,i1,j1))/pow(L,2)
;
elem[8][4]+=
((4*L*Qr1(i,j) + (complex(0,8)*L*w*Qr1(i,j))/(-12 + pow(mu,2)) + 
       (2*(L*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*Qr1(i,j) - 
            (Q22.diff(d01,i,j)*(1 + Q11(i,j))*(1 + Qtt(i,j)))/
             (-12 + pow(mu,2)) + 
            (Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j)))/
             (-12 + pow(mu,2))))/(pow(1 + Q11(i,j),2)*(1 + Q22(i,j))))*
     D[-1 + d01](i,j,i1,j1))/(4.*pow(L,2)) - 
  ((1 + Qtt(i,j))*D[-1 + d02](i,j,i1,j1))/
   (pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[8][5]+=
(complex(0,-1)*w*(2*(-12 + pow(mu,2))*
       (2 - Q11.diff(d10,i,j) + Q11(i,j))*(1 + Qtt(i,j)) + 
      (1 + Q11(i,j))*(pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*
          pow(Qr1(i,j),2) + 2*
          (36 - 5*pow(mu,2) + Qrr.diff(d10,i,j)*(-12 + pow(mu,2)) - 
            4*(-6 + pow(mu,2))*Qtt(i,j))))*D[-1 + d01](i,j,i1,j1))/
  (L*pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2))
;
elem[8][6]+=
(2*mu*exp(B1*mu*h(i,j))*(L*(-12 + pow(mu,2))*a0(i,j)*(1 + Q11(i,j))*
       Qr1(i,j) - 2*a0.diff(d01,i,j)*(1 + Qtt(i,j)))*
    D[-1 + d01](i,j,i1,j1))/
  (pow(L,2)*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*(1 + Qtt(i,j)))
;
elem[9][0]+=
0
;
elem[9][1]+=
(complex(0,2)*w*D[-1 + d01](i,j,i1,j1))/pow(-12 + pow(mu,2),2)
;
elem[9][2]+=
0
;
elem[9][3]+=
((-12 + pow(mu,2) + complex(0,4)*w)*D[-1 + d01](i,j,i1,j1))/
  (2.*pow(-12 + pow(mu,2),2))
;
elem[9][4]+=
0
;
elem[9][5]+=
0
;
elem[9][6]+=
0
;
elem[10][0]+=
0
;
elem[10][1]+=
(((complex(0,-12)*(-4 + pow(mu,2))*w)/pow(-12 + pow(mu,2),3) + 
       ((complex(0,4)*w)/(-12 + pow(mu,2)) + 
          (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
             (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
          (-2 + Q22.diff(d10,i,j) - Q22(i,j) - 
             (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)))/
        (-12 + pow(mu,2)))*D[-1 + d01](i,j,i1,j1))/2. + 
  (Qr1(i,j)*D[-1 + d02](i,j,i1,j1))/(12*L - L*pow(mu,2))
;
elem[10][2]+=
(-2*mu*(h.diff(d10,i,j)*L + L*h(i,j) - 2*h.diff(d01,i,j)*Qr1(i,j))*
    D[-1 + d01](i,j,i1,j1))/(L*(-12 + pow(mu,2)))
;
elem[10][3]+=
((-3*(-4 + pow(mu,2))*(1 + (complex(0,4)*w)/(-12 + pow(mu,2))) + 
       (-12 + pow(mu,2))*((complex(0,4)*w)/(-12 + pow(mu,2)) + 
          (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
             (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) - 
          (Qtt.diff(d01,i,j)*Qr1(i,j))/(L + L*Qtt(i,j)) + 
          (12 + pow(mu,2) + Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
             2*pow(mu,2)*Qtt(i,j))/((-12 + pow(mu,2))*(1 + Qtt(i,j))))\
)*D[-1 + d01](i,j,i1,j1))/(2.*pow(-12 + pow(mu,2),2)) + 
  (Qr1(i,j)*D[-1 + d02](i,j,i1,j1))/(12*L - L*pow(mu,2))
;
elem[10][4]+=
((-12 + pow(mu,2) + complex(0,2)*w)*D[-1 + d01](i,j,i1,j1))/
  pow(-12 + pow(mu,2),2)
;
elem[10][5]+=
(complex(0,-1)*w*Qr1(i,j)*D[-1 + d01](i,j,i1,j1))/(-12 + pow(mu,2))
;
elem[10][6]+=
(4*mu*a0(i,j)*exp(B1*mu*h(i,j))*D[-1 + d01](i,j,i1,j1))/
  (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j)))
;
}
if(j==j1 && i==i1){
elem[0][0]+=
((Q22.diff(d01,i,j)*(1 + Q11(i,j)) - 
      Q11.diff(d01,i,j)*(1 + Q22(i,j)) + 
      2*h.diff(d01,i,j)*B1*mu*(1 + Q11(i,j))*(1 + Q22(i,j)))*
    (1 + Qtt(i,j)))/(2.*L*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[0][1]+=
0
;
elem[0][2]+=
-(B1*mu*a0(i,j))
;
elem[0][3]+=
-(mu*a0(i,j))/2.
;
elem[0][4]+=
0
;
elem[0][5]+=
0
;
elem[0][6]+=
-((-12 + pow(mu,2) + complex(0,2)*w)/(-12 + pow(mu,2)))
;
elem[1][0]+=
(2*L*(-(Q11.diff(d01,i,j)*(1 + Q22(i,j))) + 
       (1 + Q11(i,j))*(Q22.diff(d01,i,j) + 
          2*h.diff(d01,i,j)*B1*mu*(1 + Q22(i,j))))*
     (((-108 + 11*pow(mu,2) - 
             2*Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
             (-84 + 9*pow(mu,2))*Q11(i,j))*(1 + Q22(i,j)) - 
          (-12 + pow(mu,2))*(1 + Q11(i,j))*
           (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j)) - 
       (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
        (Qrr.diff(d10,i,j) + Qtt(i,j))) + 
    (pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*(1 + Q22(i,j))*
       ((complex(0,1)*(-(Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*
                 (1 + Q22(i,j))) + 
              (1 + Q11(i,j))*
               (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j) + 
                 2*Qr1.diff(d02,i,j)*(1 + Q22(i,j)) + 
                 2*h.diff(d01,i,j)*Qr1.diff(d01,i,j)*B1*mu*
                  (1 + Q22(i,j))))*pow(1 + Qtt(i,j),2))/w + 
         (2*L*(-(Q11.diff(d01,i,j)*(1 + Q22(i,j))*
                 pow(1 + Qtt(i,j),2)) - 
              Q11.diff(d11,i,j)*(1 + Q22(i,j))*pow(1 + Qtt(i,j),2) + 
              (-2 + Q11.diff(d10,i,j) - Q11(i,j))*
               (Q22.diff(d01,i,j) + 
                 2*h.diff(d01,i,j)*B1*mu*(1 + Q22(i,j)))*
               pow(1 + Qtt(i,j),2) + 
              (1 + Q11(i,j))*(Q22.diff(d01,i,j) + 
                 2*h.diff(d01,i,j)*B1*mu*(1 + Q22(i,j)))*
               pow(1 + Qtt(i,j),2) + 
              Q11.diff(d01,i,j)*(1 + Qtt(i,j))*
               ((7 - Q22.diff(d10,i,j) + 6*Q22(i,j))*(1 + Qtt(i,j)) - 
                 (1 + Q22(i,j))*(Qrr.diff(d10,i,j) + Qtt(i,j)) - 
                 (1 + Q22(i,j))*(Qtt.diff(d10,i,j) + Qtt(i,j))) + 
              (1 + Q11(i,j))*(Q22.diff(d01,i,j)*
                  (-5 + Qrr.diff(d10,i,j) + Qtt.diff(d10,i,j) - 
                    3*Qtt(i,j))*(1 + Qtt(i,j)) + 
                 Q22.diff(d11,i,j)*pow(1 + Qtt(i,j),2) - 
                 2*h.diff(d01,i,j)*B1*mu*
                  (2 - Q22.diff(d10,i,j) + Q22(i,j))*
                  pow(1 + Qtt(i,j),2) + 
                 (1 + Q22(i,j))*
                  (Qtt.diff(d01,i,j)*
                     (Qrr.diff(d10,i,j) - Qtt.diff(d10,i,j)) - 
                    2*h.diff(d01,i,j)*B1*mu*(1 + Qtt(i,j))*
                     (3 - Qrr.diff(d10,i,j) - Qtt.diff(d10,i,j) + 
                       Qtt(i,j)) + 
                    (1 + Qtt(i,j))*
                     (-Qrr.diff(d11,i,j) + Qtt.diff(d11,i,j) + 
                       2*h.diff(d11,i,j)*B1*mu*(1 + Qtt(i,j)))))))/
          (-12 + pow(mu,2))))/(1 + Qtt(i,j)))/
  (4.*pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)*
    pow(1 + Q22(i,j),2))
;
elem[1][1]+=
0
;
elem[1][2]+=
(B1*mu*(-2*a0.diff(d10,i,j)*L - L*a0(i,j) + a0.diff(d01,i,j)*Qr1(i,j)))/L
;
elem[1][3]+=
-(a0.diff(d10,i,j)*mu) + (a0.diff(d01,i,j)*mu*Qr1(i,j))/(2.*L)
;
elem[1][4]+=
mu*a0(i,j)
;
elem[1][5]+=
(complex(0,0.25)*mu*(-12 + pow(mu,2))*
    (2*a0.diff(d01,i,j)*(1 + Q11(i,j))*(1 + Q22(i,j)) + 
      a0(i,j)*(Q22.diff(d01,i,j)*(1 + Q11(i,j)) - 
         Q11.diff(d01,i,j)*(1 + Q22(i,j)) + 
         2*h.diff(d01,i,j)*B1*mu*(1 + Q11(i,j))*(1 + Q22(i,j))))*
    (1 + Qtt(i,j)))/(L*w*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[1][6]+=
(complex(0,-2)*w)/(-12 + pow(mu,2))
;
elem[2][0]+=
(complex(0,4)*pow(-12 + pow(mu,2),2)*w - 
    complex(0,4)*(-12 + pow(mu,2))*(-12 + pow(mu,2) - complex(0,2)*w)*
     w + (-(Qr1.diff(d01,i,j)*pow(-12 + pow(mu,2),3)) + 
       4*L*(Qrr.diff(d10,i,j)*(-12 + pow(mu,2)) - 
          Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) - 6*(-4 + pow(mu,2))\
)*pow(w,2) + (-(Qr1.diff(d01,i,j)*pow(-12 + pow(mu,2),3)) - 
          24*L*(-4 + pow(mu,2))*pow(w,2))*Qtt(i,j))/(L*(1 + Qtt(i,j))) \
- (complex(0,1)*pow(-12 + pow(mu,2),2)*w*
       (2*L*(2 - Q11.diff(d10,i,j) + Q11(i,j))*(1 + Q22(i,j))*
          (1 + Qtt(i,j)) - (1 + Q22(i,j))*
          (-((-2 + Q11.diff(d10,i,j))*L) + L*Q11(i,j) + 
            Q11.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)) - 
         2*L*(9 - Q11.diff(d10,i,j) - Q22.diff(d10,i,j) - Qrr.diff(d\
10,i,j) - Qtt.diff(d10,i,j) + 7*Qtt(i,j) - 
            Q11.diff(d10,i,j)*Qtt(i,j) - Q22.diff(d10,i,j)*Qtt(i,j) - 
            Q22(i,j)*(-8 + Q11.diff(d10,i,j) + Qrr.diff(d10,i,j\
) + Qtt.diff(d10,i,j) + (-6 + Q11.diff(d10,i,j))*Qtt(i,j)) - 
            Q11(i,j)*(-8 + Q22.diff(d10,i,j) + Qrr.diff(d10,i,j\
) + Qtt.diff(d10,i,j) + Q22(i,j)*
                (-7 + Qrr.diff(d10,i,j) + Qtt.diff(d10,i,j) - 
                  5*Qtt(i,j)) + (-6 + Q22.diff(d10,i,j))*Qtt(i,j))) + 
         (1 + Q11(i,j))*(2*L*(2 - Q22.diff(d10,i,j) + Q22(i,j))*
             (1 + Qtt(i,j)) + (-((-2 + Q22.diff(d10,i,j))*L) + 
               L*Q22(i,j) + Q22.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)) - 
            ((1 + Q22(i,j))*(24*Qr1.diff(d01,i,j) + 72*L - 
                 12*Qrr.diff(d10,i,j)*L - 36*Qtt.diff(d10,i,j)*L - 
                 24*h.diff(d10,i,j)*B1*L*mu - 
                 2*Qr1.diff(d01,i,j)*pow(mu,2) - 2*L*pow(mu,2) + 
                 Qrr.diff(d10,i,j)*L*pow(mu,2) + 
                 3*Qtt.diff(d10,i,j)*L*pow(mu,2) + 
                 2*h.diff(d10,i,j)*B1*L*pow(mu,3) + 
                 24*h.diff(d01,i,j)*B1*mu*Qr1(i,j) - 
                 2*h.diff(d01,i,j)*B1*pow(mu,3)*Qr1(i,j) + 
                 2*(-(Qr1.diff(d01,i,j)*(-12 + pow(mu,2))) + 
                    L*(12 + pow(mu,2) + 
                       h.diff(d10,i,j)*B1*mu*(-12 + pow(mu,2))) - 
                    h.diff(d01,i,j)*B1*mu*(-12 + pow(mu,2))*Qr1(i,j)\
)*Qtt(i,j) + 2*B1*L*mu*(-12 + pow(mu,2))*h(i,j)*(1 + Qtt(i,j))))/
             (-12 + pow(mu,2)))))/
     (L*(1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j))))/
  pow(-12 + pow(mu,2),3)
;
elem[2][1]+=
(complex(0,-4)*a0.diff(d01,i,j)*mu*w)/(L*pow(-12 + pow(mu,2),2))
;
elem[2][2]+=
(complex(0,4)*a0.diff(d01,i,j)*B1*mu*w)/(L*pow(-12 + pow(mu,2),2))
;
elem[2][3]+=
(complex(0,2)*a0.diff(d01,i,j)*mu*w)/(L*pow(-12 + pow(mu,2),2))
;
elem[2][4]+=
0
;
elem[2][5]+=
mu*(-1 - (complex(0,2)*w)/(-12 + pow(mu,2)))*a0(i,j)
;
elem[2][6]+=
0
;
elem[3][0]+=
(complex(0,-2)*a0.diff(d01,i,j)*B1*mu*w*exp(B1*mu*h(i,j)))/
  (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j)))
;
elem[3][1]+=
(mu*(h.diff(d01,i,j)*Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j)) - 
      (1 + Q11(i,j))*(h.diff(d01,i,j)*Q22.diff(d01,i,j)*
          (1 + Qtt(i,j)) + 2*(1 + Q22(i,j))*
          (h.diff(d01,i,j)*Qtt.diff(d01,i,j) + 
            h.diff(d02,i,j)*(1 + Qtt(i,j))))))/
  (2.*pow(L,2)*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*(1 + Qtt(i,j)))
;
elem[3][2]+=
((16*exp((-2*mu*h(i,j))/(3.*A1)) + 24*pow(A1,2)*exp(A1*mu*h(i,j)))/
     (2 + 3*pow(A1,2)) + (-12 + pow(mu,2))/(1 + Qtt(i,j)) + 
    (complex(0,4)*w)/(1 + Qtt(i,j)) + 
    (4*pow(w,2))/((-12 + pow(mu,2))*(1 + Qtt(i,j))) - 
    (complex(0,2)*(-12 + pow(mu,2) - complex(0,2)*w)*w*
       (-(Qrr.diff(d10,i,j)*(-12 + pow(mu,2))) + 
         8*(-9 + pow(mu,2)) + (-60 + 7*pow(mu,2))*Qtt(i,j)))/
     (pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2)) + 
    (pow(B1,2)*pow(mu,2)*pow(-12 + pow(mu,2),2)*pow(a0(i,j),2)*
        exp(B1*mu*h(i,j)) - 4*pow(w,2)*
        (12 + pow(mu,2) + Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
          2*pow(mu,2)*Qtt(i,j)))/
     (pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2)) + 
    (complex(0,1)*w*(2*L*(((-(Q11.diff(d10,i,j)*
                     (-12 + pow(mu,2))) + 2*(-72 + 7*pow(mu,2)) + 
                  (-132 + 13*pow(mu,2))*Q11(i,j))*(1 + Q22(i,j)) - 
               (-12 + pow(mu,2))*(1 + Q11(i,j))*
                (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j)) - 
            2*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
             (Qrr.diff(d10,i,j) + Qtt(i,j)) - 
            (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
             (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
         (-12 + pow(mu,2))*(-2*L*(2 - Q11.diff(d10,i,j) + Q11(i,j))*
             (1 + Q22(i,j))*(1 + Qtt(i,j)) - 
            (1 + Q22(i,j))*(-((-2 + Q11.diff(d10,i,j))*L) + 
               L*Q11(i,j) + Q11.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)) + 
            (1 + Q11(i,j))*(-2*L*(2 - Q22.diff(d10,i,j) + Q22(i,j))*
                (1 + Qtt(i,j)) - 
               (-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) + 
                  Q22.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)) + 
               ((1 + Q22(i,j))*
                  (-2*Qr1.diff(d01,i,j)*(-12 + pow(mu,2)) + 
                    L*(4*pow(mu,2) + 
                       Qrr.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                       3*Qtt.diff(d10,i,j)*(-12 + pow(mu,2))) + 
                    (-2*Qr1.diff(d01,i,j)*(-12 + pow(mu,2)) + 
                       8*L*(-6 + pow(mu,2)))*Qtt(i,j)))/
                (-12 + pow(mu,2))))))/
     (L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
       pow(1 + Qtt(i,j),2)))/2.
;
elem[3][3]+=
(mu*(B1*L*mu*pow(a0(i,j),2)*exp(B1*mu*h(i,j)) - 
      complex(0,1)*w*(h.diff(d10,i,j)*L + L*h(i,j) - 
         h.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j))))/
  (2.*L*pow(1 + Qtt(i,j),2))
;
elem[3][4]+=
0
;
elem[3][5]+=
(complex(0,-1)*h.diff(d01,i,j)*mu*w)/(L + L*Q11(i,j))
;
elem[3][6]+=
(B1*mu*(-12 + pow(mu,2) + complex(0,2)*w)*a0(i,j)*exp(B1*mu*h(i,j)))/
  ((-12 + pow(mu,2))*pow(1 + Qtt(i,j),2))
;
elem[4][0]+=
(complex(0,8)*a0.diff(d01,i,j)*mu*w*exp(B1*mu*h(i,j)))/
  (L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j)))
;
elem[4][1]+=
((complex(0,16)*w)/(-12 + pow(mu,2)) - 
    (complex(0,16)*(-12 + pow(mu,2) - complex(0,2)*w)*w)/
     pow(-12 + pow(mu,2),2) - 
    (complex(0,4)*w*(2*L*(1 + Q11(i,j))*
          (2 - Q22.diff(d10,i,j) + Q22(i,j))*(1 + Qtt(i,j)) + 
         (1 + Q11(i,j))*(-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) + 
            Q22.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)) - 
         2*L*(9 - Q11.diff(d10,i,j) - Q22.diff(d10,i,j) - Qrr.diff(\
d10,i,j) - Qtt.diff(d10,i,j) + 7*Qtt(i,j) - 
            Q11.diff(d10,i,j)*Qtt(i,j) - 
            Q22.diff(d10,i,j)*Qtt(i,j) - 
            Q22(i,j)*(-8 + Q11.diff(d10,i,j) + Qrr.diff(d10,i,j\
) + Qtt.diff(d10,i,j) + (-6 + Q11.diff(d10,i,j))*Qtt(i,j)) - 
            Q11(i,j)*(-8 + Q22.diff(d10,i,j) + Qrr.diff(d10,i,j\
) + Qtt.diff(d10,i,j) + Q22(i,j)*
                (-7 + Qrr.diff(d10,i,j) + Qtt.diff(d10,i,j) - 
                  5*Qtt(i,j)) + (-6 + Q22.diff(d10,i,j))*Qtt(i,j))) + 
         (1 + Q22(i,j))*(3*L*(2 - Q11.diff(d10,i,j) + Q11(i,j))*
             (1 + Qtt(i,j)) + 
            Q11.diff(d01,i,j)*Qr1(i,j)*(1 + Qtt(i,j)) - 
            ((1 + Q11(i,j))*(-2*Qr1.diff(d01,i,j)*
                  (-12 + pow(mu,2)) + 
                 L*(72 - 2*pow(mu,2) + 
                    Qrr.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                    3*Qtt.diff(d10,i,j)*(-12 + pow(mu,2))) + 
                 2*(-(Qr1.diff(d01,i,j)*(-12 + pow(mu,2))) + 
                    L*(12 + pow(mu,2)))*Qtt(i,j)))/(-12 + pow(mu,2))))\
)/(L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j))) - 
    (-4*Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*(1 + Q22(i,j))*
        (1 + Qtt(i,j)) + 4*pow(L,2)*(-12 + pow(mu,2))*(1 + Q11(i,j))*
        (2 - Q11.diff(d10,i,j) + Q11(i,j))*(1 + Q22(i,j))*
        (1 + Qtt(i,j)) + 2*(1 + Q11(i,j))*(1 + Q22(i,j))*
        (2*pow(Qtt.diff(d01,i,j),2) + 
          4*Qtt.diff(d02,i,j)*(1 + Qtt(i,j)) - 
          (pow(L,2)*(-12 + pow(mu,2))*
             (2 - Q11.diff(d10,i,j) + Q11(i,j))*(1 + Qtt(i,j)))/2. - 
          (Q11.diff(d01,i,j)*L*(-12 + pow(mu,2))*Qr1(i,j)*
             (1 + Qtt(i,j)))/2.) - 
       4*pow(L,2)*(1 + Q11(i,j))*
        (((3*(-20 + pow(mu,2)) - 
                Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                2*(-24 + pow(mu,2))*Q11(i,j))*(1 + Q22(i,j)) - 
             (-12 + pow(mu,2))*(1 + Q11(i,j))*
              (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j)) - 
          (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
           (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
       L*pow(1 + Q11(i,j),2)*
        (-8*L*pow(mu,2)*pow(a0(i,j),2)*exp(B1*mu*h(i,j))*
           (1 + Q22(i,j)) + (-12 + pow(mu,2))*
           (-2*Qr1.diff(d01,i,j) + 6*L - 3*Q22.diff(d10,i,j)*L - 
             Q22.diff(d01,i,j)*Qr1(i,j) - 
             4*Qtt.diff(d01,i,j)*Qr1(i,j) - 
             (2*Qr1.diff(d01,i,j) + 3*(-2 + Q22.diff(d10,i,j))*L + 
                Q22.diff(d01,i,j)*Qr1(i,j))*Qtt(i,j) + 
             Q22(i,j)*(-2*Qr1.diff(d01,i,j) + 3*L - 
                4*Qtt.diff(d01,i,j)*Qr1(i,j) + 
                (-2*Qr1.diff(d01,i,j) + 3*L)*Qtt(i,j)))))/
     (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j))*
       (1 + Qtt(i,j))) + (4*pow(L,2)*(1 + Q11(i,j))*
        (2 - Q22.diff(d10,i,j) + Q22(i,j))*(1 + Qtt(i,j)) + 
       2*L*(1 + Q11(i,j))*(-((-2 + Q22.diff(d10,i,j))*L) + 
          L*Q22(i,j) + Q22.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)) - 
       4*pow(L,2)*(8 - Q11.diff(d10,i,j) - Q22.diff(d10,i,j\
) - Qrr.diff(d10,i,j) - Qtt.diff(d10,i,j) + 6*Qtt(i,j) - 
          Q11.diff(d10,i,j)*Qtt(i,j) - Q22.diff(d10,i,j)*Qtt(i,j) - 
          Q22(i,j)*(-7 + Q11.diff(d10,i,j) + Qrr.diff(d10,i,j\
) + Qtt.diff(d10,i,j) + (-5 + Q11.diff(d10,i,j))*Qtt(i,j)) - 
          Q11(i,j)*(-7 + Q22.diff(d10,i,j) + Qrr.diff(d10,i,j\
) + Qtt.diff(d10,i,j) + Q22(i,j)*
              (-6 + Qrr.diff(d10,i,j) + Qtt.diff(d10,i,j) - 
                4*Qtt(i,j)) + (-5 + Q22.diff(d10,i,j))*Qtt(i,j))) + 
       (1 + Q22(i,j))*(2*Q11.diff(d01,i,j)*L*Qr1(i,j)*(1 + Qtt(i,j)) + 
          (2*(4*pow(Qtt.diff(d01,i,j),2) + 
               L*(-12 + pow(mu,2))*
                (2*Qr1.diff(d01,i,j) - 
                  (-2 + Q11.diff(d10,i,j))*L + 
                  (2*Qr1.diff(d01,i,j) + L)*Q11(i,j))*(1 + Qtt(i,j))))/
           (-12 + pow(mu,2)) - 
          (4*pow(L,2)*(48 - 12*Q11.diff(d10,i,j) - 
               24*Qtt.diff(d10,i,j) + 
               Q11.diff(d10,i,j)*pow(mu,2) + 
               2*Qtt.diff(d10,i,j)*pow(mu,2) + 
               (Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                  2*(12 + pow(mu,2)))*Qtt(i,j) + 
               Q11(i,j)*(36 + pow(mu,2) + 
                  2*Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                  3*(4 + pow(mu,2))*Qtt(i,j))))/(-12 + pow(mu,2))))/
     (2.*pow(L,2)*(1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j))) + 
    (4*((-4*pow(w,2)*(-36 + 5*pow(mu,2) - 
              Qrr.diff(d10,i,j)*(-12 + pow(mu,2)) + 
              4*(-6 + pow(mu,2))*Qtt(i,j)))/pow(-12 + pow(mu,2),2) + 
         pow(1 + Qtt(i,j),2)*
          ((72*pow(A1,2)*exp((-2*mu*h(i,j))/(3.*A1)) + 
               48*exp(A1*mu*h(i,j)))/(2 + 3*pow(A1,2)) - 
            (2*pow(h.diff(d01,i,j),2)*pow(mu,2))/
             (pow(L,2)*(1 + Q11(i,j))) + 
            pow(Qtt.diff(d01,i,j),2)/
             (pow(L,2)*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2)) + 
            (2*((-12 + pow(mu,2))*
                  ((Q11.diff(d01,i,j)*Qtt.diff(d01,i,j)*
                       (-12 + pow(mu,2)))/2. - 
                    Qtt.diff(d02,i,j)*(-12 + pow(mu,2))*
                     (1 + Q11(i,j)) - 
                    4*pow(L,2)*pow(w,2)*(1 + Q11(i,j))*
                     (2 - Q11.diff(d10,i,j) + Q11(i,j)))*(1 + Qtt(i,j)) \
+ 2*pow(L,2)*pow(w,2)*(1 + Q11(i,j))*
                  ((3*(-20 + pow(mu,2)) - 
                       2*Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                       (-36 + pow(mu,2))*Q11(i,j))*(1 + Qtt(i,j)) - 
                    (-12 + pow(mu,2))*(1 + Q11(i,j))*
                     (Qtt.diff(d10,i,j) + Qtt(i,j)))))/
             (pow(L,2)*pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
               pow(1 + Qtt(i,j),2)))))/
     ((-12 + pow(mu,2))*(1 + Qtt(i,j))))/4.
;
elem[4][2]+=
0
;
elem[4][3]+=
(complex(0,0.5)*w*(2*Qr1.diff(d01,i,j) - Q11.diff(d10,i,j)*L + 
      Q22.diff(d10,i,j)*L + Q11.diff(d01,i,j)*Qr1(i,j) - 
      Q22.diff(d01,i,j)*Qr1(i,j) + 
      Q22(i,j)*(2*Qr1.diff(d01,i,j) + L - Q11.diff(d10,i,j)*L + 
         Q11.diff(d01,i,j)*Qr1(i,j)) + 
      Q11(i,j)*((-1 + Q22.diff(d10,i,j))*L + 
         2*Qr1.diff(d01,i,j)*(1 + Q22(i,j)) - Q22.diff(d01,i,j)*Qr1(i,j)\
)))/(L*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j)))
;
elem[4][4]+=
0
;
elem[4][5]+=
(complex(0,1)*w*(Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j)) + 
      (1 + Q11(i,j))*(-2*Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
         Q22.diff(d01,i,j)*(1 + Qtt(i,j)))))/
  (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[4][6]+=
0
;
elem[5][0]+=
(complex(0,-2)*a0.diff(d01,i,j)*mu*w*exp(B1*mu*h(i,j)))/
  (L*(-12 + pow(mu,2))*(1 + Q11(i,j)))
;
elem[5][1]+=
(complex(0,-0.25)*w*((1 + Q22(i,j))*
       (2*Qr1.diff(d01,i,j) + 2*L - Q11.diff(d10,i,j)*L + 
         (2*Qr1.diff(d01,i,j) + L)*Q11(i,j) + 
         Q11.diff(d01,i,j)*Qr1(i,j)) - 
      (1 + Q11(i,j))*(-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) + 
         Q22.diff(d01,i,j)*Qr1(i,j))))/(L*(1 + Q11(i,j))*(1 + Q22(i,j)))
;
elem[5][2]+=
(complex(0,1)*mu*w*(h.diff(d10,i,j)*L + L*h(i,j) - 
      h.diff(d01,i,j)*Qr1(i,j)))/L
;
elem[5][3]+=
0
;
elem[5][4]+=
(w*(complex(0,-2) + (8*w)/(-12 + pow(mu,2))))/4.
;
elem[5][5]+=
(complex(0,0.25)*w*(-(Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j))) + 
      (1 + Q11(i,j))*(2*Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
         Q22.diff(d01,i,j)*(1 + Qtt(i,j)))))/
  (L*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))
;
elem[5][6]+=
0
;
elem[6][0]+=
(complex(0,-1)*mu*exp(B1*mu*h(i,j))*
    (2*a0.diff(d11,i,j)*L*w*(1 + Q11(i,j)) + 
      a0.diff(d01,i,j)*(complex(0,-12)*Qr1.diff(d01,i,j) + 
         complex(0,1)*Qr1.diff(d01,i,j)*pow(mu,2) + 6*L*w - 
         2*Q11.diff(d10,i,j)*L*w + 2*h.diff(d10,i,j)*B1*L*mu*w + 
         (complex(0,1)*Qr1.diff(d01,i,j)*(-12 + pow(mu,2)) + 
            2*L*(2 + h.diff(d10,i,j)*B1*mu)*w)*Q11(i,j) + 
         2*B1*L*mu*w*h(i,j)*(1 + Q11(i,j)))))/
  (pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2))
;
elem[6][1]+=
(complex(0,-0.25)*w*(-((-4 + Q11.diff(d10,i,j) + Q22.diff(d10,i,j) + 
           Q11(i,j)*(-3 + Q22.diff(d10,i,j) - 2*Q22(i,j)) + 
           (-3 + Q11.diff(d10,i,j))*Q22(i,j))*
         ((1 + Q22(i,j))*(2*Qr1.diff(d01,i,j) + 2*L - 
              Q11.diff(d10,i,j)*L + 
              (2*Qr1.diff(d01,i,j) + L)*Q11(i,j) + 
              Q11.diff(d01,i,j)*Qr1(i,j)) - 
           (1 + Q11(i,j))*(-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) + 
              Q22.diff(d01,i,j)*Qr1(i,j)))) + 
      (1 + Q11(i,j))*(1 + Q22(i,j))*
       (-((2 - Q22.diff(d10,i,j) + Q22(i,j))*
            (2*Qr1.diff(d01,i,j) + 2*L - Q11.diff(d10,i,j)*L + 
              (2*Qr1.diff(d01,i,j) + L)*Q11(i,j) + 
              Q11.diff(d01,i,j)*Qr1(i,j))) + 
         (1 + Q22(i,j))*(2*(-1 + Q11.diff(d10,i,j))*Qr1.diff(d01,i,j\
) + Q11.diff(d01,i,j)*Qr1.diff(d10,i,j) + 2*Qr1.diff(d11,i,j) - 
            6*L + 2*Q11.diff(d10,i,j)*L - Q11.diff(d20,i,j)*L - 
            2*(-Qr1.diff(d11,i,j) + L)*Q11(i,j) + 
            Q11.diff(d11,i,j)*Qr1(i,j)) + 
         (2 - Q11.diff(d10,i,j) + Q11(i,j))*
          (-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) + 
            Q22.diff(d01,i,j)*Qr1(i,j)) + 
         (1 + Q11(i,j))*(-(Q22.diff(d01,i,j)*Qr1.diff(d10,i,j)) + 6*L - 
            2*Q22.diff(d10,i,j)*L + Q22.diff(d20,i,j)*L + 
            2*L*Q22(i,j) - Q22.diff(d11,i,j)*Qr1(i,j)))))/
  (L*pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),2))
;
elem[6][2]+=
(complex(0,1)*mu*w*(-(h.diff(d01,i,j)*Qr1.diff(d10,i,j)) + 
      3*h.diff(d10,i,j)*L + h.diff(d20,i,j)*L + L*h(i,j) - 
      (3*h.diff(d01,i,j) + h.diff(d11,i,j))*Qr1(i,j)))/L
;
elem[6][3]+=
0
;
elem[6][4]+=
(w*(8*w + (complex(0,1)*(-2*L*(-12 + pow(mu,2))*(1 + Q11(i,j))*
            (2 - Q22.diff(d10,i,j) + Q22(i,j))*(1 + Qtt(i,j)) + 
           (-12 + pow(mu,2))*(1 + Q11(i,j))*
            (-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) + 
              Q22.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j)) - 
           2*L*((-((3*(-20 + pow(mu,2)) - 
                      Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                      2*(-24 + pow(mu,2))*Q11(i,j))*(1 + Q22(i,j))) + 
                 (-12 + pow(mu,2))*(1 + Q11(i,j))*
                  (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j)) + 
              (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
               (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
           (1 + Q22(i,j))*(L*(-12 + pow(mu,2))*
               (2 - Q11.diff(d10,i,j) + Q11(i,j))*(1 + Qtt(i,j)) + 
              (-12 + pow(mu,2))*Qr1(i,j)*
               (-2*Qtt.diff(d01,i,j)*(1 + Q11(i,j)) + 
                 Q11.diff(d01,i,j)*(1 + Qtt(i,j))) + 
              4*(-(L*(-12 + pow(mu,2))*
                     (2 - Q11.diff(d10,i,j) + Q11(i,j))*(1 + Qtt(i,j)))/
                  2. + (1 + Q11(i,j))*
                  ((Qr1.diff(d01,i,j)*(-12 + pow(mu,2))*
                       (1 + Qtt(i,j)))/2. + 
                    L*(12 + pow(mu,2) + 
                       Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                       2*pow(mu,2)*Qtt(i,j)))))))/
       (L*(1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j)))))/
  (4.*(-12 + pow(mu,2)))
;
elem[6][5]+=
(complex(0,-2)*L*w*(((-108 + 11*pow(mu,2) - 
             2*Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
             (-84 + 9*pow(mu,2))*Q11(i,j))*(1 + Q22(i,j)) - 
          (-12 + pow(mu,2))*(1 + Q11(i,j))*
           (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j)) - 
       (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
        (Qrr.diff(d10,i,j) + Qtt(i,j)))*
     (-(Q22.diff(d01,i,j)*(1 + Q11(i,j))*(1 + Qtt(i,j))) + 
       (1 + Q22(i,j))*(-2*Qtt.diff(d01,i,j)*(1 + Q11(i,j)) + 
          Q11.diff(d01,i,j)*(1 + Qtt(i,j)))) + 
    pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*(1 + Q22(i,j))*
     (Q11.diff(d01,i,j)*Qr1.diff(d01,i,j)*(1 + Q22(i,j))*
        pow(1 + Qtt(i,j),2) - 
       (1 + Q11(i,j))*(1 + Qtt(i,j))*
        (Q22.diff(d01,i,j)*Qr1.diff(d01,i,j)*(1 + Qtt(i,j)) + 
          (2*(1 + Q22(i,j))*(-4*a0.diff(d01,i,j)*L*pow(mu,2)*a0(i,j)*
                exp(B1*mu*h(i,j)) + 
               (-12 + pow(mu,2))*
                (Qr1.diff(d01,i,j)*Qtt.diff(d01,i,j) + 
                  Qr1.diff(d02,i,j)*(1 + Qtt(i,j)))))/(-12 + pow(mu,2))\
) - (complex(0,2)*L*w*(-(Q22.diff(d01,i,j)*(1 + Q11(i,j))*
               pow(1 + Qtt(i,j),2)) - 
            Q22.diff(d11,i,j)*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2) + 
            (-2 + Q22.diff(d10,i,j) - Q22(i,j))*(1 + Qtt(i,j))*
             (-2*Qtt.diff(d01,i,j)*(1 + Q11(i,j)) + 
               Q11.diff(d01,i,j)*(1 + Qtt(i,j))) + 
            (1 + Q22(i,j))*(1 + Qtt(i,j))*
             (-2*Qtt.diff(d01,i,j)*(1 + Q11(i,j)) + 
               Q11.diff(d01,i,j)*(1 + Qtt(i,j))) + 
            Q22.diff(d01,i,j)*(1 + Qtt(i,j))*
             ((7 - Q11.diff(d10,i,j) + 6*Q11(i,j))*(1 + Qtt(i,j)) - 
               (1 + Q11(i,j))*(Qrr.diff(d10,i,j) + Qtt(i,j)) - 
               (1 + Q11(i,j))*(Qtt.diff(d10,i,j) + Qtt(i,j))) + 
            (1 + Q22(i,j))*(-(Q11.diff(d01,i,j)*(1 + Qtt(i,j))*
                  (5 - Qrr.diff(d10,i,j) - Qtt.diff(d10,i,j) + 
                    3*Qtt(i,j))) + 
               (1 + Qtt(i,j))*
                ((Qrr.diff(d11,i,j) - 3*Qtt.diff(d11,i,j))*
                   (1 + Q11(i,j)) + Q11.diff(d11,i,j)*(1 + Qtt(i,j))) + 
               Qtt.diff(d01,i,j)*
                (10 - 2*Q11.diff(d10,i,j) - 3*Qrr.diff(d10,i,j) + Qtt\
.diff(d10,i,j) - 2*(-4 + Q11.diff(d10,i,j))*Qtt(i,j) + 
                  Q11(i,j)*(8 - 3*Qrr.diff(d10,i,j) + Qtt.diff(d10,i,j\
) + 6*Qtt(i,j))))))/(-12 + pow(mu,2))))/
  (8.*pow(L,2)*(-12 + pow(mu,2))*pow(1 + Q11(i,j),3)*
    pow(1 + Q22(i,j),2)*(1 + Qtt(i,j)))
;
elem[6][6]+=
0
;
elem[7][0]+=
(complex(0,-4)*mu*w*exp(B1*mu*h(i,j))*
    (L*(-12 + pow(mu,2))*a0(i,j)*(1 + Q11(i,j))*Qr1(i,j) - 
      a0.diff(d01,i,j)*(1 + Qtt(i,j))))/
  (L*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*(1 + Qtt(i,j)))
;
elem[7][1]+=
(complex(0,1)*w*((2*Qr1.diff(d01,i,j))/L + 
       (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
          (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
       (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
          (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))) + 
    (Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j))*
        (Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
          Q22.diff(d01,i,j)*(1 + Qtt(i,j))) + 
       (1 + Q11(i,j))*(pow(Q22.diff(d01,i,j),2)*
           pow(1 + Qtt(i,j),2) - 
          (1 + Q22(i,j))*(1 + Qtt(i,j))*
           (Q22.diff(d01,i,j)*Qtt.diff(d01,i,j) + 
             2*Q22.diff(d02,i,j)*(1 + Qtt(i,j))) - 
          pow(1 + Q22(i,j),2)*
           (-pow(Qtt.diff(d01,i,j),2) + 
             2*Qtt.diff(d02,i,j)*(1 + Qtt(i,j)) + 
             2*pow(h.diff(d01,i,j),2)*pow(mu,2)*pow(1 + Qtt(i,j),2)\
)))/(pow(L,2)*pow(1 + Q11(i,j),2)*pow(1 + Q22(i,j),2)*(1 + Qtt(i,j))))/
  (2.*(-12 + pow(mu,2)))
;
elem[7][2]+=
(complex(0,-2)*mu*w*(h.diff(d10,i,j)*L + L*h(i,j) + 
       h.diff(d01,i,j)*Qr1(i,j)) + 
    (L*exp((-2*mu*h(i,j))/(3.*A1))*
       ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
          exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) - 
         24*A1*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*pow(1 + Qtt(i,j),2)\
))/((2 + 3*pow(A1,2))*(1 + Qtt(i,j))))/(L*(-12 + pow(mu,2)))
;
elem[7][3]+=
(complex(0,1)*w*((2*Qr1.diff(d01,i,j))/L + 
       (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
          (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
       (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
          (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j))) + 
    (2*pow(mu,2)*pow(a0(i,j),2)*exp(B1*mu*h(i,j)))/(1 + Qtt(i,j)))/
  (2.*(-12 + pow(mu,2)))
;
elem[7][4]+=
-0.5 - (complex(0,1)*w)/(-12 + pow(mu,2)) - 
  (4*pow(w,2))/pow(-12 + pow(mu,2),2)
;
elem[7][5]+=
(w*((complex(0,2) - (4*w)/(-12 + pow(mu,2)))*Qr1(i,j) - 
      (complex(0,2)*(-(Q11.diff(d01,i,j)*(1 + Q22(i,j))*
              (1 + Qtt(i,j))) + 
           (1 + Q11(i,j))*(2*Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
              Q22.diff(d01,i,j)*(1 + Qtt(i,j)))))/
       (L*(-12 + pow(mu,2))*pow(1 + Q11(i,j),2)*(1 + Q22(i,j)))))/2.
;
elem[7][6]+=
(2*mu*(-12 + pow(mu,2) + complex(0,2)*w)*a0(i,j)*exp(B1*mu*h(i,j)))/
  (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j)))
;
elem[8][0]+=
(complex(0,-2)*w*exp(B1*mu*h(i,j))*
    ((2*L*mu*a0(i,j)*(1 + Q11(i,j))*Qr1(i,j) - 
         (2*a0.diff(d01,i,j)*mu*(1 + Qtt(i,j)))/(-12 + pow(mu,2)))*
       ((-36 + 12*Q11.diff(d10,i,j) - 12*h.diff(d10,i,j)*B1*mu + 
            pow(mu,2) - Q11.diff(d10,i,j)*pow(mu,2) + 
            h.diff(d10,i,j)*B1*pow(mu,3) + 
            (-24 + h.diff(d10,i,j)*B1*mu*(-12 + pow(mu,2)))*
             Q11(i,j) + B1*mu*(-12 + pow(mu,2))*h(i,j)*(1 + Q11(i,j)))*
          (1 + Qtt(i,j)) - (-12 + pow(mu,2))*(1 + Q11(i,j))*
          (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
      2*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Qtt(i,j))*
       (L*mu*a0(i,j)*(Qr1.diff(d10,i,j)*(1 + Q11(i,j)) + 
            (-1 + Q11.diff(d10,i,j))*Qr1(i,j)) + 
         (1 + Q11(i,j))*Qr1(i,j)*
          (2*a0.diff(d10,i,j)*L*mu - a0.diff(d01,i,j)*mu*Qr1(i,j)) + 
         (mu*(-(a0.diff(d11,i,j)*(-12 + pow(mu,2))*(1 + Qtt(i,j))) + 
              a0.diff(d01,i,j)*
               (-36 + 5*pow(mu,2) - 
                 Qrr.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                 4*(-6 + pow(mu,2))*Qtt(i,j))))/pow(-12 + pow(mu,2),2)\
)))/(L*pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),2)*
    pow(1 + Qtt(i,j),2))
;
elem[8][1]+=
(complex(0,1)*(-12 + pow(mu,2))*w*
     ((2*Qr1.diff(d01,i,j))/L + 
       ((2 - Q11.diff(d10,i,j) + Q11(i,j))*
          (-((-2 + Q11.diff(d10,i,j))*L) + L*Q11(i,j) + 
            Q11.diff(d01,i,j)*Qr1(i,j)))/(L*pow(1 + Q11(i,j),2)) - 
       ((2 - Q22.diff(d10,i,j) + Q22(i,j))*
          (-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) - 
            3*Q22.diff(d01,i,j)*Qr1(i,j)))/(L*pow(1 + Q22(i,j),2)) + 
       (3*Q22.diff(d01,i,j)*Qr1.diff(d10,i,j) + 6*L - 
          2*Q22.diff(d10,i,j)*L + Q22.diff(d20,i,j)*L + 
          2*L*Q22(i,j) + 3*Q22.diff(d11,i,j)*Qr1(i,j))/(L + L*Q22(i,j)) \
+ (2 - Q11.diff(d10,i,j) + Q11(i,j) + (Q11.diff(d01,i,j)*Qr1(i,j))/L)/
        (1 + Q11(i,j)) + (-6 + 2*Q11.diff(d10,i,j) - Q11.diff(d20,i,j\
) + (Q11.diff(d01,i,j)*Qr1.diff(d10,i,j))/L - 2*Q11(i,j) + 
          (Q11.diff(d11,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
       (-2 + Q22.diff(d10,i,j) - Q22(i,j) + 
          (3*Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)) + 
       (2*(Qr1.diff(d01,i,j)*pow(1 + Qtt(i,j),2) + 
            Qr1.diff(d11,i,j)*pow(1 + Qtt(i,j),2) + 
            Qr1(i,j)*(Qtt.diff(d01,i,j)*
                (Qrr.diff(d10,i,j) - Qtt.diff(d10,i,j)) - 
               (Qrr.diff(d11,i,j) - Qtt.diff(d11,i,j))*(1 + Qtt(i,j)))\
))/(L*pow(1 + Qtt(i,j),2))) + 
    ((-(((84 - pow(mu,2) + 
                  2*Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                  (60 + pow(mu,2))*Q11(i,j))*(1 + Q22(i,j)) + 
               2*(-12 + pow(mu,2))*(1 + Q11(i,j))*
                (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j))) + 
          (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
           (Qrr.diff(d10,i,j) + Qtt(i,j)) - 
          2*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
           (Qtt.diff(d10,i,j) + Qtt(i,j)))*
        (Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j))*
           (Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
             Q22.diff(d01,i,j)*(1 + Qtt(i,j))) + 
          (1 + Q11(i,j))*(pow(Q22.diff(d01,i,j),2)*
              pow(1 + Qtt(i,j),2) - 
             (1 + Q22(i,j))*(1 + Qtt(i,j))*
              (Q22.diff(d01,i,j)*Qtt.diff(d01,i,j) + 
                2*Q22.diff(d02,i,j)*(1 + Qtt(i,j))) - 
             pow(1 + Q22(i,j),2)*
              (-pow(Qtt.diff(d01,i,j),2) + 
                2*Qtt.diff(d02,i,j)*(1 + Qtt(i,j)) + 
                2*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                 pow(1 + Qtt(i,j),2)))) + 
       (1 + Q11(i,j))*(1 + Q22(i,j))*(1 + Qtt(i,j))*
        (Q11.diff(d01,i,j)*(1 + Q22(i,j))*(1 + Qtt(i,j))*
           (Qtt.diff(d01,i,j)*
              (24 + Q22.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                (12 + pow(mu,2))*Q22(i,j)) + 
             Q22.diff(d01,i,j)*
              (24 + Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                (12 + pow(mu,2))*Qtt(i,j)) + 
             (-12 + pow(mu,2))*
              (Qtt.diff(d11,i,j)*(1 + Q22(i,j)) + 
                Q22.diff(d11,i,j)*(1 + Qtt(i,j)))) + 
          (-12 + pow(mu,2))*(-2 + Q11.diff(d10,i,j) - Q11(i,j))*
           (pow(Q22.diff(d01,i,j),2)*pow(1 + Qtt(i,j),2) - 
             (1 + Q22(i,j))*(1 + Qtt(i,j))*
              (Q22.diff(d01,i,j)*Qtt.diff(d01,i,j) + 
                2*Q22.diff(d02,i,j)*(1 + Qtt(i,j))) - 
             pow(1 + Q22(i,j),2)*
              (-pow(Qtt.diff(d01,i,j),2) + 
                2*Qtt.diff(d02,i,j)*(1 + Qtt(i,j)) + 
                2*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                 pow(1 + Qtt(i,j),2))) + 
          (Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
             Q22.diff(d01,i,j)*(1 + Qtt(i,j)))*
           (Q11.diff(d11,i,j)*(-12 + pow(mu,2))*(1 + Q22(i,j))*
              (1 + Qtt(i,j)) + 
             Q11.diff(d01,i,j)*
              ((48 - 2*pow(mu,2) + 
                   Q22.diff(d10,i,j)*(-12 + pow(mu,2)) - 
                   (-36 + pow(mu,2))*Q22(i,j))*(1 + Qtt(i,j)) + 
                (-12 + pow(mu,2))*(1 + Q22(i,j))*
                 (Qtt.diff(d10,i,j) + Qtt(i,j)))) + 
          (1 + Q11(i,j))*(2*Q22.diff(d01,i,j)*Q22.diff(d11,i,j)*
              (-12 + pow(mu,2))*pow(1 + Qtt(i,j),2) + 
             2*pow(Q22.diff(d01,i,j),2)*(1 + Qtt(i,j))*
              (24 + Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                (12 + pow(mu,2))*Qtt(i,j)) + 
             2*(-12 + pow(mu,2))*(1 + Q22(i,j))*
              (2 - Q22.diff(d10,i,j) + Q22(i,j))*
              (-pow(Qtt.diff(d01,i,j),2) + 
                2*Qtt.diff(d02,i,j)*(1 + Qtt(i,j)) + 
                2*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                 pow(1 + Qtt(i,j),2)) + 
             (Q22.diff(d01,i,j)*Qtt.diff(d01,i,j) + 
                2*Q22.diff(d02,i,j)*(1 + Qtt(i,j)))*
              (-((36 - pow(mu,2) + 
                     Q22.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                     24*Q22(i,j))*(1 + Qtt(i,j))) - 
                (-12 + pow(mu,2))*(1 + Q22(i,j))*
                 (Qtt.diff(d10,i,j) + Qtt(i,j))) - 
             (1 + Q22(i,j))*(1 + Qtt(i,j))*
              (Q22.diff(d01,i,j)*
                 (Qtt.diff(d11,i,j)*(-12 + pow(mu,2)) + 
                   Qtt.diff(d01,i,j)*(12 + pow(mu,2))) + 
                2*Q22.diff(d02,i,j)*
                 (24 + Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                   (12 + pow(mu,2))*Qtt(i,j)) + 
                (-12 + pow(mu,2))*
                 (Q22.diff(d11,i,j)*Qtt.diff(d01,i,j) + 
                   2*Q22.diff(d12,i,j)*(1 + Qtt(i,j)))) + 
             2*pow(1 + Q22(i,j),2)*
              (2*pow(Qtt.diff(d01,i,j),2)*pow(mu,2) + 
                Qtt.diff(d01,i,j)*Qtt.diff(d11,i,j)*
                 (-12 + pow(mu,2)) + 
                (-2*Qtt.diff(d02,i,j)*pow(mu,2) - 
                   Qtt.diff(d12,i,j)*(-12 + pow(mu,2)) + 
                   2*pow(a0.diff(d01,i,j),2)*pow(mu,2)*
                    exp(B1*mu*h(i,j)))*(1 + Qtt(i,j)) - 
                Qtt.diff(d02,i,j)*
                 (12 + pow(mu,2) + 
                   Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                   2*pow(mu,2)*Qtt(i,j)) - 
                2*h.diff(d01,i,j)*pow(mu,2)*(1 + Qtt(i,j))*
                 (h.diff(d11,i,j)*(-12 + pow(mu,2))*(1 + Qtt(i,j)) + 
                   h.diff(d01,i,j)*
                    (2*pow(mu,2) + 
                      Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                      3*(-4 + pow(mu,2))*Qtt(i,j)))))))/
     (pow(L,2)*pow(1 + Q11(i,j),3)*pow(1 + Q22(i,j),3)*
       pow(1 + Qtt(i,j),2)))/(2.*pow(-12 + pow(mu,2),2))
;
elem[8][2]+=
-(2*h.diff(d10,i,j)*pow(L,2)*mu + 2*pow(L,2)*mu*h(i,j) + 
     2*h.diff(d01,i,j)*L*mu*Qr1(i,j) + 
     (complex(0,4)*L*mu*w*(h.diff(d01,i,j)*Qr1.diff(d10,i,j) + 
          4*h.diff(d10,i,j)*L + h.diff(d20,i,j)*L + 2*L*h(i,j) + 
          (4*h.diff(d01,i,j) + h.diff(d11,i,j))*Qr1(i,j)))/
      (-12 + pow(mu,2)) - (2*exp((-2*mu*h(i,j))/(3.*A1))*
        (pow(L,2)*((2 + 3*pow(A1,2))*B1*pow(mu,2)*
              pow(a0(i,j),2)*exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) - 
             24*A1*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
              pow(1 + Qtt(i,j),2))*
           (((-144*A1 + 36*Q11.diff(d10,i,j)*A1 + 
                  24*h.diff(d10,i,j)*mu + 6*A1*pow(mu,2) - 
                  3*Q11.diff(d10,i,j)*A1*pow(mu,2) - 
                  2*h.diff(d10,i,j)*pow(mu,3) + 
                  (3*A1*(-36 + pow(mu,2)) - 
                     2*h.diff(d10,i,j)*mu*(-12 + pow(mu,2)))*
                   Q11(i,j) - 
                  2*mu*(-12 + pow(mu,2))*h(i,j)*(1 + Q11(i,j)))*
                (1 + Qtt(i,j)))/A1 - 
             3*(-12 + pow(mu,2))*(1 + Q11(i,j))*
              (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
          3*(-12 + pow(mu,2))*(1 + Qtt(i,j))*
           ((-2*pow(a0.diff(d01,i,j),2)*(2 + 3*pow(A1,2))*B1*
                pow(mu,2)*exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1))*
                (1 + Qtt(i,j)))/(-12 + pow(mu,2)) + 
             pow(L,2)*(-2 + Q11.diff(d10,i,j) - Q11(i,j))*
              ((2 + 3*pow(A1,2))*B1*pow(mu,2)*pow(a0(i,j),2)*
                 exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1)) - 
                24*A1*(-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
                 pow(1 + Qtt(i,j),2)) + 
             (1 + Q11(i,j))*(((2 + 3*pow(A1,2))*B1*L*pow(mu,2)*
                   a0(i,j)*exp(((2 + 3*A1*B1)*mu*h(i,j))/(3.*A1))*
                   ((2 + 3*A1*B1)*L*mu*a0(i,j)*
                      (h.diff(d10,i,j) + h(i,j)) + 
                     6*A1*(2*a0.diff(d10,i,j)*L + 
                        a0.diff(d01,i,j)*Qr1(i,j))))/(3.*A1) - 
                24*A1*pow(L,2)*(1 + Qtt(i,j))*
                 ((4 - 4*exp((2/(3.*A1) + A1)*mu*h(i,j)) + 
                      (2/(3.*A1) + A1)*mu*
                       exp((2/(3.*A1) + A1)*mu*h(i,j))*
                       (h.diff(d10,i,j) + h(i,j)))*(1 + Qtt(i,j)) + 
                   (-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
                    (Qrr.diff(d10,i,j) + Qtt(i,j)) + 
                   (-1 + exp((2/(3.*A1) + A1)*mu*h(i,j)))*
                    (Qtt.diff(d10,i,j) + Qtt(i,j)))))))/
      (3.*(2 + 3*pow(A1,2))*pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*
        pow(1 + Qtt(i,j),2)))/(2.*pow(L,2))
;
elem[8][3]+=
((2*pow(mu,2)*exp(B1*mu*h(i,j))*
       ((-12 + pow(mu,2))*(1 + Qtt(i,j))*
          (pow(L,2)*pow(a0(i,j),2)*
             (-2 + Q11.diff(d10,i,j) - Q11(i,j)) + 
            2*L*a0(i,j)*(1 + Q11(i,j))*
             (2*a0.diff(d10,i,j)*L + a0.diff(d01,i,j)*Qr1(i,j)) - 
            (2*pow(a0.diff(d01,i,j),2)*(1 + Qtt(i,j)))/
             (-12 + pow(mu,2))) + 
         pow(L,2)*pow(a0(i,j),2)*
          ((-36 + 12*Q11.diff(d10,i,j) - 12*h.diff(d10,i,j)*B1*mu + 
               pow(mu,2) - Q11.diff(d10,i,j)*pow(mu,2) + 
               h.diff(d10,i,j)*B1*pow(mu,3) + 
               (-24 + h.diff(d10,i,j)*B1*mu*(-12 + pow(mu,2)))*
                Q11(i,j) + B1*mu*(-12 + pow(mu,2))*h(i,j)*
                (1 + Q11(i,j)))*(1 + Qtt(i,j)) - 
            (-12 + pow(mu,2))*(1 + Q11(i,j))*
             (Qtt.diff(d10,i,j) + Qtt(i,j)))))/
     (pow(L,2)*(1 + Q11(i,j))*pow(1 + Qtt(i,j),2)) + 
    complex(0,1)*(-12 + pow(mu,2))*w*
     ((2*Qr1.diff(d01,i,j))/L + 
       ((2 - Q11.diff(d10,i,j) + Q11(i,j))*
          (-((-2 + Q11.diff(d10,i,j))*L) + L*Q11(i,j) + 
            Q11.diff(d01,i,j)*Qr1(i,j)))/(L*pow(1 + Q11(i,j),2)) + 
       ((2 - Q22.diff(d10,i,j) + Q22(i,j))*
          (-((-2 + Q22.diff(d10,i,j))*L) + L*Q22(i,j) + 
            Q22.diff(d01,i,j)*Qr1(i,j)))/(L*pow(1 + Q22(i,j),2)) + 
       (2 - Q11.diff(d10,i,j) + Q11(i,j) + 
          (Q11.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
       (-6 + 2*Q11.diff(d10,i,j) - Q11.diff(d20,i,j) + 
          (Q11.diff(d01,i,j)*Qr1.diff(d10,i,j))/L - 2*Q11(i,j) + 
          (Q11.diff(d11,i,j)*Qr1(i,j))/L)/(1 + Q11(i,j)) + 
       (2 - Q22.diff(d10,i,j) + Q22(i,j) + 
          (Q22.diff(d01,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)) + 
       (-6 + 2*Q22.diff(d10,i,j) - Q22.diff(d20,i,j) + 
          (Q22.diff(d01,i,j)*Qr1.diff(d10,i,j))/L - 2*Q22(i,j) + 
          (Q22.diff(d11,i,j)*Qr1(i,j))/L)/(1 + Q22(i,j)) + 
       (2*(Qr1.diff(d01,i,j)*pow(1 + Qtt(i,j),2) + 
            Qr1.diff(d11,i,j)*pow(1 + Qtt(i,j),2) + 
            Qr1(i,j)*(Qtt.diff(d01,i,j)*
                (Qrr.diff(d10,i,j) - Qtt.diff(d10,i,j)) - 
               (Qrr.diff(d11,i,j) - Qtt.diff(d11,i,j))*(1 + Qtt(i,j))))\
)/(L*pow(1 + Qtt(i,j),2))))/(2.*pow(-12 + pow(mu,2),2))
;
elem[8][4]+=
((2*Qr1.diff(d01,i,j))/L + (2*L - Q11.diff(d10,i,j)*L + L*Q11(i,j) + 
       Q11.diff(d01,i,j)*Qr1(i,j))/(L + L*Q11(i,j)) + 
    (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
       Q22.diff(d01,i,j)*Qr1(i,j))/(L + L*Q22(i,j)) - 
    (2*Qtt.diff(d01,i,j)*Qr1(i,j))/(L + L*Qtt(i,j)) - 
    (2*(12 + pow(mu,2) + Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
         2*pow(mu,2)*Qtt(i,j)))/((-12 + pow(mu,2))*(1 + Qtt(i,j))) + 
    (complex(0,2)*w*(-2 + (2*Qr1.diff(d01,i,j))/L + 
         (2*L - Q11.diff(d10,i,j)*L + L*Q11(i,j) + 
            Q11.diff(d01,i,j)*Qr1(i,j))/(L + L*Q11(i,j)) + 
         (2*L - Q22.diff(d10,i,j)*L + L*Q22(i,j) + 
            Q22.diff(d01,i,j)*Qr1(i,j))/(L + L*Q22(i,j)) - 
         (2*Qtt.diff(d01,i,j)*Qr1(i,j))/(L + L*Qtt(i,j)) - 
         (2*(12 + pow(mu,2) + 
              Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
              2*pow(mu,2)*Qtt(i,j)))/
          ((-12 + pow(mu,2))*(1 + Qtt(i,j)))))/(-12 + pow(mu,2)) - 
    (2*(Q11.diff(d01,i,j)*pow(-12 + pow(mu,2),2)*(1 + Q22(i,j))*
          (1 + Qtt(i,j))*(Qtt.diff(d01,i,j)*(1 + Q22(i,j)) + 
            Q22.diff(d01,i,j)*(1 + Qtt(i,j))) + 
         8*pow(L,2)*pow(w,2)*(1 + Q11(i,j))*(1 + Q22(i,j))*
          ((-((-108 + 7*pow(mu,2) - 
                    2*Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                    (-84 + 5*pow(mu,2))*Q11(i,j))*(1 + Q22(i,j))) + 
               2*(-12 + pow(mu,2))*(1 + Q11(i,j))*
                (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j)) + 
            (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
             (Qtt.diff(d10,i,j) + Qtt(i,j))) - 
         8*pow(L,2)*pow(w,2)*(1 + Q11(i,j))*(1 + Q22(i,j))*
          (((84 - pow(mu,2) + 
                  2*Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                  (60 + pow(mu,2))*Q11(i,j))*(1 + Q22(i,j)) + 
               2*(-12 + pow(mu,2))*(1 + Q11(i,j))*
                (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j)) - 
            (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
             (Qrr.diff(d10,i,j) + Qtt(i,j)) + 
            2*(-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
             (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
         pow(-12 + pow(mu,2),2)*(1 + Q11(i,j))*
          (pow(Q22.diff(d01,i,j),2)*pow(1 + Qtt(i,j),2) - 
            (1 + Q22(i,j))*(1 + Qtt(i,j))*
             (Q22.diff(d01,i,j)*Qtt.diff(d01,i,j) + 
               2*Q22.diff(d02,i,j)*(1 + Qtt(i,j))) - 
            pow(1 + Q22(i,j),2)*
             (-pow(Qtt.diff(d01,i,j),2) + 
               2*Qtt.diff(d02,i,j)*(1 + Qtt(i,j)) + 
               2*pow(h.diff(d01,i,j),2)*pow(mu,2)*
                pow(1 + Qtt(i,j),2)))))/
     (pow(L,2)*pow(-12 + pow(mu,2),3)*pow(1 + Q11(i,j),2)*
       pow(1 + Q22(i,j),2)*(1 + Qtt(i,j))))/4.
;
elem[8][5]+=
(w*((complex(0,2) - (4*w)/(-12 + pow(mu,2)))*
       (Qr1.diff(d10,i,j) + Qr1(i,j)) + 
      (complex(0,2)*pow(Qr1(i,j),2)*
         (Q11.diff(d01,i,j)/(1 + Q11(i,j)) - 
           Qtt.diff(d01,i,j)/(1 + Qtt(i,j))))/L + 
      Qr1(i,j)*((complex(0,2)*Qr1.diff(d01,i,j))/L - 
         (4*w)/(-12 + pow(mu,2)) + 
         (complex(0,2)*(2 - Q11.diff(d10,i,j) + Q11(i,j)))/
          (1 + Q11(i,j)) + (complex(0,2)*
            (12 + pow(mu,2) + Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
              2*pow(mu,2)*Qtt(i,j)))/((-12 + pow(mu,2))*(1 + Qtt(i,j)))\
) + (complex(0,2)*(((-((2*(24 + 
                       Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                       (12 + pow(mu,2))*Q11(i,j))*(1 + Q22(i,j)) + 
                     (-12 + pow(mu,2))*(1 + Q11(i,j))*
                      (Q22.diff(d10,i,j) + Q22(i,j)))*(1 + Qtt(i,j))) \
+ (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
                 (Qrr.diff(d10,i,j) + Qtt(i,j)) - 
                (-12 + pow(mu,2))*(1 + Q11(i,j))*(1 + Q22(i,j))*
                 (Qtt.diff(d10,i,j) + Qtt(i,j)))*
              (-(Q22.diff(d01,i,j)*(1 + Q11(i,j))*(1 + Qtt(i,j))) + 
                (1 + Q22(i,j))*
                 (-2*Qtt.diff(d01,i,j)*(1 + Q11(i,j)) + 
                   Q11.diff(d01,i,j)*(1 + Qtt(i,j)))))/(1 + Qtt(i,j)) + 
           (1 + Q11(i,j))*(1 + Q22(i,j))*
            (-(Q22.diff(d11,i,j)*(-12 + pow(mu,2))*(1 + Q11(i,j))*
                 (1 + Qtt(i,j))) - 
              (-12 + pow(mu,2))*(2 - Q22.diff(d10,i,j) + Q22(i,j))*
               (-2*Qtt.diff(d01,i,j)*(1 + Q11(i,j)) + 
                 Q11.diff(d01,i,j)*(1 + Qtt(i,j))) + 
              Q22.diff(d01,i,j)*
               ((2*(-24 + pow(mu,2)) - 
                    Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                    (-36 + pow(mu,2))*Q11(i,j))*(1 + Qtt(i,j)) - 
                 (-12 + pow(mu,2))*(1 + Q11(i,j))*
                  (Qtt.diff(d10,i,j) + Qtt(i,j))) + 
              (1 + Q22(i,j))*(-2*Qtt.diff(d01,i,j)*
                  (24 + Q11.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                    (12 + pow(mu,2))*Q11(i,j)) + 
                 Q11.diff(d01,i,j)*
                  (24 + Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
                    (12 + pow(mu,2))*Qtt(i,j)) + 
                 (-12 + pow(mu,2))*
                  (-2*Qtt.diff(d11,i,j)*(1 + Q11(i,j)) + 
                    Q11.diff(d11,i,j)*(1 + Qtt(i,j)))))))/
       (L*pow(-12 + pow(mu,2),2)*pow(1 + Q11(i,j),3)*
         pow(1 + Q22(i,j),2))))/2.
;
elem[8][6]+=
(-2*mu*exp(B1*mu*h(i,j))*(-((-12 + pow(mu,2))*
         (-12 + pow(mu,2) + complex(0,2)*w)*
         (2*a0.diff(d10,i,j)*L + a0.diff(d01,i,j)*Qr1(i,j))*
         (1 + Qtt(i,j))) - L*a0(i,j)*
       ((144 - pow(mu,4) + h.diff(d10,i,j)*B1*mu*(-12 + pow(mu,2))*
             (-12 + pow(mu,2) + complex(0,2)*w) - complex(0,48)*w + 
            B1*mu*(-12 + pow(mu,2))*
             (-12 + pow(mu,2) + complex(0,2)*w)*h(i,j))*(1 + Qtt(i,j)) - 
         (-12 + pow(mu,2))*(-12 + pow(mu,2) + complex(0,2)*w)*
          (Qtt.diff(d10,i,j) + Qtt(i,j)))))/
  (L*pow(-12 + pow(mu,2),3)*pow(1 + Qtt(i,j),2))
;
elem[9][0]+=
(complex(0,-4)*L*mu*w*a0(i,j)*exp(B1*mu*h(i,j)))/
  (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j)))
;
elem[9][1]+=
(complex(0,2)*Q22.diff(d01,i,j)*w)/
  (pow(-12 + pow(mu,2),2)*(1 + Q22(i,j)))
;
elem[9][2]+=
(complex(0,-4)*h.diff(d01,i,j)*mu*w)/pow(-12 + pow(mu,2),2)
;
elem[9][3]+=
0
;
elem[9][4]+=
0
;
elem[9][5]+=
(complex(0,1)*L*(-12 + pow(mu,2) + complex(0,2)*w)*w)/
  pow(-12 + pow(mu,2),2)
;
elem[9][6]+=
0
;
elem[10][0]+=
(complex(0,-4)*mu*w*exp(B1*mu*h(i,j))*
    (-((-12 + pow(mu,2))*(-2*a0.diff(d10,i,j)*L + 
           a0.diff(d01,i,j)*Qr1(i,j))*(1 + Qtt(i,j))) + 
      L*a0(i,j)*(mu*(-4*mu + h.diff(d10,i,j)*B1*(-12 + pow(mu,2)) + 
            B1*(-12 + pow(mu,2))*h(i,j))*(1 + Qtt(i,j)) - 
         (-12 + pow(mu,2))*(Qtt.diff(d10,i,j) + Qtt(i,j)))))/
  (pow(-12 + pow(mu,2),3)*pow(1 + Qtt(i,j),2))
;
elem[10][1]+=
(complex(0,1)*w*(-4*Q22.diff(d01,i,j)*pow(mu,2)*(1 + Q22(i,j)) + 
      (-12 + pow(mu,2))*(-2*Q22.diff(d01,i,j)*
          (-1 + Q22.diff(d10,i,j)) + 
         ((1 + Q22(i,j))*(Qtt.diff(d01,i,j)*
               (Qrr.diff(d10,i,j) - Qtt.diff(d10,i,j))*(1 + Q22(i,j)) \
+ (1 + Qtt(i,j))*(-((Qrr.diff(d11,i,j) - Qtt.diff(d11,i,j))*
                    (1 + Q22(i,j))) + 2*Q22.diff(d11,i,j)*(1 + Qtt(i,j)))\
))/pow(1 + Qtt(i,j),2))))/
  (pow(-12 + pow(mu,2),3)*pow(1 + Q22(i,j),2))
;
elem[10][2]+=
(-2*mu*(-2*a0.diff(d01,i,j)*B1*mu*(-12 + pow(mu,2))*a0(i,j)*
       exp(B1*mu*h(i,j)) + (h.diff(d01,i,j)*
          (144 - 24*pow(mu,2) + pow(mu,4) - complex(0,48)*w) + 
         complex(0,2)*h.diff(d11,i,j)*(-12 + pow(mu,2))*w)*
       (1 + Qtt(i,j))))/(pow(-12 + pow(mu,2),3)*(1 + Qtt(i,j)))
;
elem[10][3]+=
(4*a0.diff(d01,i,j)*pow(mu,2)*a0(i,j)*exp(B1*mu*h(i,j))*
     (1 + Qtt(i,j)) - complex(0,1)*w*
     (Qtt.diff(d01,i,j)*(-Qrr.diff(d10,i,j) + Qtt.diff(d10,i,j)) + 
       (Qrr.diff(d11,i,j) - Qtt.diff(d11,i,j))*(1 + Qtt(i,j))))/
  (pow(-12 + pow(mu,2),2)*pow(1 + Qtt(i,j),2))
;
elem[10][4]+=
-((Qtt.diff(d01,i,j)*(-12 + pow(mu,2) + complex(0,2)*w))/
    (pow(-12 + pow(mu,2),2)*(1 + Qtt(i,j))))
;
elem[10][5]+=
(w*(-3*L*(-4 + pow(mu,2))*(complex(0,1) - (2*w)/(-12 + pow(mu,2))) + 
      (-12 + pow(mu,2))*(complex(0,1)*Qr1.diff(d01,i,j) - 
         (2*L*w)/(-12 + pow(mu,2)) + 
         (complex(0,1)*(-((-2 + Q11.diff(d10,i,j))*L) + L*Q11(i,j) + 
              Q11.diff(d01,i,j)*Qr1(i,j)))/(1 + Q11(i,j)) - 
         (complex(0,1)*Qtt.diff(d01,i,j)*Qr1(i,j))/(1 + Qtt(i,j)) + 
         (complex(0,1)*L*(12 + pow(mu,2) + 
              Qtt.diff(d10,i,j)*(-12 + pow(mu,2)) + 
              2*pow(mu,2)*Qtt(i,j)))/((-12 + pow(mu,2))*(1 + Qtt(i,j))))\
))/pow(-12 + pow(mu,2),2)
;
elem[10][6]+=
(4*a0.diff(d01,i,j)*mu*(-12 + pow(mu,2) + complex(0,2)*w)*
    exp(B1*mu*h(i,j)))/(pow(-12 + pow(mu,2),3)*(1 + Qtt(i,j)))
;
}
