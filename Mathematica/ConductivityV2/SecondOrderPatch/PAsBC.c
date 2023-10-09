
elem[11]=
h11.diff(d02,0,j)/pow(L,2) - h22.diff(d02,0,j)/pow(L,2) - 
  htt.diff(d02,0,j)/pow(L,2) - (complex(0,2)*ht1.diff(d01,0,j)*w)/L + 
  pow(w,2)*h11(0,j)
;
elem[12]=
(complex(0,-1)*a0.diff(d02,0,j)*h11.diff(d01,0,j)*mu)/(pow(L,2)*w) + 
  (complex(0,1)*a0.diff(d02,0,j)*h22.diff(d01,0,j)*mu)/(pow(L,2)*w) + 
  (complex(0,1)*a0.diff(d02,0,j)*htt.diff(d01,0,j)*mu)/(pow(L,2)*w) - 
  complex(0,1)*A0.diff(d01,0,j)*w + L*pow(w,2)*Ax(0,j) + 
  complex(0,1.5)*a0.diff(d01,0,j)*mu*w*h11(0,j) - 
  complex(0,0.5)*a0.diff(d01,0,j)*mu*w*h22(0,j) - 
  (2*a0.diff(d02,0,j)*mu*ht1(0,j))/L - 
  complex(0,0.5)*a0.diff(d01,0,j)*mu*w*htt(0,j)
;
elem[13]=
(h11.diff(d01,0,j)*h.diff(d01,0,j)*mu)/(pow(L,2)*pow(w,2)) - 
  (h22.diff(d01,0,j)*h.diff(d01,0,j)*mu)/(pow(L,2)*pow(w,2)) - 
  (h.diff(d01,0,j)*htt.diff(d01,0,j)*mu)/(pow(L,2)*pow(w,2)) + 
  Dh(0,j) - (mu*h(0,j)*h11(0,j))/2. + (mu*h(0,j)*h22(0,j))/2. - 
  (complex(0,2)*h.diff(d01,0,j)*mu*ht1(0,j))/(L*w)
;
