
elem[11][0]=
0
;
elem[11][1]=
0
;
elem[11][2]=
0
;
elem[11][3]=
0
;
elem[11][4]=
0
;
elem[11][5]=
0
;
elem[11][6]=
0
;
elem[12][0]=
0
;
elem[12][1]=
0
;
elem[12][2]=
0
;
elem[12][3]=
0
;
elem[12][4]=
0
;
elem[12][5]=
0
;
elem[12][6]=
0
;
elem[13][0]=
0
;
elem[13][1]=
0
;
elem[13][2]=
0
;
elem[13][3]=
0
;
elem[13][4]=
0
;
elem[13][5]=
0
;
elem[13][6]=
0
;
if(j==j1){
elem[11][0]+=
0
;
elem[11][1]+=
0
;
elem[11][2]+=
0
;
elem[11][3]+=
0
;
elem[11][4]+=
0
;
elem[11][5]+=
0
;
elem[11][6]+=
0
;
elem[12][0]+=
0
;
elem[12][1]+=
0
;
elem[12][2]+=
0
;
elem[12][3]+=
0
;
elem[12][4]+=
0
;
elem[12][5]+=
0
;
elem[12][6]+=
0
;
elem[13][0]+=
0
;
elem[13][1]+=
0
;
elem[13][2]+=
0
;
elem[13][3]+=
0
;
elem[13][4]+=
0
;
elem[13][5]+=
0
;
elem[13][6]+=
0
;
}
if(0==i1){
elem[11][0]+=
0
;
elem[11][1]+=
D[-1 + d02](0,j,i1,j1)/pow(L,2)
;
elem[11][2]+=
0
;
elem[11][3]+=
-(D[-1 + d02](0,j,i1,j1)/pow(L,2))
;
elem[11][4]+=
-(D[-1 + d02](0,j,i1,j1)/pow(L,2))
;
elem[11][5]+=
(complex(0,-2)*w*D[-1 + d01](0,j,i1,j1))/L
;
elem[11][6]+=
0
;
elem[12][0]+=
0
;
elem[12][1]+=
(complex(0,-1)*a0.diff(d02,0,j)*mu*D[-1 + d01](0,j,i1,j1))/(pow(L,2)*w)
;
elem[12][2]+=
0
;
elem[12][3]+=
(complex(0,1)*a0.diff(d02,0,j)*mu*D[-1 + d01](0,j,i1,j1))/(pow(L,2)*w)
;
elem[12][4]+=
(complex(0,1)*a0.diff(d02,0,j)*mu*D[-1 + d01](0,j,i1,j1))/(pow(L,2)*w)
;
elem[12][5]+=
0
;
elem[12][6]+=
complex(0,-1)*w*D[-1 + d01](0,j,i1,j1)
;
elem[13][0]+=
0
;
elem[13][1]+=
(h.diff(d01,0,j)*mu*D[-1 + d01](0,j,i1,j1))/(pow(L,2)*pow(w,2))
;
elem[13][2]+=
0
;
elem[13][3]+=
-((h.diff(d01,0,j)*mu*D[-1 + d01](0,j,i1,j1))/(pow(L,2)*pow(w,2)))
;
elem[13][4]+=
-((h.diff(d01,0,j)*mu*D[-1 + d01](0,j,i1,j1))/(pow(L,2)*pow(w,2)))
;
elem[13][5]+=
0
;
elem[13][6]+=
0
;
}
if(j==j1 && 0==i1){
elem[11][0]+=
0
;
elem[11][1]+=
pow(w,2)
;
elem[11][2]+=
0
;
elem[11][3]+=
0
;
elem[11][4]+=
0
;
elem[11][5]+=
0
;
elem[11][6]+=
0
;
elem[12][0]+=
L*pow(w,2)
;
elem[12][1]+=
complex(0,1.5)*a0.diff(d01,0,j)*mu*w
;
elem[12][2]+=
0
;
elem[12][3]+=
complex(0,-0.5)*a0.diff(d01,0,j)*mu*w
;
elem[12][4]+=
complex(0,-0.5)*a0.diff(d01,0,j)*mu*w
;
elem[12][5]+=
(-2*a0.diff(d02,0,j)*mu)/L
;
elem[12][6]+=
0
;
elem[13][0]+=
0
;
elem[13][1]+=
-(mu*h(0,j))/2.
;
elem[13][2]+=
1
;
elem[13][3]+=
0
;
elem[13][4]+=
(mu*h(0,j))/2.
;
elem[13][5]+=
(complex(0,-2)*h.diff(d01,0,j)*mu)/(L*w)
;
elem[13][6]+=
0
;
}
