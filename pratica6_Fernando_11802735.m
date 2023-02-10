%% Pratica 6 - Projeto de Cames
clear all, close all, clc,
th=[0:.1:360]/180*pi; t=th*3/pi;

%%METODO 1: POLINOMIO 3-4-5

Ex=1; h=1; betas=pi/3; betad=pi/3; betar2=pi/3; betar1=2*pi-betas-betar2-betad;
%Ex=2; h=1; betas=2*pi/3; betad=2*pi/3; betar2=pi/240; betar1=2*pi-betas-betar2-betad;
%Ex=3; h=1; betas=2*pi/3; betad=2*pi/24; betar2=pi/3; betar1=2*pi-betas-betar2-betad;

% Polinomio de ordem 5
% x=[c0,c1,c2,c3,c4,c5];
% b=[s1,v1,a1,s2,v2,a2];
% Subida 
th1=0; th2=1;
A=[[1 1*th1 1*th1^2 1*th1^3 1*th1^4 1*th1^5];
   [0 1 2*th1 3*th1^2 4*th1^3 5*th1^4]/betas;
   [0 0 2 6*th1 12*th1^2 20*th1^3]/betas^2;
   [1 1*th2 1*th2^2 1*th2^3 1*th2^4 1*th2^5];
   [0 1 2*th2 3*th2^2 4*th2^3 5*th2^4]/betas;
   [0 0 2 6*th2 12*th2^2 20*th2^3]/betas^2];
b=[0,0,0,h,0,0]';
x=A\b;
c0s=x(1); c1s=x(2); c2s=x(3); c3s=x(4); c4s=x(5); c5s=x(6); 
% Descida 
th1=0; th2=1;
A=[[1 1*th1 1*th1^2 1*th1^3 1*th1^4 1*th1^5];
   [0 1 2*th1 3*th1^2 4*th1^3 5*th1^4]/betad;
   [0 0 2 6*th1 12*th1^2 20*th1^3]/betad^2;
   [1 1*th2 1*th2^2 1*th2^3 1*th2^4 1*th2^5];
   [0 1 2*th2 3*th2^2 4*th2^3 5*th2^4]/betad;
   [0 0 2 6*th2 12*th2^2 20*th2^3]/betad^2];
b=[h,0,0,0,0,0]';
x=A\b;
c0d=x(1); c1d=x(2); c2d=x(3); c3d=x(4); c4d=x(5); c5d=x(6); 

% Monta os graficos SVAJ
% Deslocamento S
s=zeros(size(th));
th1i=0; th1f=th1i+betar1;    id1=find(th<=th1f); s(id1)=0;
th2i=th1f; th2f=th2i+betas;  id2=find(th>=th2i&th<=th2f); thb=(th(id2)-th2i)/betas; 
s(id2)=c0s+c1s*thb+c2s*thb.^2+c3s*thb.^3+c4s*thb.^4+c5s*thb.^5;
th3i=th2f; th3f=th3i+betar2; id3=find(th>=th3i&th<=th3f); s(id3)=h;
th4i=th3f; th4f=th4i+betad;  id4=find(th>=th4i&th<=th4f); thb=(th(id4)-th4i)/betad; 
s(id4)=c0d+c1d*thb+c2d*thb.^2+c3d*thb.^3+c4d*thb.^4+c5d*thb.^5;
% Velocidade V
v=zeros(size(th));
th1i=0; th1f=th1i+betar1;    id1=find(th<=th1f); v(id1)=0;
th2i=th1f; th2f=th2i+betas;  id2=find(th>=th2i&th<=th2f); thb=(th(id2)-th2i)/betas; 
v(id2)=(1/betas)*(c1s+2*c2s*thb+3*c3s*thb.^2+4*c4s*thb.^3+5*c5s*thb.^4);
th3i=th2f; th3f=th3i+betar2; id3=find(th>=th3i&th<=th3f); v(id3)=0;
th4i=th3f; th4f=th4i+betad;  id4=find(th>=th4i&th<=th4f); thb=(th(id4)-th4i)/betad; 
v(id4)=(1/betad)*(c1d+2*c2d*thb+3*c3d*thb.^2+4*c4d*thb.^3+5*c5d*thb.^4);
% Aceleracao A
a=zeros(size(th));
th1i=0; th1f=th1i+betar1;    id1=find(th<=th1f); a(id1)=0;
th2i=th1f; th2f=th2i+betas;  id2=find(th>=th2i&th<=th2f); thb=(th(id2)-th2i)/betas; 
a(id2)=(1/betas^2)*(2*c2s+6*c3s*thb+12*c4s*thb.^2+20*c5s*thb.^3);
th3i=th2f; th3f=th3i+betar2; id3=find(th>=th3i&th<=th3f); a(id3)=0;
th4i=th3f; th4f=th4i+betad;  id4=find(th>=th4i&th<=th4f); thb=(th(id4)-th4i)/betad; 
a(id4)=(1/betad^2)*(2*c2d+6*c3d*thb+12*c4d*thb.^2+20*c5d*thb.^3);
% Jerk J
j=zeros(size(th));
th1i=0; th1f=th1i+betar1;    id1=find(th<=th1f); j(id1)=0;
th2i=th1f; th2f=th2i+betas;  id2=find(th>=th2i&th<=th2f); thb=(th(id2)-th2i)/betas; 
j(id2)=(1/betas^3)*(6*c3s+24*c4s*thb+60*c5s*thb.^2);
th3i=th2f; th3f=th3i+betar2; id3=find(th>=th3i&th<=th3f); j(id3)=0;
th4i=th3f; th4f=th4i+betad;  id4=find(th>=th4i&th<=th4f); thb=(th(id4)-th4i)/betad; 
j(id4)=(1/betad^3)*(6*c3d+24*c4d*thb+60*c5d*thb.^2);

% Plot SVAJ
figure(1), set(1,'position',[680 565 644 420]),
subplot(411), plot(t,s), ylabel('S (cm)'), axis tight,
subplot(412), plot(t,v), ylabel('V (cm/s)'), axis tight,
subplot(413), plot(t,a), ylabel('A (cm/s^2)'), axis tight,
subplot(414), plot(t,j), ylabel('J (cm/s^3)'), axis tight,
xlabel('Tempo (s)')

% Draw cam
figure(2), set(2,'position',[1325 565 560 420])
subplot(221), Rb=h/2; polar(th,Rb+s), title('R_b=h/2'),
subplot(222), Rb=h/1; polar(th,Rb+s), title('R_b=h'),
subplot(223), Rb=2*h; polar(th,Rb+s), title('R_b=2h'),
subplot(224), Rb=3*h; polar(th,Rb+s), title('R_b=3h'),
display("Deslocamento, Velocidade, Aceleração e Jerk máximos e mínimos para o método 1:")
s1_max=max(s)
v1_max=max(v)
a1_max=max(a)
j1_max=max(j)
s1_min=min(s)
v1_min=min(v)
a1_min=min(a)
j1_min=min(j)

%%METODO 2: POLINOMIO 4-5-6-7

% Polinomio de ordem 7
% x=[c0,c1,c2,c3,c4,c5,c6,c7];
% b=[s1,v1,a1,j1,s2,v2,a2,j2];
% Subida 
th1=0; th2=1;

A=[[1 1*th1 1*th1^2 1*th1^3 1*th1^4 1*th1^5 1*th1^6 1*th1^7];
   [0 1 2*th1 3*th1^2 4*th1^3 5*th1^4 6*th1^5 7*th1^6 ]/betas;
   [0 0 2 6*th1 12*th1^2 20*th1^3 30*th1^4 42*th1^5]/betas^2;
   [0 0 0 6 24*th1  60*th1^2 120*th1^3 210*th1^4]/betas^3;
   [1 1*th2 1*th2^2 1*th2^3 1*th2^4 1*th2^5 1*th2^6 1*th2^7];
   [0 1 2*th2 3*th2^2 4*th2^3 5*th2^4 6*th2^5 7*th2^6 ]/betas;
   [0 0 2 6*th2 12*th2^2 20*th2^3 30*th2^4 42*th2^5]/betas^2;
   [0 0 0 6 24*th2  60*th2^2 120*th2^3 210*th2^4]/betas^3;];
b=[0,0,0,0,h,0,0,0]';
x=A\b;
c0s=x(1); c1s=x(2); c2s=x(3); c3s=x(4); c4s=x(5); c5s=x(6);  c6s=x(7); c7s=x(8);
% Descida
th1=0; th2=1;
A=[[1 1*th1 1*th1^2 1*th1^3 1*th1^4 1*th1^5 1*th1^6 1*th1^7];
   [0 1 2*th1 3*th1^2 4*th1^3 5*th1^4 6*th1^5 7*th1^6 ]/betad;
   [0 0 2 6*th1 12*th1^2 20*th1^3 30*th1^4 42*th1^5]/betad^2;
   [0 0 0 6 24*th1  60*th1^2 120*th1^3 210*th1^4]/betad^3;
   [1 1*th2 1*th2^2 1*th2^3 1*th2^4 1*th2^5 1*th2^6 1*th2^7];
   [0 1 2*th2 3*th2^2 4*th2^3 5*th2^4 6*th2^5 7*th2^6 ]/betad;
   [0 0 2 6*th2 12*th2^2 20*th2^3 30*th2^4 42*th2^5]/betad^2;
   [0 0 0 6 24*th2  60*th2^2 120*th2^3 210*th2^4]/betad^3;];
b=[h,0,0,0,0,0,0,0]';
x=A\b;
c0d=x(1); c1d=x(2); c2d=x(3); c3d=x(4); c4d=x(5); c5d=x(6);  c6d=x(7); c7d=x(8);

% Monta os graficos SVAJ
% Deslocamento S
s=zeros(size(th));
th1i=0; th1f=th1i+betar1;    id1=find(th<=th1f); s(id1)=0;
th2i=th1f; th2f=th2i+betas;  id2=find(th>=th2i&th<=th2f); thb=(th(id2)-th2i)/betas; 
s(id2)=c0s+c1s*thb+c2s*thb.^2+c3s*thb.^3+c4s*thb.^4+c5s*thb.^5+c6s*thb.^6+c7s*thb.^7;
th3i=th2f; th3f=th3i+betar2; id3=find(th>=th3i&th<=th3f); s(id3)=h;
th4i=th3f; th4f=th4i+betad;  id4=find(th>=th4i&th<=th4f); thb=(th(id4)-th4i)/betad; 
s(id4)=c0d+c1d*thb+c2d*thb.^2+c3d*thb.^3+c4d*thb.^4+c5d*thb.^5+c6d*thb.^6+c7d*thb.^7;
% Velocidade V
v=zeros(size(th));
th1i=0; th1f=th1i+betar1;    id1=find(th<=th1f); v(id1)=0;
th2i=th1f; th2f=th2i+betas;  id2=find(th>=th2i&th<=th2f); thb=(th(id2)-th2i)/betas; 
v(id2)=(1/betas)*(c1s+2*c2s*thb+3*c3s*thb.^2+4*c4s*thb.^3+5*c5s*thb.^4+6*c6s*thb.^5+7*c7s*thb.^6);
th3i=th2f; th3f=th3i+betar2; id3=find(th>=th3i&th<=th3f); v(id3)=0;
th4i=th3f; th4f=th4i+betad;  id4=find(th>=th4i&th<=th4f); thb=(th(id4)-th4i)/betad; 
v(id4)=(1/betad)*(c1d+2*c2d*thb+3*c3d*thb.^2+4*c4d*thb.^3+5*c5d*thb.^4+6*c6d*thb.^5+7*c7d*thb.^6);
% Aceleracao A
a=zeros(size(th));
th1i=0; th1f=th1i+betar1;    id1=find(th<=th1f); a(id1)=0;
th2i=th1f; th2f=th2i+betas;  id2=find(th>=th2i&th<=th2f); thb=(th(id2)-th2i)/betas; 
a(id2)=(1/betas^2)*(2*c2s+6*c3s*thb+12*c4s*thb.^2+20*c5s*thb.^3+30*c6s*thb.^4+42*c7s*thb.^5);
th3i=th2f; th3f=th3i+betar2; id3=find(th>=th3i&th<=th3f); a(id3)=0;
th4i=th3f; th4f=th4i+betad;  id4=find(th>=th4i&th<=th4f); thb=(th(id4)-th4i)/betad; 
a(id4)=(1/betad^2)*(2*c2d+6*c3d*thb+12*c4d*thb.^2+20*c5d*thb.^3+30*c6d*thb.^4+42*c7d*thb.^5);
% Jerk J
j=zeros(size(th));
th1i=0; th1f=th1i+betar1;    id1=find(th<=th1f); j(id1)=0;
th2i=th1f; th2f=th2i+betas;  id2=find(th>=th2i&th<=th2f); thb=(th(id2)-th2i)/betas; 
j(id2)=(1/betas^3)*(6*c3s+24*c4s*thb+60*c5s*thb.^2+120*c6s*thb.^3+210*c7s*thb.^4);
th3i=th2f; th3f=th3i+betar2; id3=find(th>=th3i&th<=th3f); j(id3)=0;
th4i=th3f; th4f=th4i+betad;  id4=find(th>=th4i&th<=th4f); thb=(th(id4)-th4i)/betad; 
j(id4)=(1/betad^3)*(6*c3d+24*c4d*thb+60*c5d*thb.^2+120*c6d*thb.^3+210*c7d*thb.^4);

% Plot SVAJ
figure(3), set(3,'position',[680 565 644 420]),
subplot(411), plot(t,s), ylabel('S (cm)'), axis tight,
subplot(412), plot(t,v), ylabel('V (cm/s)'), axis tight,
subplot(413), plot(t,a), ylabel('A (cm/s^2)'), axis tight,
subplot(414), plot(t,j), ylabel('J (cm/s^3)'), axis tight,
xlabel('Tempo (s)')

% Draw cam
figure(4), set(4,'position',[1325 565 560 420])
subplot(221), Rb=h/2; polar(th,Rb+s), title('R_b=h/2'),
subplot(222), Rb=h/1; polar(th,Rb+s), title('R_b=h'),
subplot(223), Rb=2*h; polar(th,Rb+s), title('R_b=2h'),
subplot(224), Rb=3*h; polar(th,Rb+s), title('R_b=3h'), 

display("Deslocamento, Velocidade, Aceleração e Jerk máximos e mínimos para o método 2:")
s2_max=max(s)
v2_max=max(v)
a2_max=max(a)
j2_max=max(j)
s2_min=min(s)
v2_min=min(v)
a2_min=min(a)
j2_min=min(j)
