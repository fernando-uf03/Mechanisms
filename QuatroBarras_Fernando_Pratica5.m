clear all, close all, clc,
%% Simulacao de mecanismo de quatro barras
outvid=0;

p21=(2.138^2+1.236^2)^(1/2);
p31=(2.500^2+2.931^2)^(1/2);

beta2=pi/6;
beta3=pi/3;
alfa2=-1.091;
alfa3=-1.742;
delta2=2.095;
delta3=2.277;
gama2=-0.1745;
gama3=0.43633;

A1=cos(beta2)-1;
B1=sin(beta2);
C=cos(alfa2)-1;
D=sin(alfa2);
E=p21*cos(delta2);
F1=cos(beta3)-1;
G1=sin(beta3);
H=cos(alfa3)-1;
K=sin(alfa3);
L=p31*cos(delta3);
M=p21*sin(delta2);
N=p31*sin(delta3);

A2=cos(gama2)-1;
B2=sin(gama2);
F2=cos(gama3)-1;
G2=sin(gama3);

matriz_A1=[A1 -B1 C -D;F1 -G1 H -K;B1 A1 D C;G1 F1 K H];
matriz_b=[E L M N]';
matriz_x1=matriz_A1\matriz_b;

matriz_A2=[A2 -B2 C -D;F2 -G2 H -K;B2 A2 D C;G2 F2 K H];
matriz_x2=matriz_A2\matriz_b;

W=((matriz_x1(1,1)^2)+(matriz_x1(2,1)^2))^(1/2);%norma de W
Z=((matriz_x1(3,1)^2)+(matriz_x1(4,1)^2))^(1/2);%norma de Z
U=((matriz_x2(1,1)^2)+(matriz_x2(2,1)^2))^(1/2);%norma de U
S=((matriz_x2(3,1)^2)+(matriz_x2(4,1)^2))^(1/2);%norma de S

W_1x=matriz_x1(1,1); W_1y=matriz_x1(2,1);%Componentes x e y de W 
Z_1x=matriz_x1(3,1); Z_1y=matriz_x1(4,1);%Componentes x e y de Z
U_1x=matriz_x2(1,1); U_1y=matriz_x2(2,1);%Componentes x e y de U
S_1x=matriz_x2(3,1); S_1y=matriz_x2(4,1);%Componentes x e y de S

vetorW=[W_1x W_1y];
vetorZ=[Z_1x Z_1y];
vetorU=[U_1x U_1y];
vetorS=[S_1x S_1y];

%vetor do elo d (ou elo 1)
Vx=Z_1x-S_1x;
Vy=Z_1y-S_1y;
V=[Vx Vy];
Gx=W_1x+Vx-U_1x;
Gy=W_1y+Vy-U_1y;
G=norm([Gx Gy]);

%%Valores dos angulos theta
%atan2() ? uma fun??o do matlab que calcula o arctan do vetor com rela??o ?
%horizontal
theta0=atan2(W_1y,W_1x); %theta0: ?ngulo do elo 2 com horizontal
theta1=atan2(Gy,Gx);     %theta1: ?ngulo do elo 1 com horizontal
theta3=atan2(Vy,Vx);     %theta3: ?ngulo do elo 3 com horizontal
display('Theta2 inicial:')
theta2i=theta0-theta1
display('Theta2 final:')
theta2f=theta2i+beta3

%Valores de outros angulos
sigma=atan2(U_1y,U_1x);
psi=atan2(S_1y,S_1x);
phi=atan2(Z_1y,Z_1x);

display('Posi??o dos pivos fixos')
%pivo fixo O1
O1x=-Z*cos(phi)-W*cos(theta0)
O1y=-Z*sin(phi)-W*sin(theta0)
%pivo fixo O2
O2x=-S*cos(psi)-U*cos(sigma) 
O2y=-S*sin(psi)-U*sin(sigma)

display('Tamanho dos elos:')
% Lengths of links
 a=W
 b=norm(V)
 c=U
 d=G % Grashof Crank-Rocker
 e=norm(vetorZ)
 
  

% Grashof condition
lAll=[a b c d];
sL=min(lAll); lL=max(lAll); sLi=find(lAll==sL); lLi=find(lAll==lL); lOth=lAll; lOth([sLi lLi])=[];
if sL+lL>sum(lOth) disp('non-Grashof driver'), nGf=1; else nGf=0; end
% End-effector P
   BP=Z; tBP=S;
%   BP=0.5; tBP=0.5;
%   BP=0.7; tBP=0.5;
% Number of steps
N=1000;
% Vector of time instants
t=linspace(0,1,N)';

%Estudo cinem?tico
w2=4*pi; Alfa2=0;

t2v = [linspace(theta2i,theta2f,N/4) linspace(theta2f,theta2i,N/4) linspace(theta2i,theta2f,N/4) linspace(theta2f,theta2i,N/4)]';

%% Newton-Raphson algorithm for the evaluation of t3 and t4 given t2
tol=0.001;
t3=theta3;t4=acos((S.^2+V.^2-Z.^2)/(2*S*V))   
for it2=1:length(t2v)
   t2=t2v(it2); B=tol+1; iconv=0;
   while norm(B)>tol
       iconv=iconv+1;
       A=[-b*sin(t3) c*sin(t4);b*cos(t3) -c*cos(t4)];
       B=[a*cos(t2)+b*cos(t3)-c*cos(t4)-d; a*sin(t2)+b*sin(t3)-c*sin(t4)];
       Dt=-A\B;
       t3=t3+Dt(1); t4=t4+Dt(2);
   end
   t3v(it2,1)=t3; t4v(it2,1)=t4;   
end

%Declara??o dos vetores vel angular w3 e w4
w3=[]; 
w4=[];
%Declara??o dos vetores acel angular a3 e a4
a3=[];
a4=[];
%Declara??o dos vetores velocidade linear Va e Vb
Va=[];
Vb=[];
%Declara??o dos vetores acelera??o linear Aa e Ab
Aa=[];
Ab=[];

for i=1:length(t)
    Theta2=t2v(i);
    Theta3=t3v(i);
    Theta4=t4v(i);
%%C?lculo e armazenamento dos valores de vel angulares a cada passo 
%%de tempo nos vetores w3 e w4
    w3(i)=(a*w2/b)*(sin(Theta4-Theta2)/sin(Theta3-Theta4)); %omega 3
    w4(i)=(a*w2/c)*(sin(Theta2-Theta3)/sin(Theta4-Theta3)); %omega 4
   
%%C?lculo e armazenamento dos valores de aceleracoes angulares a cada passo 
%%de tempo nos vetores a3 e a4
    A=c*sin(Theta4); %A,B,C,D,E,F s?o as equa??es usadas para calcular alfa 3 e 4 nas linhas 83 e 84
    B=b*sin(Theta3);
    C=a*0*sin(Theta2)+a*(w2^2)*cos(Theta2)+b*(w3(i)^2)*cos(Theta3)-c*(w4(i)^2)*cos(Theta4); 
    D=c*cos(Theta4);
    E=b*cos(Theta3);
    F=a*0*cos(Theta2)-a*(w2^2)*sin(Theta2)-b*(w3(i)^2)*sin(Theta3)+c*(w4(i)^2)*sin(Theta4);
    a3(i)=(C*D-A*F)/(A*E-B*D); %alfa 3
    a4(i)=(C*E-B*F)/(A*E-B*D); %alfa 4
    
%%C?lculo e armazenamento dos valores de vel lineares a cada passo 
%%de tempo nos vetores Va e Vb
    Va_i=a*w2*(-sin(Theta2)); %Dire??o do versor i
    Va_j=a*w2*(cos(Theta2));  %Dire??o do versor j
    Va(i)=((Va_i^2)+(Va_j^2))^(1/2); %M?dulo da velocidade
    Vb_i=c*w4(i)*(-sin(Theta4)); %Dire??o do versor i
    Vb_j=c*w4(i)*(cos(Theta4));  %Dire??o do versor j
    Vb(i)=((Vb_i^2)+(Vb_j^2))^(1/2);  %M?dulo da velocidade
    
%%C?lculo e armazenamento dos valores de acelera??o lineares a cada passo 
%%de tempo nos vetores Aa e Ab
    Aa_i=-(a*0*sin(Theta2)+a*(w2^2)*cos(Theta2)); %%OBS: alfa 2 ? nula %%Dire??o do versor i
    Aa_j=a*0*cos(Theta2)-a*(w2^2)*sin(Theta2);    %Dire??o do versor j
    Aa(i)=((Aa_i^2)+(Aa_j^2))^(1/2);              %M?dulo da acelera??o
    Ab_i=-(c*a4(i)*sin(Theta4)+c*(w4(i)^2)*cos(Theta4));  %Dire??o do versor i
    Ab_j=c*a4(i)*cos(Theta4)-c*(w4(i)^2)*sin(Theta4);     %Dire??o do versor j
    Ab(i)=((Ab_i^2)+(Ab_j^2))^(1/2);                      %M?dulo da acelera??o
end

%Plotagem dos gr?ficos das velocidades angulares w3 e w4
figure(3)
plot(t,w3,'k-');
xlabel('Tempo (s)'); ylabel('Velocidade angular (Rad/s)');
title('Velocidade angular vs. Tempo'); grid on;
hold on;
plot(t,w4,'b-');
legend('Elo 3','Elo 4','location','northwest');
axis([0 1 -60 30]);
hold off;

%Plotagem dos gr?ficos das acelera??es angulares a3 e a4
figure(4)
plot(t,a3,'r-');
xlabel('Tempo (s)'); ylabel('Acelera??o angular (Rad/s^2)');
title('Acelera??o angular vs. Tempo'); grid on;
hold on;
plot(t,a4,'b-');
legend('Elo 3','Elo 4','location','northwest');
axis([0 1 -3000 4000]);
hold off;

%Plotagem dos gr?ficos das velocidades lineares Va e Vb
figure(5)
plot(t,Va,'b-');
xlabel('Tempo (s)'); ylabel('Velocidade linear (m/s)');
title('Velocidade linear vs. Tempo'); grid on;
hold on;
plot(t,Vb,'r');
legend('Ponto A','Ponto B','location','northwest');
axis([0 1 0 150]);
hold off;

%Plotagem dos gr?ficos das acelera??es lineares Aa e Ab
figure(6)
plot(t,Aa,'r');
xlabel('Tempo (s)'); ylabel('Acelera??o linear (m/s^2)');
title('Acelera??o linear vs. Tempo'); grid on;
hold on;
plot(t,Ab,'b');
legend('Ponto A','Ponto B','location','northwest');
axis([0 1 0 12000]);
hold off;


%% Evaluate positions of points A,B,C,D,P
rA=zeros(length(t2v),2);
rB=a*[cos(t2v) sin(t2v)];
rC=rB+b*[cos(t3v) sin(t3v)];
rD=[rA(:,1)+d rA(:,2)];
rP=rB+Z*[cos(t3v+(phi-theta3)) sin(t3v+(phi-theta3))];
tol=0.003;
posP2 = find(t2v>theta2i+beta2-tol & t2v>theta2i+beta2+tol);



%% Show simulation
figure(1), clf, set(1,'position',[1169 335 690 650])
if outvid==1
    filename='QuatroBarras_Fernando_Pratica5.gif';
    frame=getframe(gcf); im=frame2im(frame); [A,map]=rgb2ind(im,256);
    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1e-3);
end
rT=[rA; rB; rC; rD; rP]; mx=max(max(abs(rT)));
rT=[rA; rB; rC; rD; rP]; mx=max(max(abs(rT))); axis([-mx mx -mx mx]), 
axis equal,
hAB=line([rA(1,1) rB(1,1)],[rA(1,2) rB(1,2)]); set(hAB,'Color',.7*[1 1 1],'LineStyle','-'),
hBC=line([rB(1,1) rC(1,1)],[rB(1,2) rC(1,2)]); set(hBC,'Color',.7*[1 1 1],'LineStyle','-'),
hCD=line([rC(1,1) rD(1,1)],[rC(1,2) rD(1,2)]); set(hCD,'Color',.7*[1 1 1],'LineStyle','-'),
hDA=line([rD(1,1) rA(1,1)],[rD(1,2) rA(1,2)]); set(hDA,'Color',.7*[1 1 1],'LineStyle','-'),
hBP1=line([rB(1,1) rP(1,1)],[rB(1,2) rP(1,2)]); set(hBP1,'Color',.7*[1 1 1],'LineStyle','-'),
hBP2=line([rB((length(rB)/4),1) rP((length(rP)/4),1)],[rB((length(rB)/4),2) rP((length(rP)/4),2)]); set(hBP2,'Color',.7*[1 1 1],'LineStyle','-'),
hBP3=line([rB(posP2(1,1),1) rP(posP2(1,1),1)],[rB(posP2(1,1),2) rP(posP2(1,1),2)]); set(hBP3,'Color',.7*[1 1 1],'LineStyle','-'),
hBP1=line([rB(1,1) rP(1,1)],[rB(1,2) rP(1,2)]); set(hBP1,'Color',.7*[1 1 1],'LineStyle','-'),

text(rA(1,1)-mx/100,rA(1,2)-mx/20,'A')
text(rB(1,1)-mx/100,rB(1,2)-mx/20,'B')
text(rC(1,1)-mx/100,rC(1,2)-mx/20,'C')
text(rD(1,1)-mx/100,rD(1,2)-mx/20,'D')

if outvid==1 dts=5; else dts=1; end
for n=2:dts:length(t)
    set(hAB,'xdata',[rA(n,1) rB(n,1)],'ydata',[rA(n,2) rB(n,2)]);
    set(hBC,'xdata',[rB(n,1) rC(n,1)],'ydata',[rB(n,2) rC(n,2)]);
    set(hCD,'xdata',[rC(n,1) rD(n,1)],'ydata',[rC(n,2) rD(n,2)]);
    set(hBP1,'xdata',[rB(n,1) rP(n,1)],'ydata',[rB(n,2) rP(n,2)]);
    hP=line([rP(n-1,1) rP(n,1)],[rP(n-1,2) rP(n,2)]); set(hP,'Color','r','LineStyle',':');
    if outvid==1 
        frame=getframe(gcf); im=frame2im(frame); [A,map]=rgb2ind(im,256);
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1e-3);
    else pause(1e-12);
    end
end