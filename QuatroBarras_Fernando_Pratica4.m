clear all, close all, clc,
%% Simulacao de mecanismo de quatro barras
outvid=0;
% Lengths of links
  a=0.03125; b=0.12875; c=0.07125; d=0.10875; % Grashof Crank-Rocker
  m2=0.09375; m3=0.38625; m4=0.21375; % Massas dos elos em kg
  Ig2=7.629e-6; Ig3=5.335e-4; Ig4=9.043e-5; %Momentos de inércia dos elos em kg*m^2
%  a=0.14; b=0.17; c=0.32; d=0.50; % non-Grashof Triple-Rocker
%  a=0.37; b=0.32; c=0.50; d=0.14; % Grashof Double-Crank Inversion
% Grashof condition
lAll=[a b c d];
sL=min(lAll); lL=max(lAll); sLi=find(lAll==sL); lLi=find(lAll==lL); lOth=lAll; lOth([sLi lLi])=[];
if sL+lL>sum(lOth) disp('non-Grashof driver'), nGf=1; else nGf=0; end
% End-effector P
   BP=0.1; tBP=0.5;
%   BP=0.5; tBP=0.5;
%   BP=0.7; tBP=0.5;
% Number of steps
N=1000;
% Vector of time instants
t=linspace(0,1,N)';

%% Setup of input
% Initial t2
t2i=40/180*pi;
% Vector of values for t2 (for imposed t2 in constant steps)
% t2v=t2i+linspace(0,4*pi,N)';
  % Vector of values for t2 (for imposed w2 constant)
   w2=4*pi; t2v=t2i+w2*t; w2v=w2*ones(size(t)); a2v=zeros(size(t));
    % Vector of values for t2 (for imposed a2 constant)
    % a2=8*pi; t2v=t2i+a2/2*t.^2; a2v=a2*ones(size(t)); w2v=a2*t;
% Vector of values for t2 (for imposed t2 oscillatory for non-Grashof drivers)
if nGf==1 
    t2max=acos((-(b+c)^2+a^2+d^2)/(2*a*d));
    t2v=[linspace(0,t2max,N/4) linspace(t2max,-t2max,N/2) linspace(-t2max,0,N/4)]';
end

%% Newton-Raphson algorithm for the evaluation of t3 and t4 given t2
tol=0.001;
% Initial guesses for t3 and t4
  t3=pi/3; t4=2*pi/3; % Circuito aberto
  %t3=-pi/3; t4=-2*pi/3; % Circuito cruzado
disp('    it2   Iterations')
    
for it2=1:length(t2v)
   t2=t2v(it2); B=tol+1; iconv=0;
   while norm(B)>tol
       iconv=iconv+1;
       A=[-b*sin(t3) c*sin(t4);b*cos(t3) -c*cos(t4)];
       B=[a*cos(t2)+b*cos(t3)-c*cos(t4)-d; a*sin(t2)+b*sin(t3)-c*sin(t4)];
       Dt=-A\B;
       t3=t3+Dt(1); t4=t4+Dt(2);
   end
   if iconv>2 disp([it2 iconv]), end % Show number of iterations required to converge
   t3v(it2,1)=t3; t4v(it2,1)=t4;   
end

%Declaração dos vetores vel angular w3 e w4
w3=[]; 
w4=[];
%Declaração dos vetores acel angular a3 e a4
a3=[];
a4=[];
%Declaração dos vetores velocidade linear Va e Vb
Va=[];
Vb=[];
%Declaração dos vetores aceleração linear Aa e Ab
Aa=[];
Ab=[];
%Declaração dos vetores força e torque que agem nos elos
F12x=[];
F12y=[];
F32x=[];
F32y=[];
F43x=[];
F43y=[];
F14x=[];
F14y=[];
T12=[];

for i=1:length(t)
    theta2=t2v(i);
    theta3=t3v(i);
    theta4=t4v(i);
%%Cálculo e armazenamento dos valores de vel angulares a cada passo 
%%de tempo nos vetores w3 e w4
    w3(i)=(a*w2/b)*(sin(theta4-theta2)/sin(theta3-theta4)); %omega 3
    w4(i)=(a*w2/c)*(sin(theta2-theta3)/sin(theta4-theta3)); %omega 4
   
%%Cálculo e armazenamento dos valores de aceleracoes angulares a cada passo 
%%de tempo nos vetores a3 e a4
    A=c*sin(theta4); %A,B,C,D,E,F são as equações usadas para calcular alfa 3 e 4 nas linhas 83 e 84
    B=b*sin(theta3);
    C=a*0*sin(theta2)+a*(w2^2)*cos(theta2)+b*(w3(i)^2)*cos(theta3)-c*(w4(i)^2)*cos(theta4); 
    D=c*cos(theta4);
    E=b*cos(theta3);
    F=a*0*cos(theta2)-a*(w2^2)*sin(theta2)-b*(w3(i)^2)*sin(theta3)+c*(w4(i)^2)*sin(theta4);
    a3(i)=(C*D-A*F)/(A*E-B*D); %alfa 3
    a4(i)=(C*E-B*F)/(A*E-B*D); %alfa 4
    
%%Cálculo e armazenamento dos valores de vel lineares a cada passo 
%%de tempo nos vetores Va e Vb
    Va_i=a*w2*(-sin(theta2)); %Direção do versor i
    Va_j=a*w2*(cos(theta2));  %Direção do versor j
    Va(i)=((Va_i^2)+(Va_j^2))^(1/2); %Módulo da velocidade
    Vb_i=c*w4(i)*(-sin(theta4)); %Direção do versor i
    Vb_j=c*w4(i)*(cos(theta4));  %Direção do versor j
    Vb(i)=((Vb_i^2)+(Vb_j^2))^(1/2);  %Módulo da velocidade
    
%%Cálculo e armazenamento dos valores de aceleração lineares a cada passo 
%%de tempo nos vetores Aa e Ab
    Aa_i=-(a*0*sin(theta2)+a*(w2^2)*cos(theta2)); %%OBS: alfa 2 é nula %%Direção do versor i
    Aa_j=a*0*cos(theta2)-a*(w2^2)*sin(theta2);    %Direção do versor j
    Aa(i)=((Aa_i^2)+(Aa_j^2))^(1/2);              %Módulo da aceleração
    Ab_i=-(c*a4(i)*sin(theta4)+c*(w4(i)^2)*cos(theta4));  %Direção do versor i
    Ab_j=c*a4(i)*cos(theta4)-c*(w4(i)^2)*sin(theta4);     %Direção do versor j
    Ab(i)=((Ab_i^2)+(Ab_j^2))^(1/2);                      %Módulo da aceleração

%%Cálculo e armazenamento dos valores de aceleração lineares nos CGs dos elos a cada passo 
%%de tempo nos vetores Aag, Abg e Acg
    Aag_i=-((a/2)*0*sin(theta2)+(a/2)*(w2^2)*cos(theta2)); %%OBS: alfa 2 é nula %%Direção do versor i
    Aag_j=(a/2)*0*cos(theta2)-(a/2)*(w2^2)*sin(theta2);    %Direção do versor j
    Aag(i)=((Aag_i^2)+(Aag_j^2))^(1/2);              %Módulo da aceleração
    Abg_i=Aa_i-((b/2)*a3(i)*sin(theta3)+(b/2)*(w3(i)^2)*cos(theta3));  %Direção do versor i
    Abg_j=Aa_j+(b/2)*a3(i)*cos(theta3)-(b/2)*(w3(i)^2)*sin(theta3);     %Direção do versor j
    Abg(i)=(((Abg_i^2)+(Abg_j^2))^(1/2));                      %Módulo da aceleração
    Acg_i=-((c/2)*a4(i)*sin(theta4)+(c/2)*(w4(i)^2)*cos(theta4));  %Direção do versor i
    Acg_j=(c/2)*a4(i)*cos(theta4)-(c/2)*(w4(i)^2)*sin(theta4);     %Direção do versor j
    Acg(i)=((Acg_i^2)+(Acg_j^2))^(1/2);                      %Módulo da aceleração

%%Cálculo dos valores das componentes dos vetores posição R de cada elo 
%%em função dos ângulos para cada instante de tempo do movimento de rotação
    R_12x=-(a/2)*cos(theta2);
    R_12y=-(a/2)*sin(theta2);
    R_32x=-R_12x;
    R_32y=-R_12y;
    R_23x=-(b/2)*cos(theta3);
    R_23y=-(b/2)*sin(theta3);
    R_43x=-R_23x;
    R_43y=-R_23y;
    R_14x=-(c/2)*cos(theta4);
    R_14y=-(c/2)*sin(theta4);
    R_34x=-R_14x;
    R_34y=-R_14y;
%%Matriz A
    A=[1 0 1 0 0 0 0 0 0;
        0 1 0 1 0 0 0 0 0;
        -R_12y R_12x -R_32y R_32x 0 0 0 0 1;
        0 0 -1 0 1 0 0 0 0;
        0 0 0 -1 0 1 0 0 0;
        0 0 R_23y -R_23x -R_43y R_43x 0 0 0;
        0 0 0 0 -1 0 1 0 0;
        0 0 0 0 0 -1 0 1 0;
        0 0 0 0 R_34y -R_34x -R_14y R_14x 0];
%%Matriz b composta pelos valores conhecidos do sistema de equações
    matriz_b=[m2*Aag_i; m2*Aag_j; Ig2*0; m3*Abg_i; m3*Abg_j; Ig3*a3(i); m4*Acg_i; m4*Acg_j; Ig4*a4(i)];

%%Resolução do sistema de equações para cada instante de tempo da simulação
    x=A\matriz_b;
    F12x(i)=x(1,1);
    F12y(i)=x(2,1);
    F32x(i)=x(3,1);
    F32y(i)=x(4,1);
    F43x(i)=x(5,1);
    F43y(i)=x(6,1);
    F14x(i)=x(7,1);
    F14y(i)=x(8,1);
    T12(i)=x(9,1);
%%Módulo das forças aplicadas pelo mecanismo no suporte
%%OBS: Essas forças são iguais em módulo às forças aplicadas pelo suporte 
%%no mecanismo
    F12(i)=((F12x(i)^2)+(F12y(i)^2))^(1/2);
    F14(i)=((F14x(i)^2)+(F14y(i)^2))^(1/2);
end

%Plotagem do gráfico do torque T12 no elo 2
figure(3)
plot(t,T12,'b-');
xlabel('Tempo (s)'); ylabel('Torque (N.m)');
title('Torque no elo 2 vs. Tempo'); grid on;
axis([0 1 -0.1 0.06]);

%Plotagem do gráfico da força F21 aplicada pelo elo 2 aos suportes
F21=F12; %Forças iguais em módulo e opostas vetorialmente
figure(4)
plot(t,F21,'k-');
xlabel('Tempo (s)'); ylabel('Força F21 (N)');
title('Força F21 vs. Tempo'); grid on;
axis([0 1 0 7]);

%Plotagem do gráfico da força F21 aplicada pelo elo 2 aos suportes
F41=F14; %Forças iguais em módulo e opostas vetorialmente
figure(5)
plot(t,F41,'r-');
xlabel('Tempo (s)'); ylabel('Força F41 (N)');
title('Força F41 vs. Tempo'); grid on;
axis([0 1 0 4.5]);

%Plotagem dos gráficos das componentes das forças aplicadas pelo mecanismo
%nos suportes
figure(6)
F21x=-F12x;
F21y=-F12y;
plot(t,F21x,'m-');
xlabel('Tempo (s)'); ylabel('Força F21 (N)');
title('Componentes da força F21 vs. Tempo'); 
grid on; 
hold on;
plot(t,F21y,'g-');
legend('F21x','F21y','location','northwest');
axis([0 1 -2 7]);
hold off;

figure(7)
F41x=-F14x;
F41y=-F14y;
plot(t,F41x,'m-');
xlabel('Tempo (s)'); ylabel('Força F41 (N)');
title('Componentes da força F41 vs. Tempo'); 
grid on; 
hold on;
plot(t,F41y,'g-');
legend('F41x','F41y','location','northwest');
axis([0 1 -4 2]);
hold off;

%Valores máximos e mínimos das forças e torques
T12_max=max(T12)
T12_min=min(T12)
F21_max=max(F21)
F21_min=min(F21)
F41_max=max(F41)
F41_min=min(F41)
F21x_max=max(-F12x)
F21x_min=min(-F12x)
F21y_max=max(-F12y)
F21y_min=min(-F12y)
F41x_max=max(-F14x)
F41x_min=min(-F14x)
F41y_max=max(-F14y)
F41y_min=min(-F14y)

% % %Plotagem dos gráficos das velocidades angulares w3 e w4
% % figure(3)
% % plot(t,w3,'k-');
% % xlabel('Tempo (s)'); ylabel('Velocidade angular (Rad/s)');
% % title('Velocidade angular vs. Tempo'); grid on;
% % hold on;
% % plot(t,w4,'b-');
% % legend('Elo 3','Elo 4','location','northwest');
% % axis([0 1 -10 8]);
% % hold off;
% % 
% % %Plotagem dos gráficos das acelerações angulares a3 e a4
% % figure(4)
% % plot(t,a3,'r-');
% % xlabel('Tempo (s)'); ylabel('Aceleração angular (Rad/s^2)');
% % title('Aceleração angular vs. Tempo'); grid on;
% % hold on;
% % plot(t,a4,'b-');
% % legend('Elo 3','Elo 4','location','northwest');
% % axis([0 1 -200 200]);
% % hold off;
% % 
% % %Plotagem dos gráficos das velocidades lineares Va e Vb
% % figure(5)
% % plot(t,Va,'b-');
% % xlabel('Tempo (s)'); ylabel('Velocidade linear (m/s)');
% % title('Velocidade linear vs. Tempo'); grid on;
% % hold on;
% % plot(t,Vb,'r');
% % legend('Ponto A','Ponto B','location','northwest');
% % axis([0 1 0 0.7]);
% % hold off;
% % 
% % %Plotagem dos gráficos das acelerações lineares Aa e Ab
% % figure(6)
% % plot(t,Aa,'r');
% % xlabel('Tempo (s)'); ylabel('Aceleração linear (m/s^2)');
% % title('Aceleração linear vs. Tempo'); grid on;
% % hold on;
% % plot(t,Ab,'b');
% % legend('Ponto A','Ponto B','location','northwest');
% % axis([0 1 0 14]);
% % hold off;
% % 
% % %%Velocidades e acelerações angulares e lineares maximas e minimas
% % w3_e_w4_maximos=[max(w3),max(w4)]
% % Velocidade_angular_maxima=max(w3_e_w4_maximos)
% % w3_e_w4_minimos=[min(w3),min(w4)]
% % Velocidade_angular_minima=min(w3_e_w4_minimos)
% % alfa3_e_alfa4_maximos=[max(a3),max(a4)]
% % Aceleracao_angular_maxima=max(alfa3_e_alfa4_maximos)
% % alfa3_e_alfa4_minimos=[min(a3),min(a4)]
% % Aceleracao_angular_minima=min(alfa3_e_alfa4_minimos)
% % Va_e_Vb_maximos=[max(Va),max(Vb)]
% % Velocidade_linear_maxima=max(Va_e_Vb_maximos)
% % Va_e_Vb_minimos=[min(Va),min(Vb)]
% % Velocidade_linear_minima=min(Va_e_Vb_minimos)
% % Aa_e_Ab_maximos=[max(Aa),max(Ab)]
% % Aceleracao_linear_maxima=max(Aa_e_Ab_maximos)
% % Aa_e_Ab_minimos=[min(Aa),min(Ab)]
% % Aceleracao_linear_minima=min(Aa_e_Ab_minimos)

%% Evaluate positions of points A,B,C,D,P
rA=zeros(length(t2v),2);
rB=a*[cos(t2v) sin(t2v)];
rC=rB+b*[cos(t3v) sin(t3v)];
rD=[rA(:,1)+d rA(:,2)];
rP=rB+BP*[cos(t3v+tBP) sin(t3v+tBP)];


%% Show simulation
figure(1), clf, set(1,'position',[1169 335 690 650])
if outvid==1 
    filename='qbarras_nr.gif';
    frame=getframe(gcf); im=frame2im(frame); [A,map]=rgb2ind(im,256);
    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1e-3);
end
rT=[rA; rB; rC; rD; rP]; mx=max(max(abs(rT)));
rT=[rA; rB; rC; rD]; mx=max(max(abs(rT))); axis([-mx mx -mx mx]), 
axis tight, axis equal, axis off,
hBP=line([rB(1,1) rP(1,1)],[rB(1,2) rP(1,2)]); set(hBP,'Color','g','LineStyle','-','Marker','o'),
hPC=line([rP(1,1) rC(1,1)],[rP(1,2) rC(1,2)]); set(hPC,'Color','g','LineStyle','-','Marker','o'),
hAB=line([rA(1,1) rB(1,1)],[rA(1,2) rB(1,2)]); set(hAB,'Color','k','LineStyle','-','Marker','o'),
hBC=line([rB(1,1) rC(1,1)],[rB(1,2) rC(1,2)]); set(hBC,'Color','k','LineStyle','-','Marker','o'),
hCD=line([rC(1,1) rD(1,1)],[rC(1,2) rD(1,2)]); set(hCD,'Color','k','LineStyle','-','Marker','o'),
hDA=line([rD(1,1) rA(1,1)],[rD(1,2) rA(1,2)]); set(hDA,'Color','k','LineStyle','-','Marker','o'),

hABo=line([rA(1,1) rB(1,1)],[rA(1,2) rB(1,2)]); set(hABo,'Color',.7*[1 1 1],'LineStyle',':'),
hBCo=line([rB(1,1) rC(1,1)],[rB(1,2) rC(1,2)]); set(hBCo,'Color',.7*[1 1 1],'LineStyle',':'),
hCDo=line([rC(1,1) rD(1,1)],[rC(1,2) rD(1,2)]); set(hCDo,'Color',.7*[1 1 1],'LineStyle',':'),
hDAo=line([rD(1,1) rA(1,1)],[rD(1,2) rA(1,2)]); set(hDAo,'Color',.7*[1 1 1],'LineStyle',':'),
text(rA(1,1)-mx/100,rA(1,2)-mx/20,'A')
text(rB(1,1)-mx/100,rB(1,2)-mx/20,'B')
text(rC(1,1)-mx/100,rC(1,2)-mx/20,'C')
text(rD(1,1)-mx/100,rD(1,2)-mx/20,'D')

if outvid==1 dts=5; else dts=1; end
for n=2:dts:length(t)
    set(hAB,'xdata',[rA(n,1) rB(n,1)],'ydata',[rA(n,2) rB(n,2)]);
    set(hBC,'xdata',[rB(n,1) rC(n,1)],'ydata',[rB(n,2) rC(n,2)]);
    set(hCD,'xdata',[rC(n,1) rD(n,1)],'ydata',[rC(n,2) rD(n,2)]);
    set(hBP,'xdata',[rB(n,1) rP(n,1)],'ydata',[rB(n,2) rP(n,2)]);
    set(hPC,'xdata',[rP(n,1) rC(n,1)],'ydata',[rP(n,2) rC(n,2)]);
    hB=line([rB(n-1,1) rB(n,1)],[rB(n-1,2) rB(n,2)]); set(hB,'Color','r','LineStyle',':');
    hC=line([rC(n-1,1) rC(n,1)],[rC(n-1,2) rC(n,2)]); set(hC,'Color','r','LineStyle',':');
    hP=line([rP(n-1,1) rP(n,1)],[rP(n-1,2) rP(n,2)]); set(hP,'Color','m','LineStyle',':');
    if outvid==1 
        frame=getframe(gcf); im=frame2im(frame); [A,map]=rgb2ind(im,256);
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1e-3);
    else pause(1e-12);
    end
end,