%Analog elektronik - Matlab hj�lp 
%Matlab exemplel f�r att kolla fasmarginal och slutna f�rst�rkningen, samt
%ber�kna slingpoler och titta p� stegsvar mm.
%2-stegs f�rst�rkare (ASGE-GE), f�re och efter kompensering
%Med Egna v�rden

clear all; close all;

%% Definiera Data
VT = 25e-3;
Bf1 = 250; %Kolla datablad!
Bf2 = Bf1;
C1 = 100e-9; %Ers�tter Cpi1_prim
C2 = 2.2e-6; %Ers�tter Cpi2
C2ny = 40e-6; %ers�tter C2 vid capacitive narrowbandning

%% F�rst�rkare 
Rs = 50; %K�llan �r inte ideal
R1 = 100;
%R2 = 5e3;
RL = 1000; 
AtINF = -1/R1;%Asymptotiska f�rst�rkningen

Ic1ab = (9e-3)/2; %Str�mmen i ing�ngssteget
rpi_p = 2*Bf1*VT/Ic1ab; %rpi_p = 2*rpi
Ic2 = 5e-3; %Krav: max 100mV peak in --> max 1.1V peak ut, RL*Ic2>1,1V Obs!! H�r ska vi ha Krav: max 100mV peak in --> max 0.001 mA ut.  
rpi2 = Bf2*VT/Ic2;
gm2 = Ic2/VT;

Rbias = 150; %D�lig biasering med Rbias! (anv�nd str�mspegel) Rbias = 0.7/Ic1ab ,Biaserar upp GE steget.
rpi2_new = rpi2*Rbias/(rpi2+Rbias)

%DC slingf�rst�rkning och slingpoler:
ABnoll = -(Bf1*Bf2*R1)/(Rs+R1+rpi_p);
P1 = -(R1+Rs+rpi_p)/((R1+Rs)*C1*rpi_p)
P2 = -1/(rpi2*C2)

%% �r alla poler dominanta?:
w0_2p = (abs( (1-ABnoll)*P1*P2 ))^(1/2)
SummaP = P1+P2 %Summa slingpoler
SummaP_p = -sqrt(2)*w0_2p %Summa av systempoler (2st)
%Kolla att summaP > SummaP_p --> bara dominanta poler som kan flyttas till
%�nskad position
n = -(w0_2p^2)/(sqrt(2)*w0_2p+P1+P2)


%% FREKVENSKOMPENSERING

s = zpk('s'); %Definiera s

%F�re kompensering: (Betraktas som ett system med tv� poler)
ABs = ABnoll/((1-s/P1)*(1-s/P2)) 
At = AtINF*(-1)*ABs/(1-ABs); %Slutna f�rst�rkningen, icke kompenserad.
% 
%Kompensering med capacitive narrow banding
P2ny = -1/(rpi2*C2ny);
ABsny = ABnoll/((1-s/P1)*(1-s/P2ny));
Atny = AtINF*(-1)*ABsny/(1-ABsny); %Slutna f�rst�rkningen, kompenserad Capacitive Narrow Banding.
%
%%Unders�k fasmarginal f�re och efter kompensering:

[gain_margin_before, phase_margin_before] = margin((-1)*ABs) %Ignorera gain margin (mer om den i reglerteknik). 
[gain_margin_after, phase_margin_after] = margin((-1)*ABsny)
%Matlab verkar inte gilla negativa system, s� -1 beh�vs f�r att f� r�tt fasmarginal.


%% Implementera fantomnollan, unders�k alla fall:
% Lph i serie med R1
delta_Lph = 1.5; %Effektivt om delta > 7
Lph = -R1/n;
AtINF_Lph = 1/(R1 + s*Lph); 

%Efter kompensering f�r MFM med fantomnolla:
%Lph serie med R1:
ABs_n_Lph = ABnoll*(1-s/n)/...
    ((1-s/P1)*(1-s/P2)*(1-s/(delta_Lph*n)));
Atn_Lph = AtINF_Lph*(-1)*ABs_n_Lph/(1-ABs_n_Lph);

%Cph parallel med Rs
delta_Cph = 1.5;
Cph = -1/(Rs*n);
AtINF_Cph = AtINF;
%Efter kompensering f�r MFM med fantomnolla:
%Cph || Rs:
ABs_n_Cph = ABnoll*(1-s/n)/...
    ((1-s/P1)*(1-s/P2)*(1-s/(delta_Lph*n)))
Atn_Cph = AtINF_Cph*(-1)*ABs_n_Cph/(1-ABs_n_Cph) 


%% FIGURER
% Fasmarginal kollas "open loop", dvs frekvensen w0,
% d�r |AB(w0)| = 1 = 0dB, f�re = ABs och efter = ABs_n kompensering
%(Bode-funktionen beh�ver ibland ett (-1).* pga 'Phase unwrap')
figure(1)
bode((-1).*ABs,'b',(-1).*ABs_n_Lph,'k--',...
    (-1).*ABs_n_Cph, 'r--', (-1).*ABsny, 'y--')
grid on
title('Slingf�rst�rkning:');
legend('$A\beta(s)$','$A\beta\_{n,Lph}(s)$',...
    '$A\beta\_{n,Cph}(s)$','$A\beta\_{ny}(s)$',...
    'Interpreter','latex', 'Location','Best')

figure(2)
bode(At,'b',Atn_Lph,'k--',Atn_Cph,'r--', Atny,'y--')
grid on
title('Den slutna f�rst�rkningen, At');
legend('$A\_t$','$A\_{t,n,Lph}$',...
    '$A\_{t,n,Cph}$','$A\_{t_ny}$',...
    'Interpreter','latex','Location','Best')

figure(3)
stepplot((-1)*At, Atn_Lph, (-1)*Atn_Cph, (-1)*Atny)
title('Stegsvaren f�re och efter kompensering');
legend('$A\_t$','$A\_{t,n,Lph}$',...
    '$A\_{t,n,Cph}$','$A\_{t_ny}$',...
    'Interpreter','latex','Location','Best')

%% HJ�LP FÖR ATT PLOTTA M�TRESULTAT TILLSAMMANS (SAMMA FIGUR) MED SIMULERAD PRESTANDA (Kompenserat
% och okompenserat)
% W_labbet = [frekvensvektor fr�n labbet].*(2*pi);
% At_labbet_kompenserat_dB = [m�tresultat]
% At_labbet_kompenserat_fas = [m�tresultat]
%p� samma s�tt l�gger f�r At_okompenserat
W = [1:100:1e6].*(2*pi);
[MAG_At, PHASE_At] = bode(At,W);
for k = 1:length(W)
    dB_MAG_At(k) = 20*log10(MAG_At(1,1,k));
    phase_At(k) = PHASE_At(1,1,k);
end
% semilogx(W,dB_MAG_At,'b', W_labbet, At_labbet_kompenserat_dB,'r', W_labbet,At_labbet_okompenserat_dB,'k');
% semilogx(W,phase_At,'b', W_labbet, At_labbet_kompenserat_fas,'r', W_labbet,At_labbet_okompenserat_fas,'k');
% %xlabel och ylabel f�r axlarna
figure(4);
semilogx(W,dB_MAG_At,'b');%hold on; ... l�gg till m�tresultat
figure(5);
semilogx(W,phase_At,'b');%hold on; ...l�gg till m�tresultat

%% Simulated data (not finished!)
load R9_data.mat
R9_freq = R9amplitude(:,1);
R9_amp = R9amplitude(:,2);
R9_deg = R9amplitude(:,3);
subplot(1,2,1)
semilogx(R9_freq, R9_deg)
axis([1 1e6 -180 0])
grid on
xlabel('frequency (Hz)')
ylabel('Polarity (deg)')

subplot(1,2,2)
semilogx(R9_freq, R9_amp)
axis([1 1e6 -80 -35])
grid on
xlabel('frequency (Hz)')
ylabel('Amplitude dB ')


%% Save graphs in following format and settings
for k = 1:5
    figname = ['figure', num2str(k)];
    figure(k)
    title('')
    set(gca,...
    'XGrid','on',...
    'YGrid', 'on',...
    'Fontsize', 10,...
    'linewidth', 1,...
    'FontName', 'Arial')
    saveas(gcf,figname,'epsc')
end    
    