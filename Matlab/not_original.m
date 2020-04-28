%Analog elektronik - Matlab hjälp 
%Matlab exemplel för att kolla fasmarginal och slutna förstärkningen, samt
%beräkna slingpoler och titta på stegsvar mm.
%2-stegs förstärkare (ASGE-GE), före och efter kompensering

clear all; close all;

%% Definiera Data
VT = 25e-3;   % (V)
Bf1 = 250; % Kolla datablad!
Bf2 = Bf1;
C1 = 100e-9; % (F) - Ersätter Cpi1_prim
C2 = 2.2e-6; % (F) - Ersätter Cpi2

%% Förstärkare 
Rs = 50; % (Ohm) %Källan är inte ideal
R1 = 100; % (Ohm)
RL = 1000;    % (Ohm)
AtINF = -1/R1;  %Asymptotiska förstärkningen

Ic1ab = (10e-3)/2; %Strömmen i ingångssteget
rpi_p = 2*Bf1*VT/Ic1ab; %rpi_p = 2*rpi
Ic2 = 5e-3; % (A)
rpi2 = Bf2*VT/Ic2;
gm2 = Ic2/VT;

%Rbias = 150; %Dålig biasering med Rbias! (använd strömspegel) Rbias = 0.7/Ic1ab ,Biaserar upp GE steget.
%rpi2_new = rpi2*Rbias/(rpi2+Rbias)

%DC slingförstärkning och slingpoler:
ABnoll = (-100);
P1 = -(R1+Rs+rpi_p)/((R1+Rs)*C1*rpi_p)
P2 = -1/(rpi2*C2)

%% Är alla poler dominanta?:
w0_2p = sqrt((abs( (1-ABnoll)*P1*P2 )))
SummaP = P1+P2 %Summa slingpoler
SummaP_p = -sqrt(2)*w0_2p %Summa av systempoler (2st)
%Kolla att summaP > SummaP_p --> bara dominanta poler som kan flyttas till
%önskad position
n = -(w0_2p^2)/(sqrt(2)*w0_2p+P1+P2)


%%%%%%FREKVENSKOMPENSERING
%% Undersök fasmarginal före och efter kompensering:
s = zpk('s') %Definiera s

%Före kompensering: (Betraktas som ett system med två poler)
ABs = ABnoll/((1-s/P1)*(1-s/P2))
[gain_margin, phase_margin] = margin((-1)*ABs) %Ignorera gain margin (mer om den i reglerteknik).
%Matlab verkar inte gilla negativa system, s� -1 beh�vs f�r att f� r�tt fasmarginal.

At = AtINF*(-1)*ABs/(1-ABs); %Slutna förstärkningen, icke kompenserad.

%% %%%%FIGURER
%Fasmarginal kollas "open loop", dvs frekvensen w0, där |AB(w0)| = 1 = 0dB, före = ABs och efter = ABs_n kompensering
%(Bode-funktionen behöver ibland ett (-1).* pga 'Phase unwrap')
figure(1)
bode((-1).*ABs,'b'); 
title('Slingf�rst�rkning');
grid on

figure(2)
bode(At,'b'); 
title('Den slutna förstärkningen, At');
grid on
    
figure(3)
step(At);
title('Stegsvaren före och efter kompensering')



%% %%%% HJÄLP FÖR ATT PLOTTA MÄTRESULTAT TILLSAMMANS (SAMMA FIGUR) MED SIMULERAD PRESTANDA (Kompenserat
%%%%%% och okompenserat)
% W_labbet = [frekvensvektor från labbet].*(2*pi);
% At_labbet_kompenserat_dB = [mätresultat]
% At_labbet_kompenserat_fas = [mätresultat]
%på samma sätt lägger för At_okompenserat
W = [1:100:1e6].*(2*pi);
[MAG_At, PHASE_At]  =  bode(At,W);
for k = 1:length(W)
    dB_MAG_At(k) = 20*log10(MAG_At(1,1,k));
    phase_At(k) = PHASE_At(1,1,k);
end
% semilogx(W,dB_MAG_At,'b', W_labbet, At_labbet_kompenserat_dB,'r', W_labbet,At_labbet_okompenserat_dB,'k');
% semilogx(W,phase_At,'b', W_labbet, At_labbet_kompenserat_fas,'r', W_labbet,At_labbet_okompenserat_fas,'k');
% %xlabel och ylabel för axlarna
figure(4)
semilogx(W,dB_MAG_At,'b');%hold on; ... lägg till mätresultat
figure(5)
semilogx(W,phase_At,'b');%hold on; ...lägg till mätresultat
