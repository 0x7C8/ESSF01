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
Rs = 100; % (Ohm) %Källan är inte ideal
R1 = 500; % (Ohm)
RL = 1000;    % (Ohm)
AtINF = 1+R2/R1;  %Asymptotiska förstärkningen

Ic1ab = (9.4e-3)/2; %Strömmen i ingångssteget
rpi_p = 2*Bf1*VT/Ic1ab; %rpi_p = 2*rpi
Ic2 = 2e-3; %Krav: max 100mV peak in --> max 1.1V peak ut, RL*Ic2>1,1V
rpi2 = Bf2*VT/Ic2;
gm2 = Ic2/VT;

Rbias = 150; %Dålig biasering med Rbias! (använd strömspegel) Rbias = 0.7/Ic1ab ,Biaserar upp GE steget.
rpi2_new = rpi2*Rbias/(rpi2+Rbias)


%DC slingförstärkning och slingpoler:
ABnoll = -100 
P1 = -1.45e4;
P2 = -3.2e3

%% Är alla poler dominanta?:
w0_2p = sqrt((abs( (1-ABnoll)*P1*P2 ))
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
At = AtINF*(-1)*ABs/(1-ABs); %Slutna förstärkningen, icke kompenserad.

[gain_margin, phase_margin] = margin((-1)*ABs) %Ignorera gain margin (mer om den i reglerteknik).
%Matlab verkar inte gilla negativa system, s� -1 beh�vs f�r att f� r�tt fasmarginal.

%%Implementera fantomnollan, undersök alla fall:
%Här tittar vi bara på Cph || R2
delta_Cph = 10.17; %Effektivt om delta > 7
Cph = -1/(R2*n)
AtINF_Cph = AtINF*( 1 - s/(AtINF*n) ) / (1 - s/n) 

%%Efter kompensering för MFM med fantomnolla:
%%Cph || med R2:
ABs_n_Cph = ABnoll*(1-s/n)/((1-s/P1)*(1-s/P2)*(1-s/(delta_Cph*n)))
Atn_Cph = AtINF_Cph*(-1)*ABs_n_Cph/(1-ABs_n_Cph)


%% %%%%FIGURER
%Fasmarginal kollas "open loop", dvs frekvensen w0, där |AB(w0)| = 1 = 0dB, före = ABs och efter = ABs_n kompensering
%(Bode-funktionen behöver ibland ett (-1).* pga 'Phase unwrap')
figure(1)
bode((-1).*ABs,'b',(-1).*ABs_n_Cph,'k--'); 
title('Slingf�rst�rkning: före och efter fantomnolla'); legend('AB(s)','AB_n Cph(s)','Location','Best')

figure(2)
bode(At,'b',Atn_Cph,'k--'); hold on; 
title('Den slutna förstärkningen, At'); legend('A_t','A_{tn,Cph}','Location','Best')
    
figure(3)
step(At);hold on;step(Atn_Cph); 
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
