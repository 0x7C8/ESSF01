%Analog elektronik - Matlab hjälp 
%Matlab exemplel för att kolla fasmarginal och slutna förstärkningen, samt
%beräkna slingpoler och titta på stegsvar mm.
%2-stegs förstärkare (ASGE-GE), före och efter kompensering
%Med Egna v�rden

clear all; close all;

%%Definiera Data
VT=25*1e-3;
Bf1=250; %Kolla datablad!
Bf2=Bf1;
C1=100*1e-9; %Ersätter Cpi1_prim
C2=2.2*1e-6; %Ersätter Cpi2
C2ny=40*1e-6; %ers�tter C2 vid capacitive narrowbandning
%%Förstärkare 
Rs=50; %Källan är inte ideal
R1=100;
%R2=5*1e3;
RL=1000; 
AtINF=-1/R1;%Asymptotiska förstärkningen

Ic1ab=(9*1e-3)/2; %Strömmen i ingångssteget
rpi_p=2*Bf1*VT/Ic1ab; %rpi_p=2*rpi
Ic2=5*1e-3; %Krav: max 100mV peak in --> max 1.1V peak ut, RL*Ic2>1,1V Obs!! H�r ska vi ha Krav: max 100mV peak in --> max 0.001 mA ut.  
rpi2=Bf2*VT/Ic2;
gm2=Ic2/VT;

Rbias=150; %Dålig biasering med Rbias! (använd strömspegel) Rbias=0.7/Ic1ab ,Biaserar upp GE steget.
rpi2_new=rpi2*Rbias/(rpi2+Rbias)


%DC slingförstärkning och slingpoler:
ABnoll= -(Bf1*Bf2*R1)/(Rs+R1+rpi_p);
P1=-(R1+Rs+rpi_p)/((R1+Rs)*C1*rpi_p)
P2=-1/(rpi2*C2)

%%Är alla poler dominanta?:
w0_2p=(abs( (1-ABnoll)*P1*P2 ))^(1/2)
SummaP=P1+P2 %Summa slingpoler
SummaP_p=-sqrt(2)*w0_2p %Summa av systempoler (2st)
%Kolla att summaP > SummaP_p --> bara dominanta poler som kan flyttas till
%önskad position
n=-(w0_2p^2)/(sqrt(2)*w0_2p+P1+P2)


%%%%%%FREKVENSKOMPENSERING

s=zpk('s') %Definiera s

% %Före kompensering: (Betraktas som ett system med två poler)
 ABs=ABnoll/((1-s/P1)*(1-s/P2)) 
 At=AtINF*(-1)*ABs/(1-ABs); %Slutna förstärkningen, icke kompenserad.
% 
%Kompensering med capacitive narrow banding
P2ny=-1/(rpi2*C2ny);
ABsny=ABnoll/((1-s/P1)*(1-s/P2ny));
Atny=AtINF*(-1)*ABsny/(1-ABsny); %Slutna förstärkningen, kompenserad Capacitive Narrow Banding.
%
%%Undersök fasmarginal före och efter kompensering:


[gain_margin_before, phase_margin_before] = margin((-1)*ABs) %Ignorera gain margin (mer om den i reglerteknik). 
%Matlab verkar inte gilla negativa system, s� -1 beh�vs f�r att f� r�tt fasmarginal.

[gain_margin_after, phase_margin_after] = margin((-1)*ABsny)
%
%
% %%Implementera fantomnollan, undersök alla fall:
 % Lph i serie med R1
 delta_Lph=1.5; %Effektivt om delta > 7
 Lph=-R1/n;
 AtINF_Lph= 1/(R1 + s*Lph); 
% 
% %%Efter kompensering för MFM med fantomnolla:
% %%Lph serie med R1:
 ABs_n_Lph=ABnoll*(1-s/n)/((1-s/P1)*(1-s/P2)*(1-s/(delta_Lph*n)));
 Atn_Lph=AtINF_Lph*(-1)*ABs_n_Lph/(1-ABs_n_Lph);

% Cph parallel med Rs
 delta_Cph=1.5;
 Cph=-1/(Rs*n);
 AtINF_Cph= AtINF;
% %%Efter kompensering för MFM med fantomnolla:
% %%Cph || Rs:
 ABs_n_Cph=ABnoll*(1-s/n)/((1-s/P1)*(1-s/P2)*(1-s/(delta_Lph*n)))
 Atn_Cph=AtINF_Cph*(-1)*ABs_n_Cph/(1-ABs_n_Cph) 


%%%%%%FIGURER
%Fasmarginal kollas "open loop", dvs frekvensen w0, där |AB(w0)|=1=0dB, före=ABs och efter=ABs_n kompensering
%(Bode-funktionen behöver ibland ett (-1).* pga 'Phase unwrap')
figure(1);bode((-1).*ABs,'b',(-1).*ABs_n_Lph,'k--', (-1).*ABs_n_Cph, 'r--', (-1).*ABsny, 'y--'); 
title('Slingf�rst�rkning:'); legend('AB(s)','AB_n Lph(s)','AB_nCph(s)','ABny(s)','Location','Best')

figure(2);bode(At,'b',Atn_Lph,'k--',Atn_Cph,'r--', Atny,'y--'); hold on; %
title('Den slutna förstärkningen, At'); legend('A_t','A_{tn,Lph}','A_{tn,Cph}','A_t_ny','Location','Best')
    
figure(3); step((-1)*At);hold on;step(Atn_Lph);step((-1)*Atn_Cph);step((-1)*Atny); 
title('Stegsvaren före och efter kompensering'); legend('A_t','A_{tn,Lph}','A_{tn,Cph}','A_t_ny','Location','Best')



%%%%%% HJÄLP FÖR ATT PLOTTA MÄTRESULTAT TILLSAMMANS (SAMMA FIGUR) MED SIMULERAD PRESTANDA (Kompenserat
%%%%%% och okompenserat)
% W_labbet=[frekvensvektor från labbet].*(2*pi);
% At_labbet_kompenserat_dB=[mätresultat]
% At_labbet_kompenserat_fas=[mätresultat]
%på samma sätt lägger för At_okompenserat
W=[1:100:1e6].*(2*pi);
[MAG_At, PHASE_At] = bode(At,W);
for k=1:length(W)
    dB_MAG_At(k)=20*log10(MAG_At(1,1,k));
    phase_At(k)=PHASE_At(1,1,k);
end
% semilogx(W,dB_MAG_At,'b', W_labbet, At_labbet_kompenserat_dB,'r', W_labbet,At_labbet_okompenserat_dB,'k');
% semilogx(W,phase_At,'b', W_labbet, At_labbet_kompenserat_fas,'r', W_labbet,At_labbet_okompenserat_fas,'k');
% %xlabel och ylabel för axlarna
figure(4);
semilogx(W,dB_MAG_At,'b');%hold on; ... lägg till mätresultat
figure(5);
semilogx(W,phase_At,'b');%hold on; ...lägg till mätresultat