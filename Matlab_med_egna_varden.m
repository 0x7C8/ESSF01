%Analog elektronik - Matlab hjälp 
%Matlab exemplel för att kolla fasmarginal och slutna förstärkningen, samt
%beräkna slingpoler och titta på stegsvar mm.
%2-stegs förstärkare (ASGE-GE), före och efter kompensering
%Med Egna värden

clear; close all

%% Definiera Data
VT = 25e-3;
Bf1 = 250; %Kolla datablad!
Bf2 = Bf1;
C1 = 100e-9; %Ersätter Cpi1_prim
C2 = 2.2e-6; %Ersätter Cpi2
C2ny = 40e-6; %ersätter C2 vid capacitive narrowbandning

%% Förstärkare 
Rs = 50; %Källan är inte ideal
R1 = 100;
%R2 = 5e3;
RL = 1000; 
AtINF = -1/R1;%Asymptotiska förstärkningen

Ic1ab = (9e-3)/2; %Strömmen i ingångssteget
rpi_p = 2*Bf1*VT/Ic1ab; %rpi_p = 2*rpi
Ic2 = 5e-3; %Krav: max 100mV peak in --> max 1.1V peak ut, RL*Ic2>1,1V Obs!! Här ska vi ha Krav: max 100mV peak in --> max 0.001 mA ut.  
rpi2 = Bf2*VT/Ic2;
gm2 = Ic2/VT;

Rbias = 150; %Dålig biasering med Rbias! (använd strömspegel) Rbias = 0.7/Ic1ab ,Biaserar upp GE steget.
rpi2_new = rpi2*Rbias/(rpi2+Rbias)

%DC slingförstärkning och slingpoler:
ABnoll = -(Bf1*Bf2*R1)/(Rs+R1+rpi_p);
P1 = -(R1+Rs+rpi_p)/((R1+Rs)*C1*rpi_p)
P2 = -1/(rpi2*C2)

%% Är alla poler dominanta?:
w0_2p = (abs( (1-ABnoll)*P1*P2 ))^(1/2)
SummaP = P1+P2 %Summa slingpoler
SummaP_p = -sqrt(2)*w0_2p %Summa av systempoler (2st)
%Kolla att summaP > SummaP_p --> bara dominanta poler som kan flyttas till
%önskad position
n = -(w0_2p^2)/(sqrt(2)*w0_2p+P1+P2)


%% FREKVENSKOMPENSERING

s = zpk('s'); %Definiera s

%Före kompensering: (Betraktas som ett system med två poler)
ABs = ABnoll/((1-s/P1)*(1-s/P2)) 
At = AtINF*(-1)*ABs/(1-ABs); %Slutna förstärkningen, icke kompenserad.
% 
%Kompensering med capacitive narrow banding
P2ny = -1/(rpi2*C2ny);
ABsny = ABnoll/((1-s/P1)*(1-s/P2ny));
Atny = AtINF*(-1)*ABsny/(1-ABsny); %Slutna förstärkningen, kompenserad Capacitive Narrow Banding.
%
%%Undersök fasmarginal före och efter kompensering:

[gain_margin_before, phase_margin_before] = margin((-1)*ABs) %Ignorera gain margin (mer om den i reglerteknik). 
[gain_margin_after, phase_margin_after] = margin((-1)*ABsny)
%Matlab verkar inte gilla negativa system, så -1 behövs för att få rätt fasmarginal.


%% Implementera fantomnollan, undersök alla fall:
% Lph i serie med R1
delta_Lph = 1.5; %Effektivt om delta > 7
Lph = -R1/n;
AtINF_Lph = 1/(R1 + s*Lph); 

%Efter kompensering för MFM med fantomnolla:
%Lph serie med R1:
ABs_n_Lph = ABnoll*(1-s/n)/...
    ((1-s/P1)*(1-s/P2)*(1-s/(delta_Lph*n)));
Atn_Lph = AtINF_Lph*ABs_n_Lph/(1-ABs_n_Lph);

%Cph parallel med Rs
delta_Cph = 1.5;
Cph = -1/(Rs*n);
AtINF_Cph = AtINF;
%Efter kompensering för MFM med fantomnolla:
%Cph || Rs:
ABs_n_Cph = ABnoll*(1-s/n)/...
    ((1-s/P1)*(1-s/P2)*(1-s/(delta_Lph*n)))
Atn_Cph = AtINF_Cph*(-1)*ABs_n_Cph/(1-ABs_n_Cph) 


%% FIGURER
% Fasmarginal kollas "open loop", dvs frekvensen w0,
% där |AB(w0)| = 1 = 0dB, före = ABs och efter = ABs_n kompensering
%(Bode-funktionen behöver ibland ett (-1).* pga 'Phase unwrap')
figure(1)
bode((-1).*ABs,'b',(-1).*ABs_n_Lph,'k--',...
    (-1).*ABs_n_Cph, 'r--', (-1).*ABsny, 'y--')
grid on
h = gcr;
setoptions(h,'FreqUnits','Hz')

title('Slingförstärkning:');
legend('$A\beta(s)$','$A\beta\_{n,Lph}(s)$',...
    '$A\beta\_{n,Cph}(s)$','$A\beta\_{ny}(s)$',...
    'Interpreter','latex', 'Location','Best')

figure(2)
bode(At,'b',Atn_Lph,'k--',Atn_Cph,'r--', Atny,'y--')
grid on
h = gcr;
setoptions(h,'FreqUnits','Hz')
title('Den slutna förstärkningen, At');
legend('$A\_t$','$A\_{t,n,Lph}$',...
    '$A\_{t,n,Cph}$','$A\_{t_ny}$',...
    'Interpreter','latex','Location','Best')

figure(3)
stepplot((-1)*At, (-1)*Atn_Lph, (-1)*Atn_Cph, (-1)*Atny)
title('Stegsvaren före och efter kompensering');
legend('$A\_t$','$A\_{t,n,Lph}$',...
    '$A\_{t,n,Cph}$','$A\_{t_ny}$',...
    'Interpreter','latex','Location','Best')

%% HJÄLP FÃ–R ATT PLOTTA MÄTRESULTAT TILLSAMMANS (SAMMA FIGUR) MED SIMULERAD PRESTANDA (Kompenserat
% och okompenserat)
% W_labbet = [frekvensvektor från labbet].*(2*pi);
% At_labbet_kompenserat_dB = [mätresultat]
% At_labbet_kompenserat_fas = [mätresultat]
%på samma sätt lägger för At_okompenserat
W = [1:100:1e6].*(2*pi);
[MAG_At, PHASE_At] = bode(At,W);
for k = 1:length(W)
    dB_MAG_At(k) = 20*log10(MAG_At(1,1,k));
    phase_At(k) = PHASE_At(1,1,k);
end
% semilogx(W,dB_MAG_At,'b', W_labbet, At_labbet_kompenserat_dB,'r', W_labbet,At_labbet_okompenserat_dB,'k');
% semilogx(W,phase_At,'b', W_labbet, At_labbet_kompenserat_fas,'r', W_labbet,At_labbet_okompenserat_fas,'k');
% %xlabel och ylabel för axlarna

%figure(4);
%semilogx(W,dB_MAG_At,'b');%hold on; ... lägg till mätresultat
%figure(5);
%semilogx(W,phase_At,'b');%hold on; ...lägg till mätresultat

%% Simulated data
load Sim_R9_I.mat  
load Sim_R9_ejkomp.mat
load Sim_EjKomp.mat
load Sim_EjKompTran.mat
load Sim_Klippning.mat
load Sim_NarrowSvar.mat
load Sim_NarrowTran.mat

%% EJ KOMPENSERAT
% delete?
figure(11) 
subplot(2,1,1)
semilogx(LtSpiceEjKomp(:,1),LtSpiceEjKomp(:,2))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude')

subplot(2,1,2)
semilogx(LtSpiceEjKomp(:,1),LtSpiceEjKomp(:,3))
grid on
xlabel('Frequency (Hz)')
ylabel('??')

%% EJ KOMPENSERAT TRANSIENT &  NARROWBANDING TRANSIENT
figure(12)
plot(LtSpiceEjKompTran(:,1),LtSpiceEjKompTran(:,2), 'b-')
hold on
plot(LtSpiceNarrowTran(:,1),LtSpiceNarrowTran(:,2), 'r-')
grid on
axis([0 2.5e-4 4.5e-3 6.5e-3])
legend('Utan kompensering', 'Med kompensering')
xlabel('Time (s)')
ylabel('Amplitude')

%% KLIPPNINGEN
figure(13)
plot(LtSpiceKlippning(:,1),1e3*LtSpiceKlippning(:,2))
grid on
axis([0 1.5e-2 3.5 5.7])
xlabel('Time (s)')
ylabel('Amplitude (mA)')

%% NARROWBANDING SVAR
% TODO delete?
figure(14)
subplot(2,1,1)
semilogx(LtSpiceSvarNarrow(:,1),LtSpiceSvarNarrow(:,2))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude')

subplot(2,1,2)
semilogx(LtSpiceSvarNarrow(:,1),LtSpiceSvarNarrow(:,3))
grid on
xlabel('Frequency (Hz)')
ylabel('??')


%%

figure(6)
subplot(2,1,1)
semilogx(R9I2(:,1), R9I2(:,2), 'r-')
hold on
semilogx(R9Iejkomp(:,1), R9Iejkomp(:,2), 'b-')

axis([1e2 1e6 -80 -30])
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB) ')

subplot(2,1,2)
semilogx(R9I2(:,1), R9I2(:,3), 'r-')
hold on
semilogx(R9Iejkomp(:,1), R9Iejkomp(:,3), 'b-')

axis([1e2 1e6 -180 0])
grid on
xlabel('Frequency (Hz)')
ylabel('Polarity (deg)')


%% Save graphs in following format and settings
for k = 1:3
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
    