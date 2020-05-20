clear; close all

% Definiera Data 
VT = 25e-3;
Bf1 = 250; %Kolla datablad!
Bf2 = Bf1;
C1 = 100e-9; %Ers�tter Cpi1_prim
C2 = 2.2e-6; %Ers�tter Cpi2
C2ny = 40e-6; %ers�tter C2 vid capacitive narrowbandning

% F�rst�rkare 
Rs = 50; %K�llan �r inte ideal
R1 = 100;
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
Atn_Lph = AtINF_Lph*ABs_n_Lph/(1-ABs_n_Lph);

%Cph parallel med Rs
delta_Cph = 1.5;
Cph = -1/(Rs*n);
AtINF_Cph = AtINF;
%Efter kompensering f�r MFM med fantomnolla:
%Cph || Rs:
ABs_n_Cph = ABnoll*(1-s/n)/...
    ((1-s/P1)*(1-s/P2)*(1-s/(delta_Lph*n)))
Atn_Cph = AtINF_Cph*(-1)*ABs_n_Cph/(1-ABs_n_Cph) 


%% Figures (Matlab simulation)
% Fasmarginal kollas "open loop", dvs frekvensen w0,
% d�r |AB(w0)| = 1 = 0dB, f�re = ABs och efter = ABs_n kompensering
%(Bode-funktionen beh�ver ibland ett (-1).* pga 'Phase unwrap')
figure(1)
bode((-1).*ABs,'b',...
    (-1).*ABs_n_Lph,'k--',...
    (-1).*ABs_n_Cph, 'r--',...
    (-1).*ABsny, 'y--')
grid on
setoptions(gcr,'FreqUnits','Hz')
set(findall(gcf, 'Type', 'Line'),'LineWidth',2);

title('Slingf�rst�rkning:');
[~,legObj] = legend('$A\beta(s)$','$A\beta\_{n,Lph}(s)$',...
    '$A\beta\_{n,Cph}(s)$','$A\beta\_{ny}(s)$',...
    'Interpreter','latex',...
    'Location','Best',...
    'Fontsize', 11);
set(findobj(legObj,'type','line'),'linewidth',2)
% ---------------------------------------------
figure(2)
bode(At,'b',Atn_Lph,'k--',Atn_Cph,'r--', Atny,'y--')
grid on
setoptions(gcr,...
    'FreqUnits','Hz',...
    'Xlim', [1e3 1e6])
title('Den slutna f�rst�rkningen, At');
[~,legObj] = legend('$A\_t$','$A\_{t,n,Lph}$',...
    '$A\_{t,n,Cph}$','$A\_{t_ny}$',...
    'Interpreter','latex',...
    'Location','Best',...
    'Fontsize', 11);
set(findobj(legObj,'type','line'),'linewidth',2)
% ---------------------------------------------
figure(3)
stepplot((-1)*At, (-1)*Atn_Lph, (-1)*Atn_Cph, (-1)*Atny)
title('Stegsvaren f�re och efter kompensering');
legend('$A\_t$','$A\_{t,n,Lph}$',...
    '$A\_{t,n,Cph}$','$A\_{t_ny}$',...
    'Interpreter','latex','Location','Best')

%% Figures (LTspice + Matlab)
load ltspice_data_old.mat % Old data from LTspice

figure(4)
subplot(2,1,1)
semilogx(R9Iejkomp(:,1), R9Iejkomp(:,2), 'b-')
hold on
semilogx(R9I2(:,1), R9I2(:,2), 'r-')

axis([1e2 1e6 -80 -30])
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB) ')

subplot(2,1,2)
semilogx(R9Iejkomp(:,1), R9Iejkomp(:,3), 'b-')
hold on
semilogx(R9I2(:,1), R9I2(:,3), 'r-')

axis([1e2 1e6 -180 0])
legend('Utan kompensering','Med kompensering',...
    'Interpreter','latex', 'Location','Best')
xlabel('Frequency (Hz)')
ylabel('Polarity (deg)')

% EJ KOMPENSERAT TRANSIENT &  NARROWBANDING TRANSIENT
figure(5)
plot(1e6*LtSpiceEjKompTran(:,1),LtSpiceEjKompTran(:,2), 'b-')
hold on
plot(1e6*LtSpiceNarrowTran(:,1),LtSpiceNarrowTran(:,2), 'r-')
axis([0 2.5e2 4.5e-3 6.5e-3])
legend('Utan kompensering', 'Med kompensering',...
    'Interpreter','latex', 'Location','Best')
xlabel('Time (us)')
ylabel('Amplitude')

% INGEN KLIPPNING
figure(6)
plot(1e3*LtSpiceKlippning(:,1),1e3*LtSpiceKlippning(:,2))
axis([0 1.5e1 3.5 5.7])
xlabel('Time (ms)')
ylabel('Amplitude (mA)')

%% Save graphs in following format and settings
for k = 1:6
    figname = ['figure', num2str(k)];
    figure(k)
    title('')   % Remove Title before saving
%    saveas(gcf,figname,'epsc')
end    

% ------ Shared graph settings -------
function customPlotSettings()
    set(gca,...
        'Fontsize', 12,...
        'linewidth', 1,...
        'FontName', 'Arial');
end
 