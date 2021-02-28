clear all;
close all;
clc;

%%%Implantation  de  la  chaine  de  transmission  OFDM  avec  canalmultitrajets, sans bruit

%Durée symbole en nombre d’échantillons (Ts)
Ts =1;
Ns = 4;
Te=Ts/Ns;

%Nombre de bits générés
N = 14;   % le nombre des porteuses OFDM utilises
nb_bits = 1000*N; % Multipe de N 

%Génération de l’information binaire
bits=randi([0,1],1,nb_bits);

%Mapping binaire à moyenne nulle : 0->-1, 1->1
Symboles=2*bits-1;

%%Modulation  
signal_ofdm = [];
for k=1:N:nb_bits-N+1
    bloc_symboles=Symboles(k:k+N-1);
    vecteur_echantillon_signal = ifft(bloc_symboles);
    signal_ofdm = [signal_ofdm vecteur_echantillon_signal];
end


%%Le canal
alpha_0 = 0.227;
alpha_1 = 0.46;
alpha_2 = 0.688;
alpha_3 = 0.460;
alpha_4 = 0.227;
h = [alpha_0 alpha_1 alpha_2 alpha_3 alpha_4];
H = fft(h,4096);

figure;
subplot(2,1,1);
plot(linspace(1,N,length(H)),abs(H));
title("le Trace la réponse en fréquence de module du canal de propagation.");
xlabel("Les porteuses.");
ylabel("le module de H");
subplot(2,1,2);
plot(angle(H));
title("le Trace la réponse en fréquence de phase du canal de propagation.");
xlabel("fréquences normalisés");
ylabel("la phase de H");


%%le signal à la sortie du canal sans intervalle de garde

SignalRecu=filter(h,1,signal_ofdm) ;% A la sortie de canal 
figure;
plot((real(SignalRecu)));
title("le passage du signal OFDM dans le canal de propagation");
xlabel("fréquences normalisés");
ylabel("l'amplitude");

figure;
subplot(2,1,1);
pwelch(signal_ofdm)%Affichage de la DSP du signal généré
title("la DSP avant passage dans le canal");

subplot(2,1,2);
pwelch(SignalRecu)%Affichage de la DSP du signal généré
title("la DSP aprés passage dans le canal");

%%les constellations obtenues en réception sur deux porteuses sans intervalle de garde
%Démodulation 
Symboles_recu = []; % a la sortie du récepteur
for k = 1:N:nb_bits-N+1
    bloc_symboles_recu = fft(SignalRecu(k:k+N-1));
    Symboles_recu = [Symboles_recu bloc_symboles_recu];
end

%%les constellations en sortie de echantillonnage
np_1 = 3;
np_2 = 14;
bloc_1 = [];
bloc_2 = [];
for k = 1:N:nb_bits-N+1
    bloc_1 = [bloc_1 Symboles_recu(k+np_1-1)];
    bloc_2 = [bloc_2 Symboles_recu(k+np_2-1)];
end
figure;plot(real(bloc_1),imag(bloc_1),'ko');
title("les constellations obtenues en réception sur la porteuse 3 sans IG");
xlabel("phase");
ylabel("quadrature");
xlim([-4 4]);
ylim([-4 4]);
grid on

figure;plot(real(bloc_2),imag(bloc_2),'ko');
title("les constellations obtenues en réception sur la porteuse 14 sans IG");
xlabel("phase");
ylabel("quadrature");
xlim([-4 4]);
ylim([-4 4]);
grid on

%Demmaping
bits_recu_IG = (sign(real(Symboles_recu)) + 1)/2;

%TEB
er = 0;
    for k=1:nb_bits-N+1
        if bits_recu_IG(k) ~= bits(k)
            er=er+1;
        end
    end
 TEB = er/nb_bits
 
 
 %%le signal à la sortie du canal avec l'ajout de  l'intervalle de garde
 
signal_ofdm_IG = []; 
for k=1:N:nb_bits-N+1
    signal_recu_IG=signal_ofdm(k:k+N-1);
    vectsignal_recu_IG = [zeros(1,length(h)-1) signal_recu_IG];
    signal_ofdm_IG = [signal_ofdm_IG vectsignal_recu_IG];
end

SignalRecu_IG = filter(h,1,signal_ofdm_IG) ;% A la sortie de canal 


%%les constellations obtenues en réception sur deux porteuses sans intervalle de garde
%Démodulation 
Symboles_recu_IG = []; % a la sortie du récepteur
for k = 1:nb_bits/N
    bloc_OFDM_recu_IG = SignalRecu_IG(k*length(h)+(k-1)*(N-1):k*length(h)+k*(N-1));
    bloc_symboles_recu_IG = fft(bloc_OFDM_recu_IG);
    Symboles_recu_IG = [Symboles_recu_IG bloc_symboles_recu_IG];
end

%%les constellations en sortie de echantillonnage
np_1 = 3;
np_2 = 14;
bloc_1_IG = [];
bloc_2_IG = [];
for k = 1:N:nb_bits-N+1
    bloc_1_IG = [bloc_1_IG Symboles_recu_IG(k+np_1-1)];
    bloc_2_IG = [bloc_2_IG Symboles_recu_IG(k+np_2-1)];
end
figure;plot(real(bloc_1_IG),imag(bloc_1_IG),'ko');
title("les constellations obtenues en réception sur la porteuse 3 aprés l'ajout d'intervalle de garde");
xlim([-4 4]);
ylim([-4 4]);
grid on

figure;plot(real(bloc_2_IG),imag(bloc_2_IG),'ko');
title("les constellations obtenues en réception sur la porteuse 14 aprés l'ajout d'intervalle de garde");
xlim([-4 4]);
ylim([-4 4]);
grid on

%Demmaping
bits_recu_IG = (sign(real(Symboles_recu_IG)) + 1)/2;

%TEB
er = 0;
    for k=1:nb_bits-N+1
        if bits_recu_IG(k) ~= bits(k)
            er=er+1;
        end
    end
 TEB = er/nb_bits
 
 
 %%le signal à la sortie du canal avec l'ajout d'un préfixe cyclique

signal_ofdm_PC = []; 
for k=1:N:nb_bits-N+1
    bloc_symboles_PC=signal_ofdm(k:k+N-1);
    vect_signal_recu_PC = [signal_ofdm(k+N-(length(h)-1):k+N-1) bloc_symboles_PC];
    signal_ofdm_PC = [signal_ofdm_PC vect_signal_recu_PC];
end

SignalRecu_PC = filter(h,1,signal_ofdm_PC) ;% A la sortie de canal 

%Démodulation 
Symboles_recu_PC = []; % a la sortie du récepteur
for k = 1:nb_bits/N
    bloc_OFDM_recu_PC = SignalRecu_PC(k*length(h)+(k-1)*(N-1):k*length(h)+k*(N-1));
    bloc_symboles_recu_PC = fft(bloc_OFDM_recu_PC);
    Symboles_recu_PC = [Symboles_recu_PC bloc_symboles_recu_PC];
end

%%les constellations en sortie de echantillonnage
np_1 = 3;
np_2 = 14;
bloc_1_PC = [];
bloc_2_PC = [];
for k = 1:N:nb_bits-N+1
    bloc_1_PC = [bloc_1_PC Symboles_recu_PC(k+np_1-1)];
    bloc_2_PC = [bloc_2_PC Symboles_recu_PC(k+np_2-1)];
end
figure;plot(real(bloc_1_PC),imag(bloc_1_PC),'ko');
title("les constellations obtenues en réception sur la porteuse 3 avec l'ajout de préfixe cyclique");
xlim([-4 4]);
ylim([-4 4]);
grid on

figure;plot(real(bloc_2_PC),imag(bloc_2_PC),'ko');
title("les constellations obtenues en réception sur la porteuse 14 avec l'ajout de préfixe cyclique");
xlim([-4 4]);
ylim([-4 4]);
grid on

%Demmaping
bits_recu_PC = (sign(real(Symboles_recu_PC)) + 1)/2;

%TEB
er = 0;
    for k=1:nb_bits-N+1
        if bits_recu_PC(k) ~= bits(k)
            er=er+1;
        end
    end
 TEB = er/nb_bits

 
 
 %%le signal à la sortie du canal avec l'ajout d'un préfixe cyclique et
 %%Egalisation
 C =[];
 for k=0:N-1
     c = alpha_0+alpha_1*exp(-(2*1i*pi*k)/N) + alpha_2*exp(-(2*1i*pi*2*k)/N)+ alpha_3*exp(-(2*1i*pi*3*k)/N)+ alpha_4*exp(-(2*1i*pi*4*k)/N);
     C =[C c];
 end
 H = 1./C;
 
%Symboles_recu_PC_ZF  = filter(H,1,Symboles_recu_PC);

Symboles_recu_PC_ZF = []; % a la sortie du récepteur
for k = 1:nb_bits/N
    bloc_OFDM_recu_PC = SignalRecu_PC(k*length(h)+(k-1)*(N-1):k*length(h)+k*(N-1));
    bloc_symboles_recu_PC = fft(bloc_OFDM_recu_PC);
    bloc_symboles_recu_PC = bloc_symboles_recu_PC.*H;
    Symboles_recu_PC_ZF = [Symboles_recu_PC_ZF bloc_symboles_recu_PC];
end

%%les constellations en sortie de echantillonnage
np_1 = 3;
np_2 = 14;
bloc_1_PC_ZF = [];
bloc_2_PC_ZF = [];
for k = 1:N:nb_bits-N+1
    bloc_1_PC_ZF = [bloc_1_PC_ZF Symboles_recu_PC_ZF(k+np_1-1)];
    bloc_2_PC_ZF = [bloc_2_PC_ZF Symboles_recu_PC_ZF(k+np_2-1)];
end
figure;plot(real(bloc_1_PC_ZF),imag(bloc_1_PC_ZF),'ko');
title("les constellations obtenues en réception sur la porteuse 3 avec l'ajout de PC et ZFE");
xlim([-4 4]);
ylim([-4 4]);
grid on

figure;plot(real(bloc_2_PC_ZF),imag(bloc_2_PC_ZF),'ko');
title("les constellations obtenues en réception sur la porteuse 14 avec l'ajout de PC et ZFE");
xlim([-4 4]);
ylim([-4 4]);
grid on

%Demmaping
bits_recu_PC_ZF = (sign(real(Symboles_recu_PC_ZF)) + 1)/2;

%TEB
er = 0;
    for k=1:nb_bits-N+1
        if bits_recu_PC_ZF(k) ~= bits(k)
            er=er+1;
        end
    end
 TEB = er/nb_bits