clear all;
close all;
clc;

%%Implantation de la chaine de transmission OFDM sans canal

%Durée symbole en nombre d’échantillons (Ts)

M = 4 ;%QPSK

%Nombre de bits générés
N = 16;
nb_bits = 1000*N; % Multipe de N 


Es_N0_dB = [0:30]; % Eb/N0 values

%Preallocations
nErr_zfinf=zeros(1,length(Es_N0_dB));
TEB_pratic=[]
for ii = 1:length(Es_N0_dB)
   
    % QPSK symbol generations
   bits = rand(2,nb_bits)>0.5; % generating 0,1 with equal probability
   Symboles = ((1-2*bits(1,:))+1j*(1-2*bits(2,:))); % QPSK modulation following the BPSK rule for each quadatrure component: {0 -> +1; 1 -> -1} 
   %%lorsque toutes les porteuses sont utilisées
    signal_ofdm = [];
    for k=1:N:length(Symboles)-N+1
       bloc_symboles = Symboles(k:k+N-1);
       vecteur_echantillon_signal = ifft(bloc_symboles);
       signal_ofdm = [signal_ofdm vecteur_echantillon_signal];
    end
    %generer le bruit
   Pxe = mean(abs(signal_ofdm).^2);
   EbN0=10^(Es_N0_dB(ii)/10);
   sigI2=Pxe/(2*log2(M)*EbN0);
   sigQ2=Pxe/(2*log2(M)*EbN0);
   
   %n = sqrt(sig2b)*randn(1,N+Lc-1); % white gaussian noise, BPSK Case
    n = sqrt(sigI2/2)*randn(1,nb_bits)+1j*sqrt(sigQ2/2)*randn(1,nb_bits); % white gaussian noise, QPSK case
    
   % Adding Noise
   y = signal_ofdm + n; % additive white gaussian noise
   
   Symboles_recu = [];
   for k = 1:N:length(Symboles)-N+1
         vecteur_echantillon_signal_recu = y;
         bloc_symboles_recu = fft(vecteur_echantillon_signal_recu(k:k+N-1));
         Symboles_recu = [Symboles_recu bloc_symboles_recu];
   end
   bits_recu = zeros(1,2*nb_bits);
    ind1 = find(real(Symboles_recu)>0);
    ind2 = find(imag(Symboles_recu)>0);
    bits_recu(1:2:end) = real(Symboles_recu)>0;
    bits_recu(2:2:end) = imag(Symboles_recu)>0;
    
    bits_recu = ones(1,2*nb_bits) - bits_recu;

   %compte du nombre d'erreurs
    er = 0;
    for k=1:nb_bits
        if bits_recu(k) ~= bits(k)
            er=er+1;
        end
    end
    TEB_pratic(ii) = er/nb_bits;
end
figure;semilogy(Es_N0_dB,TEB_pratic);
title("le taux d'erreur binaire");
xlabel("les valeurs proposés de (Eb/N0)dB");
ylabel("l'échelle logarithmique");