clear all;
close all;
clc;

%%Implantation de la chaine de transmission OFDM sans canal

%Durée symbole en nombre d’échantillons (Ts)
Ts =1;
Ns = 4;
Te=Ts/Ns;

%Nombre de bits générés
N = 16;
nb_bits = 320*N; % Multipe de N 

%Génération de l’information binaire
bits=randi([0,1],1,nb_bits);

%Mapping binaire à moyenne nulle : 0->-1, 1->1
Symboles=2*bits-1;

%%Emission 
%une seule porteuse parmi 16 est utilisée
num_porteuse = 5;
signal_ofdm1 = [];
for k=1:length(Symboles)
    bloc_symboles1=zeros(1,N);
    bloc_symboles1(num_porteuse) = Symboles(k);
    vecteur_echantillon_signal1 = ifft(bloc_symboles1);
    signal_ofdm1 = [signal_ofdm1 vecteur_echantillon_signal1];
end

figure
pwelch(signal_ofdm1)


% %deux porteuses parmi 16 sont utilisées
num_porteuse1 = 5;
num_porteuse2 = 10;
signal_ofdm2 = [];
for k=1:length(Symboles)/2
    bloc_symboles2=zeros(1,N);
    bloc_symboles2(num_porteuse1) = Symboles(k);
    vecteur_echantillon_signal2 = ifft(bloc_symboles2);
    signal_ofdm2 = [signal_ofdm2 vecteur_echantillon_signal2];
end
figure
pwelch(signal_ofdm2)
for k=length(Symboles)/2+1:length(Symboles)
    bloc_symboles2=zeros(1,N);
    bloc_symboles2(num_porteuse2) = Symboles(k);
    vecteur_echantillon_signal2 = ifft(bloc_symboles2);
    signal_ofdm2 = [signal_ofdm2 vecteur_echantillon_signal2];
end


figure
pwelch(signal_ofdm2)
%les 8 porteuses centrales sont utilisées
num_porteuse3 = [5:12];
signal_ofdm3 = [];
k = 1;
x =length(Symboles)/length(num_porteuse3);
while k <= length(num_porteuse3);
   for j=(k-1)*x+1:k*x ;
        bloc_symboles3=zeros(1,N);
        bloc_symboles3(num_porteuse3(k)) = Symboles(j);
        vecteur_echantillon_signal3 = ifft(bloc_symboles3);
        signal_ofdm3 = [signal_ofdm3 vecteur_echantillon_signal3];
        
        
    end;
    k=k+1;
end;

figure
pwelch(signal_ofdm3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%lorsque toutes les porteuses sont utilisées
signal_ofdm = [];
for k=1:N:length(Symboles)-N+1
    bloc_symboles = Symboles(k:k+N-1);
    vecteur_echantillon_signal = ifft(bloc_symboles);
    signal_ofdm = [signal_ofdm vecteur_echantillon_signal];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Reception sans canal
Symboles_recu = [];
for k = 1:N:length(Symboles)-N+1
    vecteur_echantillon_signal_recu = signal_ofdm;
    bloc_symboles_recu = fft(vecteur_echantillon_signal_recu(k:k+N-1));
    Symboles_recu = [Symboles_recu bloc_symboles_recu];
end
%Vérifier si le TEB est nul
 erreur = 0;
    for k=1:length(Symboles)
        if sign(real(Symboles_recu(k))) ~= Symboles(k)
            erreur=erreur+1;
        end
    end
assert(erreur == 0);