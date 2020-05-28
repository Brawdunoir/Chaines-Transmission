clear;
close all;

%% Chaine de référence
Fe = 1000; % Fréquence d'échantillonage
Te=1/Fe;
Nb=50;      % Nombre de bits émis
Ns = 20;    % Nombre d'echantillons par symbole
Ts = Ns*Te;
Fs = 1/Ts;
Fmax = 4/Ts;

%Génération aléatoire 
bits = randi([0,1],1,Nb);

% Mapping binaire à Moyenne nulle et génération de la suite de Diracs pondérés par les symboles
suite_diracs = kron(2*bits - 1, [1 zeros(1,Ns-1)]);

%Génération de la réponse impulsionnelle du filtre de mise en forme (NRZ)
h = ones(1,Ns);

%Filtrage de mise en forme
x = filter(h, 1, suite_diracs);

%Affichage du signal généré
figure(1); plot(x)
title('1ere chaine - Signal généré')
xlabel('t (s)')
ylabel('Amplitude (V)');

% Calcul de la Dsp du signal 

X=fft(x);
Dsp = abs(X.^2);

%Affichage de la DSP du signal généré

frequence= linspace (-Fe/2,Fe/2,length(x));

figure(2);  
semilogy(frequence,fftshift(Dsp))
title('1ere chaine - DSP du Signal généré')
xlabel('f (hz)')
ylabel('Amplitude');

% Signal en sortie du filtre de réception
hr_ref = ones(1,Ns);
z=filter(hr_ref,1,x)/Ns;
ze = z(Ns:Ns:end);
%Affichage du signal reçu
figure(3); plot(z)
title('1ere chaine - Signal en sortie de reception')
xlabel('t (s)')
ylabel('Amplitude (V)');

%Diagramme de l'oeil

Oeil = reshape(z, Ns, length(z)/Ns);

%Tracé du Diagramme de l'oeil
figure(4); plot(Oeil);
title("1ere chaine - Diagramme de l'oeil")
xlabel('t (s)')
ylabel('Amplitude (V)');

% Demapping binaire
seuil = 0
bits_r =(ze>seuil); % Cela réalise bien le démapping et la détection par rapport au seuil

% Calcul et affichage du TEB

TEB_1 = sum(abs(bits - bits_r))/(length(bits));

fprintf('1ere chaine TEB = %d\n', TEB_1)

%% 1ere chaine avec bruit qui marche mieux

% Initialisation des constantes
Nb=100; % Nombre de bits emis
Ns=10;      % Nombre d'echantillons par symbole
M=2;        % Taille de la constellation
h=ones(1,Ns);   % reponse impulsionnelle du filtre de mise en forme
t0=Ns;
hr=ones(1,Ns);   % reponse impulsionelle du filtre de reception
Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);

TEB1=zeros(1,length(Eb_sur_N0_dB)); % Initialisation du TEB que l'on veut simuler
NerrLim=100;    % Le nombre d'erreur à partir du quel on s'arrête
Nerrcomptees=zeros(1,length(Eb_sur_N0_dB)); % Le nombre d'erreur que l'on compte
Nsimu=zeros(1,length(Eb_sur_N0_dB));    % Le nombre de simulation


for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

        % Emission
        bits = randi([0 1],1,Nb);
        symboles = 2*bits-1;
        suite_diracs = kron(symboles,[1 zeros(1,Ns-1)]);
        x=filter(h,1,suite_diracs);

        % Canal

        Pr=mean(abs(x).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        r=x+sqrt(sigma)*randn(1,Nb*Ns);

        % Reception
        z=filter(hr,1,r)/Ns;
        ze=z(t0:Ns:Nb*Ns);
        bits_estimes=(ze>seuil);
        Nerrcomptees(k)=Nerrcomptees(k)+sum(bits~=bits_estimes);
        Nsimu(k)=Nsimu(k)+1;

    end

    TEB1(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % On calcule le TEB

end
Pb=qfunc(sqrt(2*Eb_sur_N0_lin));

% Tracé de la omparaison entre le TEB simulé et le TEB théorique

figure(5)
title('1ere chaine - Comparaison du TEB simulé avec le TEB théorique')
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB1,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique','TEB simulé')
grid
set(gca,'FontSize',12)

%% Deuxieme chaine

Nb=100; % Nombre de bits emis

Ns=10;      % Nombre d'echantillons par symbole
M=2;        % Taille de la constellation
h=ones(1,Ns);   % reponse impulsionnelle du filtre de mise en forme
t0=Ns;
hr=[ones(1,Ns/2) zeros(1,Ns/2)];   % reponse impulsionelle du filtre de reception

% Génération aléatoire
bits = randi([0 1],1,Nb);

% Mapping binaire à moyenne nulle et génération de la suite de Diracs pondérés par les symboles
symboles = 2*bits-1;
suite_diracs = kron(symboles,[1 zeros(1,Ns-1)]);

% Passage dans le filtre de mise en forme
x=filter(h,1,suite_diracs);

% Reception
z=filter(hr,1,x)/(Ns/2);
z=[0 z(1:end-1)];

% Tracé du signal en sortie de réception
figure(6); plot(z)
title('2eme chaine - Signal en sortie de réception')
xlabel('t (s)')
ylabel('Amplitude (V)')
 
% Diagramme de l'oeil
Oeil = reshape(z,Ns, length(z)/(Ns));

% Tracé du diagramme de l'oeil
figure(7); plot(Oeil);
title("2eme chaine - Diagramme de l'oeil")
xlabel('t (s)')
ylabel('Amplitude (V)');
axis([0.5 Ns+0.5 -1 1])

% Demapping binaire
seuil=0;
ze = z(t0:Ns:Nb*Ns);
bits_r =(ze>seuil); % Cela réalise le démmaping et la détection par rapport au seuil

% Calcul et affichage du TEB
TEB_2 = sum(abs(bits - bits_r)/(length(bits)));
fprintf('2eme chaine TEB = %d\n', TEB_2)

%% 2eme chaine avec bruit

%Initialisation des constantes
Nb=100; % Nombre de bits emis
Ns=10;      % Nombre d'echantillons par symbole
M=2;        % Taille de la constellation
h=ones(1,Ns);   % reponse impulsionnelle du filtre de mise en forme
t0=Ns;
hr=[ones(1,Ns/2) zeros(1,Ns/2)];   % reponse impulsionelle du filtre de reception
Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);

TEB2=zeros(1,length(Eb_sur_N0_dB)); % Initialisation du TEB que l'on veut simuler
NerrLim=100;     % Le nombre d'erreur à partir du quel on s'arrête
Nerrcomptees=zeros(1,length(Eb_sur_N0_dB)); % Le nombre d'erreurs que l'on compte
Nsimu=zeros(1,length(Eb_sur_N0_dB));    % Le nombre de simulation


for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

        % Emission
        bits = randi([0 1],1,Nb);
        symboles = 2*bits-1;
        suite_diracs = kron(symboles,[1 zeros(1,Ns-1)]);
        x=filter(h,1,suite_diracs);

        % Canal

        Pr=mean(abs(x).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        r=x+sqrt(sigma)*randn(1,Nb*Ns);

        % Reception
        z=filter(hr,1,r)/(Ns/2);
        z=[0 z(1:end-1)];
        ze=z(t0:Ns:Nb*Ns);
        bits_estimes=(ze>seuil);
        Nerrcomptees(k)=Nerrcomptees(k)+sum(bits~=bits_estimes);
        Nsimu(k)=Nsimu(k)+1;

    end

    TEB2(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % Calcul du TEB

end


% Tracé de la comparaison entre le TEB simulé et le TEB théorique
Pb=qfunc(sqrt(Eb_sur_N0_lin));

figure(8)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB2,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique','TEB simulé')
title("2eme chaine - comparaison entre le TEB simulé et le TEB théorique")
grid
set(gca,'FontSize',12)

% Comparaison entre le TEB simulé et le TEB de la chaine de référence

Pb=qfunc(sqrt(2*Eb_sur_N0_lin));

figure(9)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB2,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique de référence','TEB simulé')
title("2eme chaine - comparaison entre le TEB simulé et le TEB théorique de référence")
grid
set(gca,'FontSize',12)

%% Comparaison de l'efficacité spectrale de la haine étudiée avec la chaine de référence

% Emission
bits = randi([0 1],1,Nb);
symboles = 2*bits-1;
suite_diracs = kron(symboles,[1 zeros(1,Ns-1)]);
x=filter(h,1,suite_diracs);

% Canal
Pr=mean(abs(x).^2);
sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
r=x+sqrt(sigma)*randn(1,Nb*Ns);

% Reception chaine étudiée
z=filter(hr,1,r)/(Ns/2);
z=[0 z(1:end-1)];
ze=z(t0:Ns:Nb*Ns);
bits_r = (ze +1)/2;
Dsp_r = fft(bits_r);

% Reception chaine de référence
z_ref=filter(hr_ref,1,r)/(Ns/2);
z_ref=[0 z_ref(1:end-1)];
ze_ref=z_ref(t0:Ns:Nb*Ns);
bits_r_ref = (ze_ref +1)/2;
Dsp_r_ref = fft(bits_r_ref);

% Tracé des DSP
frequence= linspace(-Fe/2,Fe/2,length(bits_r));

figure(10)
semilogy(frequence,fftshift(abs(Dsp_r.^2)),'r-','LineWidth',1.2)
hold on
semilogy(frequence,fftshift(abs(Dsp_r_ref.^2)),'b-','LineWidth',1.2)
xlabel('f (Hz)')
ylabel('DSP')
legend('DSP chaine étudiée','DSP chaine de référence')
title("2eme chaine - comparaison des DSP")
grid
set(gca,'FontSize',12)

%% 3eme chaine

Nb=100;     % Nombre de bits emis
Fe = 12000; % Fréquence d'échantillonnage en Hz
Te = 1/Fe;
Rs = 3000; % rythme symbole (en symboles par seconde)
Ts = 1/Rs;
Ns = Ts/Te;
alpha = 0.5;
h = rcosdesign(alpha,Nb,Ns); % reponse impulsionnelle du filtre de mise en forme
h_ref = ones(1,Ns); % reponse impulsionnelle du filtre de mise en forme de référence
t0 = 1;
hr = rcosdesign(alpha,Nb,Ns);   % reponse impulsionelle du filtre de reception
hr_ref = ones(1,Ns);        % reponse impulsionelle du filtre de reception  de référence

bits = randi([0 1],1,Nb);  
symboles = 2*bits-1;
suite_diracs = kron(symboles,[1 zeros(1,Ns-1)]);

retard = Nb*Ns/2;

% Passage dans le filtre de mise en forme
x = filter(h,1,[suite_diracs zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
x = x(retard+1:end); % On supprime le retard (les "retards" premières valeurs)

% Passage dans le filtre de réception
z = filter(hr,1,[x zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)


figure(11); plot(fftshift(z))
title('3eme chaine - Signal en sortie de réception')
xlabel('t (s)')
ylabel('Amplitude (V)')
 
%Diagramme de l'oeil
Oeil = reshape(z,Ns, length(z)/Ns);

figure(12); plot(Oeil);
title("3eme chaine - Diagramme de l'oeil")
xlabel('t (s)')
ylabel('Amplitude (V)');


% Demapping et détection par rapport au seuil
seuil=0; % seuil optimal
ze = z(t0:Ns:Nb*Ns);
bits_r = (ze>seuil);% Cela réalise le démmaping et la détection par rapport au seuil

% Calcul et affichage du TEB

TEB_3 = sum(abs(bits - bits_r)/(length(bits)));
fprintf('3eme chaine TEB = %d\n', TEB_3)

%% 3eme chaine avec bruit

%Initialisation des constantes
Fe = 12000; % Fréquence d'échantillonnage en Hz
Te = 1/Fe;
Rs = 3000; % rythme symbole (en symboles par seconde)
Ts = 1/Rs;
Ns = Ts/Te;      % Taille de la constellation
alpha = 0.5;
h=rcosdesign(alpha,Nb,Ns);   % reponse impulsionnelle du filtre de mise en forme
t0=1;

hr=rcosdesign(alpha,Nb,Ns);   % reponse impulsionelle du filtre de reception
Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);

TEB3=zeros(1,length(Eb_sur_N0_dB)); % Initialisation du TEB que l'on veut simuler
NerrLim=100;     % Le nombre d'erreur à partir du quel on s'arrête
Nerrcomptees=zeros(1,length(Eb_sur_N0_dB)); % Le nombre d'erreur que l'on compte
Nsimu=zeros(1,length(Eb_sur_N0_dB));    % Le nombre de simulation


for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

        % Emission
        bits = randi([0 1],1,Nb);
        symboles = 2*bits-1;
        suite_diracs = kron(symboles,[1 zeros(1,Ns-1)]);
        x = filter(h,1,[suite_diracs zeros(1,retard)]);
        x = x(retard+1:end);

        % Canal

        Pr=mean(abs(x).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        r=x+sqrt(sigma)*randn(1,Nb*Ns);

        % Reception
        z = filter(hr,1,[r zeros(1,retard)]);
        z = z(retard+1:end);
        ze=z(t0:Ns:Nb*Ns);
        bits_estimes=(ze>seuil);
        Nerrcomptees(k)=Nerrcomptees(k)+sum(bits~=bits_estimes);
        Nsimu(k)=Nsimu(k)+1;

    end

    TEB3(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % Calcul du TEB

end


% Tracé de la comparaison entre le TEB simulé et le TEB théorique
Pb=qfunc(sqrt(2*Eb_sur_N0_lin));

figure(13)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB3,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique','TEB simulé')
title("3eme chaine - comparaison entre le TEB simulé et le TEB théorique")
grid
set(gca,'FontSize',12)

% Comparaison entre le TEB simulé et le TEB de la chaine de référence

Pb=qfunc(sqrt(2*Eb_sur_N0_lin));

figure(14)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB3,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique de référence','TEB simulé')
title("3eme chaine - comparaison entre le TEB simulé et le TEB théorique de référence")
grid
set(gca,'FontSize',12)


% Comparaison de l'efficacité spectrale de la chaine étudiée avec la chaine
% de référence

% Emission
bits = randi([0 1],1,Nb);
symboles = 2*bits-1;
suite_diracs = kron(symboles,[1 zeros(1,Ns-1)]);
retard = Nb*Ns/2;

x_ref = filter(h_ref,1,suite_diracs);

% Passage dans le filtre de mise en forme pour la chaine étudié
x = filter(h,1,[suite_diracs zeros(1,retard)]);
x = x(retard+1:end);

% Canal
Pr=mean(abs(x).^2);
sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
r=x+sqrt(sigma)*randn(1,Nb*Ns);
r_ref = x_ref+sqrt(sigma)*randn(1,Nb*Ns);

% Passage dans le filtre de réception de la chaine étudiée
z = filter(hr,1,[r zeros(1,retard)]);
z = z(retard+1:end);
z=[0 z(1:end-1)];
ze=z(t0:Ns:Nb*Ns);
bits_r = (ze +1)/2;
Dsp_r = fft(bits_r);

% Reception chaine de référence
z_ref=filter(hr_ref,1,r_ref)/(Ns/2);
z_ref=[0 z_ref(1:end-1)];
ze=z_ref(t0:Ns:Nb*Ns);
bits_r_ref = (ze +1)/2;
Dsp_r_ref = fft(bits_r_ref);
Fe = 12000;

% Tracé des DSP
frequence= linspace(-Fe/2,Fe/2,length(bits_r));

figure(15)
semilogy(frequence,fftshift(abs(Dsp_r.^2)),'r-','LineWidth',1.2)
hold on
semilogy(frequence,fftshift(abs(Dsp_r_ref.^2)),'b-','LineWidth',1.2)
xlabel('f (Hz)')
ylabel('DSP')
legend('DSP chaine étudiée','DSP chaine de référence')
title("3eme chaine - comparaison de l'efficacité spectrale de la chaine étudiée et de la chaine de référence")
grid
set(gca,'FontSize',12)


% Avec un passage dans un canal de transmission

bits = randi([0 1],1,Nb);
symboles = 2*bits-1;
suite_diracs = kron(symboles,[1 zeros(1,Ns-1)]);
retard = Nb*Ns/2;

x = filter(h,1,[suite_diracs zeros(1,retard)]); 
x = x(retard+1:end); 

% Canal de transmission

N0 = 21; %Ordre du filtre
fc1 = 1500; 
fc2 = 3000;

%Génération de la réponse impulsionnelle

k = -(N0-1)/2 : 1 : (N0-1)/2;
passebas1 = 2*fc1/Fe*sinc(2*fc1/Fe*k);
passebas2 = 2*fc2/Fe*sinc(2*fc2/Fe*k);

x_a = filter(passebas1,1, [x zeros(1,retard)]);
x_b = filter(passebas2,1, [x zeros(1,retard)]);

%Tracé de la réponse impulsionnelle
figure(16);
plot(k*Fe/N0,fftshift(abs(fft(passebas1))),k*Fe/N0,fftshift(abs(fft(passebas2))))
legend('Passe-bas (fc=1500Hz)','Passe-bas (fc=3000Hz')
xlabel('frequence (Hz)');
ylabel('Amplitude (V)');



% Passage dans le filtre de réception
z_a = filter(hr,1,[x_a zeros(1,retard)]);
z_a = z_a(retard+1:end);
z_b = filter(hr,1,[x_b zeros(1,retard)]); 
z_b = z_b(retard+1:end);

figure(17);
subplot(1,2,1)
plot(z_a)
title('3eme chaine - Signal (a) en sortie de réception (fc = 1500 Hz)')
xlabel('t (s)')
ylabel('Amplitude (V)')

subplot(1,2,2)
plot(z_b)
title('3eme chaine - Signal (b) en sortie de réception (fc = 3000 Hz)')
xlabel('t (s)')
ylabel('Amplitude (V)')


%Diagramme de l'oeil

Oeil_a = reshape(z_a,Ns, length(z_a)/(Ns));
Oeil_b = reshape(z_b,Ns, length(z_b)/(Ns));

figure(18); 
subplot(1,2,1)
plot(Oeil_a);
title("3eme chaine - Diagramme de l'oeil (a) avec canal de transmission")
xlabel('t (s)')
ylabel('Amplitude (V)');

subplot(1,2,2)
plot(Oeil_b);
title("3eme chaine - Diagramme de l'oeil (b) avec canal de transmission")
xlabel('t (s)')
ylabel('Amplitude (V)');


%% 4eme chaine

Nb=100;     % Nombre de bits emis
Fe = 12000; % Fréquence d'échantillonnage en Hz
Te = 1/Fe;
Rs = 3000; % rythme symbole (en symboles par seconde)
Ts = 1/Rs;
Ns = Ts/Te; % Nombre d'echantillons par symbole
M=2;        % Taille de la constellation


h = ones(1,Ns);
h_ref = ones(1,Ns); % reponse impulsionnelle du filtre de mise en forme
t0 = Ns;
hr = ones(1,Ns);   % reponse impulsionelle du filtre de reception
hr_ref = ones(1,Ns);

% Emission chaine étudiée

bits = randi([0 1],1,Nb);
symboles=(2*bi2de(reshape(bits,2,length(bits)/2).')-3).';
suite_diracs = kron(symboles, [1 zeros(1,Ns-1)]);

% Passage dans le filtre de mise en forme
x = filter(h,1,suite_diracs);


%Affichage du signal généré
figure(19); plot(x)
title('4eme chaine - Signal généré')
xlabel('t (s)')
ylabel('Amplitude (V)');

%Calcul de la Dsp du signal par periodogramme
Dsp_x =(1/length(x))*abs(fft(x,2^nextpow2(length(x)))).^2;

%Affichage de la DSP du signal généré

frequence= linspace (-Fe/2,Fe/2,length(Dsp_x));

figure(20);  
semilogy(frequence,fftshift(Dsp_x))
title('4eme chaine - DSP du Signal généré')
xlabel('f (hz)')
ylabel('Amplitude');


%% Comparaison de la DSP avec celle de la chaîne de référence

% Suréchantillonage chaine étudiée
symboles_double = reshape([symboles ; symboles],1,2*length(symboles));
suite_diracs_double=kron(symboles_double, [1 zeros(1,Ns-1)]);
x_double = filter(h, 1, suite_diracs_double);
Dsp_x_double=(1/length(x_double))*abs(fft(x_double,2^nextpow2(length(x_double)))).^2;

% Emission chaine de référence
symboles_ref = 2*bits-1;
suite_diracs_ref = kron(symboles_ref,[1 zeros(1,Ns-1)]);
x_ref = filter(h_ref,1,suite_diracs_ref);
Dsp_x_ref=(1/length(x_ref))*abs(fft(x_ref,2^nextpow2(length(x_ref)))).^2;


% Tracé des DSP
frequence= linspace(-Fe/2,Fe/2,length(Dsp_x_double));
frequence_ref = linspace(-Fe/2,Fe/2,length(Dsp_x_ref));
figure(21)
semilogy(frequence,fftshift(Dsp_x_double),'r-','LineWidth',1)
hold on
semilogy(frequence_ref,fftshift(Dsp_x_ref),'b-','LineWidth',1)
xlabel('f (Hz)')
ylabel('DSP')
legend('DSP chaine étudiée','DSP chaine de référence')
title("4eme chaine - comparaison des DSP")
grid
set(gca,'FontSize',12)


% Reception chaine étudiée
z=filter(hr,1,x);

figure(22); plot(fftshift(z))
title('4eme chaine - Signal en sortie de réception')
xlabel('t (s)')
ylabel('Amplitude (V)')

%Diagramme de l'oeil
Oeil = reshape(z,Ns, length(z)/Ns);

figure(23); plot(Oeil);
title("4eme chaine - Diagramme de l'oeil")
xlabel('t (s)')
ylabel('Amplitude (V)');


seuil= 0;
seuil_sup = 2*Ns;
seuil_inf = -2*Ns;
ze = zeros(1,Nb/2);
for i = 1:Nb/2
    if z(Ns*i) >= seuil
        if z(Ns*i) >= seuil_sup
            ze(i) = 3;
        else
            ze(i) = 1;  
        end
    else
        if z(Ns*i) <= seuil_inf
            ze(i) = -3;
        else
            ze(i) = -1;  
        end
    end
end

% Demapping
bits_r = reshape(de2bi((ze + 3)/2).',1,length(bits));

% Calcul du nb d'erreurs
nb_erreurs = 0;
for i = 1:length(bits)
    if bits_r(i) ~= bits(i)
        nb_erreurs = nb_erreurs + 1;
    end
end
TEB_4 = nb_erreurs/Nb


fprintf('4eme chaine TEB = %d\n', TEB_4)

%% 4eme chaine avec bruit

%Initialisation des constantes
Fe = 12000; % Fréquence d'échantillonnage en Hz
Te = 1/Fe;
Rs = 3000; % rythme symbole (en symboles par seconde)
Ts = 1/Rs;
Ns = Ts/Te;        % Taille de la constellation
h=ones(1,Ns);  % reponse impulsionnelle du filtre de mise en forme
t0=Ns;
hr=ones(1,Ns);   % reponse impulsionelle du filtre de reception
Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);

TES4=zeros(1,length(Eb_sur_N0_dB)); % Initialisation du TEB que l'on veut simuler
NerrLim=100;     % Le nombre d'erreur à partir du quel on s'arrête
Nerrcomptees=zeros(1,length(Eb_sur_N0_dB)); % Le nombre d'erreur que l'on compte
Nsimu=zeros(1,length(Eb_sur_N0_dB));    % Le nombre de simulation
 


for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

        % Emission
        bits = randi([0 1],1,Nb);
        symboles=(2*bi2de(reshape(bits,2,length(bits)/2).')-3).';
        suite_diracs = kron(symboles, [1 zeros(1,Ns-1)]);
        x = filter(h,1, suite_diracs);
       

        % Canal

        Pr=mean(abs(x).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        r=x+sqrt(sigma)*randn(1,length(x));

        % Reception
        z=filter(hr,1,r);
        seuil= 0;
        seuil_sup = 2*Ns;
        seuil_inf = -2*Ns;
        ze = zeros(1,Nb/2);
        for i = 1:Nb/2
            if z(Ns*i) >= seuil
                if z(Ns*i) >= seuil_sup
                    ze(i) = 3;
                else
                    ze(i) = 1;  
                end
            else
                if z(Ns*i) <= seuil_inf
                    ze(i) = -3;
                else
                    ze(i) = -1;  
                end
            end
        end
        
        Nerrcomptees(k)=Nerrcomptees(k)+sum(symboles~=ze);
        Nsimu(k)=Nsimu(k)+1;
        
    end

    TES4(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % Calcul du TEB
    
end


% Tracé de la comparaison entre le TES simulé et le TES théorique
Pb=qfunc((3/2)*sqrt((5/4)*Eb_sur_N0_lin));

figure(24)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TES4,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TES')
legend('Coube théorique','TES simulé')
title("4eme chaine - comparaison entre le TES simulé et le TES théorique")

grid
set(gca,'FontSize',12)

% Réinitialisation des constantes
Fe = 12000; % Fréquence d'échantillonnage en Hz
Te = 1/Fe;
Rs = 3000; % rythme symbole (en symboles par seconde)
Ts = 1/Rs;
Ns = Ts/Te;        % Taille de la constellation

h=ones(1,Ns);  % reponse impulsionnelle du filtre de mise en forme
t0=Ns;
hr=ones(1,Ns);   % reponse impulsionelle du filtre de reception
Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);

TEB4=zeros(1,length(Eb_sur_N0_dB)); % Initialisation du TEB que l'on veut simuler
NerrLim=100;    % Le nombre d'erreur à partir du quel on s'arrête
Nerrcomptees=zeros(1,length(Eb_sur_N0_dB)); % Le nombre d'erreur que l'on compte
Nsimu=zeros(1,length(Eb_sur_N0_dB));    % Le nombre de simulation


for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

        % Emission
        bits = randi([0 1],1,Nb);
        symboles=(2*bi2de(reshape(bits,2,length(bits)/2).')-3).';
        suite_diracs = kron(symboles, [1 zeros(1,Ns-1)]);
        x = filter(h,1, suite_diracs);
       

        % Canal

        Pr=mean(abs(x).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        r=x+sqrt(sigma)*randn(1,length(x));

        % Reception
        z=filter(hr,1,r);
        seuil= 0;
        seuil_sup = 2*Ns;
        seuil_inf = -2*Ns;
        ze = zeros(1,Nb/2);
        for i = 1:Nb/2
            if z(Ns*i) >= seuil
                if z(Ns*i) >= seuil_sup
                    ze(i) = 3;
                else
                    ze(i) = 1;  
                end
            else
                if z(Ns*i) <= seuil_inf
                    ze(i) = -3;
                else
                    ze(i) = -1;  
                end
            end
        end
        
        
        bits_r = reshape(de2bi((ze + 3)/2).',1,length(bits));
        
        Nerrcomptees(k)=Nerrcomptees(k)+sum(bits~=bits_r);
        Nsimu(k)=Nsimu(k)+1;
        
    end

    TEB4(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % Calcul du TEB

end

% Tracé de la comparaison entre le TEB simulé et le TEB théorique
Pb=qfunc((3/2)*sqrt((5/4)*Eb_sur_N0_lin))/log2(Ns);

figure(25)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB4,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique','TEB simulé')
title("4eme chaine - comparaison entre le TEB simulé et le TEB théorique")
grid
set(gca,'FontSize',12)