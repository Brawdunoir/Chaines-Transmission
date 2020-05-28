clc
clear
close all
%% Chaine 4ASK
Nb = 3000;
fe = 10000;
Te = 1/fe;
fp = 2000;
Tp = 1/fp;
Rs = 1000;
Ts = 1/Rs;
Ns = Ts/Te;

M = 4;
%t = linspace(2*Ts, 2*Ts, Nb/2);
t = 0:1/fe:Nb/(2*Rs)-1/fe;

%Génération du signal
bits = randi([0 1],1,Nb);
x=[];
for i=1:2:Nb
    if bits(i)==0 && bits(i+1)==0
        x=[x 3];
    elseif bits(i)==0 && bits(i+1)==1
        x=[x 1];
    elseif bits(i)==1 && bits(i+1)==0
        x=[x -1];
    elseif bits(i)==1 && bits(i+1)==1
        x=[x -3];
    end
end


% %Génération de la suite de dirac
suite_diracs=kron(x, [1 zeros(1,Ns-1)]);

%Filtre de mise en forme et de réception
alpha = 0.5;
h = rcosdesign(alpha,Nb/2,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,Nb/2,Ns);  % reponse impulsionnelle du filtre de réception
retard = Nb*Ns/4;

% Passage dans le filtre de mise en forme
y = filter(h,1,[suite_diracs zeros(1,retard)]);
y = y(retard+1:end);


%Calcul de la Dsp du signal par periodogramme
DSP4ASK = pwelch(y);
frequence4ASK = linspace(-fe/2,fe/2,length(DSP4ASK));
% Passage dans le filtre de réception
z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

z_r = real(z);
z_r(1) = z_r(1)+sign(z_r(1));

%% Chaine 8PSK
Nb = 3000;
fe = 10000;
Te = 1/fe;
fp = 2000;
Tp = 1/fp;
Rs = 1000;
Ts = 1/Rs;
Ns = Ts/Te;

M = 8;
%t = linspace(2*Ts, 2*Ts, Nb/2);
t = 0:1/fe:Nb/(3*Rs)-1/fe;

%Génération du signal

n=log2(M);
% groupement des bits par n
blocs_bits=reshape(bits,n,length(bits)/n).';
% conversion des groupes de bits en une valeur décimale (de 0 à M-1)
x_decimaux=bi2de(blocs_bits);
% conversion des valeurs décimales en symboles complexes 8-PSK
x=pskmod(x_decimaux,M,0,'gray');
x = x';


%Génération de la suite de dirac
suite_diracs=kron(x, [1 zeros(1,Ns-1)]);

%Filtre de mise en forme et de réception
alpha = 0.5;
Nb3 = Nb/3;
h = rcosdesign(alpha,Nb3,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,Nb3,Ns);  % reponse impulsionnelle du filtre de réception
retard = Nb*Ns/6;

% Passage dans le filtre de mise en forme
y = filter(h,1,[suite_diracs zeros(1,retard)]);
y = y(retard+1:end);


%Calcul de la Dsp du signal par periodogramme
DSP8PSK = pwelch(y);
frequence8PSK = linspace(-fe/2,fe/2,length(DSP8PSK));
% Passage dans le filtre de réception
z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

%Echantillonage 
z_r = real(z);
z_i = imag(z);


%% Chaine QPSK
Nb = 3000;
fe = 10000;
Te = 1/fe;
fp = 2000;
Tp = 1/fp;
Rs = 1000;
Ts = 1/Rs;
Ns = Ts/Te;

M = 4;
%t = linspace(2*Ts, 2*Ts, Nb/2);
t = 0:1/fe:Nb/(2*Rs)-1/fe;

%Génération du signal

x=[];
for i=1:2:Nb
    if bits(i)==0 && bits(i+1)==0
        x=[x 1+j];
    elseif bits(i)==0 && bits(i+1)==1
        x=[x 1-j];
    elseif bits(i)==1 && bits(i+1)==0
        x=[x -1+j];
    elseif bits(i)==1 && bits(i+1)==1
        x=[x -1-j];
    end
end


%Génération de la suite de dirac
suite_diracs=kron(x, [1 zeros(1,Ns-1)]);


%Filtre de mise en forme et de réception
alpha = 0.5;
h = rcosdesign(alpha,Nb/2,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,Nb/2,Ns);  % reponse impulsionnelle du filtre de réception
retard = Nb*Ns/4;

% Passage dans le filtre de mise en forme
y = filter(h,1,[suite_diracs zeros(1,retard)]);
y = y(retard+1:end);

%Calcul de la Dsp du signal par periodogramme
DSPQPSK = pwelch(y);
frequenceQPSK = linspace(-fe/2,fe/2,length(DSPQPSK));
% Passage dans le filtre de réception
z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)


%Echantillonage 
z_r = real(z);
z_i = imag(z);
bits_recup = zeros(1,Nb);
ak = zeros(1,Nb/2);
bk = zeros(1,Nb/2);
dk = zeros(1,Nb/2);

%% Chaine 16QAM
Nb = 3000;
fe = 10000;
Te = 1/fe;
fp = 2000;
Tp = 1/fp;
Rs = 1000;
Ts = 1/Rs;
Ns = Ts/Te;

M = 16;
%t = linspace(2*Ts, 2*Ts, Nb/2);
t = 0:1/fe:Nb/(4*Rs)-1/fe;

%Génération du signal


n=log2(M);
% groupement des bits par n
blocs_bits=reshape(bits,n,length(bits)/n).';
% conversion des groupes de bits en une valeur décimale (de 0 à M-1)
x_decimaux=bi2de(blocs_bits);
% conversion des valeurs décimales en symboles complexes 8-PSK
x=qammod(x_decimaux,M,'gray');
x = x';


%Génération de la suite de dirac
suite_diracs=kron(x, [1 zeros(1,Ns-1)]);

%Filtre de mise en forme et de réception
alpha = 0.5;
Nb4 = Nb/4;
h = rcosdesign(alpha,Nb4,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,Nb4,Ns);  % reponse impulsionnelle du filtre de réception
retard = Nb*Ns/8;

% Passage dans le filtre de mise en forme
y = filter(h,1,[suite_diracs zeros(1,retard)]);
y = y(retard+1:end);


%Calcul de la Dsp du signal par periodogramme
DSP16QAM = pwelch(y);
frequence16QAM = linspace(-fe/2,fe/2,length(DSP16QAM));
% Passage dans le filtre de réception
z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

%Echantillonage 
z_r = real(z);
z_i = imag(z);
z_r(1) = z_r(1)+sign(z_r(1));
z_i(1) = z_i(1)+sign(z_i(1));

ak = zeros(1,Nb/4);
bk = zeros(1,Nb/4);
dk = zeros(1,Nb/4);

%%Comparaisons
figure(3)
semilogy(frequence4ASK,fftshift(DSP4ASK),'r-','LineWidth',1)
hold on
semilogy(frequence8PSK,fftshift(DSP8PSK),'y-','LineWidth',1)
hold on
semilogy(frequenceQPSK,fftshift(DSPQPSK),'g-','LineWidth',1)
hold on
semilogy(frequence16QAM,fftshift(DSP16QAM),'b-','LineWidth',1)
xlabel('frequence')
ylabel('DSP')
legend('DSP chaine 4ASK', 'DSP chaine 8PSK', 'DSP chaine QPSK', 'DSP chaine 16QAM' )
%title("Comparaison entre le TEB avec passe-bas équivalent et le TEB de la chaine sur fréquence porteuse")
grid
set(gca,'FontSize',12)
