clc
clear
close all
%% Initialisation 
Nb = 1000;
fe = 10000;
Te = 1/fe;
fp = 2000;
Tp = 1/fp;
Rs = 1000;
Ts = 1/Rs;
Ns = Ts/Te;

M = 4;

t = 0:1/fe:Nb/(2*Rs)-1/fe;

%%Génération du signal
bits = randi([0 1],1,Nb);
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
alpha = 0,35;
h = rcosdesign(alpha,10,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,10,Ns);  % reponse impulsionnelle du filtre de réception
retard = 10*Ns/2;

%Affichage di signal généré
figure(1); 
subplot(1,2,1)
plot(real(x))
title('chaine porteuse - partie réelle signal généré')
xlabel('t (s)')
ylabel('Amplitude (V)');

subplot(1,2,2)
plot(imag(x))
title('chaine porteuse - partie imaginaire signal généré')
xlabel('t (s)')
ylabel('Amplitude (V)');

% Passage dans le filtre de mise en forme
y = filter(h,1,[suite_diracs zeros(1,retard)]);
y = y(retard+1:end);

% Transposition sur fréquence porteuse
y = y.*exp(j*pi*fp*t);

y = real(y);

%Affichage du signal transmit surfréquence porteuse

figure(2)
plot(y)
xlabel('t (s)')
ylabel('Amplitude (V)');

%Calcul de la Dsp du signal par periodogramme
DSP = pwelch(y);

%Affichage de la DSP du signal généré
frequence= linspace(0,fe,length(DSP));
figure(3);  
semilogy(frequence,DSP)
xlabel('frequence')
ylabel('DSP')
title('chaine - DSP du signal transmit sur fréquence porteuse')

% Retour en bande de base
y_r = y.*cos(pi*fp*t);
y_i = y.*sin(pi*fp*t);

y = y_r -j*y_i;

% Passage dans le filtre de réception
z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

%Echantillonage 

z_r = real(z);
z_i = imag(z);
bits_recup = zeros(1,Nb);
ak = zeros(1,Nb);
bk = zeros(1,Nb);
dk = zeros(1,Nb);

for i = 0:(Nb/2)-1
    if z_r(1+Ns*i) >= 0  % Echantillonnage de Re(z)
        ak(1+i) = 1;
    else
        ak(1+i) = -1;
    end
    if z_i(1+Ns*i) >= 0  % Echantillonnage de Im(z)
        bk(1+i) = 1;
    else
        bk(1+i) = -1;
    end
    dk(1+i) = ak(1+i) + j*bk(1+i);  % dk = ak + j*bk
end

% Demapping 
for i = 1:Nb/2
    if dk(i) == -1+j
        bits_recup(2*i-1) = 1;
        bits_recup(2*i) = 0;
    elseif dk(i) == 1-j
        bits_recup(2*i-1) = 0;
        bits_recup(2*i) = 1; 
    elseif dk(i) == 1+j
        bits_recup(2*i-1) = 0;
        bits_recup(2*i) = 0;          
    elseif dk(i) == -1-j
        bits_recup(2*i-1) = 1;
        bits_recup(2*i) = 1;         
    end
end



TEB = sum(abs(bits_recup-bits))/Nb;
disp("TEB = " + TEB);


%% Implantation de la chaine avec buit

Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);
TEB3=zeros(1,length(Eb_sur_N0_dB));

NerrLim=100;
Nerrcomptees=zeros(1,length(Eb_sur_N0_dB));
Nsimu=zeros(1,length(Eb_sur_N0_dB));


for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

        % Emission
        bits = randi([0 1],1,Nb);
        x = [];
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


        suite_diracs = kron(x,[1 zeros(1,Ns-1)]);
        y = filter(h,1,[suite_diracs zeros(1,retard)]);
        y = y(retard+1:end);

        % Transposition sur fréquence porteuse
        y = y.*exp(2*j*pi*fp*t);
        y = real(y);
        
        % Canal

        Pr=mean(abs(y).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        y=y+sqrt(sigma)*randn(1,Nb*Ns/2);

        % Retour en bande de base
        y_r = y.*cos(2*pi*fp*t);
        y_i = y.*sin(2*pi*fp*t);

        y = y_r -j*y_i;

        % Passage dans le filtre de réception
        z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
        z = z(retard+1:end);
        z_r = real(z);
        z_i = imag(z);
        bits_recup = zeros(1,Nb);
        ak = zeros(1,Nb);
        bk = zeros(1,Nb);
        dk = zeros(1,Nb);

        for i = 0:(Nb/2)-1
            if z_r(1+Ns*i) >= 0  % Echantillonnage de Re(z)
                ak(1+i) = 1;
            else
                ak(1+i) = -1;
            end
            if z_i(1+Ns*i) >= 0  % Echantillonnage de Im(z)
                bk(1+i) = 1;
            else
                bk(1+i) = -1;
            end
            dk(1+i) = ak(1+i) + j*bk(1+i);  % dk = ak + j*bk
        end

        % Demapping 
        for i = 1:Nb/2
            if dk(i) == -1+j
                bits_recup(2*i-1) = 1;
                bits_recup(2*i) = 0;
            elseif dk(i) == 1-j
                bits_recup(2*i-1) = 0;
                bits_recup(2*i) = 1; 
            elseif dk(i) == 1+j
                bits_recup(2*i-1) = 0;
                bits_recup(2*i) = 0;          
            elseif dk(i) == -1-j
                bits_recup(2*i-1) = 1;
                bits_recup(2*i) = 1;         
            end
        end
        
        
        Nerrcomptees(k)=Nerrcomptees(k)+sum(bits~=bits_recup);
        Nsimu(k)=Nsimu(k)+1;
        
        
    end

    TEB(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end


% Tracé de la comparaison entre le TEB simulé et le TEB théorique
Pb=qfunc(sqrt(2*Eb_sur_N0_lin));

figure(4)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique','TEB simulé')
title("3eme chaine - comparaison entre le TEB simulé et le TEB théorique")
grid
set(gca,'FontSize',12)