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

M = 16;
t = 0:1/fe:Nb/(4*Rs)-1/fe;

%Génération du signal
bits = randi([0 1],1,Nb);

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
DSP = pwelch(y);

%Affichage de la DSP du signal généré
frequence= linspace(-fe/2,fe/2,length(DSP));
figure(1);  
semilogy(frequence,fftshift(DSP))
title('chaine - DSP du signal transmit sur fréquence porteuse')

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


for i = 0:(Nb/4)-1
    if z_r(1+Ns*i) >= 2 % Echantillonnage de Re(z)
        ak(1+i) = 3;
    elseif z_r(1+Ns*i) >= 0
        ak(1+i) = 1;
    elseif z_r(1+Ns*i) <= -2
        ak(1+i) = -3;
    else
        ak(1+i) = -1;
    end
    if z_i(1+Ns*i) >= 2 % Echantillonnage de Im(z)
        bk(1+i) = 3;
    elseif z_i(1+Ns*i) >= 0
        bk(1+i) = 1;
    elseif z_i(1+Ns*i) <= -2
        bk(1+i) = -3;
    else
        bk(1+i) = -1;
    end
    dk(1+i) = ak(1+i) + j*bk(1+i);  % dk = ak + j*bk
end

% Demapping 
bits_recup = [];
for i = 1:Nb/4
   if dk(i) == -3-j*3
       bits_recup = [bits_recup 0 0 0 0];
   elseif dk(i) == -3-j
       bits_recup = [bits_recup 1 0 0 0];
   elseif dk(i) == -3+j
       bits_recup = [bits_recup 1 1 0 0];
   elseif dk(i) == -3+j*3
       bits_recup = [bits_recup 0 1 0 0];
   elseif dk(i) == -1-j*3
       bits_recup = [bits_recup 0 0 1 0];
   elseif dk(i) == -1-j
       bits_recup = [bits_recup 1 0 1 0];
   elseif dk(i) == -1+j
       bits_recup = [bits_recup 1 1 1 0];
   elseif dk(i) == -1+j*3
       bits_recup = [bits_recup 0 1 1 0];
   elseif dk(i) == 1-3*j
       bits_recup = [bits_recup 0 0 1 1];
   elseif dk(i) == 1-j
       bits_recup = [bits_recup 1 0 1 1];
   elseif dk(i) == 1+j
       bits_recup = [bits_recup 1 1 1 1];
    elseif dk(i) == 1+j*3
        bits_recup = [bits_recup 0 1 1 1];
    elseif dk(i) == 3-j*3
        bits_recup = [bits_recup 0 0 0 1];
   elseif dk(i) == 3-j
       bits_recup = [bits_recup 1 0 0 1];
   elseif dk(i) == 3+j
       bits_recup = [bits_recup 1 1 0 1];
   elseif dk(i) == 3+j*3
       bits_recup = [bits_recup 0 1 0 1];
       
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
Constellations = zeros(3,Nb/2);

for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

       %Génération du signal
        bits = randi([0 1],1,Nb);

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
        
        
        % Passage dans le filtre de mise en forme
        y = filter(h,1,[suite_diracs zeros(1,retard)]);
        y = y(retard+1:end);
        
        % Canal
        Pr=mean(abs(y).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));

        y=y+sqrt(sigma)*randn(1,Nb*Ns/4)+ j*sqrt(sigma)*randn(1,Nb*Ns/4);

       
        % Passage dans le filtre de réception
        z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
        z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

        z_r = real(z);
        z_i = imag(z);
        z_r(1) = z_r(1)+sign(z_r(1));
        z_i(1) = z_i(1)+sign(z_i(1));

        ak = zeros(1,Nb/4);
        bk = zeros(1,Nb/4);
        dk = zeros(1,Nb/4);


        for i = 0:(Nb/4)-1
            if z_r(1+Ns*i) >= 2 % Echantillonnage de Re(z)
                ak(1+i) = 3;
            elseif z_r(1+Ns*i) >= 0
                ak(1+i) = 1;
            elseif z_r(1+Ns*i) <= -2
                ak(1+i) = -3;
            else
                ak(1+i) = -1;
            end
            if z_i(1+Ns*i) >= 2 % Echantillonnage de Im(z)
                bk(1+i) = 3;
            elseif z_i(1+Ns*i) >= 0
                bk(1+i) = 1;
            elseif z_i(1+Ns*i) <= -2
                bk(1+i) = -3;
            else
                bk(1+i) = -1;
            end
            dk(1+i) = ak(1+i) + j*bk(1+i);  % dk = ak + j*bk
            
            if k == 3
                ak2(1,i+1) = z_r(1+Ns*i);
                bk2(1,i+1) = z_i(1+Ns*i);
            elseif k == 5
                ak2(2,i+1) = z_r(1+Ns*i);
                bk2(2,i+1) = z_i(1+Ns*i);
            elseif k == 7
                ak2(3,i+1) = z_r(1+Ns*i);
                bk2(3,i+1) = z_i(1+Ns*i);
            end
        end

        % Demapping 
        bits_recup = [];
        for i = 1:Nb/4
           if dk(i) == -3-j*3
               bits_recup = [bits_recup 0 0 0 0];
           elseif dk(i) == -3-j
               bits_recup = [bits_recup 1 0 0 0];
           elseif dk(i) == -3+j
               bits_recup = [bits_recup 1 1 0 0];
           elseif dk(i) == -3+j*3
               bits_recup = [bits_recup 0 1 0 0];
           elseif dk(i) == -1-j*3
               bits_recup = [bits_recup 0 0 1 0];
           elseif dk(i) == -1-j
               bits_recup = [bits_recup 1 0 1 0];
           elseif dk(i) == -1+j
               bits_recup = [bits_recup 1 1 1 0];
           elseif dk(i) == -1+j*3
               bits_recup = [bits_recup 0 1 1 0];
           elseif dk(i) == 1-3*j
               bits_recup = [bits_recup 0 0 1 1];
           elseif dk(i) == 1-j
               bits_recup = [bits_recup 1 0 1 1];
           elseif dk(i) == 1+j
               bits_recup = [bits_recup 1 1 1 1];
            elseif dk(i) == 1+j*3
                bits_recup = [bits_recup 0 1 1 1];
            elseif dk(i) == 3-j*3
                bits_recup = [bits_recup 0 0 0 1];
           elseif dk(i) == 3-j
               bits_recup = [bits_recup 1 0 0 1];
           elseif dk(i) == 3+j
               bits_recup = [bits_recup 1 1 0 1];
           elseif dk(i) == 3+j*3
               bits_recup = [bits_recup 0 1 0 1];

           end
        end


        
        
        Nerrcomptees(k)=Nerrcomptees(k)+sum(bits~=bits_recup);
        Nsimu(k)=Nsimu(k)+1;
        
        
    end

    TEB(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end


% % Tracé des constellations
figure(2);
plot(ak2(1,:), bk2(1,:),'*');
hold on;
plot(ak2(2,:), bk2(2,:), '*');
hold on;
plot(ak2(3,:), bk2(3,:), '*');
hold on;
plot([3 3 3 3 1 1 1 1 -1 -1 -1 -1 -3 -3 -3 -3], [3 1 -1 -3 3 1 -1 -3 3 1 -1 -3 3 1 -1 -3],'o');
title("Constellations  en sortie du mapping et de l'échantillonnage pour différentes valeurs de Eb/N0 (8-PSK)");
xlabel('ak');
ylabel('bk');
legend("Eb/N0 = 2", "Eb/N0 = 4", "Eb/N0 = 6", 'En sortie du mapping');



% Tracé de la comparaison entre le TEB simulé et le TEB théorique

Pb=(3/4)*qfunc(sqrt((4/5)*Eb_sur_N0_lin));

figure(3)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique','TEB simulé')
title("chaine 16QAM - comparaison entre le TEB simulé et le TEB théorique")
grid
set(gca,'FontSize',12)