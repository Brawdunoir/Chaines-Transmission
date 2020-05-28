clc
clear
close all
%% Initialisation 
Nb = 1500;
fe = 10000;
Te = 1/fe;
fp = 2000;
Tp = 1/fp;
Rs = 1000;
Ts = 1/Rs;
Ns = Ts/Te;

M = 8;
t = 0:1/fe:Nb/(3*Rs)-1/fe;

%Génération du signal
bits = randi([0 1],1,Nb);

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
DSP = pwelch(y);

%Affichage de la DSP du signal généré
frequence= linspace(-fe/2,fe/2,length(DSP));
figure;  
semilogy(frequence,fftshift(DSP))
title('chaine - DSP du signal transmit sur fréquence porteuse')

% Passage dans le filtre de réception
z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

%Echantillonage 
z_r = real(z);
z_i = imag(z);

for i = 0:(Nb/3)-1 
    dk(i+1) = z_r(1+Ns*i) + j*z_i(1+Ns*i);
end


for i = 0:(Nb/3)-1
   if angle(dk(i+1)) <= pi/8 && angle(dk(i+1)) >= -pi/8
        dk(i+1) = 1.0 + j*0.0;
    elseif angle(dk(i+1)) <= 3*pi/8 && angle(dk(i+1)) >= pi/8
        dk(i+1) = sqrt(2)/2 + j*sqrt(2)/2;
    elseif angle(dk(i+1)) <= 5*pi/8 && angle(dk(i+1)) >= 3*pi/8
        dk(i+1) = 0.0 + j;
    elseif angle(dk(i+1)) <= 7*pi/8 && angle(dk(i+1)) >= 5*pi/8
        dk(i+1) = -sqrt(2)/2 + j*sqrt(2)/2;
    elseif angle(dk(i+1)) >= 7*pi/8 || angle(dk(i+1)) <= -7*pi/8
        dk(i+1) = -1.0 + j*0.0;
    elseif angle(dk(i+1)) <= -5*pi/8 && angle(dk(i+1)) >= -7*pi/8
        dk(i+1) = -sqrt(2)/2 - j*sqrt(2)/2;
    elseif angle(dk(i+1)) <= -3*pi/8 && angle(dk(i+1)) >= -5*pi/8
        dk(i+1) = 0.0 - j;
    elseif angle(dk(i+1)) <= -pi/8 && angle(dk(i+1)) >= -3*pi/8
        dk(i+1) = sqrt(2)/2 - j*sqrt(2)/2;
    end
end


% Demapping 
bits_recup = [];
for i = 1:Nb/3
   if dk(i) == j
       bits_recup = [bits_recup 1 0 1];
   elseif dk(i) == -j
       bits_recup = [bits_recup 1 1 0];
   elseif dk(i) == -sqrt(2)/2 - j*sqrt(2)/2
       bits_recup = [bits_recup 0 1 0];
   elseif dk(i) == -1
       bits_recup = [bits_recup 0 1 1];
   elseif dk(i) == sqrt(2)/2 - j*sqrt(2)/2
       bits_recup = [bits_recup 1 0 0];
   elseif dk(i) == -sqrt(2)/2 + j*sqrt(2)/2
       bits_recup = [bits_recup 1 1 1];
   elseif dk(i) == sqrt(2)/2 + j*sqrt(2)/2
       bits_recup = [bits_recup 0 0 1];
   elseif dk(i) == 1
       bits_recup = [bits_recup 0 0 0];    
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
        x=pskmod(x_decimaux,M,0,'gray');
        x = x';

        %Génération de la suite de dirac
        suite_diracs=kron(x, [1 zeros(1,Ns-1)]);
        
       
        % Passage dans le filtre de mise en forme
        y = filter(h,1,[suite_diracs zeros(1,retard)]);
        y = y(retard+1:end);
        
        % Canal
        Pr=mean(abs(y).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        ez = randn(1,Nb*Ns/2);
        y=y+sqrt(sigma)*randn(1,Nb*Ns/3)+ j*sqrt(sigma)*randn(1,Nb*Ns/3);

        
       
        % Passage dans le filtre de réception
        z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
        z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

        %Echantillonage 
        z_r = real(z);
        z_i = imag(z);

        for i = 0:(Nb/3)-1 
        dk(i+1) = z_r(1+Ns*i) + j*z_i(1+Ns*i);
            if k == 3
                ak(1,i+1) = z_r(1+Ns*i);
                bk(1,i+1) = z_i(1+Ns*i);
            elseif k == 5
                ak(2,i+1) = z_r(1+Ns*i);
                bk(2,i+1) = z_i(1+Ns*i);
            elseif k == 7
                ak(3,i+1) = z_r(1+Ns*i);
                bk(3,i+1) = z_i(1+Ns*i);
            end
        end


        for i = 0:(Nb/3)-1
           if angle(dk(i+1)) <= pi/8 && angle(dk(i+1)) >= -pi/8
                dk(i+1) = 1.0 + j*0.0;
            elseif angle(dk(i+1)) <= 3*pi/8 && angle(dk(i+1)) >= pi/8
                dk(i+1) = sqrt(2)/2 + j*sqrt(2)/2;
            elseif angle(dk(i+1)) <= 5*pi/8 && angle(dk(i+1)) >= 3*pi/8
                dk(i+1) = 0.0 + j;
            elseif angle(dk(i+1)) <= 7*pi/8 && angle(dk(i+1)) >= 5*pi/8
                dk(i+1) = -sqrt(2)/2 + j*sqrt(2)/2;
            elseif angle(dk(i+1)) >= 7*pi/8 || angle(dk(i+1)) <= -7*pi/8
                dk(i+1) = -1.0 + j*0.0;
            elseif angle(dk(i+1)) <= -5*pi/8 && angle(dk(i+1)) >= -7*pi/8
                dk(i+1) = -sqrt(2)/2 - j*sqrt(2)/2;
            elseif angle(dk(i+1)) <= -3*pi/8 && angle(dk(i+1)) >= -5*pi/8
                dk(i+1) = 0.0 - j;
            elseif angle(dk(i+1)) <= -pi/8 && angle(dk(i+1)) >= -3*pi/8
                dk(i+1) = sqrt(2)/2 - j*sqrt(2)/2;
            end
        end


        % Demapping 
        bits_recup = [];
        for i = 1:Nb/3
           if dk(i) == j
               bits_recup = [bits_recup 1 0 1];
           elseif dk(i) == -j
               bits_recup = [bits_recup 1 1 0];
           elseif dk(i) == -sqrt(2)/2 - j*sqrt(2)/2
               bits_recup = [bits_recup 0 1 0];
           elseif dk(i) == -1
               bits_recup = [bits_recup 0 1 1];
           elseif dk(i) == sqrt(2)/2 - j*sqrt(2)/2
               bits_recup = [bits_recup 1 0 0];
           elseif dk(i) == -sqrt(2)/2 + j*sqrt(2)/2
               bits_recup = [bits_recup 1 1 1];
           elseif dk(i) == sqrt(2)/2 + j*sqrt(2)/2
               bits_recup = [bits_recup 0 0 1];
           elseif dk(i) == 1
               bits_recup = [bits_recup 0 0 0];    
           end
        end


        
        
        Nerrcomptees(k)=Nerrcomptees(k)+sum(bits~=bits_recup);
        Nsimu(k)=Nsimu(k)+1;
        
        
    end

    TEB(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end


% % Tracé des constellations
figure(2);
plot(ak(1,:), bk(1,:),'*');
hold on;
plot(ak(2,:), bk(2,:), '*');
hold on;
plot(ak(3,:), bk(3,:), '*');
hold on;
plot([0 -0.7071 0.7071 1 0 0.7071 -0.7071 -1], [-1 0.7071 0.7071 0 1 -0.7071 -0.7071 0], 'o');
title("Constellations  en sortie du mapping et de l'échantillonnage pour différentes valeurs de Eb/N0 (8-PSK)");
xlabel('ak');
ylabel('bk');
legend("Eb/N0 = 2", "Eb/N0 = 4", "Eb/N0 = 6", 'En sortie du mapping');

% Tracé de la comparaison entre le TEB simulé et le TEB théorique
Pb=2*qfunc(sqrt(2*Eb_sur_N0_lin));


figure(3)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique','TEB simulé')
title("Chaine 8PSK - comparaison entre le TEB simulé et le TEB théorique")
grid
set(gca,'FontSize',12)
