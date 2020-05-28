clc
clear
close all
%% Initialisation 
Nb = 500;
fe = 10000;
Te = 1/fe;
fp = 2000;
Tp = 1/fp;
Rs = 1000;
Ts = 1/Rs;
Ns = Ts/Te;

M = 4;
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
DSP = pwelch(y);

%Affichage de la DSP du signal généré
frequence= linspace(-fe/2,fe/2,length(DSP));
figure(1);  
semilogy(frequence,fftshift(DSP))
title('chaine - DSP du signal transmit sur fréquence porteuse')

% Passage dans le filtre de réception
z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

z_r = real(z);
z_r(1) = z_r(1)+sign(z_r(1));

%Echantillonage 
bits_recup = zeros(1,Nb);
dk = zeros(1,Nb/2);

for i = 0:(Nb/2)-1
    if z_r(1+Ns*i) >=2  
        dk(1+i) = 3;
    elseif z_r(1+Ns*i) >= 0
        dk(1+i) = 1;
    
    elseif z_r(1+Ns*i) >= -2  
        dk(1+i) = -1;
    else
        dk(1+i) = -3;
    end
 
end

% Demapping 
for i = 1:Nb/2
    if dk(i) == 3
        bits_recup(2*i-1) = 0;
        bits_recup(2*i) = 0;
    elseif dk(i) == 1
        bits_recup(2*i-1) = 0;
        bits_recup(2*i) = 1; 
    elseif dk(i) == -1
        bits_recup(2*i-1) = 1;
        bits_recup(2*i) = 0;          
    elseif dk(i) == -3
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
Constellations = zeros(3,Nb/2);

for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

        % Emission
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



        suite_diracs = kron(x,[1 zeros(1,Ns-1)]);
        y = filter(h,1,[suite_diracs zeros(1,retard)]);
        y = y(retard+1:end);

        % Canal

        Pr=mean(abs(y).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        y=y+sqrt(sigma)*randn(1,Nb*Ns/2)+ j*sqrt(sigma)*randn(1,Nb*Ns/2);

        z = filter(hr,1,[y zeros(1,retard)]); % Avant de filtrer, on rajoute un nombre de 0 égal au retard
        z = z(retard+1:end); % On supprime le retard (les "retard" premières valeurs)

        z_r = real(z);
        z_r(1) = z_r(1)+sign(z_r(1));

        %Echantillonage 
        bits_recup = zeros(1,Nb);
        dk = zeros(1,Nb/2);

        for i = 0:(Nb/2)-1
            if z_r(1+Ns*i) >=2  
                dk(1+i) = 3;
            elseif z_r(1+Ns*i) >= 0
                dk(1+i) = 1;
    
            elseif z_r(1+Ns*i) >= -2  
                dk(1+i) = -1;
            else
                dk(1+i) = -3;
            end
            if k == 5
                Constellations(1,i+1) = z_r(1+Ns*i);
            elseif k == 6
                Constellations(2,i+1) = z_r(1+Ns*i);
            elseif k == 7
                Constellations(3,i+1) = z_r(1+Ns*i);
            end
 
        end

        % Demapping 
        for i = 1:Nb/2
            if dk(i) == 3
                bits_recup(2*i-1) = 0;
                bits_recup(2*i) = 0;
            elseif dk(i) == 1
                bits_recup(2*i-1) = 0;
                bits_recup(2*i) = 1; 
            elseif dk(i) == -1
                bits_recup(2*i-1) = 1;
                bits_recup(2*i) = 0;          
            elseif dk(i) == -3
                bits_recup(2*i-1) = 1;
                bits_recup(2*i) = 1;         
            end
        end


        
        
        Nerrcomptees(k)=Nerrcomptees(k)+sum(bits~=bits_recup);
        Nsimu(k)=Nsimu(k)+1;
        
        
    end

    TEB(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end


%% Tracé des constellations
figure(2);
for i = 1:3
   plot(Constellations(i,:), zeros(1,Nb/2),'*');
   hold on;
end
plot([-3,-1,1,3], [0 0 0 0], 'o');
title("Constellations en sortie du mapping et de l'échantillonnage pour différentes valeurs de Eb/N0 (4-ASK)");
xlabel('ak');
ylabel('bk');
legend('Eb/N0 = 4','Eb/N0 = 5', 'Eb/N0 = 6','en sortie du mapping');

%% Tracé de la comparaison entre le TEB simulé et le TEB théorique
Pb=(3/4)*qfunc(sqrt((4/5)*Eb_sur_N0_lin));

figure(3)
semilogy(Eb_sur_N0_dB,Pb,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB,'cp','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('Coube théorique','TEB simulé')
title("Chaine 4ASK - comparaison entre le TEB simulé et le TEB théorique")
grid
set(gca,'FontSize',12)