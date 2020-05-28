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


%Filtre de mise en forme et de réception
alpha = 0,35;
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
figure(2);  
semilogy(frequence,fftshift(DSP))
xlabel('frequence')
ylabel('DSP')
title('chaine - DSP du signal transmit')


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

        % Canal

        Pr=mean(abs(y).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        y=y+sqrt(sigma)*randn(1,Nb*Ns/2)+ j*sqrt(sigma)*randn(1,Nb*Ns/2);

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

figure(4)
semilogy(Eb_sur_N0_dB,TEB,'r-.','LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('TEB chaine avec passe-bas équivalent')
title("Chaine passe-bas équivalent - TEB simulé")
grid
set(gca,'FontSize',12)


% Tracé des constellations

ak_mapping = [1,1,-1,-1];       % ak et bk obtenus après mapping (les théoriques)
bk_mapping = [1,-1,1,-1];
ak_obtenus = zeros(1,Nb/2); % ak et bk réellement obtenus après échantillonnage
bk_obtenus = zeros(1,Nb/2);

Pr=mean(abs(y).^2);
sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
y=y+sqrt(sigma)*randn(1,Nb*Ns/2)+ j*sqrt(sigma)*randn(1,Nb*Ns/2);

z = filter(h,1,[y zeros(1,retard)]);
z = z(retard+1:end);
z_r = real(z);
z_i = imag(z);

for i = 0:(Nb/2)-1
    ak_obtenus(1+i) = z_r(1+Ns*i);
    bk_obtenus(1+i) = z_i(1+Ns*i);
end

figure;
plot(ak_obtenus, bk_obtenus, 'r*', ak_mapping, bk_mapping, 'b*');
title('Constellations pour Eb/N0 = 5');
xlabel('ak');
ylabel('bk');
legend("En sortie de l'échantillonneur", 'En sortie du mapping');

%% Comparaison des TEBs (fréquence porteuse / équivalent passe-bas)

Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);
TEB2=zeros(1,length(Eb_sur_N0_dB));

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
                x=[x 1-1i];
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

    TEB2(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end

% Tracé de la comparaison entre le TEB simulé et le TEB théorique

figure(5)
semilogy(Eb_sur_N0_dB,TEB,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB2,'b-.','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('TEB chaine avec passe-bas équivalent','TEB chaine sur fréquence porteuse')
title("Comparaison entre le TEB avec passe-bas équivalent et le TEB de la chaine sur fréquence porteuse")
grid
set(gca,'FontSize',12)
