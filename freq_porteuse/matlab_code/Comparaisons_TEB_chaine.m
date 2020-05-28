clc
clear
close all

%% Implantation de la chaine 16QAM avec buit
Nb = 3000;
fe = 10000;
Te = 1/fe;
fp = 2000;
Tp = 1/fp;
Rs = 1000;
Ts = 1/Rs;
Ns = Ts/Te;
M = 16;

alpha = 0.5;
Nb4 = Nb/4;
h = rcosdesign(alpha,Nb4,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,Nb4,Ns);  % reponse impulsionnelle du filtre de réception
retard = Nb*Ns/8;

Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);
TEBQAM=zeros(1,length(Eb_sur_N0_dB));

NerrLim=100;
Nerrcomptees=zeros(1,length(Eb_sur_N0_dB));
Nsimu=zeros(1,length(Eb_sur_N0_dB));


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

    TEBQAM(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end

%% Implantation de la chaine 4ASK avec buit
M = 4;
alpha = 0.5;
h = rcosdesign(alpha,Nb/2,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,Nb/2,Ns);  % reponse impulsionnelle du filtre de réception
retard = Nb*Ns/4;

Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);
TEBASK=zeros(1,length(Eb_sur_N0_dB));

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

    TEBASK(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end

%% Implantation de la chaine 8PSK avec buit
M = 8;
alpha = 0.5;
Nb3 = Nb/3;
h = rcosdesign(alpha,Nb3,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,Nb3,Ns);  % reponse impulsionnelle du filtre de réception
retard = Nb*Ns/6;

Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);
TEB8PSK=zeros(1,length(Eb_sur_N0_dB));

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
            if k == 5
                ak(1,i+1) = z_r(1+Ns*i);
                bk(1,i+1) = z_i(1+Ns*i);
            elseif k == 6
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

    TEB8PSK(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end

%% Implantation de la chaine QPSK avec buit
M = 4;
alpha = 0.5;
h = rcosdesign(alpha,Nb/2,Ns);  % reponse impulsionnelle du filtre de mise en forme
hr = rcosdesign(alpha,Nb/2,Ns);  % reponse impulsionnelle du filtre de réception
retard = Nb*Ns/4;

Eb_sur_N0_dB=[0:6];
Eb_sur_N0_lin=10.^(Eb_sur_N0_dB./10);
TEB3=zeros(1,length(Eb_sur_N0_dB));

NerrLim=100;
Nerrcomptees=zeros(1,length(Eb_sur_N0_dB));
Nsimu=zeros(1,length(Eb_sur_N0_dB));
Constellations = zeros(3,Nb/2);

ak2 = zeros(3,Nb/2); % ak et bk réellement obtenus après échantillonnage
bk2 = zeros(3,Nb/2);

for k=1:length(Eb_sur_N0_lin)

    while (Nerrcomptees(k)<NerrLim)

        %Génération du signal
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

        % Passage dans le filtre de mise en forme
        y = filter(h,1,[suite_diracs zeros(1,retard)]);
        y = y(retard+1:end);
        
         % Canal
        Pr=mean(abs(y).^2);
        sigma=Pr*Ns./(2*log2(M)*Eb_sur_N0_lin(k));
        ez = randn(1,Nb*Ns/2);
        y=y+sqrt(sigma)*randn(1,Nb*Ns/2)+ j*sqrt(sigma)*randn(1,Nb*Ns/2);

        
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
            if k == 5
                ak2(1,i+1) = z_r(1+Ns*i);
                bk2(1,i+1) = z_i(1+Ns*i);
            elseif k == 6
                ak2(2,i+1) = z_r(1+Ns*i);
                bk2(2,i+1) = z_i(1+Ns*i);
            elseif k == 7
                ak2(3,i+1) = z_r(1+Ns*i);
                bk2(3,i+1) = z_i(1+Ns*i);
            end
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

    TEBQPSK(k) = Nerrcomptees(k)/(Nb*Nsimu(k)); % estimation de Pb par calcul du TEB

end

figure(3)
semilogy(Eb_sur_N0_dB,TEBQAM,'r-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEBASK,'y-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEB8PSK,'g-.','LineWidth',3)
hold on
semilogy(Eb_sur_N0_dB,TEBQPSK,'b-.','MarkerSize',10,'LineWidth',3)
xlabel('SNR (dB)')
ylabel('TEB')
legend('TEB chaine 16QAM','TEB chaine 4ASK', 'TEB chaine 8PSK', 'TEB chaine QPSK' )
%title("Comparaison entre le TEB avec passe-bas équivalent et le TEB de la chaine sur fréquence porteuse")
grid
set(gca,'FontSize',12)
