function costo = stima_Neuro_correlate_pop1e2_tre_prove(p)

global P_baseline_media_norm P_affected_media_norm P_unaffected_media_norm C_baseline_media C Wp Wf np nf a 
global W13 W14 W15 W16 W23 W24 W25 W26 zp_prova1 zp_prova2 zp_prova3

% C(2,1) = p(1);
% C(2,2) = p(2);
% C(2,3) = p(3);
% C(2,4) = p(4);
% C(2,5) = p(5);
% C(2,6) = p(6);
% C(2,7) = p(7);
% C(2,8) = p(8);
% Wf(2,1) = p(9);
% W23 = p(10);
% W24 = p(11);
% W25 = p(12);
% W26 = p(13);

C(1,1) = p(1);
C(2,1) = p(2);
C(1,2) = p(3);
C(2,2) = p(4);
C(1,3) = p(5);
C(2,3) = p(6);
C(1,4) = p(7);
C(2,4) = p(8);
C(1,5) = p(9);
C(2,5) = p(10);
C(1,6) = p(11);
C(2,6) = p(12);
C(1,7) = p(13);
C(2,7) = p(14);
C(1,8) = p(15);
C(2,8) = p(16);
Wf(1,2) = p(17);
Wf(2,1) = p(18);
W13 = p(19);
W14 = p(20);
W15 = p(21);
W16 = p(22);
W23 = p(23);
W24 = p(24);
W25 = p(25);
W26 = p(26);
                 
% Wf(1,2) = p(1);
% Wf(2,1) = p(2);
% W13 = p(3);
% W14 = p(4);
% W15 = p(5);
% W16 = p(6);
% W23 = p(7);
% W24 = p(8);
% W25 = p(9);
% W26 = p(10);
% a(1) = p(11);
% a(2) = p(12);
% a(3) = p(13);

window = 50;  
zeropadding = 1000; 
Npop=2; % Number of ROIs

% time definition
dt=0.0001;
f_eulero = 1/dt;
tend = 1 + 10; %4*14;
t=(0:dt:tend);
N=length(t);
% np=sqrt(3/dt)*randn(Npop,N); %matrice
% nf=sqrt(3/dt)*randn(Npop,N);

% Parameter setting: not too small
if min(p) < 0 
    costo = 1E9;
else

e0 = 2.5; % Saturation value of the sigmoid
r = 0.56; % Slope of the sigmoid(1/mV) 


% Delays between regions (16.6 ms)

D=[0.0166; 0.0166; 0.0166; 0.0166; 0.0166; 0.0166]; 

                

% Poli (rad/s) delle sinapsi (uguali per tutte le regioni) (\omega)
%a=[75 30 150 ];   %ae = 75;
                %as = 30;
                %af = 75;

% Synaptic gains (mV)
G=[5.17 4.45 57.1]; %Ge = 5.17;
                    %Gs = 4.45;
                    %Gf = 57.1;

%% Final simulation


 v_m_prova1(1,:) = W13*zp_prova1(1,:) + W14*zp_prova1(2,:) +W15*zp_prova1(3,:) + W16*zp_prova1(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova1(2,:) = W23*zp_prova1(1,:) + W24*zp_prova1(2,:) +W25*zp_prova1(3,:) + W26*zp_prova1(4,:);     %l'ingresso al secondp neurone nella prova1
 v_m_prova2(1,:) = W13*zp_prova2(1,:) + W14*zp_prova2(2,:) +W15*zp_prova2(3,:) + W16*zp_prova2(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova2(2,:) = W23*zp_prova2(1,:) + W24*zp_prova2(2,:) +W25*zp_prova2(3,:) + W26*zp_prova2(4,:);      %l'ingresso al secondp neurone nella prova1
 v_m_prova3(1,:) = W13*zp_prova3(1,:) + W14*zp_prova3(2,:) +W15*zp_prova3(3,:) + W16*zp_prova3(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova3(2,:) = W23*zp_prova3(1,:) + W24*zp_prova3(2,:) +W25*zp_prova3(3,:) + W26*zp_prova3(4,:);      %l'ingresso al secondp neurone nella prova1

for prova = 1: 3
            
            yp=zeros(Npop,N);
            xp=zeros(Npop,N);
            vp=zeros(Npop,1);
            zp=zeros(Npop,N);
            ye=zeros(Npop,N);
            xe=zeros(Npop,N);
            ve=zeros(Npop,1);
            ze=zeros(Npop,N);
            ys=zeros(Npop,N);
            xs=zeros(Npop,N);
            vs=zeros(Npop,1);
            zs=zeros(Npop,N);
            yf=zeros(Npop,N);
            xf=zeros(Npop,N);
            zf=zeros(Npop,N);
            vf=zeros(Npop,1);
            xl=zeros(Npop,N);
            yl=zeros(Npop,N);
            
            riduzione_passo = 100;  % step reduction from 10000 to 100 Hz
            fs = f_eulero/riduzione_passo;
            eeg=zeros(Npop,(N-1-10000)/riduzione_passo);  % exclusion of the first second due to a possible transitory
            Coeher = cell(Npop,Npop);
            for j1 = 1:Npop
                for j2 = 1: Npop
                    Coher{j1,j2} = zeros(501,1);
                end
            end
            
            
            kmax=round(max(D)/dt);
            
            switch prova
                case 1 
                    m = v_m_prova1;
                case 2
                    m = v_m_prova2;   
                case 3
                    m = v_m_prova3;   %  Nx2;
            end
            
            for k=kmax + 1:N-1
                up=np(:,k);
                uf=nf(:,k);
                
                if(k>kmax)
                    for i=1:Npop
                        up(i)=up(i)+ m(i,round(k-D(i)/dt))+ Wp(i,:)*zp(:,round(k-D(i)/dt));
                        uf(i)=uf(i)+Wf(i,:)*zp(:,round(k-D(i)/dt));
                    end
                end
                
                vp(:)=C(:,2).*ye(:,k)-C(:,4).*ys(:,k)-C(:,7).*yf(:,k);
                ve(:)=C(:,1).*yp(:,k);
                vs(:)=C(:,3).*yp(:,k);
                vf(:)=C(:,6).*yp(:,k)-C(:,5).*ys(:,k)-C(:,8).*yf(:,k)+yl(:,k);
                zp(:,k)=2*e0./(1+exp(-r*(vp(:))))-e0;
                ze(:,k)=2*e0./(1+exp(-r*(ve(:))))-e0;
                zs(:,k)=2*e0./(1+exp(-r*(vs(:))))-e0;
                zf(:,k)=2*e0./(1+exp(-r*(vf(:))))-e0;
                
                xp(:,k+1)=xp(:,k)+(G(1)*a(1)*zp(:,k)-2*a(1)*xp(:,k)-a(1)*a(1)*yp(:,k))*dt;  %eulero
                yp(:,k+1)=yp(:,k)+xp(:,k)*dt; %eulero
                xe(:,k+1)=xe(:,k)+(G(1)*a(1)*(ze(:,k)+up(:)./C(:,2))-2*a(1)*xe(:,k)-a(1)*a(1)*ye(:,k))*dt;
                ye(:,k+1)=ye(:,k)+xe(:,k)*dt;
                xs(:,k+1)=xs(:,k)+(G(2)*a(2)*zs(:,k)-2*a(2)*xs(:,k)-a(2)*a(2)*ys(:,k))*dt;
                ys(:,k+1)=ys(:,k)+xs(:,k)*dt;
                xl(:,k+1)=xl(:,k)+(G(1)*a(1)*uf(:)-2*a(1)*xl(:,k)-a(1)*a(1)*yl(:,k))*dt;
                yl(:,k+1)=yl(:,k)+xl(:,k)*dt;
                xf(:,k+1)=xf(:,k)+(G(3)*a(3)*zf(:,k)-2*a(3)*xf(:,k)-a(3)*a(3)*yf(:,k))*dt;
                yf(:,k+1)=yf(:,k)+xf(:,k)*dt;
                
                
            end
            
            inizio = 10000; % exclusion of the first second due to a possible transitory
            eeg=diag(C(:,2))*ye(:,inizio:riduzione_passo:end)-diag(C(:,4))*ys(:,inizio:riduzione_passo:end)-diag(C(:,7))*yf(:,inizio:riduzione_passo:end);
            
            switch prova
                case 1
                    spettro_spe = P_baseline_media_norm;
                    [Peeg,~] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
                    Max_Peeg(1) = max(Peeg(100:end));
                    Peeg = Peeg/Max_Peeg(1);
                    costo1 = sum((Peeg(100:end) - spettro_spe(100:end,1)).^2);
                    
                    [Peeg,~] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
                    Max_Peeg(2) = max(Peeg(100:end));
                    Peeg = Peeg/Max_Peeg(2);
                    costo2 = sum((Peeg(100:end) - spettro_spe(100:end,2)).^2);                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                    Cxy = mscohere(eeg(1,:),eeg(2,:),50,[],zeropadding,fs);
                    costo3 = sum((Cxy(100:end) - C_baseline_media{1,2}(100:end)).^2);
                    
                    if max(Cxy(100:end)) > 0.2
                        costo_prova1 = abs(0.2-max(Cxy(100:end)))*150+costo1 + costo2;  % add a cost based on band centre coherence
                    else
                    %  more weight to the affected area
                        costo_prova1 = costo1 + costo2 + costo3 ;
                    end
                case 2
                    spettro_spe = P_affected_media_norm;
                    [Peeg,~] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
                    Peeg = Peeg/Max_Peeg(1);
                    costo1 = sum((Peeg(100:end) - spettro_spe(100:end,1)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,1)))*100;
                    
                    [Peeg,~] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
                    Peeg = Peeg/Max_Peeg(2);
                    costo2 = sum((Peeg(100:end) - spettro_spe(100:end,2)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,2)))*100;
                    %  more weight to the affected area
                    costo_prova2 = 2*costo1 + costo2 ;
                case 3
                    spettro_spe = P_unaffected_media_norm;
                    [Peeg,~] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
                    Peeg = Peeg/Max_Peeg(1);
                   
                    costo1 = sum((Peeg(100:end) - spettro_spe(100:end,1)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,1)))*200;
                    
                    [Peeg,~] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
                    Peeg = Peeg/Max_Peeg(2);
                    costo2 = sum((Peeg(100:end) - spettro_spe(100:end,2)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,2)))*100; 
                    %  more weight to the affected area
                    costo_prova3 = costo1 + costo2 ;
            end
      
            
%disp(costo);
            end
costo = costo_prova1 + costo_prova2 + costo_prova3;

% disp(costo)
end



