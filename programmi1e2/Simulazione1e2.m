clear all
close all 
clc
global spettro_spe 
load Spettri_dati 
load uscite_4ROI
load eeg_4ROI



% Definition of Variables
window = 50;  
zeropadding = 1000; 
width = 2;
font = 16;
Npop = 2; % Number of populations


% time definition
dt=0.0001;
f_eulero = 1/dt;
tend = 1 + 4*14;
t=(0:dt:tend);
N=length(t);


ROI = 1;   % which ROI to simulate

% Noise
rng(13)   % seed for noise generation
sigma_p = sqrt(9/dt); % Standard deviation of the input noise to excitatory neurons
sigma_f = sqrt(9/dt);% Standard deviation of the input noise to inhibitory neurons
np = randn(2,N)*sigma_p; % Generation of the input noise to excitatory neurons
nf = randn(2,N)*sigma_f; % Generation of the input noise to inhibitory neurons




%% Parameter Values

% base condition;
C(:,1) = 40.*ones(1,2); %Cep
C(:,2) = 40.*ones(1,2); %Cpe
C(:,3) = 40.*ones(1,2); %Csp
C(:,4) = 50.*ones(1,2); %Cps   
C(:,5) = 20.*ones(1,2); %Cfs
C(:,6) = 40.*ones(1,2); %Cfp
C(:,7) = 60*ones(1,2); %Cpf
C(:,8) = 20.*ones(1,2); %Cff


Wp(1,2) = 0;
Wp(2,1) = 0;



%% simulation
e0 = 2.5; % Saturation value of the sigmoid
r = 0.56; % Slope of the sigmoid(1/mV) 



%  Delay between regions (16.6 ms)

D=[0.0166; 0.0166; 0.0166; 0.0166; 0.0166; 0.0166]; 

                

% Synaptic Poles (rad/s) (\omega)
a=[75 30 300 ];   %ae = 75;
                %as = 30;
                %af = 300;

% Synaptic gains (mV)
G=[5.17 4.45 57.1]; %Ge = 5.17;
                    %Gs = 4.45;
                    %Gf = 57.1;
                    
   
load prova12_5giugno2020



         
 v_m_prova1(1,:) = W13*zp_prova1(1,:) + W14*zp_prova1(2,:) +W15*zp_prova1(3,:) + W16*zp_prova1(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova1(2,:) = W23*zp_prova1(1,:) + W24*zp_prova1(2,:) +W25*zp_prova1(3,:) + W26*zp_prova1(4,:);     %l'ingresso al secondp neurone nella prova1
 v_m_prova2(1,:) = W13*zp_prova2(1,:) + W14*zp_prova2(2,:) +W15*zp_prova2(3,:) + W16*zp_prova2(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova2(2,:) = W23*zp_prova2(1,:) + W24*zp_prova2(2,:) +W25*zp_prova2(3,:) + W26*zp_prova2(4,:);      %l'ingresso al secondp neurone nella prova1
 v_m_prova3(1,:) = W13*zp_prova3(1,:) + W14*zp_prova3(2,:) +W15*zp_prova3(3,:) + W16*zp_prova3(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova3(2,:) = W23*zp_prova3(1,:) + W24*zp_prova3(2,:) +W25*zp_prova3(3,:) + W26*zp_prova3(4,:);      %l'ingresso al secondp neurone nella prova1
% v_m_prova1(1,:) = ones(1,570001)*v_m(1,1);
% v_m_prova1(2,:) = ones(1,570001)*v_m(2,1);
% v_m_prova2(1,:) = ones(1,570001)*v_m(1,2);
% v_m_prova2(2,:) = ones(1,570001)*v_m(2,2);
% v_m_prova3(1,:) = ones(1,570001)*v_m(1,3);
% v_m_prova3(2,:) = ones(1,570001)*v_m(2,3);
            
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
            
            riduzione_passo = 100;   % step reduction from 10000 to 100 Hz
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
                    m = v_m_prova3;   % 2XN;
            end
                    
            
            for k= kmax +1:N-1
                up=np(:,k);
                uf=nf(:,k);
                
                if(k>kmax)
                    for i=1:Npop
                        up(i)=up(i)+ Wp(i,:)*zp(:,round(k-D(i)/dt)) + m(i,round(k-D(i)/dt));
                        uf(i)=uf(i)+Wf(i,:)*zp(:,round(k-D(i)/dt));
                    end
                end
                
                vp(:)=C(:,2).*ye(:,k)-C(:,4).*ys(:,k)-C(:,7).*yf(:,k);
                ve(:)=C(:,1).*yp(:,k);
                vs(:)=C(:,3).*yp(:,k);
                vf(:)=C(:,6).*yp(:,k)-C(:,5).*ys(:,k)-C(:,8).*yf(:,k)+yl(:,k);  %
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
                    colore=[0 1 0];
                    Pspe = P_baseline_media;
                    Max_spe = max(Pspe(80:end,:));
                    zp1_prova1 = zp;  
                    eeg1_prova1 = eeg;
                    style = '-';
                case 2
                    colore = [1 0 0];
                    Pspe = P_affected_media;
                    zp1_prova2 = zp;
                    eeg1_prova2 = eeg;
                    style = '-';
                case 3
                    colore = [0 0 1];
                    Pspe = P_unaffected_media;
                    zp1_prova3 = zp;
                    eeg1_prova3 = eeg;
                    style = '-';
            end
            
            figure(1)
            [Peeg,f] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
            if prova == 1
                Max_mod1 = max(Peeg(80:end));
            end
            subplot(121)
            plot(f(80:end),Peeg(80:end)/Max_mod1,'color',colore,'linewidth',width,'linestyle',style)
            if prova == 3
                xlabel('Frequency (Hz)','fontsize',14)
                ylabel('normalized PSD','fontsize',14)
                legend1 = legend('basal','affected', 'unaffected');
                set(legend1,'fontsize',10)
                set(legend1,'location','northeast','box','off')
                title('Model','fontsize',14)
                set(gca,'fontsize',14)
            elseif prova == 1
                axis([0 50 0 1.1])
            end
            hold on
            subplot(122)
            plot(f(80:end),Pspe(80:end,ROI)/Max_spe(ROI),'color',colore,'linewidth', width,'linestyle',style)
            if prova == 3
                xlabel('Frequency (Hz)','fontsize',14)
                ylabel('Normalized PSD','fontsize',14)
                legend1 = legend('basal','affected', 'unaffected');
                set(legend1,'fontsize',10)
                set(legend1,'location','southwest','box','off')
                title('Experimental','fontsize',14)
                set(gca,'fontsize',14)
                assi = axes;
                t1 = title('\fontsize{18} M1h L');
                assi.Visible = 'off';
                t1.Visible = 'on';
            elseif prova == 1
                axis([0 50 0 1.1])
            end
            hold on

            figure(3)
            subplot(2,3,1+1*(prova-1))
            plot(f(80:end),Peeg(80:end)/max(Peeg(80:end)),'color',colore,'linewidth',width)            
            hold on
            plot(f(80:end),Pspe(80:end,ROI)/(max(Pspe(80:end,ROI))),'--','color',0.5*colore+0.3,'linewidth',width)
            switch prova
                case 1
                    ylabel('PSD PMCah (normalized)','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    title('affected','fontsize',font)
                case 3
                    title('unaffected','fontsize',font)
            end
            legend1 = legend('model','experimental');
            set(legend1,'fontsize',12)
            set(legend1,'location','northeast')
            set(gca,'fontsize',font)
            
                       
            figure(2)
           [Peeg,f] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
            if prova == 1
                Max_mod2 = max(Peeg(80:end));
            end
            subplot(121)
            plot(f(80:end),Peeg(80:end)/Max_mod2,'color',colore,'linewidth',width,'linestyle',style)
            if prova == 3
                xlabel('Frequency (Hz)','fontsize',14)
                ylabel('normalized PSD','fontsize',14)
                legend1 = legend('basal','affected', 'unaffected');
                set(legend1,'fontsize',10)
                set(legend1,'location','northeast','box','off')
                title('Model','fontsize',14)
                set(gca,'fontsize',14)
            elseif prova == 1
                axis([0 50 0 1.1])
            end
            hold on
            subplot(122)
            plot(f(80:end),Pspe(80:end,ROI+1)/Max_spe(ROI+1),'color',colore,'linewidth', width,'linestyle',style)
            if prova == 3
                xlabel('Frequency (Hz)','fontsize',14)
                ylabel('Normalized PSD','fontsize',14)
                legend1 = legend('basal','affected', 'unaffected');
                set(legend1,'fontsize',10)
                set(legend1,'location','southwest','box','off')
                title('Experimental','fontsize',14)
                set(gca,'fontsize',14)
                assi = axes;
                t1 = title('\fontsize{18} M1h R');
                assi.Visible = 'off';
                t1.Visible = 'on';
            elseif prova == 1
                axis([0 50 0 1.1])
            end
            hold on

            figure(3)

            subplot(2,3,4+1*(prova-1))
            plot(f(80:end),Peeg(80:end)/max(Peeg(80:end)),'color',colore,'linewidth',width)   
            hold on
            plot(f(80:end),Pspe(80:end,ROI+1)/(max(Pspe(80:end,ROI+1))),'--','color',0.5*colore+0.5,'linewidth',width)
            switch prova
                case 1
                    ylabel('PSD PMCuh (normalized)','fontsize',font)
                case 2
                    xlabel('frequency (Hz)','fontsize',font)
            end
            legend1 = legend('model','experimental');
            set(legend1,'fontsize',12)
            set(legend1,'location','northeast')
            set(gca,'fontsize',font)
          
        figure(4)
            subplot(2,3,1+1*(prova-1))
            [Cxy f] = mscohere(eeg(1,:),eeg(2,:),50,[],zeropadding,fs);
            
            switch prova
                case 1
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{1,2}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{l} M1_{r}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{1,2}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{1,2}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
         
            
            figure(5)
            subplot(2,3,1+1*(prova-1))
            
            switch prova
                case 1
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova1(1,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{1,3}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{l} SMA_{l}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova2(1,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{1,3}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova3(1,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{1,3}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
            
            
            figure(6)
            subplot(2,3,1+1*(prova-1))
            
            switch prova
                case 1
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova1(2,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{1,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{l} SMA_{r}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova2(2,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{1,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova3(2,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{1,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
            
            figure(7)
            subplot(2,3,1+1*(prova-1))
            
            switch prova
                case 1
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova1(3,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{1,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{l} PMC_{l}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova2(3,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{1,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova3(3,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{1,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
            
            figure(8)
            subplot(2,3,1+1*(prova-1))
            
            switch prova
                case 1
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova1(4,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{1,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{l} PMC_{r}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova2(4,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{1,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    [Cxy f] = mscohere(eeg(1,:),eeg_prova3(4,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{1,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
            
        
            figure(9)
            subplot(2,3,1+1*(prova-1))
            
            switch prova
                case 1
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova1(1,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{2,3}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{r} SMA_{l}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova2(1,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{2,3}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova3(1,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{2,3}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
            
            
            figure(10)
            subplot(2,3,1+1*(prova-1))
            
            switch prova
                case 1
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova1(2,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{2,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{r} SMA_{r}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova2(2,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{2,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova3(2,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{2,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
            
            figure(11)
            subplot(2,3,1+1*(prova-1))
            
            switch prova
                case 1
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova1(3,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{2,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{r} PMC_{l}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova2(3,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{2,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova3(3,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{2,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
            
            figure(12)
            subplot(2,3,1+1*(prova-1))
            
            switch prova
                case 1
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova1(4,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_baseline_media{2,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    ylabel('Coeher. M1_{r} PMC_{r}','fontsize',font)
                    title('basal','fontsize',font)
                case 2
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova2(4,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_affected_media{2,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('affected','fontsize',font)
                    xlabel('Frequency (Hz)','fontsize',font)
                case 3
                    [Cxy f] = mscohere(eeg(2,:),eeg_prova3(4,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'color',0.5*colore+0.3,'linewidth',2)
                    hold on 
                    plot(f(100:end),C_unaffected_media{2,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
                    title('unaffected','fontsize',font)
            end
            axis([0 50 0 1])
            set(gca,'fontsize',font)
            
        end
        
        figure
        subplot(321)
        plot(t,zp1_prova1(1,:))
        subplot(322)
        plot(t,zp1_prova1(2,:))
        subplot(323)
        plot(t,zp1_prova2(1,:))    
        subplot(324)
        plot(t,zp1_prova2(2,:)) 
        subplot(325)
        plot(t,zp1_prova3(1,:))    
        subplot(326)
        plot(t,zp1_prova3(2,:)) 
        
                figure
        teeg = t(inizio:riduzione_passo:end);
        subplot(321)
        plot(teeg,eeg1_prova1(1,:))
        subplot(322)
        plot(teeg,eeg1_prova1(2,:))
        subplot(323)
        plot(teeg,eeg1_prova2(1,:))    
        subplot(324)
        plot(teeg,eeg1_prova2(2,:)) 
        subplot(325)
        plot(teeg,eeg1_prova3(1,:))    
        subplot(326)
        plot(teeg,eeg1_prova3(2,:)) 
        
        mean(zp1_prova1,2)
        mean(zp1_prova2,2)
        mean(zp1_prova3,2)
% save uscite_12 zp1_prova1 zp1_prova2 zp1_prova3
% save eeg_12 eeg1_prova1 eeg1_prova2 eeg1_prova3
