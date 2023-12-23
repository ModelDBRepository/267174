% here W23 e W14 are inhibitory
clear all
close all 
clc
global history
global P_baseline_media_norm P_affected_media_norm P_unaffected_media_norm C_baseline_media C Wp Wf np nf a 
global W13 W14 W15 W16 W23 W24 W25 W26 zp_prova1 zp_prova2 zp_prova3
load Spettri_dati 
load uscite_4ROI

% Definition of Variables

window = 50;  
zeropadding = 1000; 
Npop = 2; % Number of populations

% time definition

dt=0.0001;
f_eulero = 1/dt;
tend = 1 + 4*14;
t=(0:dt:tend);
N=length(t);

% Noise
rng(13)  % seed for noise generation
sigma_p = sqrt(9/dt); % Standard deviation of the input noise to excitatory neurons
sigma_f = sqrt(9/dt);% Standard deviation of the input noise to inhibitory neurons
np = randn(2,N)*sigma_p; % Generation of the input noise to excitatory neurons
nf = randn(2,N)*sigma_f; % Generation of the input noise to inhibitory neurons

% base condition;
C(:,1) = 40.*ones(1,2); %Cep
C(:,2) = 40.*ones(1,2); %Cpe
C(:,3) = 40.*ones(1,2); %Csp
C(:,4) = 50.*ones(1,2); %Cps   
C(:,5) = 20.*ones(1,2); %Cfs
C(:,6) = 40.*ones(1,2); %Cfp
C(:,7) = 60*ones(1,2); %Cpf
C(:,8) = 20.*ones(1,2); %Cff



%% First attempt for connections
Wf=zeros(Npop);
Wp=zeros(Npop);

Wf = [  0   10
         10    0];




Wp(1,2) = 0;
Wp(2,1) = 0;


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
                    
load prova12inib_19giugno2020
       
            
%% Parameter estimation


% p0 = [ C(2,1) C(2,2) C(2,3) C(2,4) C(2,5) C(2,6) C(2,7) C(2,8)];
% p0 = [p0 Wf(2,1) W23 W24 W25 W26];
p0 = [C(1,1) C(2,1) C(1,2) C(2,2) C(1,3) C(2,3) C(1,4) C(2,4) C(1,5) C(2,5) C(1,6) C(2,6) C(1,7) C(2,7) C(1,8) C(2,8)];
p0 = [p0 Wf(1,2) Wf(2,1) W13 W14 W15 W16 W23 W24 W25 W26];
%p0 = [p0 a(1) a(2) a(3)];


    
 % normalisation with respect to the baseline
    for ROI = 1:6
        P_baseline_media_norm(:,ROI) = P_baseline_media(:,ROI)/max(P_baseline_media(100:end,ROI));
        P_affected_media_norm(:,ROI) = P_affected_media(:,ROI)/max(P_baseline_media(100:end,ROI));
        P_unaffected_media_norm(:,ROI) = P_unaffected_media(:,ROI)/max(P_baseline_media(100:end,ROI));
    end
    
fun = @costo1e2_inib;
    % Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
%%
% % uso fiminunc
% options = optimoptions('fminunc','Display','iter','FunctionTolerance',0.1,'MaxIterations',100,'OutputFcn', @outfun);
% [p, fval] = fminunc(fun,p0,options);
%%
%uso patternsearch
LP = length(p0);
LB = [zeros(1,LP)];
UB = [3*p0(1:16) 10 15 30 60 10 10 160 160 40 40];
options = optimoptions('patternsearch','Display','iter','FunctionTolerance',0.1,'OutputFcn', @outfunp);
[p, fval] = patternsearch(fun,p0,[],[],[],[],LB,UB,[],options);
%%
%uso global search
% LP = length(p0);
% LB = [zeros(1,LP)];
% UB = [50*ones(1,12) 200*ones(1,8) 30*ones(1,4)];
% problem = createOptimProblem('fmincon','objective',fun,...
%     'x0',p0,'lb',LB,'ub',UB);
% gs = GlobalSearch('Display','iter','OutputFcn',@outfung,'NumStageOnePoints',1000,'NumTrialPoints',2000);
% [xg,fg,flg,og] = run(gs,problem)
% p=xg;
%% optimal parameters


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



                 
                 
                    
%% final simulation
width = 2;
font = 16;
ROI = 1;

     

 v_m_prova1(1,:) = W13*zp_prova1(1,:) + W15*zp_prova1(3,:) + W16*zp_prova1(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova1(2,:) = W24*zp_prova1(2,:) + W25*zp_prova1(3,:) + W26*zp_prova1(4,:);     %l'ingresso al secondp neurone nella prova1
 v_m_prova2(1,:) = W13*zp_prova2(1,:) + W15*zp_prova2(3,:) + W16*zp_prova2(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova2(2,:) = W24*zp_prova2(2,:) + W25*zp_prova2(3,:) + W26*zp_prova2(4,:);      %l'ingresso al secondp neurone nella prova1
 v_m_prova3(1,:) = W13*zp_prova3(1,:) + W15*zp_prova3(3,:) + W16*zp_prova3(4,:);      %l'ingresso al primo neurone nella prova1
 v_m_prova3(2,:) = W24*zp_prova3(2,:) +W25*zp_prova3(3,:) + W26*zp_prova3(4,:);      %l'ingresso al secondp neurone nella prova1
 
 % assuming inhibitory W14 and W23 
 v_f_prova1(1,:) = W14*zp_prova1(2,:);
 v_f_prova1(2,:) = W23*zp_prova1(1,:);
 v_f_prova2(1,:) = W14*zp_prova2(2,:);
 v_f_prova2(2,:) = W23*zp_prova2(1,:);
 v_f_prova3(1,:) = W14*zp_prova3(2,:);
 v_f_prova3(2,:) = W23*zp_prova3(1,:);
 
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
                    mf = v_f_prova1;
                case 2
                    m = v_m_prova2;
                    mf = v_f_prova2;
                case 3
                    m = v_m_prova3;
                    mf = v_f_prova3;% Nx2;
            end
            
            for k=1:N-1
                up=np(:,k);
                uf=nf(:,k);
                
                if(k>kmax)
                    for i=1:Npop
                        up(i)=up(i) + m(i,round(k-D(i)/dt))+ Wp(i,:)*zp(:,round(k-D(i)/dt));
                        uf(i)=uf(i) + mf(i,round(k-D(i)/dt))+Wf(i,:)*zp(:,round(k-D(i)/dt));
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
                    colore=[0 1 0];
                    Pspe = P_baseline_media;
                    figure(5)
                    [Cxy f] = mscohere(eeg(1,:),eeg(2,:),50,[],zeropadding,fs);
                    plot(f(100:end),Cxy(100:end),'m',f(100:end),C_baseline_media{1,2}(100:end),'g','linewidth',2)
                    axis([0 50 0 1])
                case 2
                    colore = [1 0 0];
                    Pspe = P_affected_media;
                case 3
                    colore = [0 0 1];
                    Pspe = P_unaffected_media;
            end
            
            figure(1)
            [Peeg,f] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
            subplot(121)
            title('model - g: basal, r: affected, b: unaffected')
            plot(f(80:end),Peeg(80:end),'color',colore,'linewidth',width)
            hold on
            subplot(122)
            title('experimental - g: basal, r: affected, b: unaffected')
            plot(f(80:end),Pspe(80:end,ROI),'color',colore,'linewidth', width)
            hold on
            
            figure(3)
            subplot(3,2,1+2*(prova-1))
            plot(f(80:end),Peeg(80:end)/max(Peeg(80:end)),'color',colore,'linewidth',width)
            hold on
            plot(f(80:end),Pspe(80:end,ROI)/(max(Pspe(80:end,ROI))),'--','color',0.5*colore+0.3,'linewidth',width)
  
            
            figure(2)
            [Peeg,f] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
            subplot(121)
            title('model - g: basal, r: affected, b: unaffected')
            plot(f(80:end),Peeg(80:end),'color',colore,'linewidth',width)
            hold on
            subplot(122)
            title('experimental - g: basal, r: affected, b: unaffected')
            plot(f(80:end),Pspe(80:end,ROI+1),'color',colore,'linewidth',2)
            hold on
            
            figure(3)
            subplot(3,2,2+2*(prova-1))
            plot(f(80:end),Peeg(80:end)/max(Peeg(80:end)),'color',colore,'linewidth',width)
            hold on
            plot(f(80:end),Pspe(80:end,ROI+1)/(max(Pspe(80:end,ROI+1))),'--','color',0.5*colore+0.5,'linewidth',width)
            
 end
 % save prova C Wf a W13 W14 W15 W16 W23 W24 W25 W26 fval