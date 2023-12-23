clear all
close all
clc

global history
global P_baseline_media_norm P_affected_media_norm P_unaffected_media_norm C_baseline_media Cnew Wpnew Wfnew np nf anew v_mnew
load Spettri_dati 

% Definition of Variables
window = 50;  
zeropadding = 1000; 
width = 2;
font = 16;

% time definition
dt=0.0001;
f_eulero = 1/dt;
tend =57;
t=(0:dt:tend);
N=length(t);

% Number of populations
Npop = 4;

% Noise
rng(11)   % seed for noise generation
sigma_p = sqrt(9/dt);  % Standard deviation of the input noise to excitatory neurons 
sigma_f = sqrt(9/dt); % Standard deviation of the input noise to inhibitory neurons
np = randn(Npop,N)*sigma_p; % Generation of the input noise to excitatory neurons
nf = randn(Npop,N)*sigma_f; % Generation of the input noise to inhibitory neurons

e0 = 2.5; % Saturation value of the sigmoid
r = 0.56; % Slope of the sigmoid(1/mV)
%  Delay between regions (16.6 ms)
D=[0.0166; 0.0166; 0.0166; 0.0166; 0.0166; 0.0166]; 
% Synaptic gains (mV)
G=[5.17 4.45 57.1]; 
anew=zeros(Npop,3);

% Connectivity constants
Cnew=zeros(Npop,8);
Wfnew=[ 0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0];
Wpnew=[ 0 0 0 0
        0 0 0 0 
        0 0 0 0
        0 0 0 0];
% external input
v_m=[];
v_mnew=[0 200 10   %ROI 3
     0 200 60   %ROI 4
     0 200 10   %ROI 5
     0 10 200]; %ROI 6
 

%%  this part can be useful if you want to read out the parameters of ROIs separately
% scelta = input('vuoi leggere i parametri per le ROI 3-4 (sì = 1)? ');
% if scelta == 1
%     nome_uscite = input('Nome del file dei parametri? ');
%     stringa1 = ['load ' nome_uscite];
%     eval(stringa1)
%     ver=isempty(v_m);
%     if ver == 0
%         Cnew(1:2,:)=C;
%         Wpnew(1:2,1:2)=Wp;
%         Wfnew(1:2,1:2)=Wf;
%         v_mnew(1:2,:)=v_m;
%         anew(1,:)=a;
%         anew(2,:)=a;
%     else
%         Cnew(1:2,:)=C;
%         Wpnew(1:2,1:2)=Wp;
%         Wfnew(1:2,1:2)=Wf;
%         v_mnew(1:2,:)=v_mnew(1:2,:);
%         anew(1,:)=a;
%         anew(2,:)=a;    
%     end
% end
% 
% scelta = input('vuoi leggere i parametri per le ROI 5-6 da un file (sì = 1)? ');
% if scelta == 1
%     nome_uscite = input('Nome del file dei parametri? ');
%     stringa1 = ['load ' nome_uscite];
%     eval(stringa1)
%     ver=isempty(v_m);
%     if ver == 0
%         Cnew(3:4,:)=C;
%         Wpnew(3:4,3:4)=Wp;
%         Wfnew(3:4,3:4)=Wf;
%         v_mnew(3:4,:)=v_m;
%         anew(3,:)=a;
%         anew(4,:)=a;    
%     else
%         Cnew(3:4,:)=C;
%         Wpnew(3:4,3:4)=Wp;
%         Wfnew(3:4,3:4)=Wf;
%         v_mnew(3:4,:)=v_mnew(3:4,:);
%         anew(3,:)=a;
%         anew(4,:)=a;    
%     end
% end

load prova_4ROI_2giugno2020
Wpnew(1,2) = 0;
Wpnew(2,1) = 0; 
Wpnew(2,4) = 1;


%% Parameter estimation
% p0 = [ Wpnew(1,3) Wpnew(1,4) Wpnew(2,3) Wpnew(2,4) Wpnew(3,1) Wpnew(3,2) Wpnew(4,1) Wpnew(4,2) v_mnew(1,2) v_mnew(2,2) ];
% p0 = [p0 v_mnew(3,2) v_mnew(4,2) v_mnew(1,3) v_mnew(2,3) v_mnew(3,3) v_mnew(4,3) Wfnew(1,2) Wfnew(2,1) Wfnew(3,4) Wfnew(4,3) ];
% p0 = [p0 anew(1,1) anew(1,2) anew(1,3) anew(3,1) anew(3,2) anew(3,3)];
p0 = [ Wpnew(1,3) Wpnew(1,4) Wpnew(2,3) Wpnew(2,4) Wpnew(3,1) Wpnew(3,2) Wpnew(4,1) Wpnew(4,2)  ];
p0 = [p0 Cnew(1,7) Cnew(2,7) Cnew(3,7) Cnew(4,7)];
for ROI = 1:6
    P_baseline_media_norm(:,ROI) = P_baseline_media(:,ROI)/max(P_baseline_media(100:end,ROI));
    P_affected_media_norm(:,ROI) = P_affected_media(:,ROI)/max(P_baseline_media(100:end,ROI));
    P_unaffected_media_norm(:,ROI) = P_unaffected_media(:,ROI)/max(P_baseline_media(100:end,ROI));
end

fun = @costo_4ROI;

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];

%%
% uso fminsearch
% options = optimset('Display','iter','TolFun',1,'MaxIter',300); 
% [p, fval] = fminsearch(fun,p0,options);

%%
% % uso fiminunc
% options = optimoptions('fminunc','Display','iter','FunctionTolerance',0.1,'MaxIterations',100,'OutputFcn', @outfun);
% [p, fval] = fminunc(fun,p0,options);
%%
%uso patternsearch
LP = length(p0);
LB = [zeros(1,LP)];
UB = 80*ones(1,12);
%UB = [80*ones(1,8) 600*ones(1,8) 50*ones(1,4) 800*ones(1,6)];
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
%%
% optimal parameters


Wpnew(1,3)=p(1);
Wpnew(1,4)=p(2);
Wpnew(2,3)=p(3);
Wpnew(2,4)=p(4);
Wpnew(3,1)=p(5);
Wpnew(3,2)=p(6);
Wpnew(4,1)=p(7);
Wpnew(4,2)=p(8);
Cnew(1,7)=p(9);
Cnew(2,7) = p(10);
Cnew(3,7) = p(11);
Cnew(4,7) = p(12);
% v_mnew(1,2)=p(9);
% v_mnew(2,2)=p(10);
% v_mnew(3,2)=p(11);
% v_mnew(4,2)=p(12);
% v_mnew(1,3)=p(13);
% v_mnew(2,3)=p(14);
% v_mnew(3,3)=p(15);
% v_mnew(4,3)=p(16);
% Wfnew(1,2) = p(17);
% Wfnew(2,1) = p(18);
% Wfnew(3,4) = p(19);
% Wfnew(4,3) = p(20);
% anew(1,1) = p(21);
% anew(1,2) = p(22);
% anew(1,3) = p(23);
% anew(2,1) = p(21);
% anew(2,2) = p(22);
% anew(2,3) = p(23);
% anew(3,1) = p(24);
% anew(3,2) = p(25);
% anew(3,3) = p(26);
% anew(4,1) = p(24);
% anew(4,2) = p(25);
% anew(4,3) = p(26);


%% simulazione
ROI=3;

for prova = 1: 3
    inizio = 10000;  % exclusion of the first second due to a possible transitory
    riduzione_passo = 100;  % step reduction from 10000 to 100 Hz
    fs = f_eulero/riduzione_passo;
    eeg=zeros(Npop,(N-1-10000)/riduzione_passo);  
    Coeher = cell(Npop,Npop);

    for j1 = 1:Npop
        for j2 = 1: Npop
            Coher{j1,j2} = zeros(501,1);
        end
    end
    m = v_mnew(:,prova); %ingressi esterni
    
    [zp,ze,zs,zf,vp,ve,vs,vf,yp,ye,ys,yf] = modello_fitting(Npop,D,dt,N,G,np,nf,anew,e0,r,Wpnew,Wfnew,Cnew,m);
    mean_zp(:,prova)=mean(zp,2);
    eeg=diag(Cnew(:,2))*ye(:,inizio:riduzione_passo:end)-diag(Cnew(:,4))*ys(:,inizio:riduzione_passo:end)-diag(Cnew(:,7))*yf(:,inizio:riduzione_passo:end);

    % color choice
    switch prova
        case 1
            colore=[0 1 0];
            Pspe = P_baseline_media;
            zp_prova1 = zp;
            eeg_prova1 = eeg;
        case 2
            colore = [1 0 0];
            Pspe = P_affected_media;
            zp_prova2 = zp;
            eeg_prova2 = eeg;
        case 3
            colore = [0 0 1];
            Pspe = P_unaffected_media;
            zp_prova3 = zp;
            eeg_prova3 = eeg;
    end
    % PSD plots
    figure(1)
    [Peeg3,f3] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
    subplot(121)
    plot(f3(80:end),Peeg3(80:end),'color',colore,'linewidth',width)
    if prova == 3
        xlabel('Frequenza(Hz)','fontsize',10)
        ylabel('PSD','fontsize',10)
        legend1 = legend('basale','affetta', 'non affetta');
        set(legend1,'fontsize',10)
        set(legend1,'location','northeast')
        title('Modello','fontsize',10)
        set(gca,'fontsize',10) 
    elseif prova == 1
       axis([0 50 0 max(Peeg3(100:end))*1.1])
    end
    hold on
    subplot(122)
    plot(f3(80:end),Pspe(80:end,ROI),'color',colore,'linewidth', width)
    if prova == 3
        xlabel('Frequenza (Hz)','fontsize',font)
        ylabel('PSD','fontsize',font)
        legend1 = legend('basale','affetta', 'non affetta');
        set(legend1,'fontsize',10)
        set(legend1,'location','northeast')
        title('Dati sperimentali','fontsize',10)
        set(gca,'fontsize',10)
        assi = axes;
        t1 = title('\fontsize{16} Spettri ROI 3');
        assi.Visible = 'off'; 
        t1.Visible = 'on'; 
    elseif prova == 1
        axis([0 50 0 max(Pspe(100:end,ROI))*1.1])
    end
    hold on
    
    figure(2)
    [Peeg4,f4] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
    subplot(121)
    plot(f4(80:end),Peeg4(80:end),'color',colore,'linewidth',width)
    if prova == 3
        xlabel('Frequenza (Hz)','fontsize',font)
        ylabel('PSD','fontsize',font)
        legend1 = legend('Basale','Affetta', 'Non affetta');
        set(legend1,'fontsize',10)
        set(legend1,'location','northeast')
        title('Modello','fontsize',10)
        set(gca,'fontsize',10)
    elseif prova == 1
       axis([0 50 0 max(Peeg4(100:end))*1.1])
    end
    hold on
    
    subplot(122)
    title('Dati Sperimentali')
    plot(f4(80:end),Pspe(80:end,ROI+1),'color',colore,'linewidth',2)
    if prova == 3
        xlabel('Frequenza (Hz)','fontsize',font)
        ylabel('PSD','fontsize',font)
        legend1 = legend('Basale','Affetta', 'Non affetta');
        set(legend1,'fontsize',10)
        set(legend1,'location','northeast')
        title('Dati Sperimentali','fontsize',10)
        set(gca,'fontsize',10)
        assi = axes;
        t1 = title('\fontsize{16} Spettri ROI 4');
        assi.Visible = 'off'; 
        t1.Visible = 'on'; 
    elseif prova == 1
        axis([0 50 0 max(Pspe(100:end,ROI+1))*1.1])
    end
    hold on
    
   figure (3)
    [Peeg5,f5] =  pwelch(eeg(3,:),window,[],zeropadding,fs);
    subplot(121)
    plot(f5(80:end),Peeg5(80:end),'color',colore,'linewidth',width)
    if prova == 3
        xlabel('Frequenza(Hz)','fontsize',font)
        ylabel('PSD','fontsize',font)
        legend1 = legend('basale','affetta', 'non affetta');
        set(legend1,'fontsize',12)
        set(legend1,'location','northeast')
        title('Modello','fontsize',10)
        set(gca,'fontsize',10)
    elseif prova == 1
        axis([0 50 0 max(Peeg5(100:end))*1.1])
    end
    hold on
    subplot(122)
    plot(f5(80:end),Pspe(80:end,ROI+2),'color',colore,'linewidth', width)
    if prova == 3
        xlabel('Frequenza (Hz)','fontsize',font)
        ylabel('PSD','fontsize',font)
        legend1 = legend('basale','affetta', 'non affetta');
        set(legend1,'fontsize',12)
        set(legend1,'location','northeast')
        title('Dati sperimentali','fontsize',10)
        set(gca,'fontsize',10)
        assi = axes;
        t1 = title('\fontsize{16} Spettri ROI 5');
        assi.Visible = 'off'; 
        t1.Visible = 'on'; 
    elseif prova == 1
        axis([0 50 0 max(Pspe(100:end,ROI+2))*1.1])
    end
    hold on
    
    figure (4)
    [Peeg6,f6] =  pwelch(eeg(4,:),window,[],zeropadding,fs);
    subplot(121)
    plot(f6(80:end),Peeg6(80:end),'color',colore,'linewidth',width)
    if prova == 3
        xlabel('Frequenza(Hz)','fontsize',font)
        ylabel('PSD','fontsize',font)
        legend1 = legend('basale','affetta', 'non affetta');
        set(legend1,'fontsize',12)
        set(legend1,'location','northeast')
        title('Modello','fontsize',10)
        set(gca,'fontsize',10)       
    elseif prova == 1
        axis([0 50 0 max(Peeg6(100:end))*1.1])
    end
    hold on
    subplot(122)
    plot(f6(80:end),Pspe(80:end,ROI+3),'color',colore,'linewidth', width)
    if prova == 3
        xlabel('Frequenza (Hz)','fontsize',font)
        ylabel('PSD','fontsize',font)
        legend1 = legend('basale','affetta', 'non affetta');
        set(legend1,'fontsize',12)
        set(legend1,'location','northeast')
        title('Dati sperimentali','fontsize',10)
        set(gca,'fontsize',10)
        assi = axes;
        t1 = title('\fontsize{16} Spettri ROI 6');
        assi.Visible = 'off'; 
        t1.Visible = 'on'; 
        elseif prova == 1
        axis([0 50 0 max(Pspe(100:end,ROI+3))*1.1])
    end
    hold on
    
   %Coherence
    [Cxy34,f34] = mscohere(eeg(1,:),eeg(2,:),50,[],zeropadding,fs); 
    [Cxy56,f56] = mscohere(eeg(3,:),eeg(4,:),50,[],zeropadding,fs); 
    [Cxy35,f35] = mscohere(eeg(1,:),eeg(3,:),50,[],zeropadding,fs); 
    [Cxy36,f36] = mscohere(eeg(1,:),eeg(4,:),50,[],zeropadding,fs);
    [Cxy54,f54] = mscohere(eeg(3,:),eeg(2,:),50,[],zeropadding,fs);
    [Cxy64,f64] = mscohere(eeg(4,:),eeg(2,:),50,[],zeropadding,fs);
  % Coherence plots
    figure (5)
    subplot(1,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f34(100:end),Cxy34(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f34(100:end),C_baseline_media{3,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coerenza modello vs dati ROI 3-4','fontsize',font)
            xlabel('Frequenza (Hz)','fontsize',font)
            title('Basale','fontsize',10)

        case 2
            plot(f34(100:end),Cxy34(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f34(100:end),C_affected_media{3,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 3-4','fontsize',font)

        case 3
            plot(f34(100:end),Cxy34(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f34(100:end),C_unaffected_media{3,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Non affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 3-4','fontsize',font)
    end
    
    axis([0 50 0 1])
    legend1 = legend('Modello','Dati sperimentali');
    set(legend1,'fontsize',10)
    set(legend1,'location','northeast')
    set(gca,'fontsize',font)
    
    figure (6)
    subplot(1,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f56(100:end),Cxy56(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f56(100:end),C_baseline_media{5,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coerenza modello vs dati ROI 5-6','fontsize',font)
            xlabel('Frequenza (Hz)','fontsize',font)
            title('Basale','fontsize',10)

        case 2
            plot(f56(100:end),Cxy56(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f56(100:end),C_affected_media{5,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 5-6','fontsize',font)

        case 3
            plot(f56(100:end),Cxy56(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f56(100:end),C_unaffected_media{5,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Non affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 5-6','fontsize',font)
    end
    
    axis([0 50 0 1])
    legend1 = legend('Modello','Dati sperimentali');
    set(legend1,'fontsize',10)
    set(legend1,'location','northeast')
    set(gca,'fontsize',font)
   
    
    figure (7)
    subplot(1,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f35(100:end),Cxy35(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f35(100:end),C_baseline_media{3,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coerenza modello vs dati ROI 3-5','fontsize',font)
            xlabel('Frequenza (Hz)','fontsize',font)
            title('Basale','fontsize',10)

        case 2
            plot(f35(100:end),Cxy35(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f35(100:end),C_affected_media{3,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 3-5','fontsize',font)

        case 3
            plot(f35(100:end),Cxy35(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f35(100:end),C_unaffected_media{3,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Non affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 3-5','fontsize',font)
    end
    
    axis([0 50 0 1])
    legend1 = legend('Modello','Dati sperimentali');
    set(legend1,'fontsize',10)
    set(legend1,'location','northeast')
    set(gca,'fontsize',font)
    
    figure (8)
    subplot(1,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f36(100:end),Cxy36(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f36(100:end),C_baseline_media{3,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coerenza modello vs dati ROI 3-6','fontsize',font)
            xlabel('Frequenza (Hz)','fontsize',font)
            title('Basale','fontsize',10)

        case 2
            plot(f36(100:end),Cxy36(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f36(100:end),C_affected_media{3,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 3-6','fontsize',font)

        case 3
            plot(f36(100:end),Cxy36(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f36(100:end),C_unaffected_media{3,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Non affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 3-6','fontsize',font)
    end
    
    axis([0 50 0 1])
    legend1 = legend('Modello','Dati sperimentali');
    set(legend1,'fontsize',10)
    set(legend1,'location','northeast')
    set(gca,'fontsize',font)
    
    figure (9)
    subplot(1,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f54(100:end),Cxy54(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f54(100:end),C_baseline_media{5,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coerenza modello vs dati ROI 5-4','fontsize',font)
            xlabel('Frequenza (Hz)','fontsize',font)
            title('Basale','fontsize',10)

        case 2
            plot(f54(100:end),Cxy54(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f54(100:end),C_affected_media{5,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 5-4','fontsize',font)

        case 3
            plot(f54(100:end),Cxy54(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f54(100:end),C_unaffected_media{5,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Non affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 5-4','fontsize',font)
    end
    
    axis([0 50 0 1])
    legend1 = legend('Modello','Dati sperimentali');
    set(legend1,'fontsize',10)
    set(legend1,'location','northeast')
    set(gca,'fontsize',font)
    
    figure (10)
    subplot(1,3,1+1*(prova-1))
  
    switch prova
        case 1
            plot(f64(100:end),Cxy64(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f64(100:end),C_baseline_media{6,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coerenza modello vs dati ROI 6-4','fontsize',font)
            xlabel('Frequenza (Hz)','fontsize',font)
            title('Basale','fontsize',10)

        case 2
            plot(f64(100:end),Cxy64(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f64(100:end),C_affected_media{6,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 6-4','fontsize',font)

        case 3
            plot(f64(100:end),Cxy64(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f64(100:end),C_unaffected_media{6,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Non affetta','fontsize',10)
            xlabel('Frequenza (Hz)','fontsize',font)
            ylabel('Coerenza modello vs dati ROI 6-4','fontsize',font)
    end
    
    axis([0 50 0 1])
    legend1 = legend('Modello','Dati sperimentali');
    set(legend1,'fontsize',10)
    set(legend1,'location','northeast')
    set(gca,'fontsize',font)
end

%save prova_4ROI Cnew Wfnew Wpnew anew v_mnew fval