clear all
close all
clc

global spettro_spe 
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
ROI=3;

% Noise
rng(11)  % seed for noise generation
sigma_p = sqrt(9/dt); % Standard deviation of the input noise to excitatory neurons
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

load prova_4ROI_2giugno2020
Wfnew(1,2) = 30;
Wpnew(4,1) = 5;

v_mnew(:,1) = 0;


%% simulazione
for prova = 1: 3

    inizio = 10000;  % exclusion of the first second due to a possible transitory
    riduzione_passo = 100; % step reduction from 10000 to 100 Hz
    fs = f_eulero/riduzione_passo;
    eeg=zeros(Npop,(N-1-10000)/riduzione_passo);  
    Coeher = cell(Npop,Npop);
    for j1 = 1:Npop
        for j2 = 1: Npop
            Coher{j1,j2} = zeros(501,1);
        end
    end
    m = v_mnew(:,prova); % external inputs
    
    [zp,ze,zs,zf,vp,ve,vs,vf,yp,ye,ys,yf] = modello_fitting(Npop,D,dt,N,G,np,nf,anew,e0,r,Wpnew,Wfnew,Cnew,m);
    mean_zp(:,prova)=mean(zp,2);
    eeg=diag(Cnew(:,2))*ye(:,inizio:riduzione_passo:end)-diag(Cnew(:,4))*ys(:,inizio:riduzione_passo:end)-diag(Cnew(:,7))*yf(:,inizio:riduzione_passo:end);

    % color choice
    switch prova
        case 1
            colore=[0 1 0];
            style = '-';
            Pspe = P_baseline_media;
            Max_spe = max(Pspe(80:end,:));
            zp_prova1 = zp;
            eeg_prova1 = eeg;
        case 2
            colore = [1 0 0];
            Pspe = P_affected_media;
            style = '-';
            zp_prova2 = zp;
            eeg_prova2 = eeg;
        case 3
            colore = [0 0 1];
            Pspe = P_unaffected_media;
            style = '-';
            zp_prova3 = zp;
            eeg_prova3 = eeg;
    end
    % PSD plots
    figure(1)
    [Peeg3,f3] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
    if prova == 1
        Max_mod3 = max(Peeg3(80:end)); 
    end
    subplot(121)
    plot(f3(80:end),Peeg3(80:end)/Max_mod3,'color',colore,'linewidth',width,'linestyle',style)
    if prova == 3
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('normalized PSD','fontsize',14)
        legend1 = legend('basal','affected', 'unaffected');
        set(legend1,'fontsize',10)
        set(legend1,'location','southwest','box','off')
        title('Model','fontsize',14)
        set(gca,'fontsize',14) 
    elseif prova == 1
       axis([0 50 0 1.1])
    end
    hold on
    subplot(122)
    plot(f3(80:end),Pspe(80:end,ROI)/Max_spe(ROI),'color',colore,'linewidth', width,'linestyle',style)
    if prova == 3
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('Normalized PSD','fontsize',14)
        legend1 = legend('basal','affected', 'unaffected');
        set(legend1,'fontsize',10)
        set(legend1,'location','southwest','box','off')
        title('Experimental','fontsize',14)
        set(gca,'fontsize',14)
        assi = axes;
        t1 = title('\fontsize{18} SMAp L');
        assi.Visible = 'off'; 
        t1.Visible = 'on'; 
    elseif prova == 1
        axis([0 50 0 1.1])
    end
    hold on
    
    erroreSMApL = std( Peeg3(101:300)/Max_mod3 - Pspe(101:300,ROI)/Max_spe(ROI) )
    
    figure(2)
    [Peeg4,f4] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
    if prova == 1
        Max_mod4 = max(Peeg4(80:end)); 
    end
    subplot(121)
    plot(f4(80:end),Peeg4(80:end)/Max_mod4,'color',colore,'linewidth',width,'linestyle',style)
    if prova == 3
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('normalized PSD','fontsize',14)
        legend1 = legend('basal','affected', 'unaffected');
        set(legend1,'fontsize',10)
        set(legend1,'location','southwest','box','off')
        title('Model','fontsize',14)
        set(gca,'fontsize',14) 
    elseif prova == 1
       axis([0 50 0 1.1])
    end
    hold on
    
    subplot(122)
    plot(f4(80:end),Pspe(80:end,ROI+1)/Max_spe(ROI+1),'color',colore,'linewidth', width,'linestyle',style)
   if prova == 3
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('Normalized PSD','fontsize',14)
        legend1 = legend('basal','affected', 'unaffected');
        set(legend1,'fontsize',10)
        set(legend1,'location','southwest','box','off')
        title('Experimental','fontsize',14)
        set(gca,'fontsize',14)
        assi = axes;
        t1 = title('\fontsize{18} SMAp R');
        assi.Visible = 'off'; 
        t1.Visible = 'on'; 
    elseif prova == 1
        axis([0 50 0 1.1])
    end
    hold on
    
    erroreSMApR = std( Peeg4(101:300)/Max_mod4 - Pspe(101:300,ROI+1)/Max_spe(ROI+1) )
    
   figure (3)
    [Peeg5,f5] =  pwelch(eeg(3,:),window,[],zeropadding,fs);
    if prova == 1
        Max_mod5 = max(Peeg5(80:end)); 
    end
    subplot(121)
    plot(f5(100:end),Peeg5(100:end)/Max_mod5,'color',colore,'linewidth',width,'linestyle',style)
    if prova == 3
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('normalized PSD','fontsize',14)
        legend1 = legend('basal','affected', 'unaffected');
        set(legend1,'fontsize',10)
        set(legend1,'location','southwest','box','off')
        title('Model','fontsize',14)
        set(gca,'fontsize',14) 
    elseif prova == 1
       axis([0 50 0 1.1])
    end
    hold on
    subplot(122)
    plot(f5(80:end),Pspe(80:end,ROI+2)/Max_spe(ROI+2),'color',colore,'linewidth', width,'linestyle',style)
   if prova == 3
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('Normalized PSD','fontsize',14)
        legend1 = legend('basal','affected', 'unaffected');
        set(legend1,'fontsize',10)
        set(legend1,'location','southwest','box','off')
        title('Experimental','fontsize',14)
        set(gca,'fontsize',14)
        assi = axes;
        t1 = title('\fontsize{18} PMD L');
        assi.Visible = 'off'; 
        t1.Visible = 'on'; 
    elseif prova == 1
        axis([0 50 0 1.1])
    end
    hold on
    
   errorePMDL = std( Peeg5(101:300)/Max_mod5 - Pspe(101:300,ROI+2)/Max_spe(ROI+2) )
    
    figure (4)
    [Peeg6,f6] =  pwelch(eeg(4,:),window,[],zeropadding,fs);
    if prova == 1
        Max_mod6 = max(Peeg6(80:end)); 
    end
    subplot(121)
    plot(f6(80:end),Peeg6(80:end)/Max_mod6,'color',colore,'linewidth',width,'linestyle',style)
    if prova == 3
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('normalized PSD','fontsize',14)
        legend1 = legend('basal','affected', 'unaffected');
        set(legend1,'fontsize',10)
        set(legend1,'location','southwest','box','off')
        title('Model','fontsize',14)
        set(gca,'fontsize',14) 
    elseif prova == 1
       axis([0 50 0 1.1])
    end
    hold on
    subplot(122)
    plot(f6(80:end),Pspe(80:end,ROI+3)/Max_spe(ROI+3),'color',colore,'linewidth', width,'linestyle',style)
      if prova == 3
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('Normalized PSD','fontsize',14)
        legend1 = legend('basal','affected', 'unaffected');
        set(legend1,'fontsize',10)
        set(legend1,'location','southwest','box','off')
        title('Experimental','fontsize',14)
        set(gca,'fontsize',14)
        assi = axes;
        t1 = title('\fontsize{18} PMD R');
        assi.Visible = 'off'; 
        t1.Visible = 'on'; 
    elseif prova == 1
        axis([0 50 0 1.1])
    end
    hold on
    
  errorePMDR = std( Peeg6(101:300)/Max_mod6 - Pspe(101:300,ROI+3)/Max_spe(ROI+3) )
    
   % Coherence
    [Cxy34 f34] = mscohere(eeg(1,:),eeg(2,:),50,[],zeropadding,fs); %ROI 3-4
    [Cxy56 f56] = mscohere(eeg(3,:),eeg(4,:),50,[],zeropadding,fs);% ROI 5-6
    [Cxy35 f35] = mscohere(eeg(1,:),eeg(3,:),50,[],zeropadding,fs);
    [Cxy36 f36] = mscohere(eeg(1,:),eeg(4,:),50,[],zeropadding,fs);
    [Cxy54 f54] = mscohere(eeg(3,:),eeg(2,:),50,[],zeropadding,fs);
    [Cxy64 f64] = mscohere(eeg(4,:),eeg(2,:),50,[],zeropadding,fs);
    
  % Coherence plots
    figure (5)
    subplot(2,3,1+1*(prova-1))
   
    switch prova
        case 1
            plot(f34(100:end),Cxy34(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f34(100:end),C_baseline_media{3,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coeher. SMAp L - SMAp R','fontsize',14)
            title('Basal','fontsize',14)
            error_C34 = std( Cxy34(101:301) - C_baseline_media{3,4}(101:301) )

        case 2
            plot(f34(100:end),Cxy34(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f34(100:end),C_affected_media{3,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affected','fontsize',14)
            xlabel('Frequency (Hz)','fontsize',14)
            error_C34 = std( Cxy34(101:301) - C_affected_media{3,4}(101:301) )
            
        case 3
            plot(f34(100:end),Cxy34(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f34(100:end),C_unaffected_media{3,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Unaffected','fontsize',14)
            error_C34 = std( Cxy34(101:301) - C_unaffected_media{3,4}(101:301) )
    end
    
    
    axis([0 50 0 1])
    set(gca,'fontsize',14)
%             legend1 = legend('model','experimental');
%             set(legend1,'fontsize',10)
%             set(legend1,'location','northeast','box','off')
%             set(gca,'fontsize',font)
    
    figure (6)
    subplot(2,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f56(100:end),Cxy56(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f56(100:end),C_baseline_media{5,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coeher. PMD L - PMD R','fontsize',12)
            title('Basal','fontsize',14)
            error_C56 = std( Cxy56(101:301) - C_baseline_media{5,6}(101:301) )

        case 2
            plot(f56(100:end),Cxy56(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f56(100:end),C_affected_media{5,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affected','fontsize',14)
            xlabel('Frequency (Hz)','fontsize',14)
            error_C56 = std( Cxy56(101:301) - C_affected_media{5,6}(101:301) )

        case 3
            plot(f56(100:end),Cxy56(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f56(100:end),C_unaffected_media{5,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            hold on 
            title('Unaffected','fontsize',14)
            error_C56 = std( Cxy56(101:301) - C_unaffected_media{5,6}(101:301) )
    end
    
    axis([0 50 0 1])
    set(gca,'fontsize',14)
   
   
    
    figure (7)
    subplot(2,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f35(100:end),Cxy35(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f35(100:end),C_baseline_media{3,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coeher. SMAp L - PMD L','fontsize',14)
            title('Basal','fontsize',14)
            error_C35 = std( Cxy35(101:301) - C_baseline_media{3,5}(101:301) )

        case 2
            plot(f35(100:end),Cxy35(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f35(100:end),C_affected_media{3,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affected','fontsize',14)
            xlabel('Frequency (Hz)','fontsize',14)
            error_C35 = std( Cxy35(101:301) - C_affected_media{3,5}(101:301) )

        case 3
            plot(f35(100:end),Cxy35(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f35(100:end),C_unaffected_media{3,5}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            hold on 
            title('Unaffected','fontsize',14)
            error_C35 = std( Cxy35(101:301) - C_unaffected_media{3,5}(101:301) )
    end
    
    axis([0 50 0 1])
    set(gca,'fontsize',14)
    
    figure (8)
    subplot(2,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f36(100:end),Cxy36(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f36(100:end),C_baseline_media{3,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coeher. SMAp L - PMD R','fontsize',14)
            title('Basal','fontsize',14)
            error_C36 = std( Cxy36(101:301) - C_baseline_media{3,6}(101:301) )

        case 2
            plot(f36(100:end),Cxy36(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f36(100:end),C_affected_media{3,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affected','fontsize',14)
            xlabel('Frequency (Hz)','fontsize',14)
             error_C36 = std( Cxy36(101:301) - C_affected_media{3,6}(101:301) )

        case 3
            plot(f36(100:end),Cxy36(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f36(100:end),C_unaffected_media{3,6}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            hold on 
            title('Unaffected','fontsize',14)
            error_C36 = std( Cxy36(101:301) - C_unaffected_media{3,6}(101:301) )
    end
    
    axis([0 50 0 1])
    set(gca,'fontsize',font)
    
    figure (9)
    subplot(2,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f54(100:end),Cxy54(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f54(100:end),C_baseline_media{5,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coeher. SMAp R - PMD L','fontsize',12)
            title('Basal','fontsize',14)
            error_C54 = std( Cxy54(101:301) - C_baseline_media{5,4}(101:301) )

        case 2
            plot(f54(100:end),Cxy54(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f54(100:end),C_affected_media{5,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affected','fontsize',14)
            xlabel('Frequency (Hz)','fontsize',14)
            error_C54 = std( Cxy54(101:301) - C_affected_media{5,4}(101:301) )

        case 3
            plot(f54(100:end),Cxy54(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f54(100:end),C_unaffected_media{5,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            hold on 
            title('Unaffected','fontsize',14)
            error_C54 = std( Cxy54(101:301) - C_unaffected_media{5,4}(101:301) )
    end
    
    axis([0 50 0 1])
    set(gca,'fontsize',font)
    
    figure (10)
    subplot(2,3,1+1*(prova-1))

    switch prova
        case 1
            plot(f64(100:end),Cxy64(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f64(100:end),C_baseline_media{6,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            ylabel('Coeher. SMAp R - PMD R','fontsize',14)
            title('Basal','fontsize',14)
            error_C64 = std( Cxy64(101:301) - C_baseline_media{6,4}(101:301) )


        case 2
            plot(f64(100:end),Cxy64(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f64(100:end),C_affected_media{6,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            title('Affected','fontsize',14)
            xlabel('Frequency (Hz)','fontsize',14)
            error_C64 = std( Cxy64(101:301) - C_affected_media{6,4}(101:301) )

        case 3
            plot(f64(100:end),Cxy64(100:end),'color',0.5*colore+0.3,'linewidth',2)
            hold on 
            plot(f64(100:end),C_unaffected_media{6,4}(100:end),'--','color',0.5*colore+0.3,'linewidth',2)
            hold on 
            title('Unaffected','fontsize',14)
            error_C64 = std( Cxy64(101:301) - C_unaffected_media{6,4}(101:301) )
    end
    
    axis([0 50 0 1])
    set(gca,'fontsize',14)
end
% save uscite_4ROI zp_prova1 zp_prova2 zp_prova3
% save eeg_4ROI eeg_prova1 eeg_prova2 eeg_prova3