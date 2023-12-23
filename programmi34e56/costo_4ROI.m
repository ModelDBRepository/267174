function costo = costo_4ROI(p)

global P_baseline_media_norm P_affected_media_norm P_unaffected_media_norm C_baseline_media Cnew Wpnew Wfnew np nf anew v_mnew

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



window = 50;  
zeropadding = 1000; 

Npop=4;  % Number of ROIs
dt=0.0001;
f_eulero = 1/dt;
tend=57;
t=(0:dt:tend);
N=length(t);
ROI=3;
e0 = 2.5; % Saturation value of the sigmoid
r = 0.56; % Slope of the sigmoid(1/mV) 
% Delays between regions (16.6 ms)
D=[0.0166; 0.0166; 0.0166; 0.0166; 0.0166; 0.0166]; 
% Synaptic gains (mV)
G=[5.17 4.45 57.1]; 

Lp = length(p);
% Parameter setting: not too small
if min(p(1:Lp)) < 0 
    costo = 1E9;
else
    
for prova = 1: 3     
    riduzione_passo = 100;  % step reduction from 10000 to 100 Hz
    fs = f_eulero/riduzione_passo;
    eeg=zeros(Npop,(N-1-10000)/riduzione_passo);  
    Coeher = cell(Npop,Npop);
    for j1 = 1:Npop
        for j2 = 1: Npop
            Coher{j1,j2} = zeros(501,1);
        end
    end
    m = v_mnew(:,prova); % input varies depending on the test (prova)
    [zp,ze,zs,zf,vp,ve,vs,vf,yp,ye,ys,yf] = modello_fitting(Npop,D,dt,N,G,np,nf,anew,e0,r,Wpnew,Wfnew,Cnew,m);
    inizio = 10000; % exclusion of the first second due to a possible transitory
    eeg=diag(Cnew(:,2))*ye(:,inizio:riduzione_passo:end)-diag(Cnew(:,4))*ys(:,inizio:riduzione_passo:end)-diag(Cnew(:,7))*yf(:,inizio:riduzione_passo:end);

    switch prova
        case 1
            spettro_spe = P_baseline_media_norm;
            [Peeg,~] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
            Max_Peeg(1) = max(Peeg(100:end));
            Peeg = Peeg/Max_Peeg(1);
            costo1 = sum((Peeg(100:end) - spettro_spe(100:end,ROI)).^2);
            
            [Peeg,~] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
            Max_Peeg(2) = max(Peeg(100:end));
            Peeg = Peeg/Max_Peeg(2);
            costo2 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+1)).^2);
            
            [Peeg,~] =  pwelch(eeg(3,:),window,[],zeropadding,fs);
            Max_Peeg(3) = max(Peeg(100:end));
            Peeg = Peeg/Max_Peeg(3);
            costo3 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+2)).^2);
            
            [Peeg,~] =  pwelch(eeg(4,:),window,[],zeropadding,fs);
            Max_Peeg(4) = max(Peeg(100:end));
            Peeg = Peeg/Max_Peeg(4);
            costo4 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+3)).^2);

            Cxy35 = mscohere(eeg(1,:),eeg(3,:),50,[],zeropadding,fs); %coherence between ROI 3 e 5
            Cxy36 = mscohere(eeg(1,:),eeg(4,:),50,[],zeropadding,fs); %coherence between  ROI 3 e 6
            Cxy45 = mscohere(eeg(2,:),eeg(3,:),50,[],zeropadding,fs); %coherence between  ROI 4 e 5
            Cxy46 = mscohere(eeg(2,:),eeg(4,:),50,[],zeropadding,fs); %coherence between ROI 4 e 6
            costo35 = sum((Cxy35(100:300) - C_baseline_media{ROI,ROI+2}(100:300)).^2);
            costo36 = sum((Cxy36(100:300) - C_baseline_media{ROI,ROI+3}(100:300)).^2);
            costo45 = sum((Cxy45(100:300) - C_baseline_media{ROI+1,ROI+2}(100:300)).^2);
            costo46 = sum((Cxy46(100:300) - C_baseline_media{ROI+1,ROI+3}(100:300)).^2);

            if max(Cxy35(100:end)) < 0.4 % value found from coherences between ROI
                costo_prova1_35 = (0.4-max(Cxy35(100:end)))*100+costo1 + costo3;  % add a cost based on band centre coherence
            else
                costo_prova1_35 = costo1 + costo3 + costo35;
            end
            if max(Cxy36(100:end)) < 0.15 % value found from coherences between ROI
                costo_prova1_36 = (0.15-max(Cxy36(100:end)))*100+costo1 + costo4;  % add a cost based on band centre coherence
            else
                costo_prova1_36 = costo1 + costo4 + costo36;
            end
            if max(Cxy45(100:end)) < 0.15 % value found from coherences between ROI
                costo_prova1_45 = (0.15-max(Cxy45(100:end)))*100+costo2 + costo3;  % add a cost based on band centre coherence
            else
                costo_prova1_45 = costo2 + costo3 + costo45;
            end
            if max(Cxy46(100:end)) < 0.5 % value found from coherences between ROI
                costo_prova1_46 = (0.5-max(Cxy46(100:end)))*100+costo2 + costo4;  % add a cost based on band centre coherence
            else
                costo_prova1_46 = costo2 + costo4 + costo46;
            end
            
        costo_prova1=costo_prova1_35+costo_prova1_36+costo_prova1_45+costo_prova1_46;
        
        case 2
            spettro_spe = P_affected_media_norm;
            [Peeg,~] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
            Peeg = Peeg/Max_Peeg(1);
            costo1 = sum((Peeg(100:end) - spettro_spe(100:end,ROI)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,ROI)))*100;

            [Peeg,~] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
            Peeg = Peeg/Max_Peeg(2);
            costo2 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+1)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,ROI+1)))*100;
            
            [Peeg,~] =  pwelch(eeg(3,:),window,[],zeropadding,fs);
            Peeg = Peeg/Max_Peeg(3);
            costo3 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+2)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,ROI+2)))*100;
            
            [Peeg,~] =  pwelch(eeg(4,:),window,[],zeropadding,fs);
            Peeg = Peeg/Max_Peeg(4);
            costo4 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+3)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,ROI+3)))*100;
            
            costo_prova2 = costo1 + costo2+costo3+costo4 ;
            
        case 3
            spettro_spe = P_unaffected_media_norm;
            [Peeg,~] =  pwelch(eeg(1,:),window,[],zeropadding,fs);
            Peeg = Peeg/Max_Peeg(1);
            costo1 = sum((Peeg(100:end) - spettro_spe(100:end,ROI)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,ROI)))*100;

            [Peeg,~] =  pwelch(eeg(2,:),window,[],zeropadding,fs);
            Peeg = Peeg/Max_Peeg(2);
            costo2 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+1)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,ROI+1)))*100;
            
            [Peeg,~] =  pwelch(eeg(3,:),window,[],zeropadding,fs);
            Peeg = Peeg/Max_Peeg(3);
            costo3 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+2)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,ROI+2)))*100;

            [Peeg,~] =  pwelch(eeg(4,:),window,[],zeropadding,fs);
            Peeg = Peeg/Max_Peeg(4);
            costo4 = sum((Peeg(100:end) - spettro_spe(100:end,ROI+3)).^2)+abs(max(Peeg(100:300)) - max(spettro_spe(100:300,ROI+3)))*100;

            costo_prova3 = costo1 + costo2 + costo3 + costo4 ;
    end
end
costo = costo_prova1 + costo_prova2 + costo_prova3;
end



