function [zp,ze,zs,zf,vp,ve,vs,vf,yp,ye,ys,yf] = modello_fitting(Npop,D,dt,N,G,np,nf,a,e0,r,Wp,Wf,C,m)

kmax=round(max(D)/dt);
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

    for k=1:N-1
       up=np(:,k)+m; 
       uf=nf(:,k); 

        if(k>kmax)
            for i=1:Npop
                up(i)=up(i)+Wp(i,:)*zp(:,round(k-D(i)/dt));
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
        xp(:,k+1)=xp(:,k)+(G(1)*a(:,1).*zp(:,k)-2*a(:,1).*xp(:,k)-a(:,1).*a(:,1).*yp(:,k))*dt;  
        yp(:,k+1)=yp(:,k)+xp(:,k)*dt; %eulero
        xe(:,k+1)=xe(:,k)+(G(1)*a(:,1).*(ze(:,k)+up(:)./C(:,2))-2*a(:,1).*xe(:,k)-a(:,1).*a(:,1).*ye(:,k))*dt;
        ye(:,k+1)=ye(:,k)+xe(:,k)*dt;
        xs(:,k+1)=xs(:,k)+(G(2)*a(:,2).*zs(:,k)-2*a(:,2).*xs(:,k)-a(:,2).*a(:,2).*ys(:,k))*dt;
        ys(:,k+1)=ys(:,k)+xs(:,k)*dt;
        xl(:,k+1)=xl(:,k)+(G(1)*a(:,1).*uf(:)-2*a(:,1).*xl(:,k)-a(:,1).*a(:,1).*yl(:,k))*dt;
        yl(:,k+1)=yl(:,k)+xl(:,k)*dt;
        xf(:,k+1)=xf(:,k)+(G(3)*a(:,3).*zf(:,k)-2*a(:,3).*xf(:,k)-a(:,3).*a(:,3).*yf(:,k))*dt;
        yf(:,k+1)=yf(:,k)+xf(:,k)*dt;
        
    end
    
end

