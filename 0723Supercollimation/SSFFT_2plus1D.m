%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Branch Flow  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% 3D-Non-Linear Schrodinger Equation：i*dq/dz=-（q_xx+q_yy）/(2)+ V*q +c*|q|^2*q %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii=sqrt(-1);
N=128;c=0;
Nx=N;xmax=40;dx=(2*xmax)/(N-1);x=[-xmax:dx:xmax];
Ny=N;ymax=40;dy=(2*ymax)/(N-1);y=[-ymax:dy:ymax];
dz=0.5;zmax=500;NUM=zmax/dz+1;z=[0:dz:zmax];
[X,Y,Z] = meshgrid(x,y,z);
BOOL=true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Initial Conditions：plane wave %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% LINEAR_FACTOR in spectrum space  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nx
    if (i<=Nx/2+1)
        FX(i)=((i-1)/(Nx*dx))^2;
    else
        FX(i)=((i-1-Nx)/(Nx*dx))^2;
    end
end
for i=1:Ny
    if (i<=Ny/2+1)
        FY(i)=((i-1)/(Ny*dy))^2;
    else
        FY(i)=((i-1-Ny)/(Ny*dy))^2;
    end
end
FACTOR=-4*pi^2*dz/(2); % important
LINEAR_FACTOR=zeros(Ny,Nx);
for i=1:Ny   
    for j=1:Nx
        LINEAR_FACTOR(i,j)=exp(ii*FACTOR*(FY(i)+FX(j))); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Random Potential %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        u(i,j)=exp(-    ((i-Nx/2)^2+(j-Ny/2)^2)    /Nx/Ny/0.001);
    end
end

% I_SQUARE_AVE=zeros(1,NUM);
% I_AVE_SQUARE=zeros(1,NUM);
% Index=zeros(1,NUM);

% bar = waitbar(0,'Loading...');

V = -1*abs (cos(0.2*X).^2 + cos(0.2*Y).^2 + cos(0.2*Z).^2); 
V=V-min(min(min(V)));

% V=zeros(128,128,201);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Nonlinear Propagation %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:NUM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%  SSFFT  %%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%% Half step in real space %%%%%%%%%%%%
    u=u.*exp(-ii* ( V(:,:,i)')*dz/2);
    %%%%%%%%%%%%% One step in spectre space %%%%%%%%%%
    uf=fft2(u);
    uf=uf.*LINEAR_FACTOR;
    u=ifft2(uf);
    %%%%%%%%%%% Half step in real space %%%%%%%%%%%%%%
    u=u.*exp(-ii* ( V(:,:,i)')*dz/2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%  visualization  %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mod(i*dz,0.1)==0
        
%         I_SQUARE_AVE(i)=mean( mean( abs(u).^4 ));
%         I_AVE_SQUARE(i)=mean( mean( abs(u).^2 ))^2;
%         Index(i)=I_SQUARE_AVE(i)/I_AVE_SQUARE(i)-1;
%         
%         if mod(i*dz,1)==0
%             figure(1)
%             subplot(1,2,1)
%             imagesc(abs(u))
% %             colorbar
%             axis equal
%             axis off
%             subplot(1,2,2)
%             imagesc(V(:,:,i))
%             caxis([0,3])
% %             colorbar
%             axis equal
%             axis off
%             pause(0.5)
%             
%             M(i)=getframe(gcf);
%             I=frame2im(M(i));
%             [I,map]=rgb2ind(I,256);
%             if BOOL==true
%                 imwrite(I,map,['test','.gif'],'Loopcount',inf,'DelayTime',.05)
%                 BOOL=false;
%             else
%                 imwrite(I,map,['test','.gif'],'WriteMode','append','DelayTime',.05)
%             end

%         end
        U(:,:,i)=(u);
% %         str=['Calculating...',num2str(100*i/NUM),'%'];    % 百分比形式显示处理进程,不需要删掉这行代码就行
% %         waitbar(i/NUM,bar,str)  

    end


end




figure();
isosurface(real(U),-0.8)
xlim([1,128])
ylim([1,128])
zlim([200,400])
grid on
hold on
isosurface(V,2.5)
