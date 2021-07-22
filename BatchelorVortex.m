% Batchelor Vortex, created by Yiding Zhu, 2020.09.10, changed by Ganglong Yu, 2021/05/27
clear; clc; %close all;

dim.Nx=64;     % spatial point number along x-axis;
dim.Ny=64;     % spatial point number along y-axis;
dim.Nz=3;
dim.Nt=2;     % time point number
Nx = dim.Nx;
Ny = dim.Ny;
Nz = dim.Nz;
Nt = dim.Nt;

% constant properties
prop.P0=1e5;      % coming stream pressure, 1 atm, unit: Pa
prop.q=1;         % the swirl strength, given as a ratio between the maximum tangential velocity and the core velocity
prop.niu=1.01e-6; % viscous coefficient
prop.rho0=1e3;    % coming stream density, unit: kg/m3
prop.R_max = 0.2; % max radius
P0 = prop.P0;
niu = prop.niu;
rho0 = prop.rho0;

% data diversity
Re_num = 50;
R0_num = 2;
alpha_num = 32;
smp_num = 4;   % num of sample (smp)
smp_ratio = linspace(1.5,3,smp_num);
R0_range = linspace(0.001,0.03,R0_num);
Re_range = linspace(500,5000,Re_num);
alpha_range = (0:alpha_num-1)*2*pi/alpha_num;

for h =1:Re_num
	for i =1:R0_num
        for k = 1:alpha_num
            R0 = R0_range(i);
            Re = Re_range(h);                       
            U0 = Re*niu/R0;
            alpha = alpha_range(k);

            [rl,Ur,Wr,Pr,delta] = velocity_r(R0,U0,dim,prop);
            [Uc,Vc,Wc,Pc,X,Y] = polar2cart(rl,Ur,Wr,Pr,alpha,dim,delta,prop);      

            % non-dimensionalized            
            Uc = Uc/U0;
            Vc = Vc/U0;
            Wc = Wc/U0;
            Pc = Pc/(rho0*U0^2);
            dx = mean(X(1,2:end)-X(1,1:end-1));
            X = X/dx;
            Y = Y/dx;

            save(['Data/','Bdata',num2str(h-1),'_',num2str(i-1),'_',num2str(k-1),'.mat'],'P0','U0','niu','rho0','Pc','Wc','Uc','Vc','Re','X','Y')
            fprintf('Cartesian coordinates, rotation, non-dimensionalized\n');
                 
            % region sample
            U = zeros(Ny,Nx,Nz,Nt); V = U; W = U; P = U;           
            for j = 1:smp_num            
                smp_size = (size(X,1)-1)/smp_ratio(j); % size of sample region
                pos_x = round(linspace(1,size(X,1)-smp_size-1,smp_num));   % corner position of sample region
                pos_y = pos_x;
                for t = 1:smp_num
                    smp_x = round(linspace(pos_x(t),pos_x(t)+smp_size,Nx));
                    smp_y = round(linspace(pos_y(t),pos_y(t)+smp_size,Ny));
                    smpX = X(smp_x,smp_y);
                    smpY = Y(smp_x,smp_y);
                    for iz = 1:Nz
                        for it = 1:Nt
                            U(:,:,iz,it) = Uc(smp_x,smp_y,iz,it);
                            V(:,:,iz,it) = Vc(smp_x,smp_y,iz,it);
                            W(:,:,iz,it) = Wc(smp_x,smp_y,iz,it);
                            P(:,:,iz,it) = Pc(smp_x,smp_y,iz,it);
                        end
                    end
                    save(['Data/','Bdata',num2str(h-1),'_',num2str(i-1),'_',num2str(k-1),'_',num2str(j-1),'_',num2str(t-1),'.mat'],'P0','U0','niu','rho0','P','W','U','V','Re','smpX','smpY')
                end
            end
            fprintf('sample field data\n');
        end
	end
end


function [rl,Ur,Wr,Pr,delta] = velocity_r(R0,U0,dim,prop)
%VELOCITY_R velocity along the radial axis
%===================================================
%R0:    core size
%U0:    freestream velocity, velocity scale
%dim:   structure of dimension
%prop:  structure of constant properties
%===================================================
%rl:    discretization along radial axis
%Ur:    azimuthal velocity
%Wr:    axial velocity
%Pr:    pressure
%delta: structure of discretization
%===================================================

    % constant properties
    P0=prop.P0;     
    z0=100*R0;   % effective z position
    q=prop.q;         
    niu=prop.niu; 
    rho0=prop.rho0;    
    R_max = prop.R_max;
    W0 = U0;    
    Nz = dim.Nz;
    Nt = dim.Nt;
       
    % discretization
    t0 = R0/U0;
    dL = R0;
    dz = dL/2;
    dr = dL/10;
    dt = t0/100; % Courant number <= 1

    % allocate memory
    rl=(dr:dr:R_max)';
    Ur = zeros(length(rl),Nz,Nt); % azimuthal velocity
    Wr = Ur; % axial velocity
    Pr = Ur; % pressure		
       
    for it=1:Nt % moving time point
        t=it*dt;
        % calculating velocity along radial axis with equations
        for iz=1:Nz
            z=iz*dz+z0;
            R=(R0^2+4*niu*(t+z/U0))^0.5;                   % a measure of the core size
            Ur(:,iz,it)=q*W0*(1-exp(-(rl/R).^2))./(rl/R0); % azimuthal velocity components
            Wr(:,iz,it)=U0+W0/((R/R0)^2)*exp(-(rl/R).^2);  % axial velocity
            yeta=W0*rl.^2/(4*niu*z);
            C=Ur(:,iz,it).*rl;
            C0=C./(1-exp(-yeta));
            Pyeta=(1-exp(-yeta)).^2./yeta+2*(expint(yeta)-expint(2*yeta));
            Pr(:,iz,it)=P0-rho0*W0*C0.^2/(8*niu*z).*Pyeta; % pressure, P0-rho0*C0^2/(8*niu*z)*Pyeta
        end
    end
    
    delta.dr = dr;
    delta.dz = dz;
    delta.dt = dt;
end


function [Uc,Vc,Wc,Pc,X,Y] = polar2cart(rl,Ur,Wr,Pr,alpha,dim,delta,prop)
%POLAR2CART transfer polar coordinates to Cartesian coordinates
%===================================================
%rl:    discretization along radial axis
%Ur:    azimuthal velocity
%Wr:    axial velocity
%Pr:    pressure
%alpha: rotation angle
%dim:   structure of dimension
%delta: structure of discretization
%prop:  structure of constant properties
%===================================================
%Uc:    x velocity
%Vc:    y velocity
%Uc:    z velocity (axial)
%Pc:    pressure
%X&Y:   meshgrid
%===================================================

    len_pol = size(Ur,1);
    len_cart = floor(len_pol/1.5);
    dr = delta.dr;
    Nz = dim.Nz;
    Nt = dim.Nt;
    R_max = prop.R_max;
    
    % allocate memory
    Uc = zeros(len_cart,len_cart,Nz,Nt); Vc = Uc; Wc = Uc; Utc = Uc; Pc = Uc;
    
    for it=1:Nt % moving time point
        for iz=1:Nz
            % casting onto a Cartesian coordinate
            minx = dr; miny = minx;
            maxx = R_max/1.5; maxy = maxx;
            x = linspace(minx,maxx,len_cart);
            y = linspace(miny,maxy,len_cart);
            [X,Y]=meshgrid(x,y);
            
            rotate = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];
            coord = rotate*[X(:)';Y(:)'];
            X = reshape(coord(1,:),len_cart,len_cart);
            Y = reshape(coord(2,:),len_cart,len_cart);
            cx=0; % vortex center
            cy=0;
            RR=((X-cx).^2+(Y-cy).^2).^0.5;
            Pc(:,:,iz,it)=interp1(rl,Pr(:,iz,it),RR,'linear');
            Utc(:,:,iz,it)=interp1(rl,Ur(:,iz,it),RR,'linear');
            Wc(:,:,iz,it)=interp1(rl,Wr(:,iz,it),RR,'linear');
            RY=(Y-cy)./RR;  % over zero, NaN
            RX=-(X-cx)./RR;
            Uc(:,:,iz,it)=Utc(:,:,iz,it).*RX;
            Vc(:,:,iz,it)=Utc(:,:,iz,it).*RY;
        end
    end    
end