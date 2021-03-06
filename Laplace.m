% Heat equation solver
% M\'ethodes de Jacobi et de Gauss-Seidel avec surrelaxation SOR
function Laplace_solution()

  % Input parameters
  d=0.2;                           % peak
  N=160;                            % Number of intervals
  prec=1.0e-5; % pr\'ecision $p$ requise sur le residu
  alpha=1.834    % param\`etre de relaxation SOR
  nsel_gs=1;  % s\'electeur de la m\'ethode: 0: Jacobi; 1: Gauss-Seidel

L=0.5;                          % box size
Tin =35;                         % Body temperature 
Tout=-5;                         % Outside temperature
kappa = 1.0;
xc=0.2; 
yc=L/2; 
r=0.1;          % disc 
niter=30000; % nombre maximal d'it\'erations
Ptotal = 0;

xboite = 2*L/N;
yboite = 2*L/N;
tailleboite_x = L-4*L/N;
tailleboite_y = L-4*L/N;

switch nsel_gs
  case 0
    strmethod='Jacobi'
  case 1
    strmethod='GS'
end

% ==== mesh and fields
  h  = L/N;
  L= h*N;
  
  hx = h;
  hy = h;
  Nx=N+1;
  Ny=N+1;
  Lx=L;
  Ly=L;
  
  x=[0:h:L];
  y=[0:h:L];

  [X,Y]=meshgrid(x,y); % % 1D arrays --> 2D arrays
  
  temperature  = zeros(Nx,Ny);
  flag = zeros(Nx,Ny); % =0: dans le vide; =1: dans un conducteur
  
% === Point A
  xa = xc + r + d;
  ya = yc;
  
% ==== Points tangents B and D      
  xb = xc + r*r/(xa-xc);
  yb = yc + r*sqrt(1-(r/(xa-xc))^2);
  xd = xb;
  yd = yc - r*sqrt(1-(r/(xa-xc))^2);
  
  imin=Nx;
  imax=0;
  jmin=Ny;
  jmax=0;
      
  % then the body
  for i=2:(Nx-1)
    for j=2:(Ny-1)
      % check whether we are inside the circle
      diffx = i*hx - xc;
      diffy = j*hy - yc;
      d = sqrt(diffx * diffx + diffy * diffy);
      % check whether we are on the correct side of AC, CD, DA (triangle ABD)
      dAB = (i*hx-xa)*(yb-ya) - (j*hy-ya)*(xb-xa);
      dBD = (i*hx-xb)*(yd-yb) - (j*hy-yb)*(xd-xb);
      dDA = (i*hx-xd)*(ya-yd) - (j*hy-yd)*(xa-xd);
      
      %if(((dAB <= 0.0) && (dBD<=0.0) && (dDA<=0.0) ) || ( d <= r))
      %if((i*hx >= xc -r) && (i*hx <= xc +r) && (j*hy >= yc -r) && (j*hy <= xc +r))
      %if((((i*hx-xc)/r)^2) + (((j*hx-yc)/d)^2) <= 1)
      %if((0.75*(flag(i-1,j)+flag(i,j-1)) + 1)*rand() > 0.999)
      %if(((i*hx-xc)<r)&&((j*hy-yc)<d)&&(i*hx>r)&&(j*hy>r))
      %if((i*hx>=0.07 && i*hx<=0.23 && j*hy<=1/8*hx*i+0.24125 && j*hy>=-1/8*hx*i+0.25875) ||...
       % (i*hx>=0.23 && i*hx<=0.25 && j*hy<=8*hx*i-1.57 && j*hy>=-8*hx*i+2.07) ||...
       % (i*hx>=0.25 && i*hx<=0.27 && j*hy<=-8*hx*i+2.43 && j*hy>=8*hx*i-1.93) ||...
       % (i*hx>=0.27 && i*hx<=0.43 && j*hy<=-1/8*hx*i+0.30375 && j*hy>=1/8*hx*i+0.19625))
      %if((i*hx>=0.05 && i*hx<=0.25 && (i*hx-0.05)^2+(j*hy-0.05)^2>=0.04 && (i*hx-0.05)^2+(j*hy-0.45)^2>=0.04) && j*hy<=0.45 && j*hy>=0.05||... 
       % (i*hx>=0.25 && i*hx<=0.45 && (i*hx-0.45)^2+(j*hy-0.05)^2>=0.04 && (i*hx-0.45)^2+(j*hy-0.45)^2>=0.04 && j*hy<=0.45 && j*hy>=0.05))
      %if(i*hx>=0.05 && i*hx<=0.45 && (i*hx-0.05)^2+(j*hy-0.25)^2>=0.04 && (i*hx-0.45)^2+(j*hy-0.25)^2>=0.04 && j*hy<=0.25 && j*hy>=0.05) 
      %if((i*hx>=0.05 && i*hx<=0.25 && (i*hx-0.25)^2+(j*hy-0.25)^2 <=0.04 && (i*hx-0.25)^2+(j*hy-0.35)^2>=0.01) ||...
       % (i*hx>=0.25 && i*hx<=0.35 && (i*hx-0.25)^2+(j*hy-0.15)^2 <0.01))
      	temperature(i,j) = Tin;
	      flag(i,j) = 1;
	      imin=min(i,imin);
	      imax=max(i,imax);
	      jmin=min(j,jmin);
	      jmax=max(j,jmax);
      else
	      temperature(i,j) = 0;
	      flag(i,j) = 0;
      end
    end
  end
  
temperature(1,:)=Tout;
temperature(end,:)=Tout;
temperature(:,1)=Tout;
temperature(:,end)=Tout;

flag(1,:)=1;
flag(end,:)=1;
flag(:,1)=1;
flag(:,end)=1;
% ==================

nsel_plot=0; % 0/1: do not / do contour plot temperature every iteration with pauses

% Iterations ------------------------------------------------------------------------------------

corr = zeros(niter,1);    % pour stocker les corrections au cours des it\'erations
residuall=zeros(niter,1); % pour stocker les r\'esidus au cours des it\'erations
residu=2*prec;            % on veut faire au moins 1 it\'eration
j=1;                      % compteur du nombre d'it\'erations

while ((j<=niter) & (residu>prec)) %-------------------------------------------------------------
   temperatureold=temperature;
   for jx=2:Nx-1
      for jy=2:Ny-1
	 % N.B. do not update the points on the external conductors: see loops range
	 % N.B. do not update the points on the inner conductors: see flag(.,.)
         if (nsel_gs == 1) %--- Gauss-Seidel            
        	temperaturestar = flag(jx,jy)*temperatureold(jx,jy) + ...
                              (1.0-flag(jx,jy))*0.25*(temperature(jx-1,jy)...
                                                      + temperatureold(jx+1,jy)...
                                                      + temperature(jx,jy-1)...
                                                      + temperatureold(jx,jy+1));
            temperature(jx,jy) = temperatureold(jx,jy) + alpha*(temperaturestar-temperatureold(jx,jy)); % overrelaxation         
	     elseif (nsel_gs == 0) %--- Jacobi
            temperaturestar = flag(jx,jy)*temperatureold(jx,jy) + ...
                              (1.0-flag(jx,jy))*0.25*(temperatureold(jx-1,jy)...
                                                      + temperatureold(jx+1,jy)...
                                                      + temperatureold(jx,jy-1)...
                                                      + temperatureold(jx,jy+1));
            temperature(jx,jy) = temperatureold(jx,jy) + alpha*(temperaturestar-temperatureold(jx,jy)); % overrelaxation
         end
      end
   end

   corr(j)=max(max(abs(temperature-temperatureold)));

% Calcul du residu
   residu=0.0; normres=0.0;
   for jx=2:Nx-1
      for jy=2:Ny-1
	 % N.B. do not consider the points on the external conductors: see loops range
	 % N.B. do not consider the points on the inner conductors: see flag(.,.)
         residu = residu +  (1.0-flag(jx,jy)) * ( temperature(jx,jy) ...
     	   	    - 0.25*(temperature(jx-1,jy)+temperature(jx+1,jy)+temperature(jx,jy-1)+temperature(jx,jy+1)) )^2;
      end
   end
   residu = sqrt(residu);
   residuall(j)=residu; % store for plot
   j=j+1; % increment iteration counter

   
   
fs=16; lw=1;
% Contour plot temperature at each iteration
   if (nsel_plot==1)
     figure(1)
     hc=contourf(X',Y',temperature,20);
     set(gca,'fontsize',fs)
     xlabel('x [m]')
     ylabel('y [m]')
     strj=([' ',strmethod,' h=',num2str(h,3),' \alpha=',num2str(alpha,3),' it=',num2str(j-1,3)]);
     title(['Contours of T',strj])
     axis('equal')
     axis([min(x) max(x) min(y) max(y)])
     colorbar
     pause
   end
     

end %---------------------------------------------------------------------------------------------

niteractual=j % nombre d'it\'erations accomplies
nit=niteractual;
str=([' ',strmethod,' h=',num2str(h,3),' \alpha=',num2str(alpha,3),' nit=',num2str(niteractual,3)]);

%Calcul flux de chaleur

[Jx,Jy]=meshgrid(x,y);
NormeFlux=meshgrid(x,y);

for i=2:(Nx-1)
  for j=2:(Ny-1)
    Jx(i,j)=0.5*(-kappa*(temperature(i+1,j)-temperature(i,j))/hx-kappa*(temperature(i+1,j+1)-temperature(i,j+1))/hx);
    Jy(i,j)=0.5*(-kappa*(temperature(i,j+1)-temperature(i,j))/hy-kappa*(temperature(i+1,j+1)-temperature(i+1,j))/hy);
    NormeFlux(i,j) = sqrt((Jx(i,j))^2 + (Jy(i,j))^2);
  end
end

%Calcul de la puissance totale

iboitemin=Nx-1;
iboitemax=2;
jboitemin=Ny-1;
jboitemax=2;

for i=2:(Nx-1)
  for j=2:(Ny-1)      
      if((X(i,j) > xboite) && (X(i,j) < (xboite+tailleboite_x)) && (Y(i,j) > yboite) && (Y(i,j) < (yboite+tailleboite_y)))
	      iboitemin=min(i,iboitemin);
	      iboitemax=max(i,iboitemax);
	      jboitemin=min(j,jboitemin);
	      jboitemax=max(j,jboitemax);
      end
  end
end

if iboitemin > iboitemax
   swap(iboitemin, iboitemax); 
end

if jboitemin > jboitemax
   swap(jboitemin, jboitemax); 
end

for i=iboitemin:iboitemax
   Ptotal = Ptotal + hx*(-Jy(i,jboitemin)+Jy(i,jboitemax));
end

for j=jboitemin:jboitemax
   Ptotal = Ptotal + hy*(-Jx(iboitemin,j)+Jx(iboitemax,j));
end

Ptotal = Ptotal

figure % 1 contour plot of temperature
%hc=contour(X',Y',temperature,20);
hc=contourf(X',Y',temperature,20);
set(gca,'fontsize',fs)
     xlabel('x [m]')
     ylabel('y [m]')
     title(['Contours of T',str])
     axis('equal')
     axis([min(x) max(x) min(y) max(y)])
     colorbar
     
figure
hj = quiver(X',Y',Jx,Jy);
set(gca,'fontsize',fs)
     xlabel('x [m]')
     ylabel('y [m]')
     title(['Champs des flux de chaleurs',str])
     axis('equal')
     axis([min(x) max(x) min(y) max(y)])
     
figure
hc=contourf(X',Y',NormeFlux,20);
set(gca,'fontsize',fs)
     xlabel('x [m]')
     ylabel('y [m]')
     title(['Norme of flux',str])
     axis('equal')
     axis([min(x) max(x) min(y) max(y)])
     colorbar
     
figure % 2 correction (max (temperature-temperatureold)) vs iteration number
hcorr=plot(corr,'k-');
set(gca,'fontsize',fs)
     xlabel('iteration')
     ylabel('max |T(it)-T(it-1)|')
     title(['Iterations ',str])
     set(gca,'yscale','log')

figure % 3 residual vs iteration number
hres=plot(residuall,'k-');
set(gca,'fontsize',fs)
     xlabel('iteration')
     ylabel('residual')
     title(['Iterations ',str])
     set(gca,'yscale','log')
     
figure
hc=contourf(X',Y',flag,20);
set(gca,'fontsize',fs)
     xlabel('x [m]')
     ylabel('y [m]')
     title(['forme du corps',str])
     axis('equal')
     axis([min(x) max(x) min(y) max(y)])
     colorbar

