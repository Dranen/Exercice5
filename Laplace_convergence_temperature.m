% Heat equation solver
% M\'ethodes de Jacobi et de Gauss-Seidel avec surrelaxation SOR
function Laplace_solution()

var = [50:10:500,550:50:1000,1100:100:2000];

xetude = 0.09;
yetude = 0.25;

for conv = 1:max(size(var))
      % Input parameters
      d=0.05;                           % peak
      N=var(conv);                            % Number of intervals
      prec=1e-05; % pr\'ecision $p$ requise sur le residu
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

    xboite = xc-r-L*0.1;
    yboite = yc-r-L*0.1;
    tailleboite_x = 2*r+L*0.2;
    tailleboite_y = 2*r+d+L*0.2;

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

          if(((dAB <= 0.0) && (dBD<=0.0) && (dDA<=0.0) ) || ( d <= r))
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
    
    ietude = 0;
    jetude = 0;
    
    while(ietude*hx < xetude)
        ietude = ietude +1;
    end
    ietude = ietude-1;
    while(jetude*hy < yetude)
        jetude = jetude +1;
    end
    jetude = jetude-1;
    a = (1-((xetude/hx)-ietude))*temperature(ietude,jetude+1) + (1-((ietude+1)-(xetude/hx)))*temperature(ietude+1,jetude+1);
    b = (1-((yetude/hy)-jetude))*temperature(ietude+1,jetude) + (1-((jetude+1)-(yetude/hy)))*temperature(ietude+1,jetude+1);
    c = (1-((xetude/hx)-ietude))*temperature(ietude,jetude) + (1-((ietude+1)-(xetude/hx)))*temperature(ietude+1,jetude);
    d = (1-((yetude/hy)-jetude))*temperature(ietude,jetude) + (1-((jetude+1)-(yetude/hy)))*temperature(ietude,jetude+1);
    convtemps(conv) = 0.5*((1-((xetude/hx)-ietude))*d + (1-((ietude+1)-(xetude/hx)))*b + (1-((yetude/hy)-jetude))*c + (1-((jetude+1)-(yetude/hy)))*a)
    convtemps2(conv) = interp2(X,Y,temperature,xetude,yetude)
end

figure
loglog(var,convtemps)
set(gca,'fontsize',fs)
xlabel('N')
ylabel('Température [K]')
grid on

figure
loglog(var,convtemps2)
set(gca,'fontsize',fs)
xlabel('N')
ylabel('Température [K]')
grid on