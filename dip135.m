% Makefile spin density
% last edit 15 Feb 2018
function status = dip135(Rmax,filename1,N_grid,Mu,dmin,N,gpu)
pathname = pwd; clc; tic; jam = now;
fprintf('   ============================================ \n')
fprintf('            Dipole Fields Calculation \n')
fprintf('   ============================================ \n')
fprintf('>> Reading input files \n')
if isempty(filename1); filename1 = 'CHGCAR'; end
[~,Sx,Sy,Sz,geo1] = chgcar(filename1);
latt1 = geo1.lattice; vol = abs(dot(latt1(1,:),cross(latt1(2,:),latt1(3,:))));
sumx = sum(abs(Sx(:)))^2+sum(abs(Sy(:)))^2+sum(abs(Sz(:)));
if sumx > vol
    Sx = Sx/vol; Sy = Sy/vol; Sz = Sz/vol;
end
sumx = sum(abs(Sx(:)))^2+sum(abs(Sy(:)))^2+sum(abs(Sz(:)));
if sumx > vol
    Sx = Sx/vol; Sy = Sy/vol; Sz = Sz/vol;
end
sumx = [sum(Sx(:)) sum(Sy(:)) sum(Sz(:))];
fprintf('     Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
    sumx1 = [sum(abs(Sx(:))) sum(abs(Sy(:))) sum(abs(Sz(:)))];
    fprintf('     M(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
end

dt = min((latt1(1,:) + latt1(2,:) + latt1(3,:))./(size(Sx)+1))/2;
if coord == 3; Mu = geo.coords(end,:); end
if coord ~= 2; Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:); end

dmax = ceil((Rmax+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
ngrid = size(Sx).*(2*dmax+1);

x = linspace(0,1,size(Sx,1)+1); x(end) = [];
y = linspace(0,1,size(Sx,2)+1); y(end) = [];
z = linspace(0,1,size(Sx,3)+1); z(end) = [];
[x,y,z] = ndgrid(x,y,z); x = x(:); y = y(:); z = z(:);

Sx = Sx(:)*9.274E-21; Sy = Sy(:)*9.274E-21; Sz = Sz(:)*9.274E-21;
fprintf('\n     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
fprintf('     number of grids  : [%1.0f %1.0f %1.0f]\n\n',ngrid)
fprintf('>> Dipole fields calculation \n')
Nx = N_grid(1); Ny = N_grid(2); Nz = N_grid(3);

xm = linspace(0,1,Nx+1); xm(end) = [];
ym = linspace(0,1,Ny+1); ym(end) = [];
% zm = linspace(0,1,Nz+1); zm(end) = [];
zm = 0.091667;

[xm,ym,zm] = ndgrid(xm,ym,zm); xm = xm(:); ym = ym(:); zm = zm(:);
xmu = xm*latt1(1,1) + ym*latt1(2,1) + zm*latt1(3,1);
ymu = xm*latt1(1,2) + ym*latt1(2,2) + zm*latt1(3,2);
zmu = xm*latt1(1,3) + ym*latt1(2,3) + zm*latt1(3,3);

H = zeros(1,length(xmu)); 
if gpu == 1; x = gpuArray(x); y = gpuArray(y); z = gpuArray(z); 
    Sx = gpuArray(Sx); Sy = gpuArray(Sy); Sz = gpuArray(Sz);
    Rmax = gpuArray(Rmax); H = gpuArray(H); 
end
HX = H; HY = H; HZ = H;
fprintf('     number of runs: %1.0f \n',length(xmu));
for n = 1:length(xmu); Hx = 0; Hy = 0; Hz = 0;
    for i = -dmax(1):dmax(1); X = x+i;
        for j = -dmax(2):dmax(2); Y = y+j;
            for k = -dmax(3):dmax(3); Z = z+k;
                xr = (xmu(n) - X*latt1(1,1) - Y*latt1(2,1) - Z*latt1(3,1))*1E-8;
                yr = (ymu(n) - X*latt1(1,2) - Y*latt1(2,2) - Z*latt1(3,2))*1E-8;
                zr = (zmu(n) - X*latt1(1,3) - Y*latt1(2,3) - Z*latt1(3,3))*1E-8;
                R = sqrt(xr.^2 + yr.^2 + zr.^2); SR = Sx.*xr + Sy.*yr + Sz.*zr;
                aPos = (R > dt*1E-8); R = R(aPos); SR = SR(aPos);
                xr = xr(aPos); yr = yr(aPos); zr = zr(aPos);
                SX = Sx(aPos); SY = Sy(aPos); SZ = Sz(aPos);
                aPos = (R <= Rmax*1E-8); R = R(aPos);
                Hx = Hx + sum((3*xr(aPos).*SR(aPos) - SX(aPos).*R.^2)./R.^5);
                Hy = Hy + sum((3*yr(aPos).*SR(aPos) - SY(aPos).*R.^2)./R.^5);
                Hz = Hz + sum((3*zr(aPos).*SR(aPos) - SZ(aPos).*R.^2)./R.^5);
            end
        end
    end
    HX(n) = Hx; HY(n) = Hy; HZ(n) = Hz; toc
end
if gpu == 1; HX = gather(HX); HY = gather(HY); HZ = gather(HZ); Rmax = gather(Rmax); end
fprintf('\n>> Write output file \n')
Hx = reshape(HX,[Nx Ny Nz]); Hy = reshape(HY,[Nx Ny Nz]);
Hz = reshape(HZ,[Nx Ny Nz]); H = sqrt(Hx.^2+Hy.^2+Hz.^2);
fid = fopen([pathname '/LOCPOT_dip']);
if fid ~= -1
    i = 1; fid = fopen([pathname '/LOCPOT_dip(1)']);
    while fid ~= -1
        i = i+1; fclose(fid); fid = fopen([pathname '/LOCPOT_dip(' num2str(i) ')']);
    end
    fid = fopen([pathname '/LOCPOT_dip(' num2str(i) ')'],'w');
    fprintf(['     ' pathname '/LOCPOT_dip(' num2str(i) ') \n'])
    sprint = [pathname '/LOCPOT_dip(' num2str(i) ')'];
else fid = fopen([pathname '/LOCPOT_dip'],'w');
    fprintf(['     ' pathname '/LOCPOT_dip \n'])
    sprint = [pathname '/LOCPOT_dip'];
end
fprintf(fid,'%s\n',geo1.comment);
fprintf(fid,' %1.1f\n',1);
fprintf(fid,'   %4.6f  %4.6f  %4.6f\n',geo1.lattice);
if ~isempty(geo1.symbols)
    cellfun(@(x) fprintf(fid, '%s ', x), geo1.symbols);
    fprintf(fid, '\n');
end
fprintf(fid, ' %d ', geo1.atomcount);
fprintf(fid, '\nDirect \n');
fprintf(fid, ' %19.16f %19.16f %19.16f \n', geo1.coords');
fprintf(fid,'\n   %d   %d   %d \n',size(H));
fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',H(:));
fprintf(fid,'\n   %d   %d   %d \n',size(Hx));
fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hx(:));
fprintf(fid,'\n   %d   %d   %d \n',size(Hy));
fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hy(:));
fprintf(fid,'\n   %d   %d   %d \n',size(Hz));
fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hz(:));
fclose(fid); status = 1;

fprintf('\n>> Writing log file\n');
fid = fopen([pathname '/log_make']);
if fid ~= -1
    i = 1; fid = fopen([pathname '/log_make(1)']);
    while fid ~= -1
        i = i+1; fclose(fid); fid = fopen([pathname 'log_make(' num2str(i) ')']);
    end
    fid = fopen([pathname '/log_make(' num2str(i) ')'],'w');
    fprintf(['     ' pathname '/log_make(' num2str(i) ') \n'])
else fid = fopen([pathname '/log_make'],'w');
    fprintf(['     ' pathname '/log_make \n'])
end
fprintf(fid,'   ============================================ \n');
fprintf(fid,'          Spin Density Calculation test\n');
fprintf(fid,'   ============================================ \n');
fprintf(fid,['   Running on: ' datestr(jam) '\n']);
fprintf(fid,'\n   POSCAR location: \n');
fprintf(fid,['   ' filename1 '\n']);
fprintf(fid,'\n   Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
    fprintf(fid,'     M(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
end
fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',Rmax);
fprintf(fid,'   number of grids   : [%1.0f %1.0f %1.0f] \n',Nx,Ny,Nz);
fprintf(fid,'\n   Writing output file: \n');
fprintf(fid,['   ' sprint '\n']);
fprintf(fid,'\n   ============================================ \n');
fprintf(fid,'             Calculations Completed! \n');
fprintf(fid,'   ============================================ \n');
t = toc; h = floor(t/3600); m = floor((t-h*3600)/60); t = t-h*3600 - m*60;
if h >= 1; fprintf(fid,'   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n\n',h,m,t);  
elseif m >= 1; fprintf(fid,'   Elapsed time is %1.0f min %1.4f sec.\n\n',m,t);
else fprintf(fid,'   Elapsed time is %1.4f sec.\n\n',t);
end
fclose(fid);
fprintf('\n   ============================================ \n')
fprintf('             Calculations Completed! \n')
fprintf('   ============================================ \n')
if h >= 1; fprintf('   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n',h,m,t);  
elseif m >= 1; fprintf('   Elapsed time is %1.0f min %1.4f sec.\n',m,t);
else fprintf('   Elapsed time is %1.4f sec.\n',t);
end
