% makefile point dipole manual defect
% last edit 3 Mar 2016
function status = dip214(Rmax,filename1,filename2,No,SpinVal1,Spin1,SpinVal2,Spin2,N_grid,gpu)
pathname = pwd; clc; tic; status = []; jam = now;
fprintf('   ============================================ \n')
fprintf('            Dipole Fields Calculation \n')
fprintf('   ============================================ \n')
fprintf('>> Reading input files \n')
geo = poscar(filename1); latt1 = geo.lattice;
geo2 = poscar(filename2); latt2 = geo2.lattice;
Mu = 0.5*(latt1(1,:)+latt1(2,:)+latt1(3,:)); 
dmax = ceil((Rmax+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
fprintf('\n>> Expanding atomic positions \n')
fprintf('     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
if ~isnumeric(Spin1);
    fid = fopen(Spin1); Spin1 = str2num(fgetl(fid));
    while ~feof(fid); Spin1 = [Spin1 str2num(fgetl(fid))]; end
    fclose(fid);
end
spinprt1 = Spin1; 
for i = 1:3:length(Spin1)
    norm = sqrt(Spin1(i)^2+Spin1(i+1)^2+Spin1(i+2)^2);
    if norm > 1; Spin1(i) = Spin1(i)/norm;
        Spin1(i+1) = Spin1(i+1)/norm; Spin1(i+2) = Spin1(i+2)/norm;
    end
end
Spin1 = Spin1*9.274E-21; cu1 = [];
for i = 1:length(No)
    cu1 = [cu1; geo.coords(sum(geo.atomcount(1:No(i)-1))+1:sum(geo.atomcount(1:No(i))),:)];
end
if length(Spin1) == 3 % collinear
    mag = zeros(size(cu1,1),3); for i = 1:size(cu1,1); mag(i,:) = Spin1; end; Spin1 = mag;
else Spin1 = reshape(Spin1,3,size(cu1,1))';
end
a=1; b=0;
for i = 1:length(No)
    b = b+geo.atomcount(No(i));
    Spin1(a:b,:) = Spin1(a:b,:)*SpinVal1(i); a = a+b;
end
if ~isnumeric(Spin2);
    fid = fopen(Spin2); Spin2 = str2num(fgetl(fid));
    while ~feof(fid); Spin2 = [Spin2 str2num(fgetl(fid))]; end
    fclose(fid);
end
spinprt2 = Spin2; 
for i = 1:3:length(Spin2)
    norm = sqrt(Spin2(i)^2+Spin2(i+1)^2+Spin2(i+2)^2);
    if norm > 1; Spin2(i) = Spin2(i)/norm;
        Spin2(i+1) = Spin2(i+1)/norm; Spin2(i+2) = Spin2(i+2)/norm;
    end
end
Spin2 = Spin2*9.274E-21; cu2 = []; Pos3 = [];
for i = 1:length(No)
    cu2 = [cu2; geo2.coords(sum(geo2.atomcount(1:No(i)-1))+1:sum(geo2.atomcount(1:No(i))),:)];
end
if length(Spin2) == 3 % collinear
    mag = zeros(size(cu2,1),3); for i = 1:size(cu2,1); mag(i,:) = Spin2; end; Spin2 = mag;
else Spin2 = reshape(Spin2,3,size(cu2,1))';
end
a=1; b=0;
for i = 1:length(No)
    b = b+geo2.atomcount(No(i));
    Spin2(a:b,:) = Spin2(a:b,:)*SpinVal2(i); a = a+b;
end
for n = 1:size(cu1,1)
    for i = -dmax(1):dmax(1)
        for j = -dmax(2):dmax(2)
            for k = -dmax(3):dmax(3)
                if cu1(n,1)+i >= 0 && cu1(n,1)+i <= 1 ...
                        && cu1(n,2)+j >= 0 && cu1(n,2)+j <= 1 ...
                        && cu1(n,3)+k >= 0 && cu1(n,3)+k <= 1
                    C = (cu1(n,1)+i)*latt1(1,:)+(cu1(n,2)+j)*latt1(2,:)+(cu1(n,3)+k)*latt1(3,:);
                    Pos3 = [Pos3; C Spin1(n,:)];
                else C = (cu2(n,1)+i)*latt2(1,:)+(cu2(n,2)+j)*latt2(2,:)+(cu2(n,3)+k)*latt2(3,:);
                    if sqrt(sum((C-Mu).*(C-Mu))) <= Rmax
                        Pos3 = [Pos3; C Spin2(n,:)];
                    end
                end
            end
        end
    end
end
Sp3 = Pos3(:,4:6); Pos3 = Pos3(:,1:3);
fprintf('     number of ions   : %1.0f \n',size(Pos3,1))
fprintf('\n>> Dipole fields calculation \n')
Nx = N_grid(1); Ny = N_grid(2); Nz = N_grid(3);
x = linspace(0,1,Nx+1); x(end) = [];
y = linspace(0,1,Ny+1); y(end) = [];
z = linspace(0,1,Nz+1); z(end) = [];
[x,y,z] = ndgrid(x,y,z); x = x(:); y = y(:); z = z(:);
xr = x*latt1(1,1) + y*latt1(2,1) + z*latt1(3,1);
yr = x*latt1(1,2) + y*latt1(2,2) + z*latt1(3,2);
zr = x*latt1(1,3) + y*latt1(2,3) + z*latt1(3,3);
Hx = zeros(length(xr),1); 
if gpu == 1; xr = gpuArray(xr); yr = gpuArray(yr); zr = gpuArray(zr); Hx = gpuArray(Hx);
    Sp3 = gpuArray(Sp3); Rmax = gpuArray(Rmax); Pos3 = gpuArray(Pos3);
end
Hy = Hx; Hz = Hx;
fprintf('     number of runs : %1.0f \n',length(xr))
for i = 1:length(xr)
    x = (xr(i)-Pos3(:,1))*1E-8; y = (yr(i)-Pos3(:,2))*1E-8;
    z = (zr(i)-Pos3(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
    SR = Sp3(:,1).*x + Sp3(:,2).*y + Sp3(:,3).*z;
    aPos = (R > 0); SR = SR(aPos); R = R(aPos); 
    x = x(aPos); y = y(aPos); z = z(aPos); Spin = Sp3(aPos,:);
    aPos = (R <= Rmax*1E-8); SR = SR(aPos); R = R(aPos);
    Hx(i) = sum((3*x(aPos).*SR - Spin(aPos,1).*R.^2)./R.^5); 
    Hy(i) = sum((3*y(aPos).*SR - Spin(aPos,2).*R.^2)./R.^5);
    Hz(i) = sum((3*z(aPos).*SR - Spin(aPos,3).*R.^2)./R.^5);
end
fprintf('\n>> Writing output file \n')
Hx = reshape(Hx,[Nx Ny Nz]); Hy = reshape(Hy,[Nx Ny Nz]); Hz = reshape(Hz,[Nx Ny Nz]); 
if gpu == 1; Hx = gather(Hx); Hy = gather(Hy); Hz = gather(Hz); Rmax = gather(Rmax); end
H = sqrt(Hx.^2+Hy.^2+Hz.^2);
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
fprintf(fid,'%s\n',geo.comment); fprintf(fid,' %1.1f\n',1);
fprintf(fid,'   %4.6f  %4.6f  %4.6f\n',geo.lattice);
if ~isempty(geo.symbols)
    cellfun(@(x) fprintf(fid, '%s ', x), geo.symbols); fprintf(fid, '\n');
end
fprintf(fid, ' %d ', geo.atomcount); fprintf(fid, '\nDirect \n');
fprintf(fid, ' %19.16f %19.16f %19.16f \n', geo.coords');
if length(size(H)) < 3
    fprintf(fid,'\n   %d   %d   %d \n',[size(H) 1]);
    fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',H(:));
    fprintf(fid,'\n   %d   %d   %d \n',[size(Hx) 1]);
    fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hx(:)); 
    fprintf(fid,'\n   %d   %d   %d \n',[size(Hy) 1]);
    fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hy(:)); 
    fprintf(fid,'\n   %d   %d   %d \n',[size(Hz) 1]);
    fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hz(:)); 
else
    fprintf(fid,'\n   %d   %d   %d \n',size(H));
    fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',H(:));
    fprintf(fid,'\n   %d   %d   %d \n',size(Hx));
    fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hx(:)); 
    fprintf(fid,'\n   %d   %d   %d \n',size(Hy));
    fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hy(:)); 
    fprintf(fid,'\n   %d   %d   %d \n',size(Hz));
    fprintf(fid,' %14.4e %14.4e %14.4e %14.4e %14.4e \n',Hz(:)); 
end
fclose(fid);

fprintf('\n>> Writing log file\n');
fid = fopen([pathname '/log_make']);
if fid ~= -1
    i = 1; fid = fopen([pathname '/log_make(1)']);
    while fid ~= -1
        i = i+1; fclose(fid); fid = fopen([pathname '/log_make(' num2str(i) ')']);
    end
    fid = fopen([pathname '/log_make(' num2str(i) ')'],'w');
    fprintf(['     ' pathname '/log_make(' num2str(i) ') \n'])
else fid = fopen([pathname '/log_make'],'w');
    fprintf(['     ' pathname '/log_make \n'])
end
fprintf(fid,'   ============================================ \n');
fprintf(fid,'            Point Dipole Calculation \n');
fprintf(fid,'   ============================================ \n');
fprintf(fid,['   Running on: ' datestr(jam) '\n']);
fprintf(fid,'\n   POSCAR location: \n');
fprintf(fid,['   ' filename1 '\n']); fprintf(fid,['   ' filename2 '\n']);
if isempty(geo.symbols);
    if length(No) == 1; fprintf(fid,['\n   ion : '  num2str(No) '\n']);
    else fprintf(fid,['\n   ion of : ['  num2str(No) ']\n']); end;
else ion = '';
    for i = 1:length(No); ion = [ion cell2mat(geo.symbols(No(i))) ' ']; end
    fprintf(fid,['\n   ion : ' ion '\n']);
end
fprintf(fid,['   Magmom1: '  num2str(SpinVal1) ' muB \n']);
fprintf(fid,'   Spins1: \n');
fprintf(fid,'     %4.3f   %4.3f   %4.3f\n',spinprt1);
fprintf(fid,['   Magmom2: '  num2str(SpinVal2) ' muB \n']);
fprintf(fid,'   Spins2: \n');
fprintf(fid,'     %4.3f   %4.3f   %4.3f\n',spinprt2);
fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',Rmax);
fprintf(fid,'   number of grids   : [%1.0f %1.0f %1.0f] \n',Nx,Ny,Nz);
fprintf(fid,'   number of ions    : %1.0f \n',size(Pos3,1));
fprintf(fid,'\n   Writing output file: \n');
fprintf(fid,['   ' sprint '\n\n']);
fprintf(fid,'   ============================================ \n');
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

