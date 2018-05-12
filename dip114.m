% Makefile point dipole manual
% last edit 3 Mar 2016
function status = dip114(Rmax,filename1,No,SpinVal,Spin,N_grid,gpu)
pathname = pwd; tic; jam = now;
fprintf('   ============================================ \n')
fprintf('            Dipole Fields Calculation \n')
fprintf('   ============================================ \n')
fprintf('>> Reading input files \n'); status = [];
if isempty(filename1); filename1 = 'POSCAR'; end
geo = poscar(filename1); latt = geo.lattice;
Mu = 0.5*(latt(1,:)+latt(2,:)+latt(3,:)); 
dmax = ceil((Rmax+Mu)./(latt(1,:)+latt(2,:)+latt(3,:)));
fprintf('\n>> Expanding atomic positions \n')
fprintf('     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
if ~isnumeric(Spin);
    fid = fopen(Spin); Spin = str2num(fgetl(fid));
    while ~feof(fid); Spin = [Spin str2num(fgetl(fid))]; end; fclose(fid);
end
spinprt = Spin; 
for i = 1:3:length(Spin)
    norm = sqrt(Spin(i)^2+Spin(i+1)^2+Spin(i+2)^2);
    if norm > 1; Spin(i) = Spin(i)/norm;
        Spin(i+1) = Spin(i+1)/norm; Spin(i+2) = Spin(i+2)/norm;
    end
end
Spin = Spin*9.274E-21; cu1 = []; Pos = [];
for i = 1:length(No)
    cu1 = [cu1; geo.coords(sum(geo.atomcount(1:No(i)-1))+1:sum(geo.atomcount(1:No(i))),:)];
end
if length(Spin) == 3 % collinear
    mag = zeros(size(cu1,1),3); for i = 1:size(cu1,1); mag(i,:) = Spin; end; Spin = mag;
else Spin = reshape(Spin,3,size(cu1,1))';
end
a=1; b=0;
for i = 1:length(No)
    b = b+geo.atomcount(No(i));
    Spin(a:b,:) = Spin(a:b,:)*SpinVal(i); a = a+b;
end
for n = 1:size(cu1,1)
    cu = cu1(n,1)*latt(1,:)+cu1(n,2)*latt(2,:)+cu1(n,3)*latt(3,:);
    for i = -dmax(1):dmax(1)
        for j = -dmax(2):dmax(2)
            for k = -dmax(3):dmax(3); 
                C = cu + i*latt(1,:)+j*latt(2,:)+k*latt(3,:);
                if sqrt(sum((C-Mu).*(C-Mu))) <= Rmax; 
                    Pos = [Pos; C Spin(n,:)]; 
                end
            end
        end
    end
end
Sp1 = Pos(:,4:6); Pos = Pos(:,1:3);
fprintf('     number of ions   : %1.0f \n',size(Pos,1))
fprintf('\n>> Dipole fields calculation \n')
Nx = N_grid(1); Ny = N_grid(2); Nz = N_grid(3);
x = linspace(0,1,Nx+1); x(end) = [];
y = linspace(0,1,Ny+1); y(end) = [];
z = linspace(0,1,Nz+1); z(end) = [];
% z = 0.091667;

[x,y,z] = ndgrid(x,y,z); x = x(:); y = y(:); z = z(:);
xr = x*latt(1,1) + y*latt(2,1) + z*latt(3,1);
yr = x*latt(1,2) + y*latt(2,2) + z*latt(3,2);
zr = x*latt(1,3) + y*latt(2,3) + z*latt(3,3);
Hx = zeros(length(xr),1); 
if gpu == 1; xr = gpuArray(xr); yr = gpuArray(yr); zr = gpuArray(zr); 
    Sp1 = gpuArray(Sp1); Rmax = gpuArray(Rmax); Pos = gpuArray(Pos); Hx = gpuArray(Hx);
end
Hy = Hx; Hz = Hx;
fprintf('     number of runs : %1.0f \n',length(xr))
for i = 1:length(xr)
    x = (xr(i)-Pos(:,1))*1E-8; y = (yr(i)-Pos(:,2))*1E-8;
    z = (zr(i)-Pos(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
    SR = Sp1(:,1).*x + Sp1(:,2).*y + Sp1(:,3).*z;
    aPos = (R > 0); SR = SR(aPos); R = R(aPos); 
    x = x(aPos); y = y(aPos); z = z(aPos); Spin = Sp1(aPos,:);
    aPos = (R <= Rmax*1E-8); SR = SR(aPos); R = R(aPos);
    Hx(i) = sum((3*x(aPos).*SR - Spin(aPos,1).*R.^2)./R.^5); 
    Hy(i) = sum((3*y(aPos).*SR - Spin(aPos,2).*R.^2)./R.^5);
    Hz(i) = sum((3*z(aPos).*SR - Spin(aPos,3).*R.^2)./R.^5);
end
Hx = reshape(Hx,[Nx Ny Nz]); Hy = reshape(Hy,[Nx Ny Nz]); Hz = reshape(Hz,[Nx Ny Nz]); 
if gpu == 1; Hx = gather(Hx); Hy = gather(Hy); Hz = gather(Hz); Rmax = gather(Rmax); end
H = sqrt(Hx.^2+Hy.^2+Hz.^2);

fprintf('\n>> Writing output file \n')
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
        i = i+1; fclose(fid); fid = fopen([pathname 'log_make(' num2str(i) ')']);
    end
    fid = fopen([pathname '/log_make(' num2str(i) ')'],'w');
    fprintf(['     ' pathname '/log_make(' num2str(i) ') \n'])
else fid = fopen([pathname '/log_make'],'w');
    fprintf(['     ' pathname '/log_make \n'])
end
fprintf(fid,'   ============================================ \n');
fprintf(fid,'              Makefile Point Dipole \n');
fprintf(fid,'   ============================================ \n');
fprintf(fid,['   Running on: ' datestr(jam) '\n']);
fprintf(fid,'\n   POSCAR location: \n');
fprintf(fid,['   ' filename1 '\n']);
if isempty(geo.symbols);
    if length(No) == 1; fprintf(fid,['\n   ion : '  num2str(No) '\n']);
    else fprintf(fid,['\n   ion of : ['  num2str(No) ']\n']); end;
else ion = '';
    for i = 1:length(No); ion = [ion cell2mat(geo.symbols(No(i))) ' ']; end
    fprintf(fid,['\n   ion : ' ion '\n']);
end
fprintf(fid,['   Magmom: '  num2str(SpinVal) ' muB \n']);
fprintf(fid,'   Spins: \n');
fprintf(fid,'     %4.3f   %4.3f   %4.3f\n',spinprt);
fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',Rmax);
fprintf(fid,'   number of grids   : [%1.0f %1.0f %1.0f] \n',Nx,Ny,Nz);
fprintf(fid,'   number of ions    : %1.0f \n',size(Pos,1));
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