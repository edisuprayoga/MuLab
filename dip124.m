% Makefile point dipole vasp
% last edit 3 Mar 2016
function status = dip124(Rmax,filename1,N_grid,lim,gpu)
pathname = pwd; tic; jam = now;
fprintf('   ============================================ \n')
fprintf('            Dipole Fields Calculation \n')
fprintf('   ============================================ \n')
fprintf('>> Reading input files \n'); status = [];
if isempty(filename1); filename1 = pwd; end
geo1 = poscar([filename1 '/CONTCAR']);  latt = geo1.lattice;
spin1 = outcar([filename1 '/OUTCAR']); pos1 = geo1.coords;
num = (1:size(spin1,1))'; Pos = [];
spin = sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2);
if length(lim) == 1; cond = (spin(:,1) < lim);
    pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
else cond = (spin(:,1) < lim(2)); pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
     cond = (spin(:,1) > lim(1)); pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
end
pos1 = pos1(:,1)*latt(1,:)+pos1(:,2)*latt(2,:)+pos1(:,3)*latt(3,:);
spinprt = [num spin1 sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2)]'; spin1 = spin1*9.274E-21;
Mu = 0.5*(latt(1,:)+latt(2,:)+latt(3,:)); 
dmax = ceil((Rmax+Mu)./(latt(1,:)+latt(2,:)+latt(3,:)));
fprintf('\n>> Expanding atomic positions \n')
fprintf('     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
for i = -dmax(1):dmax(1)
    for j = -dmax(2):dmax(2)
        for k = -dmax(3):dmax(3);
            Pos = [Pos; pos1+repmat(i*latt(1,:)+j*latt(2,:)+k*latt(3,:),size(pos1,1),1) spin1];
        end
    end
end
if gpu == 1; Pos = gpuArray(Pos); end; spin = Pos(:,4:6); pos = Pos(:,1:3);
fprintf('     number of ions   : %1.0f \n',size(Pos,1))
fprintf('\n>> Dipole fields calculation \n')
Nx = N_grid(1); Ny = N_grid(2); Nz = N_grid(3);
x = linspace(0,1,Nx+1); x(end) = [];
y = linspace(0,1,Ny+1); y(end) = [];
z = linspace(0,1,Nz+1); z(end) = [];
z = 0.21289;

[x,y,z] = ndgrid(x,y,z); x = x(:); y = y(:); z = z(:);
if gpu == 1; x = gpuArray(x); y = gpuArray(y); z = gpuArray(z); end;
xr = x*latt(1,1) + y*latt(2,1) + z*latt(3,1);
yr = x*latt(1,2) + y*latt(2,2) + z*latt(3,2);
zr = x*latt(1,3) + y*latt(2,3) + z*latt(3,3);
Hx = zeros(length(xr),1); if gpu == 1; Hx = gpuArray(Hx); end; Hy = Hx; Hz = Hx;
fprintf('     number of runs : %1.0f \n',length(xr))
for i = 1:length(xr); Spin = spin;
    Pos = repmat([xr(i) yr(i) zr(i)],size(pos,1),1)-pos;
    R = sqrt(Pos(:,1).^2+Pos(:,2).^2+Pos(:,3).^2);
    Pos(R>Rmax,:) = []; Spin(R>Rmax,:) = []; R(R>Rmax) = [];
    Pos(R<1E-6,:) = []; Spin(R<1E-6,:) = []; R(R<1E-6) = [];
    Pos = Pos*1E-8; R = R*1E-8; % to cm
    SR = Spin(:,1).*Pos(:,1) + Spin(:,2).*Pos(:,2) + Spin(:,3).*Pos(:,3);
    Hx(i) = sum((3*Pos(:,1).*SR - Spin(:,1).*R.^2)./R.^5); 
    Hy(i) = sum((3*Pos(:,2).*SR - Spin(:,2).*R.^2)./R.^5);
    Hz(i) = sum((3*Pos(:,3).*SR - Spin(:,3).*R.^2)./R.^5);
end
fprintf('\n>> Writing output file \n')
Hx = reshape(Hx,[Nx Ny Nz]); Hy = reshape(Hy,[Nx Ny Nz]);
Hz = reshape(Hz,[Nx Ny Nz]); H = sqrt(Hx.^2+Hy.^2+Hz.^2);
if gpu == 1;  Hx = gather(Hx); Hy = gather(Hy); Hz = gather(Hz); H = gather(H); end
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
fprintf(fid,'%s\n',geo1.comment); fprintf(fid,' %1.1f\n',1);
fprintf(fid,'   %4.6f  %4.6f  %4.6f\n',geo1.lattice);
if ~isempty(geo1.symbols)
    cellfun(@(x) fprintf(fid, '%s ', x), geo1.symbols); fprintf(fid, '\n');
end
fprintf(fid, ' %d ', geo1.atomcount); fprintf(fid, '\nDirect \n');
fprintf(fid, ' %19.16f %19.16f %19.16f \n', geo1.coords');
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
fprintf(fid,'            Dipole Fields Calculation \n');
fprintf(fid,'   ============================================ \n');
fprintf(fid,['   Running on: ' datestr(jam) '\n']);
fprintf(fid,'\n   POSCAR location: \n');
fprintf(fid,['   ' filename1 '\n']);
fprintf(fid,'\n     No      Mx      My      Mz      Mtotal');
fprintf(fid,'\n   --------------------------------------------\n');
fprintf(fid,'   %4.0f    %6.3f  %6.3f  %6.3f    %6.3f \n',spinprt(:));
fprintf(fid,'   --------------------------------------------\n');
fprintf(fid,'    sum    %6.3f  %6.3f  %6.3f    %6.3f \n',sum(spinprt(2,:)), ...
    sum(spinprt(3,:)),sum(spinprt(4,:)),sum(spinprt(5,:)));
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
