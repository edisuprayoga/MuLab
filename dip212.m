% ZPE point dipole manual defect
% last edit 3 Mar 2016
function status = dip212(Mu,coord,Rmax,filename1,filename2,No,SpinVal1, ...
    Spin1,SpinVal2,Spin2,Rz,N,gpu,log)
pathname = pwd;
if log == 1; clc; tic; jam = now;
    fprintf('   ============================================ \n')
    fprintf('            Dipole Fields Calculation \n')
    fprintf('   ============================================ \n')
    fprintf('>> Reading input files \n')
end
geo = poscar(filename1); latt1 = geo.lattice;
geo2 = poscar(filename2); latt2 = geo2.lattice;
if coord == 3; Mu = geo.coords(end,:); end
if coord ~= 2; Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:); end
if length(Rz) == 1; Rz = Rz*[1 1 1]; end; 
Rmax = Rmax + sqrt(Rz(1)^2+Rz(2)^2+Rz(3)^2)/2;
dmax = ceil((Rmax+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
if log == 1
    fprintf('\n>> Expanding atomic positions \n')
    fprintf('     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
end
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
for i = 1:3:length(Spin2)
    norm = sqrt(Spin2(i)^2+Spin2(i+1)^2+Spin2(i+2)^2);
    if norm > 1; Spin2(i) = Spin2(i)/norm;
        Spin2(i+1) = Spin2(i+1)/norm; Spin2(i+2) = Spin2(i+2)/norm;
    end
end
spinprt2 = Spin2; Spin2 = Spin2*9.274E-21; cu2 = [];
for i = 1:length(No)
    cu2 = [cu2; geo2.coords(sum(geo2.atomcount(1:No(i)-1))+1:sum(geo2.atomcount(1:No(i))),:)];
end
if length(Spin2) == 3 % collinear
    mag = zeros(size(cu2,1),3); for i = 1:size(cu2,1); mag(i,:) = Spin2; end; Spin2 = mag;
else Spin2 = reshape(Spin2,3,size(cu2,1))';
end
a=1; b=0; Pos3 = [];
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
if log == 1
    fprintf('     number of ions   : %1.0f \n',size(Pos3,1))
    fprintf('\n>> ZPE calculations \n')
end
if length(N) < 3; N = N*[1 1 1]; end
fid = fopen(filename1); geo = poscar(fid); fgetl(fid);
gridsize = fscanf(fid, '%d %d %d', [3 1])'; latt = geo.lattice;
loc = fscanf(fid, '%f', [prod(gridsize,2) 1])';
loc = reshape(loc,gridsize); fclose(fid);
SX = []; NX = size(loc);
for i = -1:1; SX = [SX; loc]; end; loc = SX; SX = [];
for i = -1:1; SX = [SX loc];  end; loc = SX; SX = [];
for i = -1:1
    if size(SX,3) == 1; SX(:,:,size(SX,3):size(SX,3)+NX(3)-1) = loc;
    else SX(:,:,size(SX,3)+1:size(SX,3)+NX(3)) = loc; end
end
loc = SX; clear SX; latt = inv(latt);
loc(end+1,:,:) = loc(1,:,:); loc(:,end+1,:) = loc(:,1,:); loc(:,:,end+1) = loc(:,:,1);
x = linspace(-1,2,size(loc,1));
y = linspace(-1,2,size(loc,2));
z = linspace(-1,2,size(loc,3));
muon = Mu(1)*latt(1,:) + Mu(2)*latt(2,:) + Mu(3)*latt(3,:);
dxyz = (Rz(1)*latt(1,:) + Rz(2)*latt(2,:) + Rz(3)*latt(3,:))/2;
xa = abs(x-muon(1) + dxyz(1)); [~,xa] = min(xa);
ya = abs(y-muon(2) + dxyz(2)); [~,ya] = min(ya);
za = abs(z-muon(3) + dxyz(3)); [~,za] = min(za);
xb = abs(x-muon(1) - dxyz(1)); [~,xb] = min(xb);
yb = abs(y-muon(2) - dxyz(2)); [~,yb] = min(yb);
zb = abs(z-muon(3) - dxyz(3)); [~,zb] = min(zb);
xi1 = linspace(x(xb),x(xa),N(1)); yi1 = linspace(y(yb),y(ya),N(2));
zi1 = linspace(z(zb),z(za),N(3)); [xi,yi,zi] = ndgrid(xi1,yi1,zi1);
xi = xi(:); yi = yi(:); zi = zi(:);  % interpolate
xi1 = xi*latt1(1,1) + yi*latt1(2,1) + zi*latt1(3,1);
yi1 = xi*latt1(1,2) + yi*latt1(2,2) + zi*latt1(3,2);
zi1 = xi*latt1(1,3) + yi*latt1(2,3) + zi*latt1(3,3);
xi = reshape(xi1,N); yi = reshape(yi1,N);
zi = reshape(zi1,N); [x,y,z] = ndgrid(x,y,z);
x = x(:); y = y(:); z = z(:);  % original
xi1 = x*latt1(1,1) + y*latt1(2,1) + z*latt1(3,1);
yi1 = x*latt1(1,2) + y*latt1(2,2) + z*latt1(3,2);
zi1 = x*latt1(1,3) + y*latt1(2,3) + z*latt1(3,3);
x = reshape(xi1,size(loc)); y = reshape(yi1,size(loc));
z = reshape(zi1,size(loc)); x = x(xa:xb, ya:yb, za:zb);
y = y(xa:xb, ya:yb, za:zb); z = z(xa:xb, ya:yb, za:zb);
Vi = interpn(x,y,z,-loc(xa:xb, ya:yb, za:zb),xi,yi,zi);
xi = xi(:); yi = yi(:); zi = zi(:);
if log == 1
    fprintf('     Hamiltonian matrix size: %1.0f x %1.0f \n',length(xi),length(xi))
end
m = 1.88353148E-28; % muon mass
h_bar = (6.62606957E-34)/(2*pi);
k = h_bar^2/(2*m)*6.24150934E38;
M = 1; L = 1; Vi = Vi(:); A = zeros(length(Vi));
dx = ((max(xi(:))-min(xi(:)))/N(1))^2; dx = k/dx;
dy = ((max(yi(:))-min(yi(:)))/N(2))^2; dy = k/dy;
dz = ((max(zi(:))-min(zi(:)))/N(3))^2; dz = k/dz;
for i = 1:length(Vi)
    if i > 1+(M-1)*N(1); A(i,i-1) = A(i,i-1) + dx; end
    if i < M*N(1); A(i,i+1) = A(i,i+1) + dx; end
    if i > N(1)*(N(2)*(L-1)+1); A(i,i-N(1)) = A(i,i-N(1)) + dy; end
    if i <= N(1)*(N(2)*L-1); A(i,i+N(1)) = A(i,i+N(1)) + dy; end
    if i > N(1)*N(2); A(i,i-N(1)*N(2)) = A(i,i-N(1)*N(2)) + dz; end
    if i <= N(1)*N(2)*(N(3)-1); A(i,i+N(1)*N(2)) = A(i,i+N(1)*N(2)) + dz; end
    A(i,i) = -(2*dx + 2*dy + 2*dz + Vi(i)); if i/N(1) == M; M = M+1; end
    if (M-1)/N(2) == L; L = L+1; end
end
if log == 1; fprintf('\n>> Solving Schrodinger equation \n'); end
[phi,Ezpe] = eig(A); phi = phi(:,end).*conj(phi(:,end));
Ezpe = (max(loc(:))-Ezpe(end))*1000;
if log == 1
    fprintf('     Zero-point Energy: %4.0f meV \n',Ezpe)
    fprintf('\n>> Dipole fields calculation \n')
    fprintf('     number of runs: %1.0f \n',length(xi));
end
Hx = zeros(length(xi),1);
if gpu == 1; xi = gpuArray(xi); yi = gpuArray(yi); zi = gpuArray(zi); 
    Sp3 = gpuArray(Sp3); Rmax = gpuArray(Rmax); Hx = gpuArray(Hx); 
end
Hy = Hx; Hz = Hx;
for i = 1:length(xi)
    x = (xi(i)-Pos3(:,1))*1E-8; y = (yi(i)-Pos3(:,2))*1E-8;
    z = (zi(i)-Pos3(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
    SR = Sp3(:,1).*x + Sp3(:,2).*y + Sp3(:,3).*z;
    aPos = (R > 0); SR = SR(aPos); R = R(aPos);
    x = x(aPos); y = y(aPos); z = z(aPos); Spin = Sp3(aPos,:);
    aPos = (R <= Rmax*1E-8); SR = SR(aPos); R = R(aPos);
    Hx(i) = sum((3*x(aPos).*SR - Spin(aPos,1).*R.^2)./R.^5); 
    Hy(i) = sum((3*y(aPos).*SR - Spin(aPos,2).*R.^2)./R.^5);
    Hz(i) = sum((3*z(aPos).*SR - Spin(aPos,3).*R.^2)./R.^5);
end
if gpu == 1; Hx = gather(Hx); Hy = gather(Hy); Hz = gather(Hz); Rmax = gather(Rmax); end
Hx = sum(phi.*Hx); Hy = sum(phi.*Hy); Hz = sum(phi.*Hz); 
H = sqrt(Hx^2+Hy^2+Hz^2); status = [H Hx Hy Hz];
if log == 1
    fprintf('     Dipole fields : %4.4f Gauss \n',H);
    fprintf('     [Hx, Hy, Hz]  : [%4.4f, %4.4f, %4.4f]\n',Hx,Hy,Hz);
    fprintf('\n>> Writing log file\n');
    fid = fopen([pathname '/log_zpe']);
    if fid ~= -1
        i = 1; fid = fopen([pathname '/log_zpe(1)']);
        while fid ~= -1
            i = i+1; fclose(fid); fid = fopen([pathname '/log_zpe(' num2str(i) ')']);
        end
        fid = fopen([pathname '/log_zpe(' num2str(i) ')'],'w');
        fprintf(['     ' pathname '/log_zpe(' num2str(i) ') \n'])
    else fid = fopen([pathname '/log_zpe'],'w');
        fprintf(['     ' pathname '/log_zpe \n'])
    end
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,'            Point Dipole with ZPE area \n');
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,['   Running on: ' datestr(jam) '\n']);
    fprintf(fid,'\n   POSCAR location: \n');
    fprintf(fid,['   ' filename1 '\n']); fprintf(fid,['   ' filename2 '\n']);
    latt = inv(latt1); Mu = Mu(1)*latt(1,:)+Mu(2)*latt(2,:)+Mu(3)*latt(3,:);
    fprintf(fid,'\n   Muon position : [%4.4f %4.4f %4.4f]\n',Mu);
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
    fprintf(fid,'   size of unit cell : [%1.0f %1.0f %1.0f] \n',2*dmax+1);
    fprintf(fid,'   number of ions    : %1.0f \n\n',size(Pos3,1));
    fprintf(fid,'   ZPE Range   : %4.4f Angstrom\n',sum(Rz)/3);
    fprintf(fid,'   Resolution  : %4.4f Angstrom\n',sum(Rz)/sum(N));
    fprintf(fid,'   Matrix size : %1.0f x %1.0f \n',length(xi),length(xi));
    fprintf(fid,'   ZPE-Energy  : %1.0f meV \n\n',Ezpe);
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,'   Dipole fields : %4.4f Gauss\n',H);
    fprintf(fid,'   [Hx, Hy, Hz]  : [%4.4f, %4.4f, %4.4f]\n',Hx,Hy,Hz);
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
end
