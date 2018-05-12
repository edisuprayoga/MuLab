% spin density + ZPE
% last edit 16 Mar 2016
function status = dip232(Mu,coord,Rmax,filename1,filename2,Rz,N,gpu,log)
pathname = pwd;
if log == 1; clc; tic; jam = now;
    fprintf('   ============================================ \n')
    fprintf('            Dipole Fields Calculation \n')
    fprintf('   ============================================ \n')
    fprintf('>> Reading input files \n')
end
[~,Sx1,Sy1,Sz1,geo1] = chgcar([filename1 '/CHGCAR']);
latt1 = geo1.lattice; vol = abs(dot(latt1(1,:),cross(latt1(2,:),latt1(3,:))));
sumx = sum(abs(Sx1(:)))^2+sum(abs(Sy1(:)))^2+sum(abs(Sz1(:)));
if sumx > vol
    Sx1 = Sx1/vol; Sy1 = Sy1/vol; Sz1 = Sz1/vol;
end
sumx = sum(abs(Sx1(:)))^2+sum(abs(Sy1(:)))^2+sum(abs(Sz1(:)));
if sumx > vol
    Sx1 = Sx1/vol; Sy1 = Sy1/vol; Sz1 = Sz1/vol;
end
sumx = [sum(Sx1(:)) sum(Sy1(:)) sum(Sz1(:))];
dt = min((latt1(1,:) + latt1(2,:) + latt1(3,:))./(size(Sx1)+1))/2;
if coord == 3; Mu = geo.coords(end,:); end
if coord ~= 2; Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:); end
Rmax = Rmax + sqrt(Rz(1)^2+Rz(2)^2+Rz(3)^2)/2; 
dmax = ceil((Rmax+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
ngrid = size(Sx1).*(2*dmax+1);
xg = linspace(0,1,size(Sx1,1)+1); xg(end) = [];
yg = linspace(0,1,size(Sx1,2)+1); yg(end) = [];
zg = linspace(0,1,size(Sx1,3)+1); zg(end) = [];
[xg,yg,zg] = ndgrid(xg,yg,zg); xg = xg(:); yg = yg(:); zg = zg(:);
Sx1 = Sx1(:)*9.274E-21; Sy1 = Sy1(:)*9.274E-21; Sz1 = Sz1(:)*9.274E-21;
[~,Sx2,Sy2,Sz2,geo2] = chgcar(filename2);
latt2 = geo2.lattice; vol = abs(dot(latt2(1,:),cross(latt2(2,:),latt2(3,:))));
sumy = sum(abs(Sx2(:)))^2+sum(abs(Sy2(:)))^2+sum(abs(Sz2(:)));
if sumy > vol
    Sx2 = Sx2/vol; Sy2 = Sy2/vol; Sz2 = Sz2/vol;
end
sumy = sum(abs(Sx2(:)))^2+sum(abs(Sy2(:)))^2+sum(abs(Sz2(:)));
if sumy > vol
    Sx2 = Sx2/vol; Sy2 = Sy2/vol; Sz2 = Sz2/vol;
end
sumy = [sum(Sx2(:)) sum(Sy2(:)) sum(Sz2(:))];
Sx2 = Sx2(:)*9.274E-21; Sy2 = Sy2(:)*9.274E-21; Sz2 = Sz2(:)*9.274E-21;
if log == 1; 
    fprintf('     Magmom1 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
    if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
        sumx1 = [sum(abs(Sx1(:))) sum(abs(Sy1(:))) sum(abs(Sz1(:)))]/9.274E-21;
        fprintf('     M1(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
    end
    fprintf('     Magmom2 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2),sumy);
    if sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2) < 1
        sumy1 = [sum(abs(Sx2(:))) sum(abs(Sy2(:))) sum(abs(Sz2(:)))]/9.274E-21;
        fprintf('     M2(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy1(1)^2+sumy1(2)^2+sumy1(3)^2),sumy1);
    end
    fprintf('\n     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
    fprintf('     number of grids  : [%1.0f %1.0f %1.0f]\n\n',ngrid)
    fprintf('>> ZPE calculations \n')
end
if length(N) < 3; N = N*[1 1 1]; end
fid = fopen([filename1 '/LOCPOT']); geo1 = poscar(fid); fgetl(fid);
gridsize = fscanf(fid, '%d %d %d', [3 1])'; latt1 = geo1.lattice;
loc = fscanf(fid, '%f', [prod(gridsize,2) 1])';
loc = reshape(loc,gridsize); fclose(fid);
SX = []; NX = size(loc);
for i = -1:1; SX = [SX; loc]; end; loc = SX; SX = [];
for i = -1:1; SX = [SX loc];  end; loc = SX; SX = [];
for i = -1:1
    if size(SX,3) == 1; SX(:,:,size(SX,3):size(SX,3)+NX(3)-1) = loc;
    else SX(:,:,size(SX,3)+1:size(SX,3)+NX(3)) = loc; end
end
loc = SX; clear SX; latt = inv(latt1);
loc(end+1,:,:) = loc(1,:,:); loc(:,end+1,:) = loc(:,1,:); 
loc(:,:,end+1) = loc(:,:,1); x = linspace(-1,2,size(loc,1));
y = linspace(-1,2,size(loc,2)); z = linspace(-1,2,size(loc,3));
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
[phi,Ezpe] = eig(A); Ezpe = (max(loc(:))-Ezpe(end))*1000;
if log == 1
    fprintf('     Zero-point Energy: %4.0f meV \n',Ezpe)
    fprintf('\n>> Dipole fields calculation \n')
    fprintf('     number of runs: %1.0f \n',length(xi));
end
phi = phi(:,end).*conj(phi(:,end));
HX = zeros(length(xi),1); 
if gpu == 1; xi = gpuArray(xi); yi = gpuArray(yi); zi = gpuArray(zi); 
    xg = gpuArray(xg); yg = gpuArray(yg); zg = gpuArray(zg); 
    Sx1 = gpuArray(Sx1); Sy1 = gpuArray(Sy1); Sz1 = gpuArray(Sz1);
    Sx2 = gpuArray(Sx2); Sy2 = gpuArray(Sy2); Sz2 = gpuArray(Sz2);
    Rmax = gpuArray(Rmax); HX = gpuArray(HX); 
end
HY = HX; HZ = HX;
for n = 1:length(xi); Hx = 0; Hy = 0; Hz = 0;
    for i = -dmax(1):dmax(1); X = xg+i;
        for j = -dmax(2):dmax(2); Y = yg+j;
            for k = -dmax(3):dmax(3); Z = zg+k;
                xr = (xi(n) - X*latt1(1,1) - Y*latt1(2,1) - Z*latt1(3,1))*1E-8;
                yr = (yi(n) - X*latt1(1,2) - Y*latt1(2,2) - Z*latt1(3,2))*1E-8;
                zr = (zi(n) - X*latt1(1,3) - Y*latt1(2,3) - Z*latt1(3,3))*1E-8;
                if i == 0 && j == 0 && k == 0; Sx = Sx1; Sy = Sy1; Sz = Sz1;
                else Sx = Sx2; Sy = Sy2; Sz = Sz2; 
                end
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
    HX(n) = Hx; HY(n) = Hy; HZ(n) = Hz;
end
if gpu == 1; HX = gather(HX); HY = gather(HY); HZ = gather(HZ); Rmax = gather(Rmax); end
Hx = sum(phi.*HX); Hy = sum(phi.*HY); Hz = sum(phi.*HZ); 
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
    fprintf(fid,'            Spin Density Calculation \n');
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,['   Running on: ' datestr(jam) '\n']);
    fprintf(fid,'\n   POSCAR location: \n');
    fprintf(fid,['   ' filename1 '\n']);
    Mu = Mu(1)*latt(1,:)+Mu(2)*latt(2,:)+Mu(3)*latt(3,:);
    fprintf(fid,'\n   Magmom1 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
    if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
        fprintf(fid,'     M1(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
    end
    fprintf(fid,'   Magmom2 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2),sumy);
    if sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2) < 1
        fprintf(fid,'     M2(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy1(1)^2+sumy1(2)^2+sumy1(3)^2),sumy1);
    end
    fprintf(fid,'\n   Muon position : [%4.4f %4.4f %4.4f]\n',Mu);
    fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',Rmax);
    fprintf(fid,'   size of unit cell : [%1.0f %1.0f %1.0f] \n',2*dmax+1);
    fprintf(fid,'   number of grids   : [%1.0f %1.0f %1.0f] \n',ngrid);
    fprintf(fid,'\n   ZPE Range   : %4.4f Angstrom\n',sum(Rz)/3);
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
