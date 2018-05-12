% spin density
% last edit 16 Mar 2016
function status = dip131(Mu,coord,Rmax,filename,gpu,log)
if log == 1; pathname = pwd; tic; jam = now;
    fprintf('   ============================================ \n')
    fprintf('            Dipole Fields Calculation \n')
    fprintf('   ============================================ \n')
    fprintf('>> Reading input files \n')
end
if isempty(filename); filename = 'CHGCAR'; end
[~,Sx,Sy,Sz,geo] = chgcar(filename);
latt = geo.lattice; vol = abs(dot(latt(1,:),cross(latt(2,:),latt(3,:))));
dt = min((latt(1,:) + latt(2,:) + latt(3,:))./(size(Sx)+1))/2;
sumx = sum(abs(Sx(:)))^2+sum(abs(Sy(:)))^2+sum(abs(Sz(:)));
if sumx > vol
    Sx = Sx/vol; Sy = Sy/vol; Sz = Sz/vol;
end
sumx = sum(abs(Sx(:)))^2+sum(abs(Sy(:)))^2+sum(abs(Sz(:)));
if sumx > vol
    Sx = Sx/vol; Sy = Sy/vol; Sz = Sz/vol;
end
sumx = [sum(Sx(:)) sum(Sy(:)) sum(Sz(:))];
if log == 1
    fprintf('     Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
    if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
        sumx1 = [sum(abs(Sx(:))) sum(abs(Sy(:))) sum(abs(Sz(:)))];
        fprintf('     M(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
    end
end
if coord == 3; Mu = geo.coords(end,:); end
if coord ~= 2; Mu = Mu(1)*latt(1,:)+Mu(2)*latt(2,:)+Mu(3)*latt(3,:); end
dmax = ceil((Rmax+Mu)./(latt(1,:)+latt(2,:)+latt(3,:)));
ngrid = size(Sx).*(2*dmax+1);
if log == 1; 
    fprintf('\n     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
    fprintf('     number of grids  : [%1.0f %1.0f %1.0f]\n\n',ngrid)
    fprintf('>> Dipole fields calculation \n')
end
x = linspace(0,1,size(Sx,1)+1); x(end) = [];
y = linspace(0,1,size(Sx,2)+1); y(end) = [];
z = linspace(0,1,size(Sx,3)+1); z(end) = [];
[x,y,z] = ndgrid(x,y,z); x = x(:); y = y(:); z = z(:);
Sx = Sx(:)*9.274E-21; Sy = Sy(:)*9.274E-21; Sz = Sz(:)*9.274E-21;
Hx = 0; 
if gpu == 1; x = gpuArray(x); y = gpuArray(y); z = gpuArray(z);
    Sx = gpuArray(Sx); Sy = gpuArray(Sy); Sz = gpuArray(Sz);
    dt = gpuArray(dt); Rmax = gpuArray(Rmax); Hx = gpuArray(Hx);
end
Hy = Hx; Hz = Hx;
for i = -dmax(1):dmax(1); X = x+i;
    for j = -dmax(2):dmax(2); Y = y+j;
        for k = -dmax(3):dmax(3); Z = z+k;
            xr = (Mu(1) - X*latt(1,1) - Y*latt(2,1) - Z*latt(3,1))*1E-8;
            yr = (Mu(2) - X*latt(1,2) - Y*latt(2,2) - Z*latt(3,2))*1E-8;
            zr = (Mu(3) - X*latt(1,3) - Y*latt(2,3) - Z*latt(3,3))*1E-8;
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
if gpu == 1; Hx = gather(Hx); Hy = gather(Hy); Hz = gather(Hz); Rmax = gather(Rmax); end
H = sqrt(Hx^2+Hy^2+Hz^2); status = [H Hx Hy Hz];
if log == 1
    fprintf('     Dipole fields : %4.4f Gauss \n',H);
    fprintf('     [Hx, Hy, Hz]  : [%4.4f, %4.4f, %4.4f]\n',Hx,Hy,Hz);
    fprintf('\n>> Writing log file\n');
    fid = fopen([pathname '/log']);
    if fid ~= -1
        i = 1; fid = fopen([pathname '/log(1)']);
        while fid ~= -1
            i = i+1; fclose(fid); fid = fopen([pathname '/log(' num2str(i) ')']);
        end
        fid = fopen([pathname '/log(' num2str(i) ')'],'w');
        fprintf(['     ' pathname '/log(' num2str(i) ') \n'])
    else fid = fopen([pathname '/log'],'w');
        fprintf(['     ' pathname '/log \n'])
    end
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,'            Spin Density Calculation \n');
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,['   Running on: ' datestr(jam) '\n']);
    fprintf(fid,'\n   POSCAR location: \n');
    fprintf(fid,['   ' filename '\n']);
    latt = inv(latt); Mu = Mu(1)*latt(1,:)+Mu(2)*latt(2,:)+Mu(3)*latt(3,:);
    fprintf(fid,'\n   Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
    if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
        fprintf(fid,'     M(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
    end
    fprintf(fid,'\n   Muon position : [%4.4f %4.4f %4.4f]\n',Mu);
    fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',Rmax);
    fprintf(fid,'   size of unit cell : [%1.0f %1.0f %1.0f] \n',2*dmax+1);
    fprintf(fid,'   number of grids   : [%1.0f %1.0f %1.0f] \n\n',ngrid);
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
