% Test spin density
% last edit 15 Mar 2016
function status = dip133(Mu,coord,Rmax,filename1,dmin,del,gpu,log)
pathname = pwd; jam = now;
if log == 1; tic;
    fprintf('   ============================================ \n')
    fprintf('            Dipole Fields Calculation \n')
    fprintf('   ============================================ \n')
    fprintf('>> Reading input files \n')
end
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
dt = min((latt1(1,:) + latt1(2,:) + latt1(3,:))./(size(Sx)+1))/2;
if coord == 3; Mu = geo.coords(end,:); end
if coord ~= 2; Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:); end
dmax = ceil((Rmax+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
ngrid = size(Sx).*(2*dmax+1);

xg = linspace(0,1,size(Sx,1)+1); xg(end) = [];
yg = linspace(0,1,size(Sx,2)+1); yg(end) = [];
zg = linspace(0,1,size(Sx,3)+1); zg(end) = [];
[xg,yg,zg] = ndgrid(xg,yg,zg); xg = xg(:); yg = yg(:); zg = zg(:);

Sx = Sx(:)*9.274E-21; Sy = Sy(:)*9.274E-21; Sz = Sz(:)*9.274E-21;
if dmin < 5; dmin = 5; end; del = linspace(dmin,Rmax,del);
if log == 1;
    fprintf('     Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
    if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
        sumx1 = [sum(abs(Sx(:))) sum(abs(Sy(:))) sum(abs(Sz(:)))]/9.274E-21;
        fprintf('     M(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
    end
    fprintf('\n     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
    fprintf('     number of grids  : [%1.0f %1.0f %1.0f]\n\n',ngrid)
    fprintf('>> Dipole fields calculation \n')
    fprintf('     number of runs: %1.0f \n',length(del)); 
end
HX = zeros(length(del),1); 
if gpu == 1; xg = gpuArray(xg); yg = gpuArray(yg); zg = gpuArray(zg);
    Sx = gpuArray(Sx); Sy = gpuArray(Sy); Sz = gpuArray(Sz);
    Rmax = gpuArray(del); HX = gpuArray(HX); 
end
HY = HX; HZ = HX;
for n = 1:length(del); Hx = 0; Hy = 0; Hz = 0;
    dmax = ceil((del(n)+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
    for i = -dmax(1):dmax(1); X = xg+i;
        for j = -dmax(2):dmax(2); Y = yg+j;
            for k = -dmax(3):dmax(3); Z = zg+k;
                xr = (Mu(1) - X*latt1(1,1) - Y*latt1(2,1) - Z*latt1(3,1))*1E-8;
                yr = (Mu(2) - X*latt1(1,2) - Y*latt1(2,2) - Z*latt1(3,2))*1E-8;
                zr = (Mu(3) - X*latt1(1,3) - Y*latt1(2,3) - Z*latt1(3,3))*1E-8;
                R = sqrt(xr.^2 + yr.^2 + zr.^2); SR = Sx.*xr + Sy.*yr + Sz.*zr;
                aPos = (R > dt*1E-8); R = R(aPos); SR = SR(aPos);
                xr = xr(aPos); yr = yr(aPos); zr = zr(aPos);
                SX = Sx(aPos); SY = Sy(aPos); SZ = Sz(aPos);
                aPos = (R <= Rmax(n)*1E-8); R = R(aPos);
                Hx = Hx + sum((3*xr(aPos).*SR(aPos) - SX(aPos).*R.^2)./R.^5);
                Hy = Hy + sum((3*yr(aPos).*SR(aPos) - SY(aPos).*R.^2)./R.^5);
                Hz = Hz + sum((3*zr(aPos).*SR(aPos) - SZ(aPos).*R.^2)./R.^5);
            end
        end
    end
    HX(n) = Hx; HY(n) = Hy; HZ(n) = Hz;
end
if gpu == 1; HX = gather(HX); HY = gather(HY); HZ = gather(HZ); Rmax = gather(Rmax); end
H = sqrt(HX.^2+HY.^2+HZ.^2); out = [del' H HX HY HZ]'; 
status = [H(end) HX(end) HY(end) HZ(end)];
if log == 1
    fprintf('     Dipole fields : %4.4f Gauss \n',H(end));
    fprintf('     [Hx, Hy, Hz]  : [%4.4f, %4.4f, %4.4f]\n',HX(end),HY(end),HZ(end));
    fprintf('\n>> Writing log file\n');
    fid = fopen([pathname '/log_test']);
    if fid ~= -1
        i = 1; fid = fopen([pathname '/log_test(1)']);
        while fid ~= -1
            i = i+1; fclose(fid); fid = fopen([pathname '/log_test(' num2str(i) ')']);
        end
        fid = fopen([pathname '/log_test(' num2str(i) ')'],'w');
        if log == 1; fprintf(['     ' pathname '/log_test(' num2str(i) ') \n']); end
    else fid = fopen([pathname '/log_test'],'w');
        if log == 1; fprintf(['     ' pathname '/log_test \n']); end
    end
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,'          Spin Density Calculation test\n');
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,['   Running on: ' datestr(jam) '\n']);
    fprintf(fid,'\n   POSCAR location: \n');
    fprintf(fid,['   ' filename1 '\n']);
    latt1 = inv(latt1); Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:);
    fprintf(fid,'\n   Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
    if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
        fprintf(fid,'     M(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
    end
    fprintf(fid,'\n   Muon position : [%4.4f %4.4f %4.4f]\n',Mu);
    fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',Rmax(end));
    fprintf(fid,'   size of unit cell : [%1.0f %1.0f %1.0f] \n',2*dmax+1);
    fprintf(fid,'   number of grids   : [%1.0f %1.0f %1.0f] \n',ngrid);
    fprintf(fid,'\n   R_dip(A)  H(G)        [Hx, Hy, Hz] \n');
    fprintf(fid,'   ----------------------------------------------------\n');
    fprintf(fid,' %8.4f  %8.4f    [%8.4f  %8.4f  %8.4f] \n',out(:));
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
end
