% Test point dipole manual
% last edit 3 Mar 2016
function status = dip113(Mu,coord,Rmax,filename1,No,SpinVal,Spin,dmin,del,gpu,log)
pathname = pwd; jam = now;
if log == 1; tic; 
    fprintf('   ============================================ \n')
    fprintf('            Dipole Fields Calculation \n')
    fprintf('   ============================================ \n')
    fprintf('>> Reading input files \n')
end
if isempty(filename1); filename1 = 'POSCAR'; end
geo = poscar(filename1); latt = geo.lattice;
if coord == 3; Mu = geo.coords(end,:); end
if coord ~= 2; Mu = Mu(1)*latt(1,:)+Mu(2)*latt(2,:)+Mu(3)*latt(3,:); end
dmax = ceil((Rmax+Mu)./(latt(1,:)+latt(2,:)+latt(3,:)));
if log == 1 
    fprintf('\n>> Expanding atomic positions \n')
    fprintf('     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
end
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
            for k = -dmax(3):dmax(3); C = cu + i*latt(1,:)+j*latt(2,:)+k*latt(3,:);
                if sqrt(sum((C-Mu).*(C-Mu))) <= Rmax; Pos = [Pos; C Spin(n,:)]; end
            end
        end
    end
end
Sp1 = Pos(:,4:6); Pos = Pos(:,1:3);
if dmin < 5; dmin = 5; end; del = linspace(dmin,Rmax,del);
if log == 1 
    fprintf('     number of ions   : %1.0f \n',size(Pos,1))
    fprintf('\n>> Dipole fields calculation \n')
    fprintf('     number of runs: %1.0f \n',length(del));
end
x = (Mu(1)-Pos(:,1))*1E-8; y = (Mu(2)-Pos(:,2))*1E-8;
z = (Mu(3)-Pos(:,3))*1E-8; Hx = zeros(length(del),1); 
if gpu == 1; x = gpuArray(x); y = gpuArray(y); z = gpuArray(z); 
    Sp1 = gpuArray(Sp1); del = gpuArray(del); Hx = gpuArray(Hx);
end
Hy = Hx; Hz = Hx; r = sqrt(x.^2+y.^2+z.^2); 
sr = Sp1(:,1).*x + Sp1(:,2).*y + Sp1(:,3).*z;
for i = 1:length(del);
    aPos = (r > 0); SR = sr(aPos); R = r(aPos);
    X = x(aPos); Y = y(aPos); Z = z(aPos); Spin = Sp1(aPos,:);
    aPos = (r <= del(i)*1E-8); SR = SR(aPos); R = R(aPos);
    Hx(i) = sum((3*X(aPos).*SR - Spin(aPos,1).*R.^2)./R.^5); 
    Hy(i) = sum((3*Y(aPos).*SR - Spin(aPos,2).*R.^2)./R.^5);
    Hz(i) = sum((3*Z(aPos).*SR - Spin(aPos,3).*R.^2)./R.^5);
end
if gpu == 1; Hx = gather(Hx); Hy = gather(Hy); Hz = gather(Hz); del = gather(del); end
H = sqrt(Hx.^2+Hy.^2+Hz.^2); out = [del' H Hx Hy Hz]';
status = [H(end) Hx(end) Hy(end) Hz(end)];
if log == 1
    fprintf('     Dipole fields : %4.4f Gauss \n',H(end));
    fprintf('     [Hx, Hy, Hz]  : [%4.4f, %4.4f, %4.4f]\n',Hx(end),Hy(end),Hz(end));
    fprintf('\n>> Writing log file\n');
end
fid = fopen([pathname '/log_test']);
if fid ~= -1
    i = 1; fid = fopen([pathname '/log_test(1)']);
    while fid ~= -1; i = i+1; fclose(fid); 
        fid = fopen([pathname '/log_test(' num2str(i) ')']);
    end
    fid = fopen([pathname '/log_test(' num2str(i) ')'],'w');
    if log == 1; fprintf(['     ' pathname '/log_test(' num2str(i) ') \n']); end
else fid = fopen([pathname '/log_test'],'w');
    if log == 1; fprintf(['     ' pathname '/log_test \n']); end
end
fprintf(fid,'   ============================================ \n');
fprintf(fid,'       Point Dipole Calculation Range Test \n');
fprintf(fid,'   ============================================ \n');
fprintf(fid,['   Running on: ' datestr(jam) '\n']);
fprintf(fid,'\n   POSCAR location: \n'); fprintf(fid,['   ' filename1 '\n']);
latt = inv(latt); Mu = Mu(1)*latt(1,:)+Mu(2)*latt(2,:)+Mu(3)*latt(3,:);
fprintf(fid,'\n   Muon position : [%4.4f %4.4f %4.4f]\n',Mu);
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
fprintf(fid,'   size of unit cell : [%1.0f %1.0f %1.0f] \n',2*dmax+1);
fprintf(fid,'   number of ions    : %1.0f \n\n',size(Pos,1));
fprintf(fid,'\n   R_dip(A)  H(G)      [Hx, Hy, Hz] \n');
fprintf(fid,'   ----------------------------------------------------\n');
fprintf(fid,' %8.4f  %8.4f    [%8.4f  %8.4f  %8.4f] \n',out(:));
fprintf(fid,'\n   ============================================ \n');
fprintf(fid,'             Calculations Completed! \n');
fprintf(fid,'   ============================================ \n');
t = toc; h = floor(t/3600); m = floor((t-h*3600)/60); t = t-h*3600 - m*60;
if h >= 1; fprintf(fid,'   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n\n',h,m,t);  
elseif m >= 1; fprintf(fid,'   Elapsed time is %1.0f min %1.4f sec.\n\n',m,t);
else fprintf(fid,'   Elapsed time is %1.4f sec.\n\n',t);
end; fclose(fid); 
if log == 1
    fprintf('\n   ============================================ \n')
    fprintf('             Calculations Completed! \n')
    fprintf('   ============================================ \n')
    if h >= 1; fprintf('   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n',h,m,t);  
    elseif m >= 1; fprintf('   Elapsed time is %1.0f min %1.4f sec.\n',m,t);
    else fprintf('   Elapsed time is %1.4f sec.\n',t);
    end
end
