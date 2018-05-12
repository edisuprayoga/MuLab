% Test point dipole manual defect
% last edit 3 Mar 2016
function status = dip213(Mu,coord,Rmax,filename1,filename2,No,SpinVal1,Spin1,SpinVal2,Spin2,dmin,del,gpu,log)
pathname = pwd; jam = now;
if log == 1; tic;
    fprintf('   ============================================ \n')
    fprintf('            Dipole Fields Calculation \n')
    fprintf('   ============================================ \n')
    fprintf('>> Reading input files \n')
end
geo = poscar(filename1); latt1 = geo.lattice;
geo2 = poscar(filename2); latt2 = geo2.lattice;
if coord == 3; Mu = geo.coords(end,:); end
if coord ~= 2; Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:); end
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
spinprt2 = Spin2; 
for i = 1:3:length(Spin2)
    norm = sqrt(Spin2(i)^2+Spin2(i+1)^2+Spin2(i+2)^2);
    if norm > 1; Spin2(i) = Spin2(i)/norm;
        Spin2(i+1) = Spin2(i+1)/norm; Spin2(i+2) = Spin2(i+2)/norm;
    end
end
Spin2 = Spin2*9.274E-21; cu2 = [];
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
if dmin < 5; dmin = 5; end; del = linspace(dmin,Rmax,del);
if log == 1
    fprintf('     number of ions   : %1.0f \n',size(Pos3,1))
    fprintf('\n>> Dipole fields calculation \n')
    fprintf('     number of runs: %1.0f \n',length(del))
end
x = (Mu(1)-Pos3(:,1))*1E-8; y = (Mu(2)-Pos3(:,2))*1E-8;
z = (Mu(3)-Pos3(:,3))*1E-8; Hx = zeros(length(del),1); 
if gpu == 1; x = gpuArray(x); y = gpuArray(y); z = gpuArray(z); 
    Sp3 = gpuArray(Sp3); del = gpuArray(del); Hx = gpuArray(Hx);
end
Hy = Hx; Hz = Hx; r = sqrt(x.^2+y.^2+z.^2); 
sr = Sp3(:,1).*x + Sp3(:,2).*y + Sp3(:,3).*z;
for i = 1:length(del);
    aPos = (r > 0); SR = sr(aPos); R = r(aPos);
    X = x(aPos); Y = y(aPos); Z = z(aPos); Spin = Sp3(aPos,:);
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
fprintf(fid,'            Point Dipole Calculation \n');
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
fprintf(fid,'   R_dip(A)  H(G)      [Hx, Hy, Hz] \n');
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
