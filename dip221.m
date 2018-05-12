% point dipole vasp defect
% last edit 3 Mar 2016
function status = dip221(Mu,coord,Rmax,filename1,filename2,lim,gpu,log)
pathname = pwd;
if log == 1; clc; tic; jam = now;
    fprintf('   ============================================ \n')
    fprintf('            Dipole Fields Calculation \n')
    fprintf('   ============================================ \n')
    fprintf('>> Reading input files \n')
end
geo1 = poscar([filename1 '/CONTCAR']); latt1 = geo1.lattice;
geo2 = poscar([filename2 '/CONTCAR']); latt2 = geo2.lattice;
spin1 = outcar([filename1 '/OUTCAR']); pos1 = geo1.coords;
spin2 = outcar([filename2 '/OUTCAR']); pos2 = geo2.coords;
num = (1:size(spin1,1))'; 
spin = sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2);
if length(lim) == 1; cond = (spin(:,1) < lim);
    pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
else cond = (spin(:,1) < lim(2)); pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
     cond = (spin(:,1) > lim(1)); pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
end
pos1 = pos1(:,1)*latt1(1,:)+pos1(:,2)*latt1(2,:)+pos1(:,3)*latt1(3,:);
spinprt1 = [num spin1 sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2)]'; spin1 = spin1*9.274E-21;
num = (1:size(spin2,1))'; Pos = [];
spin = sqrt(spin2(:,1).^2+spin2(:,2).^2+spin2(:,3).^2);
if length(lim) == 1; cond = (spin(:,1) < lim);
    pos2(cond,:) = []; spin2(cond,:) = [];  num(cond,:) = [];
else cond = (spin(:,1) < lim(2)); pos2(cond,:) = []; spin2(cond,:) = []; num(cond,:) = [];
     cond = (spin(:,1) > lim(1)); pos2(cond,:) = []; spin2(cond,:) = []; num(cond,:) = [];
end
pos2 = pos2(:,1)*latt2(1,:)+pos2(:,2)*latt2(2,:)+pos2(:,3)*latt2(3,:);
spinprt2 = [num spin2 sqrt(spin2(:,1).^2+spin2(:,2).^2+spin2(:,3).^2)]'; spin2 = spin2*9.274E-21;
if coord == 3; Mu = geo1.coords(end,:); end
if coord ~= 2; Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:); end
dmax = ceil((Rmax+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
if log == 1
    fprintf('\n>> Expanding atomic positions \n')
    fprintf('     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
end
for i = -dmax(1):dmax(1)
    for j = -dmax(2):dmax(2)
        for k = -dmax(3):dmax(3)
            if i == 0 && j == 0 && k == 0; Pos = [Pos; pos1 spin1];
            else Pos = [Pos; pos2+repmat(i*latt2(1,:)+j*latt2(2,:)+k*latt2(3,:),size(pos2,1),1) spin2];
            end
        end
    end
end
spin = Pos(:,4:6); pos = Pos(:,1:3);
if log == 1
    fprintf('     number of ions   : %1.0f \n',size(Pos,1))
    fprintf('\n>> Dipole fields calculation \n')
end
Pos = repmat(Mu,size(pos,1),1)-pos;
if gpu == 1; spin = gpuArray(spin); Pos = gpuArray(Pos); Rmax = gpuArray(Rmax); end
R = sqrt(Pos(:,1).^2+Pos(:,2).^2+Pos(:,3).^2);
Pos(R>Rmax,:) = []; spin(R>Rmax,:) = []; R(R>Rmax) = [];
Pos(R<1E-6,:) = []; spin(R<1E-6,:) = []; R(R<1E-6) = [];
Pos = Pos*1E-8; R = R*1E-8; % to cm
SR = spin(:,1).*Pos(:,1) + spin(:,2).*Pos(:,2) + spin(:,3).*Pos(:,3);
Hx = sum((3*Pos(:,1).*SR - spin(:,1).*R.^2)./R.^5); 
Hy = sum((3*Pos(:,2).*SR - spin(:,2).*R.^2)./R.^5);
Hz = sum((3*Pos(:,3).*SR - spin(:,3).*R.^2)./R.^5);
if gpu == 1; Hx = gather(Hx); Hy = gather(Hy); Hz = gather(Hz); Rmax = gather(Rmax); end
H = sqrt(Hx.^2+Hy.^2+Hz.^2); status = [H Hx Hy Hz];
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
    fprintf(fid,'               Point Dipole (VASP) \n');
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,['   Running on: ' datestr(jam) '\n']);
    fprintf(fid,'\n   POSCAR location: \n');
    fprintf(fid,['   ' filename1 '\n']); fprintf(fid,['   ' filename2 '\n']);
    latt = inv(latt1); Mu = Mu(1)*latt(1,:)+Mu(2)*latt(2,:)+Mu(3)*latt(3,:);
    fprintf(fid,'\n   Muon position : [%4.4f %4.4f %4.4f]\n',Mu);
    fprintf(fid,'\n     No      Mx      My      Mz      Mtotal');
    fprintf(fid,'\n   --------------------------------------------\n');
    fprintf(fid,'   %4.0f    %6.3f  %6.3f  %6.3f    %6.3f \n',spinprt1(:));
    fprintf(fid,'   --------------------------------------------\n');
    fprintf(fid,'    sum    %6.3f  %6.3f  %6.3f    %6.3f \n',sum(spinprt1(2,:)), ...
        sum(spinprt1(3,:)),sum(spinprt1(4,:)),sum(spinprt1(5,:)));
    fprintf(fid,'\n     No2     Mx      My      Mz      Mtotal');
    fprintf(fid,'\n   --------------------------------------------\n');
    fprintf(fid,'   %4.0f    %6.3f  %6.3f  %6.3f    %6.3f \n',spinprt2(:));
    fprintf(fid,'   --------------------------------------------\n');
    fprintf(fid,'    sum    %6.3f  %6.3f  %6.3f    %6.3f \n',sum(spinprt2(2,:)), ...
        sum(spinprt2(3,:)),sum(spinprt2(4,:)),sum(spinprt2(5,:)));
    fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',Rmax);
    fprintf(fid,'   size of unit cell : [%1.0f %1.0f %1.0f] \n',2*dmax+1);
    fprintf(fid,'   number of ions    : %1.0f \n\n',size(Pos,1));
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
