function varargout = GUI_Dipole_Run(varargin)

% GUI_DIPOLE_ALL MATLAB code by Edi Suprayoga
% feel free to contact me at suprayoga.edi@gmail.com
% last edit 17 Mar 2016

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Dipole_Run_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Dipole_Run_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout; [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else gui_mainfcn(gui_State, varargin{:});
end

function GUI_Dipole_Run_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject; guidata(hObject, handles);
set(handles.figure1,'Color',[0.8 0.8 0.8])

function varargout = GUI_Dipole_Run_OutputFcn(~, ~, handles) 
varargout{1} = handles.output; set(handles.autorized,'Value',1);

function run_Callback(~, ~, handles)
pathname = pwd; jam = now; clc; tic
fprintf('   ============================================ \n')
fprintf('            Dipole Fields Calculation \n')
fprintf('   ============================================ \n')
fprintf('>> Reading input files \n')
opt1 = get(handles.opt1,'Value'); opt2 = get(handles.opt2,'Value');
opt3 = get(handles.opt3,'Value'); opt_job = get(handles.opt_job,'Value');
filename1 = get(handles.pos1,'String');
if opt1 == 1; % manual
    geo1 = poscar(filename1); latt1 = geo1.lattice; 
    SpinVal1 = str2num(get(handles.S1,'String'));
    No = str2num(get(handles.opt_in,'String'));
    fid = fopen(get(handles.mag1,'String'));
    if fid == -1; Spin1 = str2num(get(handles.mag1,'String'));
    else Spin1 = str2num(fgetl(fid));
        while ~feof(fid); Spin1 = [Spin1 str2num(fgetl(fid))]; end; fclose(fid);
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
        cu1 = [cu1; geo1.coords(sum(geo1.atomcount(1:No(i)-1))+1:sum(geo1.atomcount(1:No(i))),:)];
    end
    if length(Spin1) == 3 % collinear
        mag = zeros(size(cu1,1),3); for i = 1:size(cu1,1); mag(i,:) = Spin1; end; Spin1 = mag;
    else Spin1 = reshape(Spin1,3,size(cu1,1))';
    end
    a=1; b=0;
    for i = 1:length(No)
        b = b+geo1.atomcount(No(i));
        Spin1(a:b,:) = Spin1(a:b,:)*SpinVal1(i); a = a+b;
    end
    fprintf('\n>> Expanding atomic positions \n')
elseif opt2 == 1; % vasp
    geo1 = poscar([filename1 '/CONTCAR']);
    spin1 = outcar([filename1 '/OUTCAR']); num = (1:size(spin1,1))';
    pos1 = geo1.coords; latt1 = geo1.lattice;
    spin = sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2);
    lim = str2num(get(handles.opt_in,'String'));
    if length(lim) == 1; cond = (spin(:,1) < lim);
        pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
    else cond = (spin(:,1) < lim(2)); pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
         cond = (spin(:,1) > lim(1)); pos1(cond,:) = []; spin1(cond,:) = []; num(cond,:) = [];
    end
    pos1 = pos1(:,1)*latt1(1,:)+pos1(:,2)*latt1(2,:)+pos1(:,3)*latt1(3,:);
    spinprt1 = [num spin1 sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2)]'; spin1 = spin1*9.274E-21;
    fprintf('\n>> Expanding atomic positions \n')
elseif opt3 == 1; % chgcar
    if opt_job == 2; [~,Sx,Sy,Sz,geo1] = chgcar([filename1 '/CHGCAR']);
    else [~,Sx,Sy,Sz,geo1] = chgcar(filename1);
    end
    latt1 = geo1.lattice; vol = abs(dot(latt1(1,:),cross(latt1(2,:),latt1(3,:))));
    dt = min((latt1(1,:) + latt1(2,:) + latt1(3,:))./(size(Sx)+1))/2;
    sumx = sum(abs(Sx(:)))^2+sum(abs(Sy(:)))^2+sum(abs(Sz(:)));
    if sumx > vol
        Sx = Sx/vol; Sy = Sy/vol; Sz = Sz/vol;
    end
    sumx = sum(abs(Sx(:)))^2+sum(abs(Sy(:)))^2+sum(abs(Sz(:)));
    if sumx > vol
        Sx = Sx/vol; Sy = Sy/vol; Sz = Sz/vol;
    end
    sumx = [sum(Sx(:)) sum(Sy(:)) sum(Sz(:))]; ngrid = size(Sx);
    Sx = Sx(:)*9.274E-21; Sy = Sy(:)*9.274E-21; Sz = Sz(:)*9.274E-21;
    if get(handles.opt_pos,'Value') == 1
        fprintf('     Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
        if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
            sumx1 = [sum(abs(Sx(:))) sum(abs(Sy(:))) sum(abs(Sz(:)))]/9.274E-21;
            fprintf('     M(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
        end
    else fprintf('     Magmom1 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
        if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
            sumx1 = [sum(abs(Sx(:))) sum(abs(Sy(:))) sum(abs(Sz(:)))]/9.274E-21;
            fprintf('     M1(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
        end
        filename2 = get(handles.pos2,'String');
        [~,Sx2,Sy2,Sz2,geo2] = chgcar(filename2); Sx1 = Sx; Sy1 = Sy; Sz1 = Sz;
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
        fprintf('     Magmom2 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2),sumy);
        if sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2) < 1
            sumy1 = [sum(abs(Sx2(:))) sum(abs(Sy2(:))) sum(abs(Sz2(:)))];
            fprintf('     M2(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy1(1)^2+sumy1(2)^2+sumy1(3)^2),sumy1);
        end
        Sx2 = Sx2(:)*9.274E-21; Sy2 = Sy2(:)*9.274E-21; Sz2 = Sz2(:)*9.274E-21;
    end
    xg = linspace(0,1,ngrid(1)+1); xg(end) = [];
    yg = linspace(0,1,ngrid(2)+1); yg(end) = [];
    zg = linspace(0,1,ngrid(3)+1); zg(end) = [];
    [xg,yg,zg] = ndgrid(xg,yg,zg); xg = xg(:); yg = yg(:); zg = zg(:);
end
if opt_job == 4; Mu = [0.5 0.5 0.5];
elseif get(handles.opt_ax,'Value') == 3; Mu = geo1.coords(end,:);
else Mu = str2num(get(handles.muon,'String'));
end
if get(handles.opt_ax,'Value') ~= 2 || get(handles.opt_job,'Value') == 4
    Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:);
end
if opt_job == 3; Rmax = str2double(get(handles.dmax,'String'));
elseif opt_job == 2; Rz = str2num(get(handles.Rzpe,'String'));
    if length(Rz) == 1; Rz = Rz*[1 1 1]; end
    Rmax = str2double(get(handles.Rmax,'String')) + sqrt(Rz(1)^2+Rz(2)^2+Rz(3)^2)/2;
else Rmax = str2double(get(handles.Rmax,'String'));
end
dmax = ceil((Rmax+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
fprintf('\n     size of unit cell: [%1.0f %1.0f %1.0f] \n',2*dmax+1)
if opt1 == 1;  Pos1 = [];
    if get(handles.opt_pos,'Value') == 1;
        for n = 1:size(cu1,1)
            cu = cu1(n,1)*latt1(1,:)+cu1(n,2)*latt1(2,:)+cu1(n,3)*latt1(3,:);
            for i = -dmax(1):dmax(1)
                for j = -dmax(2):dmax(2)
                    for k = -dmax(3):dmax(3); 
                        C = cu + i*latt1(1,:)+j*latt1(2,:)+k*latt1(3,:);
                        if sqrt(sum((C-Mu).*(C-Mu))) <= Rmax; 
                            Pos1 = [Pos1; C Spin1(n,:)]; 
                        end
                    end
                end
            end
        end
        Sp1 = Pos1(:,4:6); Pos1 = Pos1(:,1:3);
        fprintf('     number of ions   : %1.0f \n',size(Pos1,1))
    elseif get(handles.opt_pos,'Value') == 2 % defect
        filename2 = get(handles.pos2,'String'); cu2 = []; Pos3 = [];
        geo2 = poscar(filename2); latt2 = geo2.lattice; a=1; b=0;
        SpinVal2 = str2num(get(handles.S2,'String'));
        fid = fopen(get(handles.mag2,'String'));
        if fid == -1; Spin2 = str2num(get(handles.mag2,'String'));
        else Spin2 = str2num(fgetl(fid));
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
        Spin2 = Spin2*9.274E-21;
        if length(Spin2) == 3 % ferro
            mag = zeros(size(cu1,1),3);
            for i = 1:size(cu1,1); mag(i,:) = Spin2; end; Spin2 = mag;
        else Spin2 = reshape(Spin2,3,size(cu1,1))';
        end
        for i = 1:length(No)
            cu2 = [cu2; geo2.coords(sum(geo2.atomcount(1:No(i)-1))+1:sum(geo2.atomcount(1:No(i))),:)];
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
    end
elseif opt2 == 1; Pos = [];
    if get(handles.opt_pos,'Value') == 1; % perfect
        for i = -dmax(1):dmax(1)
            for j = -dmax(2):dmax(2)
                for k = -dmax(3):dmax(3);
                    Pos = [Pos; pos1+repmat(i*latt1(1,:)+j*latt1(2,:)+k*latt1(3,:),size(pos1,1),1) spin1];
                end
            end
        end
    elseif get(handles.opt_pos,'Value') == 2 % defect
        filename2 = get(handles.pos2,'String');
        geo2 = poscar([filename2 '/CONTCAR']); 
        spin2 = outcar([filename2 '/OUTCAR']);
        pos2 = geo2.coords; latt2 = geo2.lattice; num = (1:size(spin2,1))';
        spin = sqrt(spin2(:,1).^2+spin2(:,2).^2+spin2(:,3).^2);
        if length(lim) == 1; cond = (spin(:,1) < lim);
            pos2(cond,:) = []; spin2(cond,:) = [];  num(cond,:) = [];
        else cond = (spin(:,1) < lim(2)); pos2(cond,:) = []; spin2(cond,:) = []; num(cond,:) = [];
             cond = (spin(:,1) > lim(1)); pos2(cond,:) = []; spin2(cond,:) = []; num(cond,:) = [];
        end
        pos2 = pos2(:,1)*latt2(1,:)+pos2(:,2)*latt2(2,:)+pos2(:,3)*latt2(3,:);
        spinprt2 = [num spin2 sqrt(spin2(:,1).^2+spin2(:,2).^2+spin2(:,3).^2)]'; spin2 = spin2*9.274E-21;
        for i = -dmax(1):dmax(1)
            for j = -dmax(2):dmax(2)
                for k = -dmax(3):dmax(3)
                    if i == 0 && j == 0 && k == 0; % inside supercell
                         Pos = [Pos; pos1 spin1];
                    else Pos = [Pos; pos2+repmat(i*latt2(1,:)+j*latt2(2,:)+k*latt2(3,:),size(pos2,1),1) spin2];
                    end
                end
            end
        end
    end
    spin = Pos(:,4:6); pos = Pos(:,1:3);
    fprintf('     number of ions   : %1.0f \n',size(Pos,1))
elseif opt3 == 1;
    ngrid = ngrid.*(2*dmax+1);
    fprintf('     number of grids  : [%1.0f %1.0f %1.0f]\n',ngrid)
end
if get(handles.opt_job,'Value') == 1 % dipole calc
    fprintf('\n>> Dipole fields calculation \n')
    if opt1 == 1;
        if get(handles.opt_pos,'Value') == 1 % perfect
            x = (Mu(1)-Pos1(:,1))*1E-8; y = (Mu(2)-Pos1(:,2))*1E-8;
            z = (Mu(3)-Pos1(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
            SR = Sp1(:,1).*x + Sp1(:,2).*y + Sp1(:,3).*z;
            aPos = (R > 0) & (R <= Rmax*1E-8);
            SR = SR(aPos); R = R(aPos); Spin = Sp1(aPos,:);
        elseif get(handles.opt_pos,'Value') == 2 % defect
            x = (Mu(1)-Pos3(:,1))*1E-8; y = (Mu(2)-Pos3(:,2))*1E-8;
            z = (Mu(3)-Pos3(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
            SR = Sp3(:,1).*x + Sp3(:,2).*y + Sp3(:,3).*z;
            aPos = (R > 0) & (R <= Rmax*1E-8);
            SR = SR(aPos); R = R(aPos); Spin = Sp3(aPos,:);
        end
        Hx = sum((3*x(aPos).*SR - Spin(:,1).*R.^2)./R.^5); 
        Hy = sum((3*y(aPos).*SR - Spin(:,2).*R.^2)./R.^5);
        Hz = sum((3*z(aPos).*SR - Spin(:,3).*R.^2)./R.^5);
        
    elseif opt2 == 1; Pos = repmat(Mu,size(pos,1),1)-pos;
        R = sqrt(Pos(:,1).^2+Pos(:,2).^2+Pos(:,3).^2);
        Pos(R>Rmax,:) = []; spin(R>Rmax,:) = []; R(R>Rmax) = [];
        Pos(R<1E-6,:) = []; spin(R<1E-6,:) = []; R(R<1E-6) = [];
        Pos = Pos*1E-8; R = R*1E-8; % to cm
        SR = spin(:,1).*Pos(:,1) + spin(:,2).*Pos(:,2) + spin(:,3).*Pos(:,3);
        Hx = sum((3*Pos(:,1).*SR - spin(:,1).*R.^2)./R.^5); 
        Hy = sum((3*Pos(:,2).*SR - spin(:,2).*R.^2)./R.^5);
        Hz = sum((3*Pos(:,3).*SR - spin(:,3).*R.^2)./R.^5);
    elseif opt3 == 1; Hx = 0; Hy = 0; Hz = 0;
        for i = -dmax(1):dmax(1); X = xg+i;
            for j = -dmax(2):dmax(2); Y = yg+j;
                for k = -dmax(3):dmax(3); Z = zg+k;
                    xr = (Mu(1) - X*latt1(1,1) - Y*latt1(2,1) - Z*latt1(3,1))*1E-8;
                    yr = (Mu(2) - X*latt1(1,2) - Y*latt1(2,2) - Z*latt1(3,2))*1E-8;
                    zr = (Mu(3) - X*latt1(1,3) - Y*latt1(2,3) - Z*latt1(3,3))*1E-8;
                    if get(handles.opt_pos,'Value') == 2
                        if i == 0 && j == 0 && k == 0; Sx = Sx1; Sy = Sy1; Sz = Sz1;
                        else Sx = Sx2; Sy = Sy2; Sz = Sz2; 
                        end
                    end
                    R = sqrt(xr.^2 + yr.^2 + zr.^2); SR = Sx.*xr + Sy.*yr + Sz.*zr;
                    aPos = (R >= dt*1E-8) & (R <= Rmax*1E-8); R = R(aPos);
                    Hx = Hx + sum((3*xr(aPos).*SR(aPos) - Sx(aPos).*R.^2)./R.^5);
                    Hy = Hy + sum((3*yr(aPos).*SR(aPos) - Sy(aPos).*R.^2)./R.^5);
                    Hz = Hz + sum((3*zr(aPos).*SR(aPos) - Sz(aPos).*R.^2)./R.^5);
                end
            end
        end
    end
    H = sqrt(Hx^2+Hy^2+Hz^2);
elseif get(handles.opt_job,'Value') == 2 % Dipole ZPE
    fprintf('\n>> ZPE calculations \n')
    N = str2num(get(handles.Nzpe,'String'));
    if length(N) < 3; N = N*[1 1 1]; end
    if get(handles.opt_grid,'Value') == 2; N = ceil(Rz./N); end
    if opt1 == 1; fid = fopen(filename1);
    else fid = fopen([filename1 '/LOCPOT']);
    end
    geo1 = poscar(fid); fgetl(fid);
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
    loc = SX; SX = []; latt = inv(latt1);
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
    fprintf('     Hamiltonian matrix size: %1.0f x %1.0f \n',length(xi),length(xi))
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
    fprintf('\n>> Solving Schrodinger equation \n')
    [phi,Ezpe] = eig(A); Ezpe = (max(loc(:))-Ezpe(end))*1000;
    fprintf('     Zero-point Energy: %4.0f meV \n',Ezpe)
    fprintf('\n>> Dipole fields calculation \n')
    phi = phi(:,end).*conj(phi(:,end));
    Hx = zeros(length(xi),1); Hy = Hx; Hz = Hx;
    fprintf('     number of runs: %1.0f \n',length(xi))
    if opt1 == 1;
        if get(handles.opt_pos,'Value') == 1 % perfect
            for i = 1:length(xi)
                x = (xi(i)-Pos1(:,1))*1E-8; y = (yi(i)-Pos1(:,2))*1E-8;
                z = (zi(i)-Pos1(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
                SR = Sp1(:,1).*x + Sp1(:,2).*y + Sp1(:,3).*z;
                aPos = (R > 0) & (R <= Rmax*1E-8); SR = SR(aPos); R = R(aPos);
                Hx(i) = sum((3*x(aPos).*SR - Sp1(aPos,1).*R.^2)./R.^5); 
                Hy(i) = sum((3*y(aPos).*SR - Sp1(aPos,2).*R.^2)./R.^5);
                Hz(i) = sum((3*z(aPos).*SR - Sp1(aPos,3).*R.^2)./R.^5);
            end
        elseif get(handles.opt_pos,'Value') == 2 % defect
            for i = 1:length(xi)
                x = (xi(i)-Pos3(:,1))*1E-8; y = (yi(i)-Pos3(:,2))*1E-8;
                z = (zi(i)-Pos3(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
                SR = Sp3(:,1).*x + Sp3(:,2).*y + Sp3(:,3).*z;
                aPos = (R > 0) & (R <= Rmax*1E-8); SR = SR(aPos); R = R(aPos);
                Hx(i) = sum((3*x(aPos).*SR - Sp3(aPos,1).*R.^2)./R.^5); 
                Hy(i) = sum((3*y(aPos).*SR - Sp3(aPos,2).*R.^2)./R.^5);
                Hz(i) = sum((3*z(aPos).*SR - Sp3(aPos,3).*R.^2)./R.^5);
            end
        end
    elseif opt2 == 1;
        for i = 1:length(xi); Spin = spin;
            Pos = repmat([xi(i) yi(i) zi(i)],size(pos,1),1)-pos;
            R = sqrt(Pos(:,1).^2+Pos(:,2).^2+Pos(:,3).^2);
            Pos(R>Rmax,:) = []; Spin(R>Rmax,:) = []; R(R>Rmax) = [];
            Pos(R<1E-6,:) = []; Spin(R<1E-6,:) = []; R(R<1E-6) = [];
            Pos = Pos*1E-8; R = R*1E-8; % to cm
            SR = Spin(:,1).*Pos(:,1) + Spin(:,2).*Pos(:,2) + Spin(:,3).*Pos(:,3);
            Hx(i) = sum((3*Pos(:,1).*SR - Spin(:,1).*R.^2)./R.^5); 
            Hy(i) = sum((3*Pos(:,2).*SR - Spin(:,2).*R.^2)./R.^5);
            Hz(i) = sum((3*Pos(:,3).*SR - Spin(:,3).*R.^2)./R.^5);
        end
    elseif opt3 == 1; HX = zeros(length(xi),1); HY = HX; HZ = HX;
        for n = 1:length(xi); Hx = 0; Hy = 0; Hz = 0;
            for i = -dmax(1):dmax(1); X = xg+i;
                for j = -dmax(2):dmax(2); Y = yg+j;
                    for k = -dmax(3):dmax(3); Z = zg+k;
                        xr = (xi(n) - X*latt1(1,1) - Y*latt1(2,1) - Z*latt1(3,1))*1E-8;
                        yr = (yi(n) - X*latt1(1,2) - Y*latt1(2,2) - Z*latt1(3,2))*1E-8;
                        zr = (zi(n) - X*latt1(1,3) - Y*latt1(2,3) - Z*latt1(3,3))*1E-8;
                        if get(handles.opt_pos,'Value') == 2
                            if i == 0 && j == 0 && k == 0; Sx = Sx1; Sy = Sy1; Sz = Sz1;
                            else Sx = Sx2; Sy = Sy2; Sz = Sz2; 
                            end
                        end
                        R = sqrt(xr.^2 + yr.^2 + zr.^2); 
                        SR = Sx.*xr + Sy.*yr + Sz.*zr;
                        aPos = (R >= dt*1E-8) & (R <= Rmax*1E-8); R = R(aPos);
                        Hx = Hx + sum((3*xr(aPos).*SR(aPos) - Sx(aPos).*R.^2)./R.^5);
                        Hy = Hy + sum((3*yr(aPos).*SR(aPos) - Sy(aPos).*R.^2)./R.^5);
                        Hz = Hz + sum((3*zr(aPos).*SR(aPos) - Sz(aPos).*R.^2)./R.^5);
                    end
                end
            end
            HX(n) = Hx; HY(n) = Hy; HZ(n) = Hz;
        end
        Hx = HX; Hy = HY; Hz = HZ;
    end
    Hx = sum(phi.*Hx); Hy = sum(phi.*Hy); Hz = sum(phi.*Hz); H = sqrt(Hx^2+Hy^2+Hz^2);
elseif get(handles.opt_job,'Value') == 3 % converge
    fprintf('\n>> Dipole fields calculation \n')
    dmin = str2double(get(handles.dmin,'String')); if dmin < 5; dmin = 5; end;
    del = str2double(get(handles.del,'String')); del = linspace(dmin,Rmax,del);
    Hx = zeros(length(del),1); Hy = Hx; Hz = Hx;
    fprintf('     number of runs: %1.0f \n',length(del))
    if opt1 == 1;
        if get(handles.opt_pos,'Value') == 1 % perfect
            x = (Mu(1)-Pos1(:,1))*1E-8; y = (Mu(2)-Pos1(:,2))*1E-8;
            z = (Mu(3)-Pos1(:,3))*1E-8; r = sqrt(x.^2+y.^2+z.^2);
            sr = Sp1(:,1).*x + Sp1(:,2).*y + Sp1(:,3).*z;
            for i = 1:length(del);
                aPos = (r > 0) & (r <= del(i)*1E-8);
                SR = sr(aPos); R = r(aPos); Spin = Sp1(aPos,:);
                Hx(i) = sum((3*x(aPos).*SR - Spin(:,1).*R.^2)./R.^5); 
                Hy(i) = sum((3*y(aPos).*SR - Spin(:,2).*R.^2)./R.^5);
                Hz(i) = sum((3*z(aPos).*SR - Spin(:,3).*R.^2)./R.^5);
            end
        elseif get(handles.opt_pos,'Value') == 2 % defect
            x = (Mu(1)-Pos3(:,1))*1E-8; y = (Mu(2)-Pos3(:,2))*1E-8;
            z = (Mu(3)-Pos3(:,3))*1E-8; r = sqrt(x.^2+y.^2+z.^2);
            sr = Sp3(:,1).*x + Sp3(:,2).*y + Sp3(:,3).*z;
            for i = 1:length(del);
                aPos = (r > 0) & (r <= del(i)*1E-8);
                SR = sr(aPos); R = r(aPos); Spin = Sp3(aPos,:);
                Hx(i) = sum((3*x(aPos).*SR - Spin(:,1).*R.^2)./R.^5);
                Hy(i) = sum((3*y(aPos).*SR - Spin(:,2).*R.^2)./R.^5);
                Hz(i) = sum((3*z(aPos).*SR - Spin(:,3).*R.^2)./R.^5);
            end
        end
    elseif opt2 == 1; Pos = repmat(Mu,size(pos,1),1)-pos;
        R = sqrt(Pos(:,1).^2+Pos(:,2).^2+Pos(:,3).^2);
        Pos(R>Rmax,:) = []; spin(R>Rmax,:) = []; R(R>Rmax) = [];
        Pos(R<1E-6,:) = []; spin(R<1E-6,:) = []; R(R<1E-6) = [];
        Pos = Pos*1E-8; R = R*1E-8; % to cm
        SR = spin(:,1).*Pos(:,1) + spin(:,2).*Pos(:,2) + spin(:,3).*Pos(:,3);
        for i = 1:length(del);
            aPos = (R <= del(i)*1E-8);
            sr = SR(aPos); r = R(aPos); spn = spin(aPos,:); Ion = Pos(aPos,:);
            Hx(i) = sum((3*Ion(:,1).*sr - spn(:,1).*r.^2)./r.^5); 
            Hy(i) = sum((3*Ion(:,2).*sr - spn(:,2).*r.^2)./r.^5);
            Hz(i) = sum((3*Ion(:,3).*sr - spn(:,3).*r.^2)./r.^5);
        end
    elseif opt3 == 1; HX = zeros(length(del),1); HY = HX; HZ = HX;
        for n = 1:length(del); Hx = 0; Hy = 0; Hz = 0;
            dmax = ceil((del(n)+Mu)./(latt1(1,:)+latt1(2,:)+latt1(3,:)));
            for i = -dmax(1):dmax(1); X = xg+i;
                for j = -dmax(2):dmax(2); Y = yg+j;
                    for k = -dmax(3):dmax(3); Z = zg+k;
                        xr = (Mu(1) - X*latt1(1,1) - Y*latt1(2,1) - Z*latt1(3,1))*1E-8;
                        yr = (Mu(2) - X*latt1(1,2) - Y*latt1(2,2) - Z*latt1(3,2))*1E-8;
                        zr = (Mu(3) - X*latt1(1,3) - Y*latt1(2,3) - Z*latt1(3,3))*1E-8;
                        if get(handles.opt_pos,'Value') == 2
                            if i == 0 && j == 0 && k == 0
                                Sx = Sx1; Sy = Sy1; Sz = Sz1;
                            else Sx = Sx2; Sy = Sy2; Sz = Sz2;
                            end
                        end
                        R = sqrt(xr.^2 + yr.^2 + zr.^2); SR = Sx.*xr + Sy.*yr + Sz.*zr;
                        aPos = (R >= dt*1E-8) & (R <= del(n)*1E-8); R = R(aPos);
                        Hx = Hx + sum((3*xr(aPos).*SR(aPos) - Sx(aPos).*R.^2)./R.^5);
                        Hy = Hy + sum((3*yr(aPos).*SR(aPos) - Sy(aPos).*R.^2)./R.^5);
                        Hz = Hz + sum((3*zr(aPos).*SR(aPos) - Sz(aPos).*R.^2)./R.^5);
                    end
                end
            end
            HX(n) = Hx; HY(n) = Hy; HZ(n) = Hz;
        end
        Hx = HX; Hy = HY; Hz = HZ;
    end
    H = sqrt(Hx.^2+Hy.^2+Hz.^2); out = [del' H Hx Hy Hz]';
    figure; hold on; box on;
    if del(1) < 5; plot([0 del(end)],H(end)*[1 1],'k--')
    else plot([del(1) del(end)],H(end)*[1 1],'k--')
    end
    plot(del,H,'Marker','o','LineWidth',2); H = H(end); Hz = Hz(end);
    legend([num2str(round(H)) ' G']); Hx = Hx(end); Hy = Hy(end);
    xlabel('Calc range (Angs)','FontSize',18);
    ylabel('Dipole fields (G)','FontSize',18);
    
elseif get(handles.opt_job,'Value') == 4 % makefile
    fprintf('\n>> Dipole fields calculation \n')
    Nz = str2double(get(handles.Nz,'String'));
    Nx = str2double(get(handles.Nx,'String'));
    Ny = str2double(get(handles.Ny,'String'));
    xr = linspace(0,1,Nx+1); xr(end) = [];
    yr = linspace(0,1,Ny+1); yr(end) = [];
    zr = linspace(0,1,Nz+1); zr(end) = [];
    [xr,yr,zr] = ndgrid(xr,yr,zr); x = xr(:); y = yr(:); z = zr(:);
    xr = x*latt1(1,1) + y*latt1(2,1) + z*latt1(3,1);
    yr = x*latt1(1,2) + y*latt1(2,2) + z*latt1(3,2);
    zr = x*latt1(1,3) + y*latt1(2,3) + z*latt1(3,3);
    Hx = zeros(length(xr),1); Hy = Hx; Hz = Hx;
    fprintf('     number of runs : %1.0f \n',length(xr))
    if opt1 == 1;
        if get(handles.opt_pos,'Value') == 1 % perfect
            for i = 1:length(xr)
                x = (xr(i)-Pos1(:,1))*1E-8; y = (yr(i)-Pos1(:,2))*1E-8;
                z = (zr(i)-Pos1(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
                SR = Sp1(:,1).*x + Sp1(:,2).*y + Sp1(:,3).*z;
                aPos = (R > 0) & (R <= Rmax*1E-8); SR = SR(aPos); R = R(aPos);
                Hx(i) = sum((3*x(aPos).*SR - Sp1(aPos,1).*R.^2)./R.^5); 
                Hy(i) = sum((3*y(aPos).*SR - Sp1(aPos,2).*R.^2)./R.^5);
                Hz(i) = sum((3*z(aPos).*SR - Sp1(aPos,3).*R.^2)./R.^5);
            end
        elseif get(handles.opt_pos,'Value') == 2 % defect
            for i = 1:length(xr)
                x = (xr(i)-Pos3(:,1))*1E-8; y = (yr(i)-Pos3(:,2))*1E-8;
                z = (zr(i)-Pos3(:,3))*1E-8; R = sqrt(x.^2+y.^2+z.^2);
                SR = Sp3(:,1).*x + Sp3(:,2).*y + Sp3(:,3).*z;
                aPos = (R > 0) & (R <= Rmax*1E-8); SR = SR(aPos); R = R(aPos);
                Hx(i) = sum((3*x(aPos).*SR - Sp3(aPos,1).*R.^2)./R.^5); 
                Hy(i) = sum((3*y(aPos).*SR - Sp3(aPos,2).*R.^2)./R.^5);
                Hz(i) = sum((3*z(aPos).*SR - Sp3(aPos,3).*R.^2)./R.^5);
            end
        end
    elseif opt2 == 1;
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
    elseif opt3 == 1; HX = Hx; HY = Hx; HZ = Hx;
        for n = 1:length(xr); Hx = 0; Hy = 0; Hz = 0;
            for i = -dmax(1):dmax(1); X = xg+i;
                for j = -dmax(2):dmax(2); Y = yg+j;
                    for k = -dmax(3):dmax(3); Z = zg+k;
                        x = (xr(n) - X*latt1(1,1) - Y*latt1(2,1) - Z*latt1(3,1))*1E-8;
                        y = (yr(n) - X*latt1(1,2) - Y*latt1(2,2) - Z*latt1(3,2))*1E-8;
                        z = (zr(n) - X*latt1(1,3) - Y*latt1(2,3) - Z*latt1(3,3))*1E-8;
                        if get(handles.opt_pos,'Value') == 2
                            if i == 0 && j == 0 && k == 0
                                Sx = Sx1; Sy = Sy1; Sz = Sz1;
                            else Sx = Sx2; Sy = Sy2; Sz = Sz2; 
                            end
                        end
                        R = sqrt(x.^2 + y.^2 + z.^2); SR = Sx.*x + Sy.*y + Sz.*z;
                        aPos = (R >= dt*1E-8) & (R <= Rmax*1E-8); R = R(aPos);
                        Hx = Hx + sum((3*x(aPos).*SR(aPos) - Sx(aPos).*R.^2)./R.^5);
                        Hy = Hy + sum((3*y(aPos).*SR(aPos) - Sy(aPos).*R.^2)./R.^5);
                        Hz = Hz + sum((3*z(aPos).*SR(aPos) - Sz(aPos).*R.^2)./R.^5);
                    end
                end
            end
            HX(n) = Hx; HY(n) = Hy; HZ(n) = Hz; 
        end
        Hx = HX; Hy = HY; Hz = HZ;
    end
    fprintf('\n>> Writing output file \n')
    Hx = reshape(Hx,[Nx Ny Nz]); Hy = reshape(Hy,[Nx Ny Nz]);
    Hz = reshape(Hz,[Nx Ny Nz]); H = sqrt(Hx.^2+Hy.^2+Hz.^2);
    fid = fopen([pathname '/LOCPOT_dip']);
    if fid ~= -1
        i = 1; fid = fopen([pathname '/LOCPOT_dip(1)']);
        while fid ~= -1
            i = i+1; fclose(fid); fid = fopen([pathname '/LOCPOT_dip(' num2str(i) ')']);
        end
        fid = fopen([pathname '/LOCPOT_dip(' num2str(i) ')'],'w');
        fprintf(['     ' pathname '/LOCPOT_dip(' num2str(i) ') \n'])
    else fid = fopen([pathname '/LOCPOT_dip'],'w');
        fprintf(['     ' pathname '/LOCPOT_dip \n'])
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
    Hx = ''; Hy = ''; Hz = ''; H = 'Finish!';
end
if get(handles.opt_job,'Value') < 4
    if abs(Hx) < 1E-5; Hx = 0; end; if abs(Hy) < 1E-5; Hy = 0; end
    if abs(Hz) < 1E-5; Hz = 0; end; if abs(H) < 1E-5; H = 0; end
    fprintf('     dipole fields : %4.4f Gauss \n',H);
    fprintf('     [Hx, Hy, Hz]  : [%4.4f, %4.4f, %4.4f]\n',Hx,Hy,Hz);
end
set(handles.Ht,'String',H); set(handles.Hx,'String',Hx)
set(handles.Hy,'String',Hy); set(handles.Hz,'String',Hz)
if get(handles.fig,'Value') == 1 && get(handles.opt_job,'Value') < 4
    figure; Cu = []; fprintf('\n>> Drawing figure \n')
    if opt2 == 1; latt = inv(latt1); Spin1 = spin1/9.274E-21;
        cu1 = pos1(:,1)*latt(1,:)+pos1(:,2)*latt(2,:)+pos1(:,3)*latt(3,:);
    end
    for n = 1:size(cu1,1)
        for i = -1:1
            for j = -1:1
                for k = -1:1
                    if cu1(n,1)+i > -0.02 && cu1(n,1)+i < 1.02 ...
                            && cu1(n,2)+j > -0.02 && cu1(n,2)+j < 1.02 ...
                            && cu1(n,3)+k > -0.02 && cu1(n,3)+k < 1.02
                        Cu = [Cu; (cu1(n,1)+i)*latt1(1,:)+(cu1(n,2)+j)*latt1(2,:)+(cu1(n,3)+k)*latt1(3,:) Spin1(n,:)];
                    end
                end
            end
        end
    end
    atom(Cu,0.4,'b',1); atom(Mu,0.3,'m',1); arrow(Cu(:,1:3),Cu(:,4:6),'k');
    arrow(Mu,[Hx Hy Hz],'r')
    if get(handles.opt_job,'Value') == 2
        drbox([max(xi(:))-min(xi(:)) 0 0; 0 max(yi(:))-min(yi(:)) 0; 0 0 ...
            max(zi(:))-min(zi(:))],[min(xi(:)) min(yi(:)) min(zi(:))],'b');
    end
    xlabel('x'); ylabel('y'); zlabel('z'); daspect([1 1 1])
    drbox(latt1); camlight right local; rotate3d on; axis tight
end
if get(handles.fig_range,'Value') == 1 && get(handles.opt_job,'Value') < 4
    if get(handles.fig,'Value') == 0; fprintf('\n>> drawing figure \n'); end
    figure;
    if opt2 == 1; Pos = repmat(Mu,size(pos,1),1)-pos;
        R = sqrt(Pos(:,1).^2+Pos(:,2).^2+Pos(:,3).^2);
        pos(R>Rmax,:) = []; R(R>Rmax) = []; pos(R<1E-6,:) = [];
        atom(pos,0.3,'b',1);
    elseif get(handles.opt_pos,'Value') == 1 % perfect
       atom(Pos1,0.3,'b',1); 
    elseif get(handles.opt_pos,'Value') == 2 % defect
        atom(Pos3,0.3,'b',1);
    end
    atom(Mu,0.3,'r',0.8); drbox(latt1); 
    xlabel('x'); ylabel('y'); zlabel('z'); daspect([1 1 1])
    camlight right local; rotate3d on; axis tight
end
if get(handles.log,'Value') == 1;
    fprintf('\n>> Writing log file \n')
    if get(handles.opt_job,'Value') == 1; job = '';
    elseif get(handles.opt_job,'Value') == 2; job = '_zpe';
    elseif get(handles.opt_job,'Value') == 3; job = '_test';
    elseif get(handles.opt_job,'Value') == 4; job = '_make';
    end
    fid = fopen([pathname '/log' job]);
    if fid ~= -1; i = 1; fid = fopen([pathname '/log' job '(1)']);
        while fid ~= -1
            i = i+1; fclose(fid); fid = fopen([pathname '/log' job '(' num2str(i) ')']);
        end
        fid = fopen([pathname '/log' job '(' num2str(i) ')'],'w');
        fprintf(['     ' pathname '/log' job '(' num2str(i) ') \n'])
    else fid = fopen([pathname '/log' job],'w');
        fprintf(['     ' pathname '/log' job ' \n'])
    end
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,'            Dipole Fields Calculation \n');
    fprintf(fid,'   ============================================ \n');
    fprintf(fid,['   Running on: ' datestr(jam) '\n']);
    fprintf(fid,'\n   POSCAR location: \n');
    fprintf(fid,['   ' filename1 '\n']);
    if get(handles.opt_pos,'Value') == 2; fprintf(fid,['   ' filename2 '\n']); end
    if opt1 == 1;
        if isempty(geo1.symbols);
            if length(No) == 1; fprintf(fid,['\n   ion : '  num2str(No) '\n']);
            else fprintf(fid,['\n   ion of : ['  num2str(ion_no) ']\n']); end;
        else ion = '';
            for i = 1:length(No); ion = [ion cell2mat(geo1.symbols(No(i))) ' ']; end
            fprintf(fid,['\n   ion : ' ion '\n']);
        end
        fprintf(fid,['   Magmom_1  : '  num2str(SpinVal1) ' muB \n']);
        fprintf(fid,'   Spins1: \n');
        fprintf(fid,'     %4.3f   %4.3f   %4.3f\n',spinprt1);
        if get(handles.opt_pos,'Value') == 1; jum = size(Pos1,1);
        else jum = size(Pos3,1);
            fprintf(fid,['   Magmom_2  : '  num2str(SpinVal2) ' muB \n']);
            fprintf(fid,'   Spins2: \n');
            fprintf(fid,'     %4.3f   %4.3f   %4.3f\n',spinprt2);
        end
    elseif opt2 == 1
        fprintf(fid,'\n     No      Mx      My      Mz      Mtotal');
        fprintf(fid,'\n   --------------------------------------------\n');
        fprintf(fid,'   %4.0f    %6.3f  %6.3f  %6.3f    %6.3f \n',spinprt1(:));
        fprintf(fid,'   --------------------------------------------\n');
        fprintf(fid,'    sum    %6.3f  %6.3f  %6.3f    %6.3f \n',sum(spinprt1(2,:)), ...
            sum(spinprt1(3,:)),sum(spinprt1(4,:)),sum(spinprt1(5,:)));
        if get(handles.opt_pos,'Value') == 2 
            fprintf(fid,'\n     No2     Mx      My      Mz      Mtotal');
            fprintf(fid,'\n   --------------------------------------------\n');
            fprintf(fid,'   %4.0f    %6.3f  %6.3f  %6.3f    %6.3f \n',spinprt2(:));
            fprintf(fid,'   --------------------------------------------\n');
            fprintf(fid,'    sum    %6.3f  %6.3f  %6.3f    %6.3f \n',sum(spinprt2(2,:)), ...
                sum(spinprt2(3,:)),sum(spinprt2(4,:)),sum(spinprt2(5,:)));
        end
    elseif opt3 == 1;
        if get(handles.opt_pos,'Value') == 1
            fprintf(fid,'\n   Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
            if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
                fprintf(fid,'     M(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
            end
        else fprintf(fid,'\n   Magmom1 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
            if sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2) < 1
                fprintf(fid,'     M1(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx1(1)^2+sumx1(2)^2+sumx1(3)^2),sumx1);
            end
            fprintf(fid,'   Magmom2 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2),sumy);
            if sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2) < 1
                fprintf(fid,'     M2(abs) : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy1(1)^2+sumy1(2)^2+sumy1(3)^2),sumy1);
            end
        end
    end
    if get(handles.opt_job,'Value') < 4; 
        latt1 = inv(latt1); Mu = Mu(1)*latt1(1,:)+Mu(2)*latt1(2,:)+Mu(3)*latt1(3,:);
        fprintf(fid,'\n   Muon position : [%4.4f %4.4f %4.4f]\n',Mu);
    end
    fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',Rmax);
    fprintf(fid,'   size of unit cell : [%1.0f %1.0f %1.0f] \n',2*dmax+1);
    if opt1 == 1; fprintf(fid,'   number of ions    : %1.0f \n',jum);
    elseif opt2 == 1; fprintf(fid,'   number of ions    : %1.0f \n',size(Pos,1)); 
    elseif opt3 == 1 && get(handles.opt_job,'Value') < 4
        fprintf(fid,'   number of grids   : [%1.0f %1.0f %1.0f] \n',ngrid); 
    end
    if get(handles.opt_job,'Value') == 2 
        fprintf(fid,'\n   ZPE Range   : %4.4f Angstrom\n',sum(Rz)/3);
        fprintf(fid,'   Resolution  : %4.4f Angstrom\n',sum(Rz)/sum(N));
        fprintf(fid,'   Matrix size : %1.0f x %1.0f \n',length(xi),length(xi));
        fprintf(fid,'   ZPE-Energy  : %1.0f meV \n',Ezpe);
    elseif get(handles.opt_job,'Value') == 3
        fprintf(fid,'\n   R_dip(A)  H(G)       [Hx, Hy, Hz] \n');
        fprintf(fid,'   ----------------------------------------------------\n');
        fprintf(fid,' %8.4f  %8.4f     [%8.4f  %8.4f  %8.4f] \n',out(:));
    elseif get(handles.opt_job,'Value') == 4
        fprintf(fid,'   number of grids   : [%1.0f %1.0f %1.0f] \n',Nx,Ny,Nz);
    end
    if get(handles.opt_job,'Value') < 3
        fprintf(fid,'\n   ============================================ ');
        fprintf(fid,'\n   Dipole fields : %4.4f Gauss\n',H);
        fprintf(fid,'   [Hx, Hy, Hz]  : [%4.4f, %4.4f, %4.4f]\n',Hx,Hy,Hz);
    else fprintf(fid,'             Calculations Completed! \n');
    end
    fprintf(fid,'   ============================================ \n');
    t = toc; h = floor(t/3600); m = floor((t-h*3600)/60); t = t-h*3600 - m*60;
    if h >= 1; fprintf(fid,'   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n\n',h,m,t);  
    elseif m >= 1; fprintf(fid,'   Elapsed time is %1.0f min %1.4f sec.\n\n',m,t);
    else fprintf(fid,'   Elapsed time is %1.4f sec.\n\n',t);
    end
    fclose(fid);
end
fprintf('   ============================================ \n')
fprintf('             Calculations Completed! \n')
fprintf('   ============================================ \n')
if get(handles.log,'Value') == 0;
    t = toc; h = floor(t/3600); m = floor((t-h*3600)/60); t = t-h*3600 - m*60;
end
if h >= 1; fprintf('   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n',h,m,t);  
elseif m >= 1; fprintf('   Elapsed time is %1.0f min %1.4f sec.\n',m,t);
else fprintf('   Elapsed time is %1.4f sec.\n',t);
end

function Rmax_Callback(~, ~, ~)
function Rmax_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ht_Callback(~, ~, ~)
function Ht_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Hz_Callback(~, ~, ~)
function Hz_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Hy_Callback(~, ~, ~)
function Hy_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Hx_Callback(~, ~, ~)
function Hx_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function log_Callback(~, ~, ~)
function fig_Callback(~, ~, handles)
if get(handles.opt3,'Value') == 1 || get(handles.opt_job,'Value') == 4;
    set(handles.fig,'Value',0);
end

function fig_range_Callback(~, ~, handles)
if get(handles.opt3,'Value') == 1 || get(handles.opt_job,'Value') == 4;
    set(handles.fig_range,'Value',0);
end

function load1_Callback(~, ~, handles)
global filename pathname; 
[filename, pathname] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname); 
    if get(handles.opt2,'Value') == 1
        set(handles.pos1,'String',pathname); 
    elseif get(handles.opt_job,'Value') == 2
        set(handles.pos1,'String',pathname);
    else set(handles.pos1,'String',[pathname filename]);
    end
end

function pos1_Callback(~, ~, ~)
function pos1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mag1_Callback(~, ~, ~)
function mag1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mag2_Callback(~, ~, ~)
function mag2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function S1_Callback(~, ~, ~)
function S1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function S2_Callback(~, ~, ~)
function S2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load2_Callback(~, ~, handles)
if get(handles.opt_pos,'Value') == 2
    [filename2, pathname2] = uigetfile({'*.*', 'All Files (*.*)'});
    if ~isnumeric(pathname2); 
        if get(handles.opt2,'Value') == 1
            set(handles.pos2,'String',pathname2); 
        else set(handles.pos2,'String',[pathname2 filename2]); 
        end
    end
end

function pos2_Callback(~, ~, ~)
function pos2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_pos_Callback(~, ~, handles)
if get(handles.opt_pos,'Value') == 1
    set(handles.pos2,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.pos2,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.mag2,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.mag2,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.S2,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.S2,'BackgroundColor',[0.8 0.8 0.8])
else set(handles.pos2,'BackgroundColor',[1 1 1])
    set(handles.pos2,'ForegroundColor',[0 0 0])
    if get(handles.opt1,'Value') == 1
        set(handles.mag2,'BackgroundColor',[1 1 1])
        set(handles.mag2,'ForegroundColor',[0 0 0])
        set(handles.S2,'BackgroundColor',[1 1 1])
        set(handles.S2,'ForegroundColor',[0 0 0])
    end
end

function opt_pos_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_job_Callback(~, ~, handles)
if get(handles.opt_job,'Value') == 2
    if get(handles.opt1,'Value') == 1
        set(handles.load1,'String','LOCPOT');
    else set(handles.load1,'String','Folder');
    end
    set(handles.uipanel5,'ForegroundColor',[0 0 1])
    set(handles.Rzpe,'ForegroundColor',[0 0 0]); set(handles.Nzpe,'ForegroundColor',[0 0 0]);
    set(handles.Rzpe,'BackgroundColor',[1 1 1]); set(handles.Nzpe,'BackgroundColor',[1 1 1]);
else set(handles.uipanel5,'ForegroundColor',[0.25 0.25 0.25])
    if get(handles.opt1,'Value') == 1
        set(handles.load1,'String','POSCAR')
    elseif get(handles.opt3,'Value') == 1
        set(handles.load1,'String','CHGCAR')
    end
    set(handles.Rzpe,'ForegroundColor',[0.8 0.8 0.8]);
    set(handles.Nzpe,'ForegroundColor',[0.8 0.8 0.8]);
    set(handles.Rzpe,'BackgroundColor',[0.8 0.8 0.8]);
    set(handles.Nzpe,'BackgroundColor',[0.8 0.8 0.8]);
end
if get(handles.opt_job,'Value') == 3
    set(handles.uipanel2,'ForegroundColor',[0 0 1])
    set(handles.dmin,'ForegroundColor',[0 0 0]);
    set(handles.dmin,'BackgroundColor',[1 1 1]);
    set(handles.dmax,'ForegroundColor',[0 0 0]);
    set(handles.dmax,'BackgroundColor',[1 1 1]);
    set(handles.del,'ForegroundColor',[0 0 0]);
    set(handles.del,'BackgroundColor',[1 1 1]);
    set(handles.Rmax,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.Rmax,'BackgroundColor',[0.8 0.8 0.8])
else set(handles.Rmax,'BackgroundColor',[1 1 1])
    set(handles.uipanel2,'ForegroundColor',[0.25 0.25 0.25])
    set(handles.dmin,'ForegroundColor',[0.8 0.8 0.8]);
    set(handles.dmin,'BackgroundColor',[0.8 0.8 0.8]);
    set(handles.dmax,'ForegroundColor',[0.8 0.8 0.8]);
    set(handles.dmax,'BackgroundColor',[0.8 0.8 0.8]);
    set(handles.del,'ForegroundColor',[0.8 0.8 0.8]);
    set(handles.del,'BackgroundColor',[0.8 0.8 0.8]);
    set(handles.Rmax,'ForegroundColor',[0 0 0])
end
if get(handles.opt_job,'Value') == 4;
    set(handles.fig,'Value',0); set(handles.fig_range,'Value',0);
    set(handles.uipanel4,'ForegroundColor',[0 0 1])
    set(handles.Nx,'ForegroundColor',[0 0 0]);
    set(handles.Nx,'BackgroundColor',[1 1 1]);
    set(handles.Ny,'ForegroundColor',[0 0 0]);
    set(handles.Ny,'BackgroundColor',[1 1 1]);
    set(handles.Nz,'ForegroundColor',[0 0 0]);
    set(handles.Nz,'BackgroundColor',[1 1 1]);
    set(handles.muon,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.muon,'BackgroundColor',[0.8 0.8 0.8])
else set(handles.Nx,'ForegroundColor',[0.8 0.8 0.8]);
    set(handles.uipanel4,'ForegroundColor',[0.25 0.25 0.25])
    set(handles.Nx,'BackgroundColor',[0.8 0.8 0.8]);
    set(handles.Ny,'ForegroundColor',[0.8 0.8 0.8]);
    set(handles.Ny,'BackgroundColor',[0.8 0.8 0.8]);
    set(handles.Nz,'ForegroundColor',[0.8 0.8 0.8]);
    set(handles.Nz,'BackgroundColor',[0.8 0.8 0.8]);
    if get(handles.opt_ax,'Value') < 3;
        set(handles.muon,'ForegroundColor',[0 0 0])
        set(handles.muon,'BackgroundColor',[1 1 1])
    end
end

function opt_job_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt1_Callback(~, ~, handles)
set(handles.opt1,'Value',1); set(handles.opt2,'Value',0); set(handles.opt3,'Value',0)
set(handles.opt_text,'String','Ion Number: '); set(handles.opt_in,'String','1');
if get(handles.opt_pos,'Value') == 2
    set(handles.mag2,'BackgroundColor',[1 1 1])
    set(handles.mag2,'ForegroundColor',[0 0 0])
    set(handles.S2,'BackgroundColor',[1 1 1])
    set(handles.S2,'ForegroundColor',[0 0 0])
end
set(handles.opt_in,'ForegroundColor',[0 0 0])
set(handles.opt_in,'BackgroundColor',[1 1 1])
set(handles.mag1,'BackgroundColor',[1 1 1])
set(handles.mag1,'ForegroundColor',[0 0 0])
set(handles.S1,'BackgroundColor',[1 1 1])
set(handles.S1,'ForegroundColor',[0 0 0])
set(handles.load2,'String','POSCAR_out');
if get(handles.opt_job,'Value') == 2
    set(handles.load1,'String','LOCPOT');
else set(handles.load1,'String','POSCAR');
end

function opt2_Callback(~, ~, handles)
set(handles.opt1,'Value',0); set(handles.opt2,'Value',1); set(handles.opt3,'Value',0)
set(handles.opt_text,'String','Cutoff SpinVal: '); set(handles.opt_in,'String','0');
set(handles.opt_in,'ForegroundColor',[0 0 0])
set(handles.opt_in,'BackgroundColor',[1 1 1])
set(handles.mag2,'ForegroundColor',[0.8 0.8 0.8])
set(handles.mag2,'BackgroundColor',[0.8 0.8 0.8])
set(handles.S2,'ForegroundColor',[0.8 0.8 0.8])
set(handles.S2,'BackgroundColor',[0.8 0.8 0.8])
set(handles.mag1,'ForegroundColor',[0.8 0.8 0.8])
set(handles.mag1,'BackgroundColor',[0.8 0.8 0.8])
set(handles.S1,'ForegroundColor',[0.8 0.8 0.8])
set(handles.S1,'BackgroundColor',[0.8 0.8 0.8])
set(handles.load1,'String','Folder');
set(handles.load2,'String','Folder_out');

function opt3_Callback(~, ~, handles)
set(handles.fig,'Value',0); set(handles.fig_range,'Value',0);
set(handles.opt1,'Value',0); set(handles.opt2,'Value',0); set(handles.opt3,'Value',1)
if get(handles.opt_job,'Value') == 2
    set(handles.load1,'String','Folder');
else set(handles.load1,'String','CHGCAR');
end
set(handles.load2,'String','CHGCAR_out');
set(handles.opt_in,'ForegroundColor',[0.8 0.8 0.8])
set(handles.opt_in,'BackgroundColor',[0.8 0.8 0.8])
set(handles.mag2,'ForegroundColor',[0.8 0.8 0.8])
set(handles.mag2,'BackgroundColor',[0.8 0.8 0.8])
set(handles.S2,'ForegroundColor',[0.8 0.8 0.8])
set(handles.S2,'BackgroundColor',[0.8 0.8 0.8])
set(handles.mag1,'ForegroundColor',[0.8 0.8 0.8])
set(handles.mag1,'BackgroundColor',[0.8 0.8 0.8])
set(handles.S1,'ForegroundColor',[0.8 0.8 0.8])
set(handles.S1,'BackgroundColor',[0.8 0.8 0.8])

function opt_in_Callback(~, ~, ~)
function opt_in_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function muon_Callback(~, ~, ~)
function muon_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_ax_Callback(~, ~, handles)
if get(handles.opt_ax,'Value') == 3
    set(handles.muon,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.muon,'BackgroundColor',[0.8 0.8 0.8])
elseif get(handles.opt_job,'Value') < 3
    set(handles.muon,'ForegroundColor',[0 0 0])
    set(handles.muon,'BackgroundColor',[1 1 1])
end

function opt_ax_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Rzpe_Callback(~, ~, ~)
function Rzpe_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nzpe_Callback(~, ~, ~)
function Nzpe_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_grid_Callback(~, ~, ~)
function opt_grid_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ny_Callback(~, ~, ~)
function Ny_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nz_Callback(~, ~, ~)
function Nz_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nx_Callback(~, ~, ~)
function Nx_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function del_Callback(~, ~, ~)
function del_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dmax_Callback(~, ~, ~)
function dmax_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dmin_Callback(~, ~, ~)
function dmin_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function geometry = poscar(filename)
if ~isnumeric(filename)
	geometry.filename = filename;
    fid = fopen(filename);
else fid = filename;
    geometry.filename = fopen(fid);
end
geometry.comment = fgetl(fid);
scale = fscanf(fid, '%f',1); cartesian = 0;
geometry.lattice = fscanf(fid, '%f %f %f', [3 3])'; 
geometry.lattice = geometry.lattice*scale; fgetl(fid);
line = fgetl(fid); has_symbols_line = false;
if sum(isstrprop(line, 'digit')) == 0
    geometry.symbols = regexp(line, '([^ ]*)', 'match');
    line = fgetl(fid); has_symbols_line = true;
else geometry.symbols = {};
end
geometry.atomcount = sscanf(line,'%d');
natoms = sum(geometry.atomcount);
line = fgetl(fid); geometry.selective = 0;
if line(1) == 's' || line(1) == 'S'
    geometry.selective = 1; line = fgetl(fid);
end
if line(1) == 'C' || line(1) == 'c' || line(1) == 'K' || line(1) =='k'
	cartesian = 1;
end
for i = 1:natoms
	line = fgetl(fid); geometry.coords(i,:) = sscanf(line, '%f %f %f');
    if ~has_symbols_line
        if numel(strfind(line,'!') == 1)
            str = regexp(line, '[^ ]*', 'match');
            str = str{end}; newelement = false;
            if numel(geometry.symbols) == 0
            	newelement = true;
            elseif strcmp(geometry.symbols{end},str) == 0
                newelement = true;
            end
            if newelement; geometry.symbols{end+1} = str; end
        end
    end
end
if cartesian == 1
	geometry.coords = geometry.coords*scale;
    geometry.coords = geometry.coords/geometry.lattice;
end 
if ~isnumeric(filename); fclose(fid); end

function [chg, mag_x, mag_y, mag_z, geo] = chgcar(filename)
fid = fopen(filename); geo = poscar(fid);
vol = abs(dot(geo.lattice(1,:),cross(geo.lattice(2,:),geo.lattice(3,:))));
natoms = sum(geo.atomcount); fgetl(fid); 
gridsize = fscanf(fid, '%d %d %d', [3 1])';
chg = fscanf(fid, '%f', [prod(gridsize,2) 1])';
chg = reshape(chg,gridsize);
mag_x = []; mag_y = []; mag_z = [];
fgetl(fid); pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'
    chg = chg/vol;
    for i = 1:natoms; line = fgetl(fid);
        nentries = sscanf(line,['augmentation occupancies ' num2str(i) ' %d']);
        fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
    end    
end
pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'; % MAGX
    for i = 1:natoms; line = fgetl(fid);
        nentries = sscanf(line,['augmentation occupancies (imaginary part) ' num2str(i) ' %d']); 
        fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
    end
    fscanf(fid, '%f', [natoms 1]);
    fgetl(fid); pos = ftell(fid); fgetl(fid);
    if ~feof(fid); fseek(fid,pos,'bof');
        gridsize = fscanf(fid, '%d %d %d', [3 1])';
        mag_x = fscanf(fid, '%f', [prod(gridsize,2) 1])';
        mag_x = reshape(mag_x,gridsize); mag_x = mag_x/vol;
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) > 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_x = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_x = reshape(mag_x,gridsize);
    if size(str2num(line),2) == 3; % CHG
        mag_x = mag_x/vol;
    end
end
fgetl(fid); pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'
    for i = 1:natoms; line = fgetl(fid);
        nentries = sscanf(line,['augmentation occupancies ' num2str(i) ' %d']); 
        fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
    end    
end
pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'; % MAGY
    for i = 1:natoms; line = fgetl(fid);
        nentries = sscanf(line,['augmentation occupancies (imaginary part) ' num2str(i) ' %d']); 
        fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
    end
    fscanf(fid, '%f', [natoms 1]); 
    fgetl(fid); pos = ftell(fid); fgetl(fid);
    if ~feof(fid); fseek(fid,pos,'bof');
        gridsize = fscanf(fid, '%d %d %d', [3 1])';
        mag_y = fscanf(fid, '%f', [prod(gridsize,2) 1])';
        mag_y = reshape(mag_y,gridsize); mag_y = mag_y/vol;
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) > 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_y = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_y = reshape(mag_y,gridsize);
    if size(str2num(line),2) == 3; % CHG
        mag_y = mag_y/vol;
    end
end
fgetl(fid); pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'
    for i = 1:natoms; line = fgetl(fid);
        nentries = sscanf(line,['augmentation occupancies ' num2str(i) ' %d']); 
        fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
    end    
end
pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'; % MAGZ
    for i = 1:natoms; line = fgetl(fid);
        nentries = sscanf(line,['augmentation occupancies (imaginary part) ' num2str(i) ' %d']); 
        fscanf(fid, '%f', [nentries 1]); fgetl(fid);
    end
    fscanf(fid, '%f', [natoms 1]);
    fgetl(fid); pos = ftell(fid); fgetl(fid);
    if ~feof(fid); fseek(fid,pos,'bof');
        gridsize = fscanf(fid, '%d %d %d', [3 1])';
        mag_z = fscanf(fid, '%f', [prod(gridsize,2) 1])';
        mag_z = reshape(mag_z,gridsize); mag_z = mag_z/vol;
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) > 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_z = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_z = reshape(mag_z,gridsize);
    if size(str2num(line),2) == 3; % CHG
        mag_z = mag_z/vol;
    end
end
fclose(fid);

function result = outcar(filename)
fid = fopen(filename);
while ~feof(fid); line = fgetl(fid);
    if numel(regexp(line,'magnetization \(x\)'))==1
        pos = ftell(fid);
    end
end
fseek(fid,pos,'bof'); fgetl(fid); line = fgetl(fid); fgetl(fid);
if length(line) > 40; sp = 6; else sp = 5; end
mx = fscanf(fid,'%f',[sp inf])'; mx = mx(:,end);
while ~feof(fid); line = fgetl(fid);
    if numel(regexp(line,'magnetization \(y\)'))==1
        pos = ftell(fid);
    end
end
fseek(fid,pos,'bof'); fgetl(fid); fgetl(fid); fgetl(fid);
my = fscanf(fid,'%f',[sp inf])'; my = my(:,end);
while ~feof(fid); line = fgetl(fid);
    if numel(regexp(line,'magnetization \(z\)'))==1
        pos = ftell(fid);
    end
end
fseek(fid,pos,'bof'); fgetl(fid); fgetl(fid); fgetl(fid);
mz = fscanf(fid,'%f',[sp inf])'; mz = mz(:,end);
fclose(fid); result = [mx my mz];

function atom(target,r,col,alp)
hold on
if size(target,1) < 100
    for i = 1:size(target,1)
        [rx,ry,rz] = sphere(10); rx = rx*r; ry = ry*r; rz = rz*r;
        h = surf((rx+target(i,1)),(ry+target(i,2)),(rz+target(i,3)));
        set(h,'edgecolor','none','facecolor',col,'facelighting','gouraud');
        alpha(h,alp)
    end
else
    for i = 1:size(target,1)
        plot3(target(i,1),target(i,2),target(i,3),[col '.'])
    end
end

function arrow(Point,Arg,col)
for i = 1:size(Point,1)
    ar = [Arg(i,1) Arg(i,2) Arg(i,3)]/sqrt(Arg(i,1)^2+Arg(i,2)^2+Arg(i,3)^2);
    if sum(ar) ~= 0
        arrow3(Point(i,:)-ar,Point(i,:)+1.75*ar,'color',col,'stemWidth',0.1,'tipWidth',0.25);
    end
end

function h = arrow3(p1,p2,varargin)                            
propertyNames = {'edgeColor'}; propertyValues = {'none'};
for argno = 1:2:nargin-2
    switch varargin{argno}
        case 'color'
            propertyNames = {propertyNames(:),'facecolor'};
            propertyValues = {propertyValues(:),varargin{argno+1}};
        case 'stemWidth'
            if isreal(varargin{argno+1}); stemWidth = varargin{argno+1};
            else warning('arrow3:stemWidth','stemWidth must be a real number');
            end
        case 'tipWidth'
            if isreal(varargin{argno+1}); tipWidth = varargin{argno+1};
            else warning('arrow3:tipWidth','tipWidth must be a real number');
            end
        otherwise
            propertyNames = {propertyNames(:),varargin{argno}};
            propertyValues = {propertyValues(:),varargin{argno+1}};
    end
end
if ~exist('stemWidth','var'); ax = axis;
    if numel(ax)==4
        stemWidth = norm(ax([2 4])-ax([1 3]))/300;
    elseif numel(ax)==6
        stemWidth = norm(ax([2 4 6])-ax([1 3 5]))/300;
    end
end
if ~exist('tipWidth','var'); tipWidth = 3*stemWidth; end
tipAngle = 22.5/180*pi; tipLength = tipWidth/tan(tipAngle/2);
ppsc = 50; ppbc = 250; p1 = p1(:); p2 = p2(:);
x = (p2-p1)/norm(p2-p1); y = cross(x,[0;0;1]);
if norm(y)<0.1; y = cross(x,[0;1;0]); end
y = y/norm(y); z = cross(x,y); z = z/norm(z);
theta = 0:2*pi/ppsc:2*pi; sintheta = sin(theta); costheta = cos(theta);
upsilon = 0:2*pi/ppbc:2*pi; sinupsilon = sin(upsilon); cosupsilon = cos(upsilon);
f = NaN([ppsc+ppbc+2 ppbc+1]);
if norm(p2-p1)>tipLength
    for idx = 1:ppsc+1
        v(idx,:) = p1 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
    end
    p3 = p2-tipLength*x;
    for idx = 1:ppsc+1
        v(ppsc+1+idx,:) = p3 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
    end
    for idx = 1:ppbc+1
        v(2*ppsc+2+idx,:) = p3 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
    end
    v(2*ppsc+ppbc+4,:) = p2; f(1,1:ppsc+1) = 1:ppsc+1;
    for idx = 1:ppsc; f(1+idx,1:4) = [idx idx+1 ppsc+1+idx+1 ppsc+1+idx]; end
    f(ppsc+2,:) = 2*ppsc+3:(2*ppsc+3)+ppbc;
    for idx = 1:ppbc
        f(ppsc+2+idx,1:3) = [2*ppsc+2+idx 2*ppsc+2+idx+1 2*ppsc+ppbc+4];
    end
else tipWidth = 2*sin(tipAngle/2)*norm(p2-p1);
    for idx = 1:ppbc+1
        v(idx,:) = p1 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
    end
    v(ppbc+2,:) = p2; f(1,:) = 1:ppbc+1;
    for idx = 1:ppbc; f(1+idx,1:3) = [idx idx+1 ppbc+2]; end
end
fv.faces = f; fv.vertices = v; h = patch(fv);
for propno = 1:numel(propertyNames)
    try set(h,propertyNames{propno},propertyValues{propno});
    catch 
        disp(lasterr)
    end
end

function drbox(latt,O,col)
if nargin == 1; O = [0 0 0]; col = 'k'; end
if nargin == 2; if isnumeric(O); col = 'k'; else col = O; O = [0 0 0]; end; end
line([0 latt(1,1)]+O(1),[0 latt(1,2)]+O(2),[0 latt(1,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(2,1) latt(1,1)+latt(2,1)]+O(1),[latt(2,2) latt(1,2)+latt(2,2)]+O(2), ...
    [latt(2,3) latt(1,3)+latt(2,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(3,1) latt(1,1)+latt(3,1)]+O(1),[latt(3,2) latt(1,2)+latt(3,2)]+O(2), ...
    [latt(3,3) latt(1,3)+latt(3,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(2,1)+latt(3,1) latt(1,1)+latt(2,1)+latt(3,1)]+O(1), ...
    [latt(2,2)+latt(3,2) latt(1,2)+latt(2,2)+latt(3,2)]+O(2), ...
    [latt(2,3)+latt(3,3) latt(1,3)+latt(2,3)+latt(3,3)]+O(3),'LineWidth',2,'Color',col)
line([0 latt(2,1)]+O(1),[0 latt(2,2)]+O(2),[0 latt(2,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(1,1) latt(1,1)+latt(2,1)]+O(1),[latt(1,2) latt(1,2)+latt(2,2)]+O(2), ...
    [latt(1,3) latt(1,3)+latt(2,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(3,1) latt(2,1)+latt(3,1)]+O(1),[latt(3,2) latt(2,2)+latt(3,2)]+O(2), ...
    [latt(3,3) latt(2,3)+latt(3,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(1,1)+latt(3,1) latt(1,1)+latt(2,1)+latt(3,1)]+O(1), ...
    [latt(1,2)+latt(3,2) latt(1,2)+latt(2,2)+latt(3,2)]+O(2), ...
    [latt(1,3)+latt(3,3) latt(1,3)+latt(2,3)+latt(3,3)]+O(3),'LineWidth',2,'Color',col)
line([0 latt(3,1)]+O(1),[0 latt(3,2)]+O(2),[0 latt(3,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(1,1) latt(1,1)+latt(3,1)]+O(1),[latt(1,2) latt(1,2)+latt(3,2)]+O(2), ...
    [latt(1,3) latt(1,3)+latt(3,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(2,1) latt(2,1)+latt(3,1)]+O(1),[latt(2,2) latt(2,2)+latt(3,2)]+O(2), ...
    [latt(2,3) latt(2,3)+latt(3,3)]+O(3),'LineWidth',2,'Color',col)
line([latt(1,1)+latt(2,1) latt(1,1)+latt(2,1)+latt(3,1)]+O(1), ...
    [latt(1,2)+latt(2,2) latt(1,2)+latt(2,2)+latt(3,2)]+O(2), ...
    [latt(1,3)+latt(2,3) latt(1,3)+latt(2,3)+latt(3,3)]+O(3),'LineWidth',2,'Color',col)
