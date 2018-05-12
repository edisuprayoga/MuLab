function varargout = GUI_Relaxation(varargin)

% GUI_Relaxation MATLAB code by Edi Suprayoga
% feel free to contact me at suprayoga.edi@gmail.com
% last edit 3 Mar 2016

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Relaxation_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Relaxation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout; [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else gui_mainfcn(gui_State, varargin{:});
end

function GUI_Relaxation_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject; guidata(hObject, handles);
set(handles.figure1,'Color',[0.8 0.8 0.8])

function varargout = GUI_Relaxation_OutputFcn(~, ~, handles) 
varargout{1} = handles.output; set(handles.autorized,'Value',1);

function ion_no_Callback(~, ~, ~)
function ion_no_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ion_rad_Callback(~, ~, ~)
function ion_rad_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ion_col_Callback(~, ~, ~)
function ion_col_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xb_Callback(~, ~, ~)
function xb_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yb_Callback(~, ~, ~)
function yb_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zb_Callback(~, ~, ~)
function zb_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bond_l_Callback(~, ~, ~)
function bond_l_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bond_Callback(~, ~, handles)
if get(handles.bond,'Value') == 1
    set(handles.bond_l,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.bond_l,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.detil_o,'Value',0); set(handles.detil_b,'Value',0)
else
    set(handles.bond_l,'ForegroundColor',[0 0 0])
    set(handles.bond_l,'BackgroundColor',[1 1 1])
    if get(handles.bond,'Value') == 2; set(handles.detil_o,'Value',0);
    elseif get(handles.bond,'Value') == 3; set(handles.detil_b,'Value',0);
    end
end

function bond_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cek_no_Callback(~, ~, handles)
if get(handles.opt,'Value') == 3
    set(handles.cek_no,'Value',1)
    set(handles.no,'BackgroundColor',[1 1 1])
    set(handles.no,'ForegroundColor',[0 0 0])
elseif get(handles.cek_no,'Value') == 1
    set(handles.no,'BackgroundColor',[1 1 1])
    set(handles.no,'ForegroundColor',[0 0 0])
else set(handles.no,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.no,'ForegroundColor',[0.8 0.8 0.8])
end

function load1_Callback(~, ~, handles)
global filename pathname
[filename, pathname] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname)
    if get(handles.opt,'Value') >= 2
        set(handles.pathname1,'String',[pathname filename])
    else set(handles.pathname1,'String',pathname)
    end
end

function opt_Callback(~, ~, handles)
if get(handles.opt,'Value') == 2
    set(handles.pathname2,'BackgroundColor',[1 1 1])
    set(handles.pathname2,'ForegroundColor',[0 0 0])
    set(handles.load2,'String','CONTCAR'); set(handles.load1,'String','POSCAR')
else set(handles.load2,'String','')
    set(handles.pathname2,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.pathname2,'ForegroundColor',[0.8 0.8 0.8])
    if get(handles.opt,'Value') == 1
        set(handles.load1,'String','Folder')
    else set(handles.load1,'String','Filename')
    end
end
if get(handles.opt,'Value') == 3; set(handles.cek_no,'Value',1)
    set(handles.cek_no,'String','Number of files being calculated:')
    set(handles.no,'ForegroundColor',[0 0 0]); set(handles.no,'BackgroundColor',[1 1 1])
else set(handles.cek_no,'String','Distance between Muon and Ion no:')
    set(handles.cek_no,'Value',0); set(handles.no,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.no,'BackgroundColor',[0.8 0.8 0.8])
end

function opt_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pathname1_Callback(~, ~, ~)
function pathname1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pathname2_Callback(~, ~, ~)
function pathname2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load2_Callback(~, ~, handles)
global filename2 pathname2
if get(handles.opt,'Value') == 2
    [filename2, pathname2] = uigetfile({'*.*', 'All Files (*.*)'});
    if ~isnumeric(pathname2)
        set(handles.pathname2,'String',[pathname2 filename2])
    end
end

function run_Callback(~, ~, handles)
pathname = get(handles.pathname1,'String'); tic
if get(handles.cls,'Value') == 1; clc; end
if get(handles.opt,'Value') == 1
    geo1 = poscar([pathname '/POSCAR']); 
    geo2 = poscar([pathname '/CONTCAR']);
elseif get(handles.opt,'Value') == 2
    pathname2 = get(handles.pathname2,'String');
    geo1 = poscar(pathname); geo2 = poscar(pathname2);
elseif get(handles.opt,'Value') == 3; geo1 = poscar([pathname '1']);
    aa = str2num(get(handles.no,'String'));
    geo2 = poscar([pathname num2str(aa(end))]);
end
latt1 = geo1.lattice; latt2 = geo2.lattice;
latt = latt1(1,:)+latt1(2,:)+latt1(3,:);
d1 = geo1.coords; d2 = geo2.coords;
for i = 1:size(d1,1)
    d1(i,:) = d1(i,1)*latt1(1,:)+d1(i,2)*latt1(2,:)+d1(i,3)*latt1(3,:);
    d2(i,:) = d2(i,1)*latt2(1,:)+d2(i,2)*latt2(2,:)+d2(i,3)*latt2(3,:);
end
d = abs(d1-d2); out = zeros(size(d,1),1);
for i = 1:size(d,1)
    for j = 1:3
        if abs(d(i,j)-latt(j)) < 2; d(i,j) = d(i,j)-latt(j); end
    end
    out(i) = sqrt(dot(d(i,:),d(i,:))); 
end
x = sqrt(mean(out.*out)); d_mu = out(end,:); indx = [1 2 3 1 2];
vol1 = abs(dot(geo1.lattice(1,:),cross(geo1.lattice(2,:),geo1.lattice(3,:))));
vol2 = abs(dot(geo2.lattice(1,:),cross(geo2.lattice(2,:),geo2.lattice(3,:))));
d_vol = (vol2/vol1-1)*100; a1 = zeros(1,3); a2 = a1; alp1 = a1; alp2 = a1;
for i = 1:3
    a1(i) = sqrt(latt1(i,1)^2+latt1(i,2)^2+latt1(i,3)^2);
    a2(i) = sqrt(latt2(i,1)^2+latt2(i,2)^2+latt2(i,3)^2);
end
for i = 1:3
    alp1(i) = acosd(dot(latt1(indx(i+1),:),latt1(indx(i+2),:))/(a1(indx(i+1))*a1(indx(i+2))));
    alp2(i) = acosd(dot(latt2(indx(i+1),:),latt2(indx(i+2),:))/(a2(indx(i+1))*a2(indx(i+2))));
end
if get(handles.cls,'Value') == 0; fprintf('\n'); end
fprintf('   ============================================================== \n')
fprintf('      a      b      c       alpha   beta    gamma     Vol(A^3) \n')
fprintf('   %4.4f %4.4f %4.4f   %4.4f %4.4f %4.4f   %4.4f \n',a1,alp1,vol1)
fprintf('   %4.4f %4.4f %4.4f   %4.4f %4.4f %4.4f   %4.4f \n',a2,alp2,vol2)
if max(out) == 0; fprintf('   Both crystals are exactly the same. \n')
else
    if d_vol == 0; fprintf('   Crystal volume doesnt change. \n')
    else fprintf('\n   Shifting Volume\t   = %4.4f %%',d_vol)
    end
    fprintf('\n   Average displacement\t   = %4.4f Angs \n',x)
    fprintf('   Range of displacement   = %4.4f - %4.4f Angs \n',min(out), max(out))
end
if geo1.atomcount(end) == 1
    fprintf('\n   Muon position 1\t   = (%4.4f %4.4f %4.4f) \n',geo1.coords(end,:))
    fprintf('   Muon position 2\t   = (%4.4f %4.4f %4.4f) \n',geo2.coords(end,:))
    fprintf('   Shifting Muon position  = %4.4f Angs \n',d_mu)
    Mu1 = geo1.coords(end,:); Mu2 = geo2.coords(end,:);
    Mu1 = Mu1(1)*latt1(1,:)+Mu1(2)*latt1(2,:)+Mu1(3)*latt1(3,:);
    Mu2 = Mu2(1)*latt2(1,:)+Mu2(2)*latt2(2,:)+Mu2(3)*latt2(3,:);
    if get(handles.cek_no,'Value') == 1 && get(handles.opt,'Value') < 3;
        no = str2double(get(handles.no,'String'));
        Cu1 = geo1.coords(sum(geo1.atomcount(1:no-1))+1:sum(geo1.atomcount(1:no)),:);
        Cu2 = geo2.coords(sum(geo2.atomcount(1:no-1))+1:sum(geo2.atomcount(1:no)),:);
        out1 = zeros(size(Cu1,1),1); out2 = zeros(size(Cu2,1),1);
        for i = 1:size(Cu1,1)
            C1 = (Cu1(i,1)*latt1(1,:)+Cu1(i,2)*latt1(2,:)+Cu1(i,3)*latt1(3,:))-Mu1;
            C2 = (Cu2(i,1)*latt2(1,:)+Cu2(i,2)*latt2(2,:)+Cu2(i,3)*latt2(3,:))-Mu2;
            out1(i) = sqrt(dot(C1,C1)); out2(i) = sqrt(dot(C2,C2));
        end
        [C1, p1] = min(out1); [C2, p2] = min(out2); 
        dNu = Cu2(p2,:) - Cu1(p1,:); dNu = sqrt(dot(dNu,dNu));
        fprintf('\n   Shift Nuclear position  = %4.4f Angs \n',dNu)
        fprintf('   Distance Muon - Nuclear = %4.4f & %4.4f Angs \n',min(out1),min(out2))
        fprintf('   Shifting Muon - Nuclear = %4.4f Angs \n',C2-C1)
    end
end
fprintf('\n>> Drawing figure \n')
No = str2num(get(handles.ion_no,'String')); col = get(handles.ion_col,'String');
xb = str2num(get(handles.xb,'String')); rd = str2num(get(handles.ion_rad,'String'));
yb = str2num(get(handles.yb,'String')); zb = str2num(get(handles.zb,'String'));
if xb(1)<0; xb(1)=0; end; if xb(2)>1; xb(2)=1; end; if xb(2)<=xb(1); xa=xb; xb=[xa(2) xa(1)]; end
if yb(1)<0; yb(1)=0; end; if yb(2)>1; yb(2)=1; end; if yb(2)<=yb(1); ya=yb; yb=[ya(2) ya(1)]; end
if zb(1)<0; zb(1)=0; end; if zb(2)>1; zb(2)=1; end; if zb(2)<=zb(1); za=zb; xb=[za(2) za(1)]; end
if get(handles.bond,'Value') > 1; d0 = str2num(get(handles.bond_l,'String')); end
if get(handles.opt,'Value') < 3; num = 1; else num = str2num(get(handles.no,'String')); end
if length(No) == 1; col_oct = col(1); else col_oct = col(2);end
if length(num) == 1; num = [1 num]; end; mn = 1;

for m = num(1):num(2)
    if get(handles.opt,'Value') == 3;
        geo2 = poscar([pathname num2str(m)]); latt2 = geo2.lattice; h = figure(m);
        Mu2 = geo2.coords(end,:); Mu2 = Mu2(1)*latt1(1,:)+Mu2(2)*latt1(2,:)+Mu2(3)*latt1(3,:);
    elseif get(handles.cls,'Value') == 1; figure(1), clf; else figure;
    end
    for no = 1:length(No); Cu1 = []; Cu2 = [];
        cu1 = geo1.coords(sum(geo1.atomcount(1:No(no)-1))+1:sum(geo1.atomcount(1:No(no))),:);
        cu2 = geo2.coords(sum(geo2.atomcount(1:No(no)-1))+1:sum(geo2.atomcount(1:No(no))),:);
        for n = 1:size(cu1,1)
            for i = -1:1
                for j = -1:1
                    for k = -1:1
                        if get(handles.bond,'Value') > 1  && no <= 2
                            latt = inv(latt1); d1 = (d0(1)-1)*(latt(1,:)+latt(2,:)+latt(3,:));
                            latt = inv(latt2); d2 = (d0(1)-1)*(latt(1,:)+latt(2,:)+latt(3,:));
                        else d1 = 0.02*[1 1 1]; d2 = 0.02*[1 1 1]; 
                        end
                        if cu1(n,1)+i >= xb(1)-d1(1) && cu1(n,1)+i <= xb(2)+d1(1) ...
                                && cu1(n,2)+j >= yb(1)-d1(2) && cu1(n,2)+j <= yb(2)+d1(2) ...
                                && cu1(n,3)+k >= zb(1)-d1(3) && cu1(n,3)+k <= zb(2)+d1(3)
                            Cu1 = [Cu1; (cu1(n,1)+i)*latt1(1,:)+(cu1(n,2)+j)*latt1(2,:)+(cu1(n,3)+k)*latt1(3,:)];
                        end
                        if cu2(n,1)+i > xb(1)-d2(1) && cu2(n,1)+i < xb(2)+d2(1) ...
                                && cu2(n,2)+j > yb(1)-d2(2) && cu2(n,2)+j < yb(2)+d2(2) ...
                                && cu2(n,3)+k > zb(1)-d2(3) && cu2(n,3)+k < zb(2)+d2(3)
                            Cu2 = [Cu2; (cu2(n,1)+i)*latt2(1,:)+(cu2(n,2)+j)*latt2(2,:)+(cu2(n,3)+k)*latt2(3,:)];
                        end
                    end
                end
            end
        end
        atom(Cu2,rd(no),col(no),1)
        if get(handles.cek_ion,'Value') == 1; atom(Cu1,rd(no),col(no),0.3); end
        if geo1.atomcount(end) == 1; atom(Mu2,0.2,'r',1);
            if get(handles.cek_mu,'Value') == 1; atom(Mu1,0.2,'r',0.3); end
        end
        if get(handles.bond,'Value') == 2 && no <= 2
            if no == 1; ion2 = Cu2; ion1 = Cu1;
            elseif no == 2; drbond(ion2,Cu2,d0,col(1:2),2);
                if get(handles.detil_b,'Value') == 1; drbond(ion1,Cu1,d0,'cc',1.5); end
            end
        elseif get(handles.bond,'Value') == 3 && no == 1; drplane(Cu2,d0,col_oct,[0.5 0.5 0.5],0.2)
            if get(handles.detil_o,'Value') == 1; drplane(Cu1,d0,[1 1 1],'m',0); end
        elseif get(handles.bond,'Value') == 4 && no <= 2
            if no == 1; drplane(Cu2,d0,col_oct,[0.5 0.5 0.5],0.2); ion2 = Cu2; ion1 = Cu1;
                if get(handles.detil_o,'Value') == 1; drplane(Cu1,d0,[1 1 1],'m',0); end
            elseif no == 2; drbond(ion2,Cu2,d0,col(1:2),2);
                if get(handles.detil_b,'Value') == 1; drbond(ion1,Cu1,d0,'cc',1.5); end
            end
        end
    end
    if get(handles.cek_orientation,'Value') == 1; ori = str2num(get(handles.ori,'String')); 
        if length(ori) == 1; ori = ori*[1 1]; end; view(ori(1), ori(2))
    end
    camlight right local; rotate3d on; mn = mn+1; daspect([1 1 1])
    if get(handles.box,'Value') == 1
        drbox(latt1,[0 0 0],'b',1); drbox(latt2,[0 0 0],'k',2);
    end
    if get(handles.def,'Value') == 1; axis tight
        xlabel('a','FontSize',12); ylabel('b','FontSize',12); zlabel('c','FontSize',12)
    else axis off
    end
    if get(handles.opt,'Value') == 3; print(h,'-dpng',[pathname num2str(m)]); 
        if get(handles.cls,'Value') == 1; close(h); end
    end
end
fprintf('   ============================================================== \n   '); toc

function def_Callback(~, ~, ~)
function no_Callback(~, ~, ~)
function no_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function geometry = poscar(filename)
if ~isnumeric(filename)
	geometry.filename = filename; fid = fopen(filename);
else fid = filename; geometry.filename = fopen(fid);
end
geometry.comment = fgetl(fid); scale = fscanf(fid, '%f',1); 
geometry.lattice = fscanf(fid, '%f %f %f', [3 3])'; 
geometry.lattice = geometry.lattice*scale;
fgetl(fid); line = fgetl(fid); has_symbols_line = false;
if sum(isstrprop(line, 'digit')) == 0
    geometry.symbols = regexp(line, '([^ ]*)', 'match');
    line = fgetl(fid); has_symbols_line = true;
else geometry.symbols = {};
end
geometry.atomcount = sscanf(line,'%d');
natoms = sum(geometry.atomcount); cartesian = 0;
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
            str = regexp(line, '[^ ]*', 'match'); str = str{end};
            newelement = false; 
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

function atom(target,r,col,alp); hold on
for i = 1:size(target,1)
    [rx,ry,rz] = sphere(10); rx = rx*r; ry = ry*r; rz = rz*r;
    h = surf((rx+target(i,1)),(ry+target(i,2)),(rz+target(i,3)));
    set(h,'edgecolor','none','facecolor',col,'facelighting','gouraud'); alpha(h,alp)
end

function drplane(atom,d,col_oct,col_ln,alp)
if length(d) == 1; d = [0 d]; end; hold on
for i = 1:size(atom,1)-2
    for j = i+1:size(atom,1)-1
        d1 = sqrt((atom(i,1)-atom(j,1))^2+(atom(i,2)-atom(j,2))^2+...
            (atom(i,3)-atom(j,3))^2);
        if d1 >= d(1) && d1 <= d(2)
            for k = j+1:size(atom,1)
                d2 = sqrt((atom(j,1)-atom(k,1))^2+(atom(j,2)-atom(k,2))^2+...
                    (atom(j,3)-atom(k,3))^2);
                if d2 >= d(1) && d2 <= d(2)
                    d3 = sqrt((atom(i,1)-atom(k,1))^2+(atom(i,2)-atom(k,2))^2+...
                        (atom(i,3)-atom(k,3))^2);
                    if d3 >= d(1) && d3 <= d(2)
                        X = [atom(i,1) atom(j,1) atom(k,1)];
                        Y = [atom(i,2) atom(j,2) atom(k,2)];
                        Z = [atom(i,3) atom(j,3) atom(k,3)];
                        h = fill3(X,Y,Z,col_oct,'EdgeColor',col_ln); alpha(h,alp)
                    end
                end 
            end 
        end  
    end
end

function drbond(atom1,atom2,d0,col,width)
if nargin < 5; col = 'kk'; width = 1; end
if length(d0) < 2; d0 = [0 d0]; end; hold on
if atom1(1,:) == atom2(1,:); a = 1; else a = 0; end
for i = 1:size(atom1,1) - a
    for j = 1:size(atom2,1)
        d = sqrt((atom1(i,1)-atom2(j,1))^2+(atom1(i,2)-atom2(j,2))^2+...
            (atom1(i,3)-atom2(j,3))^2);
        if d >= d0(1) && d <= d0(2)
            R = atom2(j,:)-atom1(i,:);
            mid = atom1(i,:)+R/2;
            line([atom1(i,1) mid(1)],[atom1(i,2) mid(2)],...
                [atom1(i,3) mid(3)],'LineWidth',width,'Color',col(1))
            line([mid(1) atom2(j,1)],[mid(2) atom2(j,2)],...
                [mid(3) atom2(j,3)],'LineWidth',width,'Color',col(2))
        end
    end
end

function drbox(latt,O,col,d)
if nargin == 1; O = [0 0 0]; col = 'k'; d = 2; end
if nargin == 2; d = 2; if isnumeric(O); col = 'k'; else col = O; O = [0 0 0]; end; end
line([0 latt(1,1)]+O(1),[0 latt(1,2)]+O(2),[0 latt(1,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(2,1) latt(1,1)+latt(2,1)]+O(1),[latt(2,2) latt(1,2)+latt(2,2)]+O(2), ...
    [latt(2,3) latt(1,3)+latt(2,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(3,1) latt(1,1)+latt(3,1)]+O(1),[latt(3,2) latt(1,2)+latt(3,2)]+O(2), ...
    [latt(3,3) latt(1,3)+latt(3,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(2,1)+latt(3,1) latt(1,1)+latt(2,1)+latt(3,1)]+O(1), ...
    [latt(2,2)+latt(3,2) latt(1,2)+latt(2,2)+latt(3,2)]+O(2), ...
    [latt(2,3)+latt(3,3) latt(1,3)+latt(2,3)+latt(3,3)]+O(3),'LineWidth',d,'Color',col)
line([0 latt(2,1)]+O(1),[0 latt(2,2)]+O(2),[0 latt(2,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(1,1) latt(1,1)+latt(2,1)]+O(1),[latt(1,2) latt(1,2)+latt(2,2)]+O(2), ...
    [latt(1,3) latt(1,3)+latt(2,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(3,1) latt(2,1)+latt(3,1)]+O(1),[latt(3,2) latt(2,2)+latt(3,2)]+O(2), ...
    [latt(3,3) latt(2,3)+latt(3,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(1,1)+latt(3,1) latt(1,1)+latt(2,1)+latt(3,1)]+O(1), ...
    [latt(1,2)+latt(3,2) latt(1,2)+latt(2,2)+latt(3,2)]+O(2), ...
    [latt(1,3)+latt(3,3) latt(1,3)+latt(2,3)+latt(3,3)]+O(3),'LineWidth',d,'Color',col)
line([0 latt(3,1)]+O(1),[0 latt(3,2)]+O(2),[0 latt(3,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(1,1) latt(1,1)+latt(3,1)]+O(1),[latt(1,2) latt(1,2)+latt(3,2)]+O(2), ...
    [latt(1,3) latt(1,3)+latt(3,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(2,1) latt(2,1)+latt(3,1)]+O(1),[latt(2,2) latt(2,2)+latt(3,2)]+O(2), ...
    [latt(2,3) latt(2,3)+latt(3,3)]+O(3),'LineWidth',d,'Color',col)
line([latt(1,1)+latt(2,1) latt(1,1)+latt(2,1)+latt(3,1)]+O(1), ...
    [latt(1,2)+latt(2,2) latt(1,2)+latt(2,2)+latt(3,2)]+O(2), ...
    [latt(1,3)+latt(2,3) latt(1,3)+latt(2,3)+latt(3,3)]+O(3),'LineWidth',d,'Color',col)

function box_Callback(~, ~, ~)
function detil_b_Callback(~, ~, handles)
if get(handles.bond,'Value') == 1 || get(handles.bond,'Value') == 3
    set(handles.detil_b,'Value',0); 
end

function cls_Callback(~, ~, ~)
function detil_o_Callback(~, ~, handles)
if get(handles.bond,'Value') <= 2; set(handles.detil_o,'Value',0); 
end

function ori_Callback(~, ~, ~)
function ori_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cek_orientation_Callback(~, ~, handles)
if get(handles.cek_orientation,'Value') == 1
    set(handles.ori,'ForegroundColor',[0 0 0])
    set(handles.ori,'BackgroundColor',[1 1 1])
else set(handles.ori,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.ori,'BackgroundColor',[0.8 0.8 0.8])
end

function cek_ion_Callback(~, ~, ~)
function cek_mu_Callback(~, ~, ~)
