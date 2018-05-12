function varargout = GUI_VESTA(varargin)

% GUI_VESTA MATLAB code by Edi Suprayoga
% feel free to contact me at suprayoga.edi@gmail.com
% last edit 3 Mar 2016

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_VESTA_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_VESTA_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout; [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else gui_mainfcn(gui_State, varargin{:});
end

function GUI_VESTA_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject; guidata(hObject, handles);
set(handles.figure1,'Color',[0.8 0.8 0.8])

function varargout = GUI_VESTA_OutputFcn(~, ~, handles) 
varargout{1} = handles.output; set(handles.autorized,'Value',1);

function filename_Callback(~, ~, ~)
function filename_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load_file_Callback(~, ~, handles)
global filename pathname
[filename, pathname] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname)
    set(handles.filename,'String',[pathname filename])
end

function opt_file_Callback(~, ~, handles)
if (get(handles.opt_file,'Value') == 1); set(handles.load_file,'String','LOCPOT')
elseif (get(handles.opt_file,'Value') == 6); set(handles.load_file,'String','ELF')
elseif (get(handles.opt_file,'Value') == 7); set(handles.load_file,'String','DIPOLE')
else set(handles.load_file,'String','CHGCAR')
end
if (get(handles.opt_file,'Value') == 5); set(handles.H,'String','S')
    set(handles.Hx,'String','Sx'); set(handles.Hy,'String','Sy'); set(handles.Hz,'String','Sz')
elseif (get(handles.opt_file,'Value') == 7); set(handles.H,'String','H')
    set(handles.Hx,'String','Hx'); set(handles.Hy,'String','Hy'); set(handles.Hz,'String','Hz')
end

function opt_file_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_plot_Callback(~, ~, handles)
if get(handles.opt_plot,'Value') == 2 || get(handles.opt_plot,'Value') == 3
    set(handles.lbl,'String','Label')
    set(handles.cek_ax,'String','Level')
    set(handles.c_axis,'String','10')
elseif get(handles.opt_plot,'Value') > 4
    set(handles.lbl,'String','Colour')
    set(handles.cek_ax,'String','Scale')
    set(handles.c_axis,'String','1')
else set(handles.cek_ax,'String','Axis (%)')
    set(handles.c_axis,'String','0 100')
    set(handles.lbl,'Value',0)
end

function opt_plot_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Atom_no_Callback(~, ~, ~)
function Atom_no_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function run_2D_Callback(~, ~, handles)
at = str2num(get(handles.level,'String')); at2 = 0; clc; tic
if get(handles.super,'Value') == 1 && length(at) > 1
    at = at - [0.5 0.5 0.5];
    for i = 1:3
        if at(i) < 0; at(i) = at(i) + 1; end
    end
end
fprintf('   ===================================== \n')
fprintf('   Read input files. \n')
filename = get(handles.filename,'String');

if (get(handles.opt_file,'Value') == 1) % locpot
    [rho, geo] = locpot(filename); rho = -rho;
    
elseif (get(handles.opt_file,'Value') == 2) % chg
    [rho, ~, geo] = chgcar(filename);
    
elseif (get(handles.opt_file,'Value') == 3) % chg gradien
    [rho, ~, geo] = chgcar(filename);
    [dx,dy,dz] = gradient(rho,0.01,0.01,0.01);
    rho = abs(sqrt(dx.^2+dy.^2+dz.^2))./rho.^(4/3);
    
elseif (get(handles.opt_file,'Value') == 4) % spin coll
    [~, rho, geo] = chgcar(filename);
    
elseif (get(handles.opt_file,'Value') == 5) % spin non-coll
    [~, mag_x, mag_y, mag_z, geo] = magcar(filename);
    rho = sqrt(mag_x.^2+mag_y.^2+mag_z.^2);
    
    if get(handles.Habs,'Value') == 1
        if get(handles.Hx,'Value') == 1; rho = abs(mag_x);
        elseif get(handles.Hy,'Value') == 1; rho = abs(mag_y);
        elseif get(handles.Hz,'Value') == 1; rho = abs(mag_z);
        end
    elseif get(handles.Hx,'Value') == 1; rho = mag_x;
    elseif get(handles.Hy,'Value') == 1; rho = mag_y;
    elseif get(handles.Hz,'Value') == 1; rho = mag_z;    
    end
    mag_x(end+1,:,:) = mag_x(1,:,:); mag_x(:,end+1,:) = mag_x(:,1,:); mag_x(:,:,end+1) = mag_x(:,:,1); 
    mag_y(end+1,:,:) = mag_y(1,:,:); mag_y(:,end+1,:) = mag_y(:,1,:); mag_y(:,:,end+1) = mag_y(:,:,1); 
    mag_z(end+1,:,:) = mag_z(1,:,:); mag_z(:,end+1,:) = mag_z(:,1,:); mag_z(:,:,end+1) = mag_z(:,:,1); 
    
elseif (get(handles.opt_file,'Value') == 6) % ELF
    [rho, geo] = locpot(filename);
    
elseif (get(handles.opt_file,'Value') == 7) % Dipole
    [rho, mag_x, mag_y, mag_z, geo] = dipcar(filename);
    
    if get(handles.Habs,'Value') == 1
        if get(handles.Hx,'Value') == 1; rho = abs(mag_x);
        elseif get(handles.Hy,'Value') == 1; rho = abs(mag_y);
        elseif get(handles.Hz,'Value') == 1; rho = abs(mag_z);
        end
    elseif get(handles.Hx,'Value') == 1; rho = mag_x;
    elseif get(handles.Hy,'Value') == 1; rho = mag_y;
    elseif get(handles.Hz,'Value') == 1; rho = mag_z;
    elseif get(handles.Hxy,'Value') == 1; rho = sqrt(mag_x.^2+mag_y.^2);
%     else rho = rho*5;
    end
    mag_x(end+1,:,:) = mag_x(1,:,:); mag_x(:,end+1,:) = mag_x(:,1,:); mag_x(:,:,end+1) = mag_x(:,:,1); 
    mag_y(end+1,:,:) = mag_y(1,:,:); mag_y(:,end+1,:) = mag_y(:,1,:); mag_y(:,:,end+1) = mag_y(:,:,1); 
    mag_z(end+1,:,:) = mag_z(1,:,:); mag_z(:,end+1,:) = mag_z(:,1,:); mag_z(:,:,end+1) = mag_z(:,:,1); 
end
rho(end+1,:,:) = rho(1,:,:); rho(:,end+1,:) = rho(:,1,:); rho(:,:,end+1) = rho(:,:,1); 
[dx,dy,dz] = gradient(rho,0.01,0.01,0.01);
latt = geo.lattice; lattice = inv(geo.lattice);
fprintf('   Drawing 2D Figure. \n')
figure; hold on;

if get(handles.bot_xy,'Value') == 1 % xy
    x_label = 'a'; y_label = 'b';
    if length(at) > 1; at2 = at; at = at(3);
        if (get(handles.opt_axis_2D,'Value') == 1); 
            at2 = at2(1)*latt(1,:)+at2(2)*latt(2,:)+at2(3)*latt(3,:); 
        end
        at2 = [at2(1), at2(2)];
    end
    
    if (get(handles.opt_axis_2D,'Value') == 2); at = at*lattice(3,3); end
    c = linspace(0,1,size(rho,3)); at = abs(c-at); [~,at] = min(at);
    latt = latt(1,:)+latt(2,:)+latt(3,:);
    nx = linspace(0,latt(1),size(rho,1)); ny = linspace(0,latt(2),size(rho,2));
    
    if (get(handles.opt_file,'Value') == 5) || (get(handles.opt_file,'Value') == 7) 
        Vx(:,:) = mag_x(:,:,at); Vy(:,:) = mag_y(:,:,at);
    else Vx(:,:) = dx(:,:,at); Vy(:,:) = dy(:,:,at);
    end
    V(:,:) = rho(:,:,at);
    
elseif get(handles.bot_xz,'Value') == 1 % xz
    x_label = 'a'; y_label = 'c';
    if length(at) > 1; at2 = at; at = at(2);
        if (get(handles.opt_axis_2D,'Value') == 1); 
            at2 = at2(1)*latt(1,:)+at2(2)*latt(2,:)+at2(3)*latt(3,:); 
        end
        at2 = [at2(1), at2(3)];
    end
    
    if (get(handles.opt_axis_2D,'Value') == 2); at = at*lattice(2,2); end
    c = linspace(0,1,size(rho,2)); at = abs(c-at); [~,at] = min(at);
    latt = latt(1,:)+latt(2,:)+latt(3,:); 
    nx = linspace(0,latt(1),size(rho,1)); ny = linspace(0,latt(3),size(rho,3));
    
    if (get(handles.opt_file,'Value') == 5) || (get(handles.opt_file,'Value') == 7) 
        Vx(:,:) = mag_x(:,at,:); Vy(:,:) = mag_z(:,at,:);
    else Vx(:,:) = dx(:,at,:); Vy(:,:) = dz(:,at,:);
    end
    V(:,:) = rho(:,at,:);
    
elseif get(handles.bot_yz,'Value') == 1 % yz
    x_label = 'b'; y_label = 'c';
    if length(at) > 1; at2 = at; at = at(1);
        if (get(handles.opt_axis_2D,'Value') == 1); 
            at2 = at2(1)*latt(1,:)+at2(2)*latt(2,:)+at2(3)*latt(3,:); 
        end
        at2 = [at2(2), at2(3)];
    end
    
    if (get(handles.opt_axis_2D,'Value') == 2); at = at*lattice(1,1); end
    c = linspace(0,1,size(rho,1)); at = abs(c-at); [~,at] = min(at);
    latt = latt(1,:)+latt(2,:)+latt(3,:);
    nx = linspace(0,latt(2),size(rho,2)); ny = linspace(0,latt(3),size(rho,3));
    
    if (get(handles.opt_file,'Value') == 5) || (get(handles.opt_file,'Value') == 7) 
        Vx(:,:) = mag_y(at,:,:); Vy(:,:) = mag_z(at,:,:);
    else Vx(:,:) = dy(at,:,:); Vy(:,:) = dz(at,:,:);
    end
    V(:,:) = rho(at,:,:);
end

if get(handles.cek_ax,'Value') == 1
    ax = str2num(get(handles.c_axis,'String'));
end

if get(handles.super,'Value') == 1
    if get(handles.bot_xy,'Value') == 1 
        nx = linspace(-0.5*latt(1),1.5*latt(1),2*size(rho,1)); ny = linspace(-0.5*latt(2),1.5*latt(2),2*size(rho,2));
    elseif get(handles.bot_xz,'Value') == 1 
        nx = linspace(-0.5*latt(1),1.5*latt(1),2*size(rho,1)); ny = linspace(-0.5*latt(3),1.5*latt(3),2*size(rho,3));
    elseif get(handles.bot_yz,'Value') == 1 
        nx = linspace(-0.5*latt(2),1.5*latt(2),2*size(rho,2)); ny = linspace(-0.5*latt(3),1.5*latt(3),2*size(rho,3));
    end
    Vc = V; sn = size(Vc);
    for i = 1:size(Vc,1); Vc(sn(1)+i,:) = V(i,:); end
    V = Vc;
    for i = 1:size(Vc,2); Vc(:,sn(2)+i) = V(:,i); end
    V = Vc;
    if (get(handles.opt_file,'Value') == 5) || (get(handles.opt_file,'Value') == 7) 
        Vcx = Vx; Vcy = Vy;
        for i = 1:size(Vcx,1); 
            Vcx(sn(1)+i,:) = Vx(i,:); Vcy(sn(1)+i,:) = Vy(i,:); 
        end
        Vx = Vcx;  Vy = Vcy;
        for i = 1:size(Vcx,2); 
            Vcx(:,sn(2)+i) = Vx(:,i); Vcy(:,sn(2)+i) = Vy(:,i); 
        end
        Vx = Vcx; Vy = Vcy;
    end
end

if get(handles.opt_plot,'Value') == 1 % surf
    surf(nx,ny,V'); view(-38,66); rotate3d on;
    if get(handles.cek_ax,'Value') == 1
        caxis([min(V(:))+ax(1)*max(V(:))/100 ax(2)*max(V(:))/100]);
    end
    
elseif get(handles.opt_plot,'Value') == 2 % contour
    
    if get(handles.cek_ax,'Value') == 1
        [~,h] = contour(nx,ny,V',ax);
    else [~,h] = contour(nx,ny,V');
    end
    if get(handles.lbl,'Value') == 1
        set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
    end
    daspect([1 1 1]); box on
    
elseif get(handles.opt_plot,'Value') == 3 % contf
    if get(handles.cek_ax,'Value') == 1
        [~,h] = contourf(nx,ny,V',ax); daspect([1 1 1]); box on
    else [~,h] = contourf(nx,ny,V'); daspect([1 1 1]); box on
    end
    if get(handles.lbl,'Value') == 1
        set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
    end
    
elseif get(handles.opt_plot,'Value') == 4 % image
    
    imagesc(nx,ny,interp2(V',2)); daspect([1 1 1]); box on
    if get(handles.cek_ax,'Value') == 1
        caxis([min(V(:))+ax(1)*max(V(:))/100 ax(2)*max(V(:))/100]);
    end
    
elseif get(handles.opt_plot,'Value') == 5 % quiver
    [nx,ny] = ndgrid(nx,ny);
    if get(handles.lbl,'Value') == 0
        if get(handles.cek_ax,'Value') == 1
            quiver(nx,ny,Vx,Vy,ax);
        else quiver(nx,ny,Vx,Vy); 
        end
    elseif get(handles.lbl,'Value') == 1
        if get(handles.cek_ax,'Value') == 1
            quiverc(nx,ny,Vx,Vy,ax);
        else quiverc(nx,ny,Vx,Vy); 
        end
    end
    daspect([1 1 1]); axis tight; box on
    
elseif get(handles.opt_plot,'Value') == 6 % quiver + contour
    contour(nx,ny,V'); [nx,ny] = ndgrid(nx,ny); hold on
    if get(handles.lbl,'Value') == 0
        if get(handles.cek_ax,'Value') == 1
            quiver(nx,ny,Vx,Vy,ax);
        else quiver(nx,ny,Vx,Vy); 
        end
    elseif get(handles.lbl,'Value') == 1
        if get(handles.cek_ax,'Value') == 1
            quiverc(nx,ny,Vx,Vy,ax);
        else quiverc(nx,ny,Vx,Vy); 
        end
    end 
    daspect([1 1 1]); axis tight; box on
end

if length(at2) > 1 && get(handles.opt_plot,'Value') > 1
    plot(at2(1), at2(2), 'k', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize',8)
end

axis tight; xlabel(x_label); ylabel(y_label);
fprintf('   ===================================== \n   '); toc

function bot_xy_Callback(~, ~, handles)
set(handles.bot_xy,'Value',1); set(handles.bot_xy,'ForegroundColor',[0 0 1])
set(handles.bot_xz,'Value',0); set(handles.bot_xz,'ForegroundColor',[0 0 0])
set(handles.bot_yz,'Value',0); set(handles.bot_yz,'ForegroundColor',[0 0 0])
set(handles.status,'String','Z position:')

function bot_xz_Callback(~, ~, handles)
set(handles.bot_xy,'Value',0); set(handles.bot_xy,'ForegroundColor',[0 0 0])
set(handles.bot_xz,'Value',1); set(handles.bot_xz,'ForegroundColor',[0 0 1])
set(handles.bot_yz,'Value',0); set(handles.bot_yz,'ForegroundColor',[0 0 0])
set(handles.status,'String','Y position:')

function bot_yz_Callback(~, ~, handles)
set(handles.bot_xy,'Value',0); set(handles.bot_xy,'ForegroundColor',[0 0 0])
set(handles.bot_xz,'Value',0); set(handles.bot_xz,'ForegroundColor',[0 0 0])
set(handles.bot_yz,'Value',1); set(handles.bot_yz,'ForegroundColor',[0 0 1])
set(handles.status,'String','X position:')

function level_Callback(~, ~, ~)
function level_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [geometry] = get_poscar(filename)
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

function [locpot, geometry] = locpot(filename)
fid = fopen(filename);
geometry = get_poscar(fid); fgetl(fid);
gridsize = fscanf(fid, '%d %d %d', [3 1])';
locpot = fscanf(fid, '%f', [prod(gridsize,2) 1])';
locpot = reshape(locpot,gridsize); fclose(fid);

function [chg, mag, geo] = chgcar(filename)
fid = fopen(filename); geo = get_poscar(fid); mag = [];
vol = abs(dot(geo.lattice(1,:),cross(geo.lattice(2,:),geo.lattice(3,:))));
natoms = sum(geo.atomcount); fgetl(fid); 
gridsize = fscanf(fid, '%d %d %d', [3 1])';
chg = fscanf(fid, '%f', [prod(gridsize,2) 1])'; chg = reshape(chg,gridsize);
fgetl(fid); pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'; chg = chg/vol;
    for i = 1:natoms; line = fgetl(fid);
        nentries = sscanf(line,['augmentation occupancies ' num2str(i) ' %d']);
        fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
    end
    pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
    if line(1) == 'a';
        for i = 1:natoms; line = fgetl(fid);
            nentries = sscanf(line,['augmentation occupancies (imaginary part) ' num2str(i) ' %d']); 
            fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
        end
    end
    fscanf(fid, '%f', [natoms 1]); fgetl(fid); pos = ftell(fid); fgetl(fid);
    if ~feof(fid) ; fseek(fid,pos,'bof');
        gridsize = fscanf(fid, '%d %d %d', [3 1])';
        mag = fscanf(fid, '%f', [prod(gridsize,2) 1])';
        mag = reshape(mag,gridsize); mag = mag/vol;
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) ~= 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag = reshape(mag,gridsize);
    if size(str2num(line),2) == 3; % CHG
        mag = mag/vol; chg = chg/vol;
    end
end
fclose(fid);

function [chg, mag_x, mag_y, mag_z, geo] = magcar(filename)
fid = fopen(filename); geo = get_poscar(fid);
vol = abs(dot(geo.lattice(1,:),cross(geo.lattice(2,:),geo.lattice(3,:))));
natoms = sum(geo.atomcount); fgetl(fid); 
gridsize = fscanf(fid, '%d %d %d', [3 1])';
chg = fscanf(fid, '%f', [prod(gridsize,2) 1])';
chg = reshape(chg,gridsize); mag_x = []; mag_y = []; mag_z = [];
fgetl(fid); pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'; chg = chg/vol;
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
    fscanf(fid, '%f', [natoms 1]); fgetl(fid); pos = ftell(fid); fgetl(fid);
    if ~feof(fid); fseek(fid,pos,'bof');
        gridsize = fscanf(fid, '%d %d %d', [3 1])';
        mag_x = fscanf(fid, '%f', [prod(gridsize,2) 1])';
        mag_x = reshape(mag_x,gridsize); mag_x = mag_x/vol;
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) ~= 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_x = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_x = reshape(mag_x,gridsize);
    if size(str2num(line),2) == 3; mag_x = mag_x/vol; end
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
    fscanf(fid, '%f', [natoms 1]); fgetl(fid); pos = ftell(fid); fgetl(fid);
    if ~feof(fid); fseek(fid,pos,'bof');
        gridsize = fscanf(fid, '%d %d %d', [3 1])';
        mag_y = fscanf(fid, '%f', [prod(gridsize,2) 1])';
        mag_y = reshape(mag_y,gridsize); mag_y = mag_y/vol;
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) ~= 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_y = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_y = reshape(mag_y,gridsize);
    if size(str2num(line),2) == 3; mag_y = mag_y/vol; end
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
    fscanf(fid, '%f', [natoms 1]); fgetl(fid); pos = ftell(fid); fgetl(fid);
    if ~feof(fid); fseek(fid,pos,'bof');
        gridsize = fscanf(fid, '%d %d %d', [3 1])';
        mag_z = fscanf(fid, '%f', [prod(gridsize,2) 1])';
        mag_z = reshape(mag_z,gridsize); mag_z = mag_z/vol;
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) ~= 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_z = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_z = reshape(mag_z,gridsize);
    if size(str2num(line),2) == 3; mag_z = mag_z/vol; end
end
fclose(fid);

function [H, Hx, Hy, Hz, geo] = dipcar(filename)
fid = fopen(filename); 
geo = poscar(fid); fgetl(fid); 
gridsize = fscanf(fid, '%d %d %d', [3 1])';
H = fscanf(fid, '%f', [prod(gridsize,2) 1])';
H = reshape(H,gridsize);
if fgetl(fid) == ''; fgetl(fid); end
gridsize = fscanf(fid, '%d %d %d', [3 1])';
Hx = fscanf(fid, '%f', [prod(gridsize,2) 1])';
Hx = reshape(Hx,gridsize);
if fgetl(fid) == ''; fgetl(fid); end
gridsize = fscanf(fid, '%d %d %d', [3 1])';
Hy = fscanf(fid, '%f', [prod(gridsize,2) 1])';
Hy = reshape(Hy,gridsize);
if fgetl(fid) == ''; fgetl(fid); end
gridsize = fscanf(fid, '%d %d %d', [3 1])';
Hz = fscanf(fid, '%f', [prod(gridsize,2) 1])';
Hz = reshape(Hz,gridsize);
fclose(fid);

function drawbox(latt)
line([0 latt(1,1)],[0 latt(1,2)],[0 latt(1,3)],'LineWidth',2,'Color','k')
line([latt(2,1) latt(1,1)+latt(2,1)],[latt(2,2) latt(1,2)+latt(2,2)], ...
    [latt(2,3) latt(1,3)+latt(2,3)],'LineWidth',2,'Color','k')
line([latt(3,1) latt(1,1)+latt(3,1)],[latt(3,2) latt(1,2)+latt(3,2)], ...
    [latt(3,3) latt(1,3)+latt(3,3)],'LineWidth',2,'Color','k')
line([latt(2,1)+latt(3,1) latt(1,1)+latt(2,1)+latt(3,1)],...
    [latt(2,2)+latt(3,2) latt(1,2)+latt(2,2)+latt(3,2)], ...
    [latt(2,3)+latt(3,3) latt(1,3)+latt(2,3)+latt(3,3)],'LineWidth',2,'Color','k')
line([0 latt(2,1)],[0 latt(2,2)],[0 latt(2,3)],'LineWidth',2,'Color','k')
line([latt(1,1) latt(1,1)+latt(2,1)],[latt(1,2) latt(1,2)+latt(2,2)],...
    [latt(1,3) latt(1,3)+latt(2,3)],'LineWidth',2,'Color','k')
line([latt(3,1) latt(2,1)+latt(3,1)],[latt(3,2) latt(2,2)+latt(3,2)],...
    [latt(3,3) latt(2,3)+latt(3,3)],'LineWidth',2,'Color','k')
line([latt(1,1)+latt(3,1) latt(1,1)+latt(2,1)+latt(3,1)],...
    [latt(1,2)+latt(3,2) latt(1,2)+latt(2,2)+latt(3,2)],...
    [latt(1,3)+latt(3,3) latt(1,3)+latt(2,3)+latt(3,3)],'LineWidth',2,'Color','k')
line([0 latt(3,1)],[0 latt(3,2)],[0 latt(3,3)],'LineWidth',2,'Color','k')
line([latt(1,1) latt(1,1)+latt(3,1)],[latt(1,2) latt(1,2)+latt(3,2)],...
    [latt(1,3) latt(1,3)+latt(3,3)],'LineWidth',2,'Color','k')
line([latt(2,1) latt(2,1)+latt(3,1)],[latt(2,2) latt(2,2)+latt(3,2)],...
    [latt(2,3) latt(2,3)+latt(3,3)],'LineWidth',2,'Color','k')
line([latt(1,1)+latt(2,1) latt(1,1)+latt(2,1)+latt(3,1)],...
    [latt(1,2)+latt(2,2) latt(1,2)+latt(2,2)+latt(3,2)],...
    [latt(1,3)+latt(2,3) latt(1,3)+latt(2,3)+latt(3,3)],'LineWidth',2,'Color','k')

function Iso_Callback(~, ~, ~)
function Iso_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_axis1D_Callback(~, ~, ~)
function opt_axis1D_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function point2_Callback(~, ~, ~)
function point2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function point1_Callback(~, ~, ~)
function point1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function run_1D_Callback(~, ~, handles)
pos1 = str2num(get(handles.point1,'String')); clc; tic
fprintf('   ===================================== \n')
fprintf('   Read input files. \n')
pos2 = str2num(get(handles.point2,'String'));
N = str2double(get(handles.Npoint,'String')); 
filename = get(handles.filename,'String');
if (get(handles.opt_file,'Value') == 1) % locpot
    [rho, geo] = locpot(filename); rho = -rho; y_label = 'Energy (eV)';
elseif (get(handles.opt_file,'Value') == 2) % chg
    [rho, ~, geo] = chgcar(filename); y_label = 'Charge';
elseif (get(handles.opt_file,'Value') == 3) % chg gradien
    [rho, ~, geo] = chgcar(filename);
    [dx,dy,dz] = gradient(rho,0.01,0.01,0.01);
    rho = abs(sqrt(dx.^2+dy.^2+dz.^2))./rho.^(4/3); y_label = 'Gradient Density';
elseif (get(handles.opt_file,'Value') == 4) % spin coll
    [~, rho, geo] = chgcar(filename); y_label = 'Spin';
elseif (get(handles.opt_file,'Value') == 5) % spin non-coll
    [~, mag_x, mag_y, mag_z, geo] = magcar(filename); y_label = 'Spin';
    if get(handles.H,'Value') == 1; rho = sqrt(mag_x.^2+mag_y.^2+mag_z.^2);
    elseif get(handles.Hx,'Value') == 1; rho = mag_x;
    elseif get(handles.Hy,'Value') == 1; rho = mag_y;
    elseif get(handles.Hz,'Value') == 1; rho = mag_z;
    end
elseif (get(handles.opt_file,'Value') == 6) % ELF
    [rho, geo] = locpot(filename); y_label = 'ELF';
elseif (get(handles.opt_file,'Value') == 7) % Dipole
    [rho, mag_x, mag_y, mag_z, geo] = dipcar(filename); y_label = 'Dipole Field (G)';
    if get(handles.Hx,'Value') == 1; rho = mag_x;
    elseif get(handles.Hy,'Value') == 1; rho = mag_y;
    elseif get(handles.Hz,'Value') == 1; rho = mag_z;
    end
end
rho(end+1,:,:) = rho(1,:,:); rho(:,end+1,:) = rho(:,1,:); rho(:,:,end+1) = rho(:,:,1); 
x = linspace(0,1,size(rho,1)); y = linspace(0,1,size(rho,2));
z = linspace(0,1,size(rho,3)); pos = [pos1; pos2]; pot = []; 
fprintf('   Drawing 1D Figure. \n')
figure; hold on; box on
if get(handles.opt_axis1D,'Value') == 2
    latt = inv(geo.lattice);
    pos1 = pos(1,1)*latt(1,:)+pos(1,2)*latt(2,:)+pos(1,3)*latt(3,:);
    pos2 = pos(2,1)*latt(1,:)+pos(2,2)*latt(2,:)+pos(2,3)*latt(3,:);
end
a = abs(pos(1,1)-x); [~,a] = min(a); b = abs(pos(1,2)-y); [~,b] = min(b);
c = abs(pos(1,3)-z); [~,c] = min(c); p1 = [a b c];
a = abs(pos(2,1)-x); [~,a] = min(a); b = abs(pos(2,2)-y); [~,b] = min(b);
c = abs(pos(2,3)-z); [~,c] = min(c); p2 = [a b c];
if length(N) > 1
    cx = ':'; cy = ':'; cz = ':';
    dx = 1:size(rho,1); dy = 1:size(rho,2);
else cz = [num2str(p1(3)) ':' num2str(p2(3))];
    cx = [num2str(p1(1)) ':' num2str(p2(1))]; dx = p1(1):p2(1);
    cy = [num2str(p1(2)) ':' num2str(p2(2))]; dy = p1(2):p2(2);
end
if p1(3) == p2(3)
    if p1(2) == p2(2)
        if p1(1) == p2(1)
            error('points are in the same position')
        else  V = rho(cx,p1(2),p1(3)); r = x;
        end
    elseif p1(1) == p2(1) % in b-axis
        V = rho(p1(1),cy,p1(3)); r = y;
    else l = 1;
        for i = dx
            j = round(p1(2)+(p2(2)-p1(2)).*(i-p1(1))./(p2(1)-p1(1)));
            if j > 0 && j <= size(rho,2)
                V(l) = rho(i,j,p1(3)); r(l) = sqrt(x(i)^2+y(j)^2); l = l+1;
            end
        end
    end
elseif p1(1) == p2(1)
    if p1(2) == p2(2); V = rho(p1(1),p1(2),cz); r = z;
    else l = 1;
        for j = dy; k = round(p1(3)+(p2(3)-p1(3)).*(j-p1(2))./(p2(2)-p1(2)));
            if k > 0 && k <= size(rho,3)
                V(l) = rho(p1(1),j,k); r(l) = sqrt(y(j)^2+z(k)^2); l = l+1;
            end
        end
    end
elseif p1(2) == p2(2); l = 1;
    for i = dx; k = round(p1(3)+(p2(3)-p1(3)).*(i-p1(1))./(p2(1)-p1(1)));
        if k > 0 && k <= size(rho,3)
            V(l) = rho(i,p1(1),k); r(l) = sqrt(x(i)^2+z(k)^2); l = l+1;
        end
    end
else l = 1;
    for i = dx; j = round(p1(2)+(p2(2)-p1(2)).*(i-p1(1))./(p2(1)-p1(1)));
        k = round(p1(3)+(p2(3)-p1(3)).*(i-p1(1))./(p2(1)-p1(1)));
        if j > 0 && j <= size(rho,2) && k > 0 && k <= size(rho,3)
            V(l) = rho(i,j,k); r(l) = sqrt(x(i)^2+y(j)^2+z(k)^2); l = l+1;
        end
    end
end
V = -V(:)'; ri = linspace(r(1),r(end),N); Vi = interp1(r,V,ri); pot = [pot, Vi];
if get(handles.opt_axis1D,'Value') == 1
    r = linspace(0,1,N);
elseif get(handles.opt_axis1D,'Value') == 2
    pos1 = pos(1,:); pos2 = pos(2,:);
    r = sqrt(dot(pos1-pos2,pos1-pos2)); r = linspace(0,r,N);
end
plot(r,pot,'b','LineWidth',2)
xlabel('Positions','FontSize',12); ylabel(y_label,'FontSize',12)
fprintf('   ===================================== \n   '); toc

function opt_axis_2D_Callback(~, ~, ~)
function opt_axis_2D_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function run_3D_Callback(~, ~, handles)
iso = str2double(get(handles.Iso,'String')); clc
fprintf('   ===================================== \n'); tic
fprintf('   Read input files. \n')
No = str2num(get(handles.Atom_no,'String'));
filename = get(handles.filename,'String');
if (get(handles.opt_file,'Value') == 1) % locpot
    [rho, geo] = locpot(filename); rho = -rho;
elseif (get(handles.opt_file,'Value') == 2) % chg
    [rho, ~, geo] = chgcar(filename);
elseif (get(handles.opt_file,'Value') == 3) % chg gradien
    [rho, ~, geo] = chgcar(filename);
    [dx,dy,dz] = gradient(rho,0.01,0.01,0.01);
    rho = abs(sqrt(dx.^2+dy.^2+dz.^2))./rho.^(4/3);
elseif (get(handles.opt_file,'Value') == 4) % spin coll
    [~, rho, geo] = chgcar(filename);
elseif (get(handles.opt_file,'Value') == 5) % spin non-coll
    [~, mag_x, mag_y, mag_z, geo] = magcar(filename);
    if get(handles.H,'Value') == 1; rho = sqrt(mag_x.^2+mag_y.^2+mag_z.^2);
    elseif get(handles.Hx,'Value') == 1; rho = mag_x;
    elseif get(handles.Hy,'Value') == 1; rho = mag_y;
    elseif get(handles.Hz,'Value') == 1; rho = mag_z;
    end
elseif (get(handles.opt_file,'Value') == 6) % ELF
    [rho, geo] = locpot(filename);
elseif (get(handles.opt_file,'Value') == 7) % Dipole
    [rho, mag_x, mag_y, mag_z, geo] = dipcar(filename);
    if get(handles.Hx,'Value') == 1; rho = mag_x;
    elseif get(handles.Hy,'Value') == 1; rho = mag_y;
    elseif get(handles.Hz,'Value') == 1; rho = mag_z;
    end
end
rho(end+1,:,:) = rho(1,:,:); rho(:,end+1,:) = rho(:,1,:); rho(:,:,end+1) = rho(:,:,1); 
if (get(handles.opt_iso,'Value') == 2)
    iso = -iso;
elseif (get(handles.opt_iso,'Value') == 4)
    iso = min(rho(:))+iso;
end
fprintf('   Drawing 3D Figure. \n')
col = get(handles.ion_col,'String'); rd = str2num(get(handles.ion_rad,'String'));
if get(handles.opt_bond,'Value') > 1; d0 = str2num(get(handles.bond_l,'String')); end
xb = str2num(get(handles.xb,'String')); figure; latt = geo.lattice;
yb = str2num(get(handles.yb,'String')); zb = str2num(get(handles.zb,'String'));
if xb(2)<=xb(1); xa=xb; xb=[xa(2) xa(1)]; end; if xb(1)<0; xb(1)=0; end; if xb(2)>1; xb(2)=1; end
if yb(2)<=yb(1); ya=yb; yb=[ya(2) ya(1)]; end; if yb(1)<0; yb(1)=0; end; if yb(2)>1; yb(2)=1; end
if zb(2)<=zb(1); za=zb; xb=[za(2) za(1)]; end; if zb(1)<0; zb(1)=0; end; if zb(2)>1; zb(2)=1; end
if length(No) == 1; col_oct = col(1); else col_oct = col(2);end
for no = 1:length(No)
    cu = geo.coords(sum(geo.atomcount(1:No(no)-1))+1:sum(geo.atomcount(1:No(no))),:); Cu = [];
    for n = 1:size(cu,1)
        for i = -1:1
            for j = -1:1
                for k = -1:1
                    if get(handles.opt_bond,'Value') > 1  && no <= 2
                        latt1 = inv(latt); d1 = (d0(1)-1)*(latt1(1,:)+latt1(2,:)+latt1(3,:));
                    else d1 = 0.02*[1 1 1];
                    end
                    if cu(n,1)+i >= xb(1)-d1(1) && cu(n,1)+i <= xb(2)+d1(1) ...
                            && cu(n,2)+j >= yb(1)-d1(2) && cu(n,2)+j <= yb(2)+d1(2) ...
                            && cu(n,3)+k >= zb(1)-d1(3) && cu(n,3)+k <= zb(2)+d1(3)
                        Cu = [Cu; (cu(n,1)+i)*latt(1,:)+(cu(n,2)+j)*latt(2,:)+(cu(n,3)+k)*latt(3,:)];
                    end
                end
            end
        end
    end
    draw_atom(Cu,rd(no),col(no),1)
    if get(handles.opt_bond,'Value') == 2 && no <= 2
        if no == 1; ion1 = Cu;
        elseif no == 2; drbond(ion1,Cu,d0,col(1:2),2); end
    elseif get(handles.opt_bond,'Value') == 3 && no == 1
        drplane(Cu,d0,col_oct,0.4*[1 1 1],0.3)
    elseif get(handles.opt_bond,'Value') == 4 && no <= 2
        if no == 1; ion1 = Cu; drplane(Cu,d0,col_oct,0.4*[1 1 1],0.3)
        elseif no == 2; drbond(ion1,Cu,d0,col(1:2),2); end
    end
end
xa = linspace(0,1,size(rho,1)); ya = linspace(0,1,size(rho,2));
za = linspace(0,1,size(rho,3)); 
[x,y,z] = ndgrid(xa,ya,za); x = x(:); y = y(:); z = z(:);
xi = x*latt(1,1) + y*latt(2,1) + z*latt(3,1);
yi = x*latt(1,2) + y*latt(2,2) + z*latt(3,2);
zi = x*latt(1,3) + y*latt(2,3) + z*latt(3,3);
x = reshape(xi,size(rho)); y = reshape(yi,size(rho)); z = reshape(zi,size(rho)); 
for i = 1:2
    a = abs(xa-xb(i)); [~,xb(i)] = min(a);
    a = abs(ya-yb(i)); [~,yb(i)] = min(a);
    a = abs(za-zb(i)); [~,zb(i)] = min(a);
end
x = x(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
y = y(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
z = z(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
rho = rho(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
p = patch(isosurface(x,y,z,rho,iso));
set(p,'FaceColor','y','EdgeColor','none','FaceAlpha',0.5)
if (get(handles.opt_iso,'Value') == 3)
    p = patch(isosurface(x,y,z,rho,-iso));
    set(p,'FaceColor','c','EdgeColor','none','FaceAlpha',0.5)
end
if get(handles.box,'Value') == 1; drawbox(geo.lattice); end
daspect([1 1 1]); camlight right local; rotate3d on;
if get(handles.fix,'Value') == 1
    xlabel('a','FontSize',12); ylabel('b','FontSize',12); 
    zlabel('c','FontSize',12); axis tight
else axis off
end
fprintf('   ===================================== \n   '); toc

function draw_atom(target,r,col,alp); hold on
if nargin == 3; alp = 1;
elseif nargin == 2
    if ischar(r); col = r; r = 0.2; alp = 1; else col = 'k'; alp = 1; end
elseif nargin == 1; r = 0.2; col = 'k'; alp = 1;
end
[rx,ry,rz] = sphere(10); rx = rx*r; ry = ry*r; rz = rz*r;
for i = 1:size(target,1)
    h = surf((rx+target(i,1)),(ry+target(i,2)),(rz+target(i,3)));
    set(h,'edgecolor','none','facecolor',col,'facelighting','gouraud'); alpha(h,alp)
end

function opt_iso_Callback(~, ~, ~)
function opt_iso_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Npoint_Callback(~, ~, ~)
function Npoint_CreateFcn(hObject, ~, ~)
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

function opt_bond_Callback(~, ~, handles)
if get(handles.opt_bond,'Value') > 1
    set(handles.bond_l,'BackgroundColor',[1 1 1])
    set(handles.bond_l,'ForegroundColor',[0 0 0])
else set(handles.bond_l,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.bond_l,'ForegroundColor',[0.8 0.8 0.8])
end

function opt_bond_CreateFcn(hObject, ~, ~)
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

function box_Callback(~, ~, ~)
function fix_Callback(~, ~, ~)

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

function cek_ax_Callback(~, ~, handles)
if get(handles.cek_ax,'Value') == 1
    set(handles.c_axis,'BackgroundColor',[1 1 1])
    set(handles.c_axis,'ForegroundColor',[0 0 0])
else set(handles.c_axis,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.c_axis,'ForegroundColor',[0.8 0.8 0.8])
end

function c_axis_Callback(~, ~, ~)
function c_axis_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function run_dis_Callback(~, ~, handles)
filename = get(handles.filename,'String'); clc; tic
fprintf('   ===================================== \n')
fprintf('   Read input files. \n')
figure, hold on, box on; xlabel('Distribution','FontSize',18); 
if (get(handles.opt_file,'Value') == 1) % locpot
    [rho, ~] = locpot(filename); rho = -rho;
    xlabel('Energy (eV)','FontSize',18); 
elseif (get(handles.opt_file,'Value') == 2) % chg
    [rho, ~, ~] = chgcar(filename);
elseif (get(handles.opt_file,'Value') == 3) % chg gradien
    [rho, ~, ~] = chgcar(filename);
    [dx,dy,dz] = gradient(rho,0.01,0.01,0.01);
    rho = abs(sqrt(dx.^2+dy.^2+dz.^2))./rho.^(4/3);
elseif (get(handles.opt_file,'Value') == 4) % spin coll
    [~, rho, ~] = chgcar(filename);
elseif (get(handles.opt_file,'Value') == 5) % spin non-coll
    [~, mag_x, mag_y, mag_z, ~] = magcar(filename);
    if get(handles.H,'Value') == 1; rho = sqrt(mag_x.^2+mag_y.^2+mag_z.^2);
    elseif get(handles.Hx,'Value') == 1; rho = mag_x;
    elseif get(handles.Hy,'Value') == 1; rho = mag_y;
    elseif get(handles.Hz,'Value') == 1; rho = mag_z;
    end
elseif (get(handles.opt_file,'Value') == 6) % ELF
    [rho, ~] = locpot(filename);
elseif (get(handles.opt_file,'Value') == 7) % Dipole
    [rho, mag_x, mag_y, mag_z, ~] = dipcar(filename);
    xlabel('Fields (G)','FontSize',18); 
    if get(handles.Hx,'Value') == 1; rho = mag_x;
    elseif get(handles.Hy,'Value') == 1; rho = mag_y;
    elseif get(handles.Hz,'Value') == 1; rho = mag_z;
    end
end
rho = rho(:); fprintf('   Plot Distribution. \n')
if (get(handles.opt_file,'Value') == 7) && (get(handles.freq,'Value') == 1)
    rho = rho*2*pi*0.013554; xlabel('Freq (MHz)','FontSize',18); 
end
dis = str2num(get(handles.dis,'String'));
if length(dis) == 1; dis = [0 dis]; end
if dis(2) < dis(1); dis = [dis(2) dis(1)]; end
del = str2double(get(handles.del,'String'));
if get(handles.opt_del,'Value') == 2
    del = round((dis(2)-dis(2))/del);
end
Chg = linspace(dis(1),dis(2),del);
Dist = zeros(1,length(Chg));
for i = 1:length(Chg)-1;
    n = (rho > Chg(i) & rho <= Chg(i+1));
    Dist(i+1) = sum(n(:));
end
aa = Chg(Dist == max(Dist)); aa = max(aa(:));
plot(aa*[1 1],[0 max(Dist)],'g','LineWidth',1)
plot(Chg,Dist,'k','linewidth',2.0)
legend(['Distribution max '  num2str(aa)])
ylabel('Probability','FontSize',18);
fprintf('   ===================================== \n   '); toc

function del_Callback(~, ~, ~)
function del_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dis_Callback(~, ~, ~)
function dis_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_del_Callback(~, ~, ~)
function opt_del_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function hh = quiverc(varargin)
alpha = 0.33; beta = 0.23; autoscale = 1; plotarrows = 1; 
filled = 0; ms = ''; col = ''; lw=1; nin = nargin;
while ischar(varargin{nin}); vv = varargin{nin};
    if ~isempty(vv) && strcmpi(vv(1),'f')
        filled = 1; nin = nin-1;
    else [~,c,m,~] = colstyle(vv); nin = nin-1;
        if ~isempty(c), col = c; end
        if ~isempty(m), ms = m; plotarrows = 0; end
        if isequal(m,'.'), ms = ''; end
    end 
end
error(nargchk(2,5,nin));
if nin < 4; [msg,x,y,u,v] = xyzchk(varargin{1:2});
else [msg,x,y,u,v] = xyzchk(varargin{1:4});
end
if ~isempty(msg), error(msg); end
if nin == 3 || nin == 5; autoscale = varargin{nin}; end
if numel(u)==1, u = u(ones(size(x))); end
if numel(v)==1, v = v(ones(size(u))); end
if autoscale
    if min(size(x))==1, n=sqrt(numel(x)); m=n; else [m,n]=size(x); end
    delx = diff([min(x(:)) max(x(:))])/n;
    dely = diff([min(y(:)) max(y(:))])/m;
    len = sqrt((u.^2 + v.^2)/(delx.^2 + dely.^2));
    autoscale = autoscale*0.9 / max(len(:));
    u = u*autoscale; v = v*autoscale;
end
vr = sqrt(u.^2+v.^2); vrn = round(vr/max(vr(:))*64);
CC = colormap; ax = newplot;
next = lower(get(ax,'NextPlot'));
hold_state = ishold; hold on
x = x(:).'; y = y(:).'; u = u(:).'; v = v(:).';
vrn = vrn(:).'; uu = [x; x+u; NaN(size(u))];
vv = [y; y+v; NaN(size(u))]; vrn1 = [vrn; NaN(size(u)); NaN(size(u))];
uui = uu(:);  vvi = vv(:);  vrn1 = vrn1(:); imax = size(uui);
for i = 1:3:imax-1
    ii = int8(round(vrn1(i)));
    if ii == 0; ii = 1; end        
    c1 = CC(ii,1); c2 = CC(ii,2); c3 = CC(ii,3);
    plot(uui(i:i+1),vvi(i:i+1),'linewidth',lw,'color',[c1 c2 c3]);
end
if plotarrows,
    hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
        x+u-alpha*(u-beta*(v+eps));NaN(size(u))];
    hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
        y+v-alpha*(v+beta*(u+eps));NaN(size(v))];
    vrn2 = [vrn;vrn;vrn;vrn];
    uui = hu(:); vvi = hv(:); vrn2 = vrn2(:); imax = size(uui);
    for i = 1:imax-1; ii = int8(round(vrn2(i)));
        if ii == 0; ii=1; end 
        c1 = CC(ii,1); c2 = CC(ii,2); c3 = CC(ii,3);
        plot(uui(i:i+1),vvi(i:i+1),'linewidth',lw,'color',[c1 c2 c3]);
    end
else h2 = [];
end
if ~isempty(ms); hu = x; hv = y; hold on
    h3 = plot(hu(:),hv(:),[col ms]);
    if filled, set(h3,'markerfacecolor',get(h1,'color')); end
else h3 = [];
end
if ~hold_state, hold off, view(2); set(ax,'NextPlot',next); end
if nargout>0, hh = [h1;h2;h3]; end

function lbl_Callback(~, ~, handles)
a = get(handles.opt_plot,'Value');
if a == 1 || a == 4; set(handles.lbl,'Value',0); end


function H_Callback(~, ~, handles)
set(handles.Hx,'Value',0); set(handles.Hx,'ForegroundColor',[0 0 0])
set(handles.Hy,'Value',0); set(handles.Hy,'ForegroundColor',[0 0 0])
set(handles.Hz,'Value',0); set(handles.Hz,'ForegroundColor',[0 0 0])
set(handles.H,'Value',1); set(handles.H,'ForegroundColor',[0 0 1])
set(handles.Hxy,'Value',0); set(handles.Hxy,'ForegroundColor',[0 0 0])

function Hx_Callback(~, ~, handles)
set(handles.Hx,'Value',1); set(handles.Hx,'ForegroundColor',[0 0 1])
set(handles.Hy,'Value',0); set(handles.Hy,'ForegroundColor',[0 0 0])
set(handles.Hz,'Value',0); set(handles.Hz,'ForegroundColor',[0 0 0])
set(handles.H,'Value',0); set(handles.H,'ForegroundColor',[0 0 0])
set(handles.Hxy,'Value',0); set(handles.Hxy,'ForegroundColor',[0 0 0])

function Hy_Callback(~, ~, handles)
set(handles.Hx,'Value',0); set(handles.Hx,'ForegroundColor',[0 0 0])
set(handles.Hy,'Value',1); set(handles.Hy,'ForegroundColor',[0 0 1])
set(handles.Hz,'Value',0); set(handles.Hz,'ForegroundColor',[0 0 0])
set(handles.H,'Value',0); set(handles.H,'ForegroundColor',[0 0 0])
set(handles.Hxy,'Value',0); set(handles.Hxy,'ForegroundColor',[0 0 0])

function Hz_Callback(~, ~, handles)
set(handles.Hx,'Value',0); set(handles.Hx,'ForegroundColor',[0 0 0])
set(handles.Hy,'Value',0); set(handles.Hy,'ForegroundColor',[0 0 0])
set(handles.Hz,'Value',1); set(handles.Hz,'ForegroundColor',[0 0 1])
set(handles.H,'Value',0); set(handles.H,'ForegroundColor',[0 0 0])
set(handles.Hxy,'Value',0); set(handles.Hxy,'ForegroundColor',[0 0 0])

function Hxy_Callback(hObject, eventdata, handles)
set(handles.Hx,'Value',0); set(handles.Hx,'ForegroundColor',[0 0 0])
set(handles.Hy,'Value',0); set(handles.Hy,'ForegroundColor',[0 0 0])
set(handles.Hz,'Value',0); set(handles.Hz,'ForegroundColor',[0 0 0])
set(handles.H,'Value',0); set(handles.H,'ForegroundColor',[0 0 0])
set(handles.Hxy,'Value',1); set(handles.Hxy,'ForegroundColor',[0 0 1])

function freq_Callback(~, ~, ~)
function Habs_Callback(hObject, eventdata, handles)
function super_Callback(hObject, eventdata, handles)

