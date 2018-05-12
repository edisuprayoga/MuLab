function varargout = GUI_Makefile(varargin)

% GUI_MAKEFILE MATLAB code by Edi Suprayoga
% feel free to contact me at suprayoga.edi@gmail.com
% last edit 3 Mar 2016

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Makefile_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Makefile_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout; [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else gui_mainfcn(gui_State, varargin{:});
end

function GUI_Makefile_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject; guidata(hObject, handles);
set(handles.figure1,'Color',[0.8 0.8 0.8])

function varargout = GUI_Makefile_OutputFcn(~, ~, handles) 
varargout{1} = handles.output; set(handles.autorized,'Value',1);

function load_sup_Callback(~, ~, handles)
global filename_sup pathname_sup
[filename_sup, pathname_sup] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname_sup)
    set(handles.disp_pos_sup,'String',[pathname_sup filename_sup])
end

function disp_pos_sup_Callback(~, ~, ~)
function disp_pos_sup_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_sup_Callback(~, ~, handles)
if get(handles.opt_sup,'Value') < 4
     set(handles.muon_sup,'BackgroundColor',[0.8 0.8 0.8])
     set(handles.muon_sup,'ForegroundColor',[0.8 0.8 0.8])
     set(handles.range,'BackgroundColor',[0.8 0.8 0.8])
     set(handles.range,'ForegroundColor',[0.8 0.8 0.8])
else set(handles.muon_sup,'BackgroundColor',[1 1 1])
     set(handles.muon_sup,'ForegroundColor',[0 0 0])
    if get(handles.opt_sup,'Value') < 8
         set(handles.range,'BackgroundColor',[0.8 0.8 0.8])
         set(handles.range,'ForegroundColor',[0.8 0.8 0.8])
    else set(handles.range,'ForegroundColor',[0 0 0])
         set(handles.range,'BackgroundColor',[1 1 1])
    end
end

function opt_sup_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nx_sup_Callback(~, ~, ~)
function Nx_sup_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function range_Callback(~, ~, ~)
function range_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function muon_sup_Callback(~, ~, ~)
function muon_sup_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function run_sup_Callback(~, ~, handles)
filename = get(handles.disp_pos_sup,'String'); clc; tic
fprintf('   ======================================= \n');
fprintf('>> Reading input files \n')
Nx = str2num(get(handles.Nx_sup,'String'));
muon = str2num(get(handles.muon_sup,'String'));
Rmax = str2num(get(handles.range,'String'));
if get(handles.center2,'Value') == 1
    if get(handles.opt_sup,'Value') == 1 % default
        make_poscar(filename, Nx)
    elseif get(handles.opt_sup,'Value') == 2 % all static
        make_poscar(filename, Nx, 'all_F')
    elseif get(handles.opt_sup,'Value') == 3 % all relax
        make_poscar(filename, Nx, 'all_T')
    elseif get(handles.opt_sup,'Value') == 4 % default with muon
        make_poscar(filename, Nx, [], muon)
    elseif get(handles.opt_sup,'Value') == 5 % all static with muon
        make_poscar(filename, Nx, 'static_mu', muon)
    elseif get(handles.opt_sup,'Value') == 6 % all relax with muon
        make_poscar(filename, Nx , 'all_T', muon)
    elseif get(handles.opt_sup,'Value') == 7 % only relax muon
        make_poscar(filename, Nx, 'all_F', muon)
    elseif get(handles.opt_sup,'Value') == 8 % relax muon range
        make_poscar(filename, Nx, 'sel_D', muon, Rmax)
    end
else
    if get(handles.opt_sup,'Value') == 1 % default
        make_poscar2(filename, Nx)
    elseif get(handles.opt_sup,'Value') == 2 % all static
        make_poscar2(filename, Nx, 'all_F')
    elseif get(handles.opt_sup,'Value') == 3 % all relax
        make_poscar2(filename, Nx, 'all_T')
    elseif get(handles.opt_sup,'Value') == 4 % default with muon
        make_poscar2(filename, Nx, [], muon)
    elseif get(handles.opt_sup,'Value') == 5 % all static with muon
        make_poscar2(filename, Nx, 'static_mu', muon)
    elseif get(handles.opt_sup,'Value') == 6 % all relax with muon
        make_poscar2(filename, Nx, 'all_T', muon)
    elseif get(handles.opt_sup,'Value') == 7 % only relax muon
        make_poscar2(filename, Nx, 'all_F', muon)
    elseif get(handles.opt_sup,'Value') == 8 % relax muon range
        make_poscar2(filename, Nx, 'sel_D', muon, Rmax)
    end
end
fprintf('>> Makefile Completed!'); 
fprintf('\n   ======================================= \n');toc

function load_xdat_Callback(~, ~, handles)
global filename_xdat pathname_xdat
[filename_xdat, pathname_xdat] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname_xdat)
    if get(handles.opt_xdat,'Value') < 3
        set(handles.disp_xdat,'String',pathname_xdat)
    else set(handles.disp_xdat,'String',[pathname_xdat filename_xdat])
    end
end

function disp_xdat_Callback(~, ~, ~)
function disp_xdat_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function make_chgcar(geo,chg,filename)
fprintf('>> Writing CHGCAR file \n')
fid = fopen(filename);
if fid ~= -1; i = 1; 
    fid = fopen([filename '(1)']);
    while fid ~= -1; i = i+1; 
        fclose(fid); fid = fopen([filename '(' num2str(i) ')']);
    end
    fid = fopen([filename '(' num2str(i) ')'],'w');
    fprintf(['     ' filename '(' num2str(i) ')\n'])
else fid = fopen(filename,'w');
    fprintf(['     ' filename '\n'])
end
fprintf(fid,'%s\n',geo.comment);
fprintf(fid,' %1.1f\n',1);
fprintf(fid,'   %4.6f  %4.6f  %4.6f\n',geo.lattice);
if ~isempty(geo.symbols)
    cellfun(@(x) fprintf(fid, '%s ', x), geo.symbols);
    fprintf(fid, '\n');
end
fprintf(fid, ' %d ', geo.atomcount);
fprintf(fid, '\nDirect \n');
fprintf(fid, ' %19.16f %19.16f %19.16f \n', geo.coords');
fprintf(fid,'\n   %d   %d   %d \n',size(chg));
fprintf(fid,' %e %e %e %e %e\n',chg(:));
fclose(fid);

function run_xdat_Callback(~, ~, handles)
pathname = get(handles.disp_xdat,'String'); clc; tic
fprintf('   ======================================= \n');
N = str2double(get(handles.N_xdat,'String'));
if get(handles.opt_xdat,'Value') < 3 % make xdatcar
    fprintf('>> Reading OUTCAR file \n')
    fid = fopen([pathname '/OUTCAR']); vasp = fscanf(fid, '%s', 1); 
    vsp = str2double(vasp(8)); vasp = str2double(vasp(6)); fclose(fid);
    if get(handles.cek_num,'Value') == 0
        fprintf('>> Reading OSZICAR file \n')
        fid = fopen([pathname '/OSZICAR']); energy = [];
        while ~feof(fid); line = fgetl(fid);
            [match, tok] = regexp(line, 'E0= ([\.-+E0-9]*)', 'match', 'tokens');
            if numel(match) > 0; energy = [energy str2double(tok{1}{1})]; end 
        end
        N = length(energy); fclose(fid);
    end
    geo1 = get_poscar([pathname '/POSCAR']); fprintf('>> Reading POSCAR file \n')
    geo2 = get_poscar([pathname '/CONTCAR']); fprintf('>> Reading CONTCAR file \n')
    fprintf('\n>> Writing POSCAR files \n');
    fprintf('     Number of files: %1.0f \n',N)
    fid = fopen([pathname '/XDAT/POSCAR_1']);
    if fid ~= -1; i = 1;
        fid = fopen([pathname '/XDAT(1)/POSCAR_1']);
        while fid ~= -1; i = i+1; 
            fclose(fid); fid = fopen([pathname '/XDAT(' num2str(i) ')/POSCAR_1']);
        end
        dir = [pathname '/XDAT(' num2str(i) ')/']; mkdir(dir); 
        fprintf(['     ' pathname '/XDAT(' num2str(i) ')\n'])
    else fprintf(['     ' pathname '/XDAT \n'])
        dir = [pathname '/XDAT/']; mkdir(dir); 
    end
    file = textread([pathname '/XDATCAR'],'%s','delimiter','\n','whitespace','');
    if get(handles.opt_xdat,'Value') == 2 || vsp == 4 || vasp == 4
        shift = geo2.lattice - geo1.lattice + repmat([0.5 0.5 0.5], 3, 1);
        shift = mod(shift,1) - repmat([0.5 0.5 0.5], 3, 1); 
        if get(handles.opt_xdat,'Value') == 2 || vsp == 4; n = 8; elseif vasp == 4; n = 6; end
        for i = 1:N
            fid = fopen([dir 'POSCAR_' num2str(i)],'w');
            fprintf(fid,[geo1.comment '\n']);
            fprintf(fid,'1.0\n');
            geo1.lattice = geo1.lattice + i/N*shift;
            fprintf(fid, '%19.16f %19.16f %19.16f\n', geo1.lattice');
            if ~isempty(geo1.symbols)
                cellfun(@(x) fprintf(fid, '%s ', x), geo1.symbols);
                fprintf(fid, '\n');
            end
            fprintf(fid, '%d ', geo1.atomcount);
            fprintf(fid, '\nDirect\n');
            for j = 1:sum(geo1.atomcount)
                geo1.coords(j,:) = cell2mat(textscan(cell2mat(file(n+j)),'%f64 %f64 %f64'));
            end
            fprintf(fid, '%19.16f %19.16f %19.16f\n', geo1.coords');
            fclose(fid); n = n+j+1;
        end
    elseif vasp == 5; n = 9; % vasp.5.3
        for i = 1:N
            fid = fopen([dir 'POSCAR_' num2str(i)],'w');
            fprintf(fid,[geo1.comment '\n']); fprintf(fid,'1.0\n');
            for j = 1:3
                lattice = cell2mat(textscan(cell2mat(file(n+j)),'%f64 %f64 %f64'));
                fprintf(fid, '%19.16f %19.16f %19.16f\n', lattice');
            end
            n = n+j+3;
            if ~isempty(geo1.symbols)
                cellfun(@(x) fprintf(fid, '%s ', x), geo1.symbols);
                fprintf(fid, '\n');
            end
            fprintf(fid, '%d ', geo1.atomcount);
            fprintf(fid, '\nDirect\n');
            for j = 1:sum(geo1.atomcount)
                geo1.coords(j,:) = cell2mat(textscan(cell2mat(file(n+j)),'%f64 %f64 %f64'));
            end
            fprintf(fid, '%19.16f %19.16f %19.16f\n', geo1.coords');
            fclose(fid); n = n+j+2;
        end
    end
elseif get(handles.opt_xdat,'Value') == 3 % make chgcar
    if get(handles.cek_num,'Value') == 0
        [~, mag, geo] = chgcar(pathname);
        make_chgcar(geo,mag,[pathname '_mag']);
    elseif get(handles.cek_num,'Value') == 1
        for n = 1:N
            [~, mag, geo] = chgcar([pathname num2str(n)]);
            make_chgcar(geo,mag,[pathname num2str(n) '_mag']);
        end
    end
elseif get(handles.opt_xdat,'Value') == 4 % make chgcar LNONCOLLINEAR
    if get(handles.cek_num,'Value') == 0
        [~, mag_x, mag_y, mag_z, geo] = magcar(pathname);
        if get(handles.Mx,'Value') == 1
            make_chgcar(geo,mag_x,[pathname '_mag_x']);
        end
        if get(handles.My,'Value') == 1
            make_chgcar(geo,mag_y,[pathname '_mag_y']);
        end
        if get(handles.Mz,'Value') == 1
            make_chgcar(geo,mag_z,[pathname '_mag_z']);
        end
        if get(handles.Mt,'Value') == 1
            mag = sqrt(mag_x.^2+mag_y.^2+mag_z.^2);
            make_chgcar(geo,mag,[pathname '_mag']);
        end
    elseif get(handles.cek_num,'Value') == 1
        for n = 1:N
            [~, mag_x, mag_y, mag_z, geo] = magcar([pathname num2str(n)]);
            if get(handles.Mx,'Value') == 1
                make_chgcar(geo,mag_x,[pathname num2str(n) '_mag_x']);
            end
            if get(handles.My,'Value') == 1
                make_chgcar(geo,mag_y,[pathname num2str(n) '_mag_y']);
            end
            if get(handles.Mz,'Value') == 1
                make_chgcar(geo,mag_z,[pathname num2str(n) '_mag_z']);
            end
            if get(handles.Mt,'Value') == 1
                mag = sqrt(mag_x.^2+mag_y.^2+mag_z.^2);
                make_chgcar(geo,mag,[pathname '_mag']);
            end
        end
    end
end
fprintf('\n>> Makefile Completed!')
fprintf('\n   ======================================= \n   '); toc

function run_trace_Callback(~, ~, handles)
filename = get(handles.disp_trace,'String'); clc; tic
fprintf('   ======================================= \n')
fprintf('>> Reading input files \n')
Nx = str2num(get(handles.Nx_trace,'String'));
Rx = str2num(get(handles.Rx,'String'));
if length(Nx) ~= 3; Nx = Nx(1)*[1 1 1]; end
if length(Rx) ~= 3; Rx = Rx(1)*[1 1 1]; end
muon = str2num(get(handles.muon_trace,'String'));
if get(handles.opt_trace,'Value') == 1
    make_trace(filename,muon,Rx,Nx)
elseif get(handles.opt_trace,'Value') == 2
    make_trace(filename,Nx)
end
fprintf('>> Makefile Completed!'); 
fprintf('\n   ======================================= \n'); toc

function muon_trace_Callback(~, ~, ~)
function muon_trace_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_trace_Callback(~, ~, handles)
if get(handles.opt_trace,'Value') == 1
    set(handles.muon_trace,'BackgroundColor',[1 1 1])
    set(handles.muon_trace,'ForegroundColor',[0 0 0])
    set(handles.Rx,'BackgroundColor',[1 1 1]);
    set(handles.Rx,'ForegroundColor',[0 0 0]); 
else 
    set(handles.muon_trace,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.muon_trace,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.Rx,'BackgroundColor',[0.8 0.8 0.8]);
    set(handles.Rx,'ForegroundColor',[0.8 0.8 0.8]);
end

function opt_trace_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function disp_trace_Callback(~, ~, ~)
function disp_trace_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load_trace_Callback(~, ~, handles)
global filename_trace pathname_trace
[filename_trace, pathname_trace] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname_trace)
    set(handles.disp_trace,'String',[pathname_trace filename_trace])
end

function Rx_Callback(~, ~, ~)
function Rx_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nx_trace_Callback(~, ~, ~)
function Nx_trace_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function geometry = get_poscar(filename)
if ~isnumeric(filename)
	geometry.filename = filename; fid = fopen(filename);
else fid = filename; geometry.filename = fopen(fid);
end
geometry.comment = fgetl(fid); scale = fscanf(fid, '%f',1);
geometry.lattice = fscanf(fid, '%f %f %f', [3 3])'; 
geometry.lattice = geometry.lattice*scale; fgetl(fid); 
line = fgetl(fid); has_symbols_line = false; cartesian = 0;
if sum(isstrprop(line, 'digit')) == 0
    geometry.symbols = regexp(line, '([^ ]*)', 'match');
    line = fgetl(fid); has_symbols_line = true;
else geometry.symbols = {};
end
geometry.atomcount = sscanf(line,'%d');
natoms = sum(geometry.atomcount); line = fgetl(fid);
geometry.selective = 0;
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
            str = regexp(line, '[^ ]*', 'match'); str = str{end}; newelement = false;
            if numel(geometry.symbols) == 0; newelement = true;
            elseif strcmp(geometry.symbols{end},str) == 0; newelement = true; end
            if newelement; geometry.symbols{end+1} = str; end
        end
    end
end
if cartesian == 1
	geometry.coords = geometry.coords*scale;
    geometry.coords = geometry.coords/geometry.lattice;
end 
if ~isnumeric(filename); fclose(fid); end

function [chg, mag, geo] = chgcar(filename)
fid = fopen(filename); fprintf('>> Reading CHGCAR file \n')
geo = get_poscar(fid); mag = [];
natoms = sum(geo.atomcount); fgetl(fid); 
gridsize = fscanf(fid, '%d %d %d', [3 1])';
chg = fscanf(fid, '%f', [prod(gridsize,2) 1])';
chg = reshape(chg,gridsize);
fgetl(fid); pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a' % CHGCAR
    for i = 1:natoms
        line = fgetl(fid);
        nentries = sscanf(line,['augmentation occupancies ' num2str(i) ' %d']);
        fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
    end
    pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
    if line(1) == 'a';
        for i = 1:natoms
            line = fgetl(fid);
            nentries = sscanf(line,['augmentation occupancies (imaginary part) ' num2str(i) ' %d']); 
            fscanf(fid, '%f', [nentries 1]); fgetl(fid); 
        end
    end
    fscanf(fid, '%f', [natoms 1]); fgetl(fid); pos = ftell(fid); fgetl(fid);
    if ~feof(fid)     
        fseek(fid,pos,'bof');
        gridsize = fscanf(fid, '%d %d %d', [3 1])';
        mag = fscanf(fid, '%f', [prod(gridsize,2) 1])';
        mag = reshape(mag,gridsize);
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) > 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag = reshape(mag,gridsize);
end
fclose(fid);

function make_trace(filename,mu,a,N)
fprintf('>> Writing Makefile \n')
geo = get_poscar(filename);
if nargin < 3 % trace all
    x = linspace(0,1,mu(1)); x(end) = [];
    y = linspace(0,1,mu(2)); y(end) = [];
    z = linspace(0,1,mu(3)); z(end) = [];
    [x,y,z] = ndgrid(x,y,z); mu = [];
    x = x(:); y = y(:); z = z(:); 
    for i = 1:length(x); mu = [mu; x(i) y(i) z(i)]; end
else  mu = coord('dir2real',mu,geo);
    x = linspace(mu(1)-a(1)/2, mu(1)+a(1)/2, N(1)); 
    y = linspace(mu(2)-a(2)/2, mu(2)+a(2)/2, N(2)); 
    z = linspace(mu(3)-a(3)/2, mu(3)+a(3)/2, N(3)); mu = [];
    for i = 1:length(x)
        for j = 1:length(y)
            for k = 1:length(z)
                mu = [mu; coord('real2dir',[x(i) y(j) z(k)], geo)];
            end
        end
    end
end
fid = fopen([filename '_TRACE/POSCAR']);
if fid ~= -1; i = 1; 
    fid = fopen([filename '_TRACE(1)/POSCAR']);
    while fid ~= -1; i = i+1; 
        fclose(fid); fid = fopen([filename '_TRACE(' num2str(i) ')/POSCAR']);
    end
    dir = [filename '_TRACE(' num2str(i) ')/'];
    fprintf(['     ' filename '_TRACE(' num2str(i) ')\n'])
else fprintf(['     ' filename '_TRACE\n'])
    dir = [filename '_TRACE/'];
end
mkdir(dir); atomcount = [geo.atomcount; 1];
if ~isempty(geo.symbols)
    geo.symbols = [geo.symbols, 'H'];
end
for i = 1:size(mu,1)
    coords = [geo.coords; mu(i,:)];
    fid = fopen([dir 'POSCAR_' num2str(i)],'w');
    fprintf(fid,[geo.comment '\n']);
    fprintf(fid,'1.0\n');
    fprintf(fid, '%19.16f %19.16f %19.16f\n', geo.lattice');
    if ~isempty(geo.symbols)
        cellfun(@(x) fprintf(fid, '%s ', x), geo.symbols);
        fprintf(fid, '\n');
    end
    fprintf(fid, '%d ', atomcount);
    fprintf(fid, '\nDirect\n');
    fprintf(fid, '%19.16f %19.16f %19.16f\n', coords');
    fclose(fid);
end
atomcount = [geo.atomcount; size(mu,1)];
coords = [geo.coords; mu];
fid = fopen([dir 'POSCAR'],'w');
fprintf(fid,[geo.comment '\n']);
fprintf(fid,'1.0\n');
fprintf(fid, '%19.16f %19.16f %19.16f\n', geo.lattice');
if ~isempty(geo.symbols)
    cellfun(@(x) fprintf(fid, '%s ', x), geo.symbols);
    fprintf(fid, '\n');
end
fprintf(fid, '%d ', atomcount);
fprintf(fid, '\nDirect\n');
fprintf(fid, '%19.16f %19.16f %19.16f\n', coords');
fclose(fid);

function make_poscar(filename, array, opt, muon, R)
fprintf('>> Writing Makefile \n')
if nargin < 3; opt = []; end
geo1 = get_poscar(filename); geo2 = geo1; opt1 = 0;
if nargin > 1
    if isnumeric(array); geo2 = supercell(geo1,array); else opt = array; end
end
if strcmp(opt,'all_T')
    opt = []; for i = 1:length(geo2.atomcount); opt = [opt; 'T']; end
elseif strcmp(opt,'all_F')
    opt = []; for i = 1:length(geo2.atomcount); opt = [opt; 'F']; end
elseif strcmp(opt,'sel_D')
    opt = []; for i = 1:sum(geo2.atomcount); opt = [opt; 'F']; end
elseif strcmp(opt,'static_mu'); opt1=1;
    opt = []; for i = 1:length(geo2.atomcount); opt = [opt; 'F']; end
end
if nargin > 3
    muon = coord('dir2sup',muon,geo1,geo2);
    geo2.atomcount = [geo2.atomcount; size(muon,1)];
    geo2.coords = [geo2.coords; muon];
    if ~isempty(geo2.symbols)
        geo2.symbols = [geo2.symbols, 'H'];
    end
    if ~isnumeric(opt)
        if opt1 == 1; opt = [opt; 'F']; else opt = [opt; 'T']; end
    end
end
fid = fopen(filename);
if fid ~= -1; i = 1; 
    fid = fopen([filename '(1)']);
    while fid ~= -1; i = i+1; 
        fclose(fid); fid = fopen([filename '(' num2str(i) ')']);
    end
    fid = fopen([filename '(' num2str(i) ')'],'w');
    fprintf(['     ' filename '(' num2str(i) ')\n'])
else fid = fopen(filename,'w');
    fprintf(['     ' filename '\n'])
end
fprintf(fid,[geo2.comment '\n']);
fprintf(fid,'1.0\n'); % scale factor
fprintf(fid, '%19.16f %19.16f %19.16f\n', geo2.lattice');
if ~isempty(geo2.symbols)
    cellfun(@(x) fprintf(fid, '%s ', x), geo2.symbols);
    fprintf(fid, '\n');
end
fprintf(fid, '%d ', geo2.atomcount); % number of each species
if ischar(opt); fprintf(fid, '\nSelected Dynamics\nDirect\n');
    if length(opt) > length(geo2.atomcount)
        ion = coord('dir2real',geo2.coords,geo2.lattice);
        if nargin < 5; R = 2; end
        for i = 1:size(geo2.coords,1)
            d = ion(i,:) - ion(end,:);
            if sqrt(sum(d.*d)) <= R; opt(i) = 'T'; end
            fprintf(fid,'%19.16f %19.16f %19.16f %s %s %s\n', ...
                geo2.coords(i,:),opt(i),opt(i),opt(i));
        end
    else N = 1;
        for i = 1:length(geo2.atomcount)
            for j = 1:geo2.atomcount(i);
                fprintf(fid,'%19.16f %19.16f %19.16f %s %s %s\n', ...
                    geo2.coords(N,:),opt(i),opt(i),opt(i)); N = N+1;
            end
        end
    end
else fprintf(fid, '\nDirect\n');
    fprintf(fid, '%19.16f %19.16f %19.16f\n', geo2.coords');
end
fclose(fid);

function make_poscar2(filename, array, opt, muon, R)
fprintf('>> Writing Makefile \n')
if nargin < 3; opt = []; end
geo1 = get_poscar(filename); geo2 = geo1; opt1 = 0;
if nargin > 1
    if isnumeric(array); geo2 = supercell2(geo1,array); else opt = array; end
end
if strcmp(opt,'all_T')
    opt = []; for i = 1:length(geo2.atomcount); opt = [opt; 'T']; end
elseif strcmp(opt,'all_F')
    opt = []; for i = 1:length(geo2.atomcount); opt = [opt; 'F']; end
elseif strcmp(opt,'sel_D')
    opt = []; for i = 1:sum(geo2.atomcount); opt = [opt; 'F']; end
elseif strcmp(opt,'static_mu'); opt1 = 1; 
    opt = []; for i = 1:length(geo2.atomcount); opt = [opt; 'F']; end
end
if nargin > 3
    muon = coord('dir2sup',muon,geo1,geo2);
    geo2.atomcount = [geo2.atomcount; size(muon,1)];
    geo2.coords = [geo2.coords; muon];
    if ~isempty(geo2.symbols); geo2.symbols = [geo2.symbols, 'H']; end
    if ~isnumeric(opt)
        if opt1 == 1; opt = [opt; 'F']; else opt = [opt; 'T']; end
    end
end
fid = fopen(filename);
if fid ~= -1; i = 1; 
    fid = fopen([filename '(1)']);
    while fid ~= -1; i = i+1; 
        fclose(fid); fid = fopen([filename '(' num2str(i) ')']);
    end
    fid = fopen([filename '(' num2str(i) ')'],'w');
    fprintf(['     ' filename '(' num2str(i) ')\n'])
else fid = fopen(filename,'w');
    fprintf(['     ' filename '\n'])
end
fprintf(fid,[geo2.comment '\n']);
fprintf(fid,'1.0\n'); % scale factor
fprintf(fid, '%19.16f %19.16f %19.16f\n', geo2.lattice');
if ~isempty(geo2.symbols)
    cellfun(@(x) fprintf(fid, '%s ', x), geo2.symbols);
    fprintf(fid, '\n');
end
fprintf(fid, '%d ', geo2.atomcount); % number of each species
if ischar(opt)
    fprintf(fid, '\nSelected Dynamics\nDirect\n');
    if length(opt) > length(geo2.atomcount)
        ion = coord('dir2real',geo2.coords,geo2.lattice);
        if nargin < 5; R = 2; end
        for i = 1:size(geo2.coords,1)
            d = ion(i,:) - ion(end,:);
            if sqrt(sum(d.*d)) <= R; opt(i) = 'T'; end
            fprintf(fid,'%19.16f %19.16f %19.16f %s %s %s\n', ...
                geo2.coords(i,:),opt(i),opt(i),opt(i));
        end
    else N = 1;
        for i = 1:length(geo2.atomcount)
            for j = 1:geo2.atomcount(i);
                fprintf(fid,'%19.16f %19.16f %19.16f %s %s %s\n', ...
                    geo2.coords(N,:),opt(i),opt(i),opt(i)); N = N+1;
            end
        end
    end
else fprintf(fid, '\nDirect\n');
    fprintf(fid, '%19.16f %19.16f %19.16f\n', geo2.coords');
end
fclose(fid);

function geo2 = supercell(geo1, array)
geo2 = geo1; geo2.coords = [];
for i = 1:numel(geo1.atomcount)
  for a1 =-(array(1)-1)/2:(array(1)-1)/2
      for a2 =-(array(2)-1)/2:(array(2)-1)/2
           for a3 =-(array(3)-1)/2:(array(3)-1)/2
              start = sum(geo1.atomcount(1:i-1))+1;
              geo2.coords = [geo2.coords; geo1.coords(start:start + ...
                  geo1.atomcount(i)-1,:)+repmat([a1 a2 a3],geo1.atomcount(i),1)];
          end
      end
  end        
end
geo2.atomcount = geo1.atomcount*prod(array);
for i = 1:3
    geo2.coords(:,i) = geo2.coords(:,i)/array(i);
    geo2.lattice(i,:) = geo2.lattice(i,:)*array(i);
end

function geo2 = supercell2(geo1, array)
geo2 = geo1; geo2.coords = [];
for i = 1:numel(geo1.atomcount)
  for a1 =0:array(1)-1
      for a2 =0:array(2)-1
           for a3 =0:array(3)-1
              start = sum(geo1.atomcount(1:i-1))+1;
              geo2.coords = [geo2.coords; geo1.coords(start:start + ...
                  geo1.atomcount(i)-1,:)+repmat([a1 a2 a3],geo1.atomcount(i),1)];
          end
      end
  end        
end
geo2.atomcount = geo1.atomcount*prod(array);
for i = 1:3
    geo2.coords(:,i) = geo2.coords(:,i)/array(i);
    geo2.lattice(i,:) = geo2.lattice(i,:)*array(i);
end

function pos = coord(mode,pos,geo,grd)
if nargin == 2; geo = []; elseif nargin < 4; grd = geo; end
if isnumeric(geo); latt = geo;
elseif ischar(geo); geo = get_poscar(geo); latt = geo.lattice;
else latt = geo.lattice;
end
switch mode
    case 'dir2real'
        for i = 1:size(pos,1)
            pos(i,1:3) = pos(i,1)*latt(1,:)+pos(i,2)*latt(2,:)+pos(i,3)*latt(3,:);
        end
    case 'dir2sup'
        if isnumeric(grd); geo2 = supercell(geo,grd); else geo2 = grd; end
        for i = 1:size(pos,1); pos(i,:) = pos(i,:)-0.5*[1 1 1]; end
        pos = coord('dir2real',pos,geo); pos = coord('real2dir',pos,geo2);
        for i = 1:size(pos,1); pos(i,:) = pos(i,:)+0.5*[1 1 1]; end
    case 'real2dir'
        latt = inv(latt);
        for i = 1:size(pos,1)
            pos(i,1:3) = pos(i,1)*latt(1,:)+pos(i,2)*latt(2,:)+pos(i,3)*latt(3,:);
        end
end

function [chg, mag_x, mag_y, mag_z, geo] = magcar(filename)
fprintf('>> Reading CHGCAR file \n')
fid = fopen(filename); geo = get_poscar(fid);
natoms = sum(geo.atomcount); fgetl(fid); 
gridsize = fscanf(fid, '%d %d %d', [3 1])';
chg = fscanf(fid, '%f', [prod(gridsize,2) 1])';
chg = reshape(chg,gridsize); mag_x = []; mag_y = []; mag_z = [];
fgetl(fid); pos = ftell(fid); line = fgetl(fid); fseek(fid,pos,'bof');
if line(1) == 'a'
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
        mag_x = reshape(mag_x,gridsize);
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) > 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_x = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_x = reshape(mag_x,gridsize);
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
        mag_y = reshape(mag_y,gridsize);
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) > 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_y = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_y = reshape(mag_y,gridsize);
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
        mag_z = reshape(mag_z,gridsize);
    end
elseif ~isnumeric(line)
    if size(str2num(line),2) > 3; % LOCPOT
        fscanf(fid, '%f', [natoms 1]);
    end
    gridsize = fscanf(fid, '%d %d %d', [3 1])';
    mag_z = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    mag_z = reshape(mag_z,gridsize);
end
fclose(fid);

function center2_Callback(~, ~, handles)
set(handles.center0,'Value',0); set(handles.center0,'ForegroundColor',[0 0 0])
set(handles.center2,'Value',1); set(handles.center2,'ForegroundColor',[0 0 1])

function center0_Callback(~, ~, handles)
set(handles.center2,'Value',0); set(handles.center2,'ForegroundColor',[0 0 0])
set(handles.center0,'Value',1); set(handles.center0,'ForegroundColor',[0 0 1])

function opt_xdat_Callback(~, ~, handles)
if get(handles.opt_xdat,'Value') < 3
    set(handles.load_xdat,'String','Folder')
    set(handles.uipanel3,'Title','XDATCAR')
else set(handles.uipanel3,'Title','SPIN')
    set(handles.load_xdat,'String','CHGCAR')
end
if get(handles.opt_xdat,'Value') < 4
    set(handles.Mx,'Value',0); set(handles.My,'Value',0)
    set(handles.Mz,'Value',0); set(handles.Mt,'Value',0)
end

function opt_xdat_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function N_xdat_Callback(~, ~, ~)
function N_xdat_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cek_num_Callback(~, ~, handles)
if get(handles.cek_num,'Value') == 1
    set(handles.N_xdat,'BackgroundColor',[1 1 1])
    set(handles.N_xdat,'ForegroundColor',[0 0 0])
else set(handles.N_xdat,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.N_xdat,'ForegroundColor',[0.8 0.8 0.8])
end

function Mt_Callback(~, ~, handles)
if get(handles.opt_xdat,'Value') < 4; set(handles.Mt,'Value',0); end

function Mx_Callback(~, ~, handles)
if get(handles.opt_xdat,'Value') < 4; set(handles.Mx,'Value',0); end

function My_Callback(~, ~, handles)
if get(handles.opt_xdat,'Value') < 4; set(handles.My,'Value',0); end

function Mz_Callback(~, ~, handles)
if get(handles.opt_xdat,'Value') < 4; set(handles.Mz,'Value',0); end
