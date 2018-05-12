function varargout = GUI_Dipole_Input(varargin)

% GUI_MAIN MATLAB code by Edi Suprayoga
% feel free to contact me at suprayoga.edi@gmail.com
% last edit 16 Mar 2016

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Dipole_Input_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Dipole_Input_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout; [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else gui_mainfcn(gui_State, varargin{:});
end

function GUI_Dipole_Input_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject; guidata(hObject, handles);
set(handles.figure1,'Color',[0.8 0.8 0.8])

function varargout = GUI_Dipole_Input_OutputFcn(~, ~, handles) 
varargout{1} = handles.output; set(handles.autorized,'Value',1);

function run_Callback(~, ~, handles)
pathname = pwd; clc
fprintf('>> Writing INPUT file \n');
opt_job = get(handles.opt_job,'Value');
if opt_job == 1; job = '';
elseif opt_job == 2; job = '_zpe';
elseif opt_job == 3; job = '_test';
elseif opt_job == 4; job = '_make';
end
fid = fopen([pathname '/INPUT' job '.m']);
if fid ~= -1; i = 1; fid = fopen([pathname '/INPUT' job '_1.m']);
    while fid ~= -1; i = i+1; fclose(fid); 
        fid = fopen([pathname '/INPUT' job '_' num2str(i) '.m']);
    end
    fid = fopen([pathname '/INPUT' job '_' num2str(i) '.m'],'w');
    fprintf(['     ' pathname '/INPUT' job '_' num2str(i) '.m \n'])
else fid = fopen([pathname '/INPUT' job '.m'],'w');
    fprintf(['     ' pathname '/INPUT' job '.m \n'])
end
if get(handles.gpu,'Value') == 1; 
    fprintf(fid,'gpu = 1; %% using GPU \n');
else fprintf(fid,'gpu = 0; %% using CPU \n');
end
opt_pos = get(handles.opt_pos,'Value');
if opt_pos == 1;     fprintf(fid,'opt_pos = 1; %% Perfect crystal \n');
elseif opt_pos == 2; fprintf(fid,'opt_pos = 2; %% Defect crystal \n');
end
opt1 = get(handles.opt1,'Value'); opt2 = get(handles.opt2,'Value'); 
opt3 = get(handles.opt3,'Value');
if opt1 == 1; file1 = 'POSCAR'; file2 = 'POSCAR_out';
    fprintf(fid,'opt_spn = 1; %% Manual input spin \n');
elseif opt2 == 1; file1 = 'folder'; file2 = 'folder_out';
    fprintf(fid,'opt_spn = 2; %% Input spin from VASP\n');
elseif opt3 == 1; file1 = 'CHGCAR'; file2 = 'CHGCAR_out';
    fprintf(fid,'opt_spn = 3; %% Spin distributed (from CHGCAR) \n');
end
if opt_job == 1;     fprintf(fid,'opt_job = 1; %% Dipole fields \n');
elseif opt_job == 2; fprintf(fid,'opt_job = 2; %% Dipole fields with ZPE \n');
    if opt1 == 1; file1 = 'LOCPOT'; else file1 = 'folder'; end
elseif opt_job == 3; fprintf(fid,'opt_job = 3; %% Converge test \n');
elseif opt_job == 4; fprintf(fid,'opt_job = 4; %% Makefile \n');
end
fprintf(fid,'\n%% Input files \n'); lim = ',[]';
fprintf(fid,[file1 ' = ''' get(handles.pos1,'String') '''; \n']);
if opt_pos == 2
    fprintf(fid,[file2 ' = ''' get(handles.pos2,'String') '''; \n']);
end
if opt_job < 4; fprintf(fid,'\n%% Muon position \n');
    coord = get(handles.opt_ax,'Value');
    fprintf(fid,['muon = [' get(handles.muon,'String') ']; \n']);
    if coord == 1;     fprintf(fid,'coord = 1; %% Direct \n');
    elseif coord == 2; fprintf(fid,'coord = 2; %% Cartesian \n');
    elseif coord == 3; fprintf(fid,'coord = 3; %% from POSCAR \n');
    end
    if opt_job < 3
        fprintf(fid,['\nR_dip = ' get(handles.Rmax,'String') '; %% Calculation range \n']);
    end
end
if opt1 == 1;
    if length(str2num(get(handles.opt_in,'String'))) == 1
        fprintf(fid,['ion_no = ' get(handles.opt_in,'String') '; \n']);
    else fprintf(fid,['ion_no = [' get(handles.opt_in,'String') ']; \n']);
    end
    fprintf(fid,'\n%% Spin structure \n');
    if length(str2num(get(handles.S1,'String'))) == 1
        fprintf(fid,['mag1  = ' get(handles.S1,'String') '; %% Bohr magneton \n']);
    else fprintf(fid,['mag1  = [' get(handles.S1,'String') ']; %% Bohr magneton \n']);
    end
    fid1 = fopen(get(handles.mag1,'String'));
    if fid1 == -1; fprintf(fid,['spin1 = [' get(handles.mag1,'String') ']; \n']);
    else fclose(fid1);
        fprintf(fid,['spin1 = ''' get(handles.mag1,'String') '''; \n']);
    end
    if opt_pos == 2
        if length(str2num(get(handles.S2,'String'))) == 1
            fprintf(fid,['mag2  = ' get(handles.S2,'String') '; %% Bohr magneton \n']);
        else fprintf(fid,['mag2  = [' get(handles.S2,'String') ']; %% Bohr magneton \n']);
        end
        fid1 = fopen(get(handles.mag2,'String'));
        if fid1 == -1; fprintf(fid,['spin2 = [' get(handles.mag2,'String') ']; \n']);
        else fclose(fid1);
            fprintf(fid,['spin2 = ''' get(handles.mag2,'String') '''; \n']);
        end
    end
elseif opt2 == 1; lim = str2num(get(handles.opt_in,'String'));
    if length(lim) == 1;
        fprintf(fid,['cut_spin = ' num2str(lim) '; %%cut off spin \n']); 
    else fprintf(fid,['cut_spin = [' num2str(lim) ']; %%cut off spin \n']);
    end
    lim = ',cut_spin';
end
if opt_job == 2 % zpe
    fprintf(fid,'\n%% ZPE option \n');
    Rzpe = str2num(get(handles.Rzpe,'String'));
    Nzpe = str2num(get(handles.Nzpe,'String'));
    if get(handles.opt_grid,'Value') == 2; Nzpe = round(Rzpe./Nzpe); end
    if length(Rzpe) == 1; Rzpe = num2str(Rzpe);
        fprintf(fid,['zpe_range = [' Rzpe ' ' Rzpe ' ' Rzpe ']; %% Angstrom\n']);
    else fprintf(fid,['zpe_range = [' num2str(Rzpe) ']; %% Angstrom \n']);
    end
    if length(Nzpe) == 1; Nzpe = num2str(Nzpe);
        fprintf(fid,['zpe_ngrid = [' Nzpe ' ' Nzpe ' ' Nzpe ']; %% Number of grids\n']);
    else fprintf(fid,['zpe_ngrid = [' num2str(Nzpe) ']; %% Number of grids \n']);
    end
elseif opt_job == 3 % test
    fprintf(fid,'\n%% Converge test option \n');
    fprintf(fid,['R_min = ' get(handles.dmin,'String') '; %% Angstrom \n']);
    fprintf(fid,['R_max = ' get(handles.dmax,'String') '; %% Angstrom \n']);
    fprintf(fid,['N_cal = ' get(handles.del,'String') '; %% Number of calculation \n']);
elseif opt_job == 4 % makefile
    fprintf(fid,'\n%% Makefile option \n');
    Nx = get(handles.Nx,'String'); Ny = get(handles.Ny,'String'); Nz = get(handles.Nz,'String');
    fprintf(fid,['N_grid = [' Nx ' ' Ny ' ' Nz ']; %% Number of grids \n']);
    fprintf(fid,['R_dip = ' get(handles.Rmax,'String') '; %% Angstrom \n']);
end
fprintf(fid,'\n%% ===== END OF INPUT FILE ===== %% \n');
if opt_job == 4; muon = ',[]'; coord = ',[]'; else muon = ',muon'; coord = ',coord'; end
if opt_job == 3; R_dip = ',R_max'; else R_dip = ',R_dip'; end
if opt_pos == 1; file2 = '[]'; end
if opt1 == 1; ion_no = ',ion_no'; mag1 = 'mag1'; spin1 = ',spin1';
    if opt_pos == 2; mag2 = ',mag2'; spin2 = ',spin2'; else mag2 = ',[]'; spin2 = ',[]'; end;
else ion_no = ',[]'; mag1 = '[]'; spin1 = ',[]'; mag2 = ',[]'; spin2 = ',[]';
end
if opt_job == 2; zpe_range = ',zpe_range'; zpe_ngrid = ',zpe_ngrid'; 
else zpe_range = ',[]'; zpe_ngrid = ',[]'; 
end
if opt_job == 3; R_min = ',R_min'; N_cal = ',N_cal'; else R_min = ',[]'; N_cal = ',[]'; end
if opt_job == 4; N_grid = ',N_grid'; else N_grid = ',[]'; end
fprintf(fid,['main(opt_pos,opt_spn,opt_job' muon coord R_dip ',' file1 ',' file2 ion_no ', ... \n']);
fprintf(fid,['    ' mag1 spin1 mag2 spin2 zpe_range zpe_ngrid R_min N_cal N_grid lim ',gpu); \n']);
fclose(fid);

function Rmax_Callback(~, ~, ~)
function Rmax_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cpu_Callback(~, ~, handles)
set(handles.cpu,'Value',1); set(handles.gpu,'Value',0);

function gpu_Callback(~, ~, handles)
set(handles.cpu,'Value',0); set(handles.gpu,'Value',1);

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
