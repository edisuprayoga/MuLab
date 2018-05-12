function varargout = GUI_DOS(varargin)

% GUI_DOS MATLAB code by Edi Suprayoga
% feel free to contact me at suprayoga.edi@gmail.com
% last edit 3 Mar 2016

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_DOS_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_DOS_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout; [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else gui_mainfcn(gui_State, varargin{:});
end

function GUI_DOS_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject; guidata(hObject, handles);
set(handles.figure1,'Color',[0.8 0.8 0.8])

function varargout = GUI_DOS_OutputFcn(~, ~, handles) 
varargout{1} = handles.output; set(handles.autorized,'Value',1);

function cek_reverse_Callback(~, ~, ~)
function load_file_Callback(~, ~, handles)
global filename pathname
[filename, pathname] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname); set(handles.disp_pathname,'String',pathname); end

function run_plot_Callback(~, ~, handles)
tic; clc; fprintf('   ============================================ \n')
no = str2num(get(handles.disp_no,'String'));
ax = str2num(get(handles.disp_axis,'String'));
lwidth = str2num(get(handles.disp_width,'String'));
Msize = str2num(get(handles.disp_msize,'String'));
pathname = get(handles.disp_pathname,'String'); pathname = [pathname '/'];
opt_dos = get(handles.opt_dos,'Value'); geo = poscar(pathname);
[en, dos, ef, pdos] = doscar(pathname); [eig,k,~] = eigenval(pathname); 
if ~isempty(dos); eig = eig-ef; en = en-ef; end
file = textread([pathname 'KPOINTS'],'%s','delimiter','\n','whitespace','');
line = cell2mat(file(3)); kpoint = linspace(0,1,size(eig,1));
if line(1) == 'l' || line(1) == 'L'
    Npt = cell2mat(textscan(cell2mat(file(2)),'%f64')); kpt=[]; d=[];
    for i = 5:length(file)
        if ~isempty(cell2mat(file(i)))
            kpt = [kpt; cell2mat(textscan(cell2mat(file(i)),'%f64 %f64 %f64'))];
        end
    end
    for i = 1:size(kpt,1)-1
        r = kpt(i,:)-kpt(i+1,:); r = sqrt(r(1)^2+r(2)^2+r(3)^2);
        j = i; while j/2 >= 1; j = j-2; end
        if j == 0 && r ~= 0; r = 0; end; d = [d; r];
    end
    d = d/sum(d); b = 0; r = []; kpoint = [];
    for i = 1:size(d,1); j = i;
        while j/2 >= 1; j = j-2; end; b = [b b(end) + d(i)];
        if j == 1 && d(i) ~=0; r = [r linspace(b(end-1), b(end), Npt)]; end
    end
    if length(r) ~= size(eig,1)
        b = linspace(0,1,size(eig,1)/Npt+1); r = [];
        for i = 1:size(b,2)-1; r = [r linspace(b(i),b(i+1),Npt)]; end; b(1) = [];
    end
    kpoint.line = b; kpoint.k = r;
end
if get(handles.plot_clc,'Value') == 1; figure(1); cla; hold on
elseif get(handles.plot_hold,'Value') == 1; figure(1); hold on
else figure; hold on
end
ion = []; nor = 0; j0 = []; j01 = []; j1 = []; j2 = []; y1 = 0; y2 = 0;
if get(handles.opt_no,'Value') == 1 % All
    no = 1:size(geo.atomcount,1); nor = sum(geo.atomcount);
    for i = 1:length(no)
        ion = [ion sum(geo.atomcount(1:no(i)-1))+1:sum(geo.atomcount(1:no(i)))];
    end
elseif get(handles.opt_no,'Value') == 2 % by Element
    for i = 1:length(no)
        ion = [ion sum(geo.atomcount(1:no(i)-1))+1:sum(geo.atomcount(1:no(i)))];
        nor = nor + geo.atomcount(no(i));
    end
elseif get(handles.opt_no,'Value') == 3 % individual ion
    for i = 1:length(no); ion = [ion no(i)]; nor = nor+1; end
end
if isempty(pdos)
    y1 = dos(:,1); y2 = []; if size(dos,2) > 1; y2 = dos(:,2); end
elseif get(handles.opt_no,'Value') == 1 && get(handles.cek_all,'Value') == 1 ...
        && get(handles.cek_pband,'Value') ~= 1
    y1 = dos(:,1); y2 = []; if size(dos,2) > 1; y2 = dos(:,2); end
    if size(pdos,2) > 9 && size(pdos,2) < 17; j0 = 1:9; 
    elseif size(pdos,2) > 18; j0 = 1:16; end
else
    if size(pdos,2) < 33 % spin collinear
        if (get(handles.orb_s,'Value') == 1);   j1 = [j1 1];  j2 = [j2 2];  j0 = [j0 1]; end
        if (get(handles.orb_px,'Value') == 1);  j1 = [j1 7];  j2 = [j2 8];  j0 = [j0 4]; end
        if (get(handles.orb_py,'Value') == 1);  j1 = [j1 3];  j2 = [j2 4];  j0 = [j0 2]; end
        if (get(handles.orb_pz,'Value') == 1);  j1 = [j1 5];  j2 = [j2 6];  j0 = [j0 3]; end
        if (get(handles.orb_dxy,'Value') == 1); j1 = [j1 9];  j2 = [j2 10]; j0 = [j0 5]; end
        if (get(handles.orb_dxz,'Value') == 1); j1 = [j1 15]; j2 = [j2 16]; j0 = [j0 8]; end
        if (get(handles.orb_dyz,'Value') == 1); j1 = [j1 11]; j2 = [j2 12]; j0 = [j0 6]; end
        if (get(handles.orb_dz2,'Value') == 1); j1 = [j1 13]; j2 = [j2 14]; j0 = [j0 7]; end
        if (get(handles.orb_dx2,'Value') == 1); j1 = [j1 17]; j2 = [j2 18]; j0 = [j0 9]; end
        if size(pdos,2) > 9 && size(pdos,2) < 17 % f-elements para
            if (get(handles.orb_f_3,'Value') == 1); j0 = [j0 10]; end
            if (get(handles.orb_f_2,'Value') == 1); j0 = [j0 11]; end
            if (get(handles.orb_f_1,'Value') == 1); j0 = [j0 12]; end
            if (get(handles.orb_f0,'Value') == 1);  j0 = [j0 13]; end
            if (get(handles.orb_f1,'Value') == 1);  j0 = [j0 14]; end
            if (get(handles.orb_f2,'Value') == 1);  j0 = [j0 15]; end
            if (get(handles.orb_f3,'Value') == 1);  j0 = [j0 16]; end
        elseif size(pdos,2) > 18 % f-elements spin
            if get(handles.orb_f_3,'Value') == 1; j1 = [j1 19]; j2 = [j2 20]; j0 = [j0 10]; end
            if get(handles.orb_f_2,'Value') == 1; j1 = [j1 21]; j2 = [j2 22]; j0 = [j0 11]; end
            if get(handles.orb_f_1,'Value') == 1; j1 = [j1 23]; j2 = [j2 24]; j0 = [j0 12]; end
            if get(handles.orb_f0,'Value') == 1;  j1 = [j1 25]; j2 = [j2 26]; j0 = [j0 13]; end
            if get(handles.orb_f1,'Value') == 1;  j1 = [j1 27]; j2 = [j2 28]; j0 = [j0 14]; end
            if get(handles.orb_f2,'Value') == 1;  j1 = [j1 29]; j2 = [j2 30]; j0 = [j0 15]; end
            if get(handles.orb_f3,'Value') == 1;  j1 = [j1 31]; j2 = [j2 32]; j0 = [j0 16]; end
        end; j01 = j0;
    else % spin non-collinear
        if get(handles.orb_s,'Value') == 1;   j01 = [j01 1];  j0 = [j0 1]; end
        if get(handles.orb_px,'Value') == 1;  j01 = [j01 13]; j0 = [j0 4]; end
        if get(handles.orb_py,'Value') == 1;  j01 = [j01 5];  j0 = [j0 2]; end
        if get(handles.orb_pz,'Value') == 1;  j01 = [j01 9];  j0 = [j0 3]; end
        if get(handles.orb_dxy,'Value') == 1; j01 = [j01 17]; j0 = [j0 5]; end
        if get(handles.orb_dxz,'Value') == 1; j01 = [j01 29]; j0 = [j0 8]; end
        if get(handles.orb_dyz,'Value') == 1; j01 = [j01 21]; j0 = [j0 6]; end
        if get(handles.orb_dz2,'Value') == 1; j01 = [j01 25]; j0 = [j0 7]; end
        if get(handles.orb_dx2,'Value') == 1; j01 = [j01 33]; j0 = [j0 9]; end
        if size(pdos,2) > 36 % f-elements
            if (get(handles.orb_f_3,'Value') == 1); j01 = [j01 37]; j0 = [j0 10]; end
            if (get(handles.orb_f_2,'Value') == 1); j01 = [j01 41]; j0 = [j0 11]; end
            if (get(handles.orb_f_1,'Value') == 1); j01 = [j01 45]; j0 = [j0 12]; end
            if (get(handles.orb_f0,'Value') == 1);  j01 = [j01 49]; j0 = [j0 13]; end
            if (get(handles.orb_f1,'Value') == 1);  j01 = [j01 53]; j0 = [j0 14]; end
            if (get(handles.orb_f2,'Value') == 1);  j01 = [j01 57]; j0 = [j0 15]; end
            if (get(handles.orb_f3,'Value') == 1);  j01 = [j01 61]; j0 = [j0 16]; end
        end
    end
    if size(dos,2) == 1; y2 = []; % para
        for i = ion; for j = j01; y1 = y1 + pdos(:,j,i); end; end
    else % Spin 
        for i = ion
            for j = j1; y1 = y1 + pdos(:,j,i); end
            for j = j2; y2 = y2 + pdos(:,j,i); end
        end
    end
end
if get(handles.cek_norm,'Value') == 1; % normalisation
    y1 = y1/nor; y2 = y2/nor; 
elseif get(handles.cek_normby,'Value') == 1; % normalisation
    nor = str2double(get(handles.disp_norm,'String'));
    y1 = y1*nor; y2 = y2*nor; 
end
if get(handles.cek_reverse,'Value') == 1 && size(eig,3) > 1; % reverse axis
    y0 = y1; y1 = y2; y2 = y0; y0 = eig(:,:,1); 
    eig(:,:,1) = eig(:,:,2); eig(:,:,2) = y0;
end
if (get(handles.opt_line,'Value') == 1);     lcol = 'k';
elseif (get(handles.opt_line,'Value') == 2); lcol = 'b';
elseif (get(handles.opt_line,'Value') == 3); lcol = 'r';
elseif (get(handles.opt_line,'Value') == 4); lcol = 'g';
elseif (get(handles.opt_line,'Value') == 5); lcol = 'y';
elseif (get(handles.opt_line,'Value') == 6); lcol = 'm';
elseif (get(handles.opt_line,'Value') == 7); lcol = 'c';
elseif (get(handles.opt_line,'Value') == 8); lcol = 0;
end
if (get(handles.opt_marker,'Value') == 1);     mcol = 1;
elseif (get(handles.opt_marker,'Value') == 2); mcol = 'k';
elseif (get(handles.opt_marker,'Value') == 3); mcol = 'b';
elseif (get(handles.opt_marker,'Value') == 4); mcol = 'r';
elseif (get(handles.opt_marker,'Value') == 5); mcol = 'g';
elseif (get(handles.opt_marker,'Value') == 6); mcol = 'y';
elseif (get(handles.opt_marker,'Value') == 7); mcol = 'm';
elseif (get(handles.opt_marker,'Value') == 8); mcol = 'c';
elseif (get(handles.opt_marker,'Value') == 9); mcol = 2;
elseif (get(handles.opt_marker,'Value') == 10); mcol = 0;
end
if get(handles.cek_axis,'Value') == 0;
    if y1 == 0; ax(2) = 1; ax(3) = min(eig(:)); ax(4) = max(eig(:));
         if ~isempty(y2); ax(1) = -1; else ax(1) = 0; end 
    else if ~isempty(y2); ax(1) = -max(y2(2:end)); else ax(1) = 0; end 
        ax(2) = max(y1(2:end)); ax(3) = min(en(2:end)); ax(4) = max(en(2:end));
    end
elseif opt_dos == 1; ax = [ax(3) ax(4) ax(1) ax(2)];
end
if opt_dos == 1; fprintf('>> Plot DOS \n'); % DOS
    if isnumeric(lcol); lcol = mcol; end
    plt_dos(en,[y1 y2],lcol,lwidth,[ax(3) ax(4) ax(1) ax(2)])
elseif opt_dos == 2; fprintf('>> Plot DOS vertical \n'); % DOS ver
    if isnumeric(lcol); lcol = mcol; end
    plt_dosver(en,[y1 y2],lcol,lwidth,ax)
    ylabel('E-E_F (eV)','FontSize',18)
elseif opt_dos == 3
    if get(handles.cek_pband,'Value') == 0 % BANDS
        fprintf('>> Plot Band Structures \n');
        if size(eig,3) == 2; subplot(1,2,1); ylabel('E-E_F (eV)','FontSize',18)
            plt_band(kpoint,eig(:,:,2),lcol,mcol,ax); subplot(1,2,2)
        else ylabel('E-E_F (eV)','FontSize',18)
        end
        plt_band(kpoint,eig(:,:,1),lcol,mcol,ax)
    elseif get(handles.cek_pband,'Value') == 1 % PBANDS
        if size(pdos,2) < 33; [pband1, pband2] = procar(pathname);
        else pband1 = procar_lnon(pathname); end
        if get(handles.cek_reverse,'Value') == 1 && size(eig,3) > 1;
            y0 = pband1; pband1 = pband2; pband2 = y0; % reverse axis
        end
        fprintf('>> Plot partial Bands \n')
        if size(eig,3) == 2; subplot(1,2,1)
            plt_pband(kpoint,eig(:,:,2),pband1,ion,j0,lcol,mcol,1,Msize,ax)
            ylabel('E-E_F (eV)','FontSize',18); subplot(1,2,2)
            plt_pband(kpoint,eig(:,:,1),pband2,ion,j0,lcol,mcol,1,Msize,ax)
        else plt_pband(kpoint,eig(:,:,1),pband1,ion,j0,lcol,mcol,1,Msize,ax)
            ylabel('E-E_F (eV)','FontSize',18)
        end
    end
elseif opt_dos == 4
    if get(handles.cek_pband,'Value') == 0 % DOS & BANDS
        fprintf('>> Plot DOS & Band Structures \n');
        if size(eig,3) == 2; subplot(1,3,1)
            plt_band(kpoint,eig(:,:,2),lcol,mcol,ax)
            ylabel('E-E_F (eV)','FontSize',18); subplot(1,3,2); 
            plt_band(kpoint,eig(:,:,1),lcol,mcol,ax); subplot(1,3,3)
        else subplot(1,2,1); plt_band(kpoint,eig(:,:,1),lcol,mcol,ax)
            ylabel('E-E_F (eV)','FontSize',18); subplot(1,2,2)
        end
    elseif get(handles.cek_pband,'Value') == 1 % PDOS & PBANDS
        if size(pdos,2) < 33; [pband1, pband2] = procar(pathname);
        else pband1 = procar_lnon(pathname); end
        if get(handles.cek_reverse,'Value') == 1 && size(eig,3) > 1;
            y0 = pband1; pband1 = pband2; pband2 = y0; % reverse axis
        end
        fprintf('>> Plot DOS & partial Bands \n')
        if size(eig,3) == 2; subplot(1,3,1)
            plt_pband(kpoint,eig(:,:,2),pband1,ion,j0,lcol,mcol,1,Msize,ax)
            ylabel('E-E_F (eV)','FontSize',18); subplot(1,3,2)
            plt_pband(kpoint,eig(:,:,1),pband2,ion,j0,lcol,mcol,1,Msize,ax)
            subplot(1,3,3)
        else subplot(1,2,1)
            plt_pband(kpoint,eig(:,:,1),pband1,ion,j0,lcol,mcol,1,Msize,ax)
            ylabel('E-E_F (eV)','FontSize',18); subplot(1,2,2)
        end
    end
    if isnumeric(lcol); lcol = mcol; end; plt_dosver(en,[y1 y2],lcol,lwidth,ax)
elseif opt_dos == 5; fprintf('>> Plot DOS & Kpoints \n');  % DOS & KPT
    latt = geo.lattice; a1 = latt(1,:); a2 = latt(2,:); a3 = latt(3,:);
    b1 = 2*pi*cross(a2,a3)/dot(a1,cross(a2,a3));
    b2 = 2*pi*cross(a3,a1)/dot(a1,cross(a2,a3));
    b3 = 2*pi*cross(a1,a2)/dot(a1,cross(a2,a3)); latt = [b1; b2; b3];
    for j = 1:size(k,1)
        k(j,1:3) = k(j,1)*latt(1,:)+k(j,2)*latt(2,:)+k(j,3)*latt(3,:);
    end
    if size(eig,3) == 2
        subplot(1,4,1); plot3(k(:,1),k(:,2),k(:,3),'o')
        xlabel('kx'); ylabel('ky'); zlabel('kz'); 
        drawbox(latt); daspect([1 1 1])
        subplot(1,4,2); plt_band(kpoint,eig(:,:,2),lcol,mcol,ax)
        ylabel('E-E_F (eV)','FontSize',18)
        subplot(1,4,3); plt_band(kpoint,eig(:,:,1),lcol,mcol,ax)
        subplot(1,4,4)
    else
        subplot(1,3,1); plot3(k(:,1),k(:,2),k(:,3),'o')
        xlabel('kx'); ylabel('ky'); zlabel('kz')
        drawbox(latt); daspect([1 1 1])
        subplot(1,3,2); plt_band(kpoint,eig(:,:,1),lcol,mcol,ax)
        ylabel('E-E_F (eV)','FontSize',18)
        subplot(1,3,3)
    end
    if isnumeric(lcol); lcol = mcol; end; plt_dosver(en,[y1 y2],lcol,lwidth,ax)
end
fprintf('   ============================================ \n   '); toc

function plt_band(kpoint,eig,lcol,mcol,ax)
hold on; box on
if mcol == 1 && isnumeric(lcol); mcol = 'k'; 
elseif mcol == 1; mcol = lcol; elseif mcol == 2; mcol = 'k'; 
end
if ~isnumeric(kpoint); b = kpoint.line; kpt = kpoint.k;
else kpt = kpoint; end
if ~isnumeric(lcol) && ~isnumeric(mcol)
    plot(kpt,eig(:,:,1),'.-','Color',lcol,'MarkerEdgeColor',mcol,'MarkerSize',5)
elseif ~isnumeric(lcol) &&  isnumeric(mcol)
    plot(kpt,eig(:,:,1),'-','Color',lcol)
elseif  isnumeric(lcol) && ~isnumeric(mcol)
    plot(kpt,eig(:,:,1),'.','MarkerEdgeColor',mcol,'MarkerSize',5)
end
c = min(eig(:)); if ax(3) < c; c = ax(3); end
d = max(eig(:)); if ax(4) > d; d = ax(4); end
if ~isnumeric(kpoint)
    for i = 1:length(b)-1; plot(b(i)*[1 1],[c d],'k--'); end
end
plot([0 1],[0 0],'k--'); xlabel('Kpoints','FontSize',18)
axis([0 1 ax(3) ax(4)]); box on

function plt_pband(kpoint,eig,pband,ion,j0,lcol,mcol,Lwidth,Msize,ax)
hold on; box on
if ~isnumeric(kpoint); kpt = kpoint.k; b = kpoint.line; else kpt = kpoint; end
y1 = zeros(size(pband,1), size(pband,2));
for i = ion; for j = j0; y1(:,:) = y1 + pband(:,:,i,j); end; end
Msize = 10*Msize/max(y1(:));
for i = 1:size(pband,2)
    x = pband(:,i,ion,j0);
    if sum(x(:)) > 0 && ~isnumeric(lcol)
        plot(kpt,eig(:,i,1), lcol, 'LineWidth', Lwidth);
    end
    for j = 1:size(pband,1)
        if y1(j,i) > 0 && ~isnumeric(mcol)
            scatter(kpt(j),eig(j,i,1),Msize*y1(j,i),'filled',mcol);
        elseif y1(j,i) > 0 && mcol == 1
            scatter(kpt(j),eig(j,i,1),Msize*y1(j,i),y1(j,i),'filled');
        elseif y1(j,i) > 0 && mcol == 2
            scatter(kpt(j),eig(j,i,1),Msize,y1(j,i),'filled');
        end
    end
end
c = min(eig(:)); if ax(3) < c; c = ax(3); end
d = max(eig(:)); if ax(4) > d; d = ax(4); end
if ~isnumeric(kpoint)
    for i = 1:length(b)-1; plot(b(i)*[1 1],[c d],'k--'); end
end
plot([0 1],[0 0],'k--'); xlabel('Kpoints','FontSize',18)
axis([0 1 ax(3) ax(4)]); box on

function plt_dos(en,dos,lcol,lwidth,ax)
hold on; box on
if isnumeric(lcol); lcol = 'k'; end;
if size(dos,2) > 1
    plot(en,dos(:,1),lcol,'LineWidth',lwidth); a = -max(dos(:,2));
    plot(en,-dos(:,2),lcol,'LineWidth',lwidth); b = max(dos(:,1)); 
else plot(en,dos(:,1),lcol,'LineWidth',lwidth)
    a = 0; b = max(dos(:,1));
end
if ax(3) < a; a = ax(3); end; if ax(4) > b; b = ax(4); end
c = min(en(2:end)); if ax(1) < c; c = ax(1); end
d = max(en(2:end)); if ax(2) > d; d = ax(2); end
plot([0 0],[a b],'k--'); plot([c d],[0 0],'k') % horizontal
xlabel('E-E_F (eV)','FontSize',18); axis(ax);
ylabel('DOS (states/eV)','FontSize',18)

function plt_dosver(en,dos,lcol,lwidth,ax)
hold on; box on
if isnumeric(lcol); lcol = 'k'; end;
if size(dos,2) > 1
    plot(dos(:,1),en,lcol,'LineWidth',lwidth); a = -max(dos(:,2)); 
    plot(-dos(:,2),en,lcol,'LineWidth',lwidth); b = max(dos(:,1));
else plot(dos(:,1),en,lcol,'LineWidth',lwidth)
    a = 0; b = max(dos(:,1));
end
if ax(1) < a; a = ax(1); end; if ax(2) > b; b = ax(2); end
c = min(en(2:end)); if ax(3) < c; c = ax(3); end
d = max(en(2:end)); if ax(4) > d; d = ax(4); end
plot([a b],[0 0],'k--'); plot([0 0],[c d],'k') % vertical
xlabel('DOS (states/eV)','FontSize',18); axis(ax); 

function disp_no_Callback(~, ~, ~)
function disp_no_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_line_Callback(~, ~, ~)
function opt_line_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function disp_pathname_Callback(~, ~, ~)
function disp_pathname_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_no_Callback(~, ~, handles)
if get(handles.opt_no,'Value') == 1
    set(handles.disp_no,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.disp_no,'BackgroundColor',[0.8 0.8 0.8])
else set(handles.disp_no,'ForegroundColor',[0 0 0])
    set(handles.disp_no,'BackgroundColor',[1 1 1])
end

function opt_no_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_dos_Callback(~, ~, handles)
if get(handles.opt_dos,'Value') < 3 || get(handles.opt_dos,'Value') > 4
    set(handles.cek_pband,'Value',0); set(handles.cek_pband,'Foreground',[0 0 0])
end
function opt_dos_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function disp_axis_Callback(~, ~, ~)
function disp_axis_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cek_axis_Callback(~, ~, handles)
if get(handles.cek_axis,'Value') == 0
    set(handles.disp_axis,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.disp_axis,'BackgroundColor',[0.8 0.8 0.8])
else set(handles.disp_axis,'ForegroundColor',[0 0 0])
    set(handles.disp_axis,'BackgroundColor',[1 1 1])
end

function opt_marker_Callback(~, ~, ~)
function opt_marker_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_new_Callback(~, ~, handles)
set(handles.plot_new,'Value',1); set(handles.plot_new,'ForegroundColor',[0 0 1])
set(handles.plot_hold,'Value',0); set(handles.plot_hold,'ForegroundColor',[0 0 0])
set(handles.plot_clc,'Value',0); set(handles.plot_clc,'ForegroundColor',[0 0 0])

function plot_hold_Callback(~, ~, handles)
set(handles.plot_new,'Value',0); set(handles.plot_new,'ForegroundColor',[0 0 0])
set(handles.plot_hold,'Value',1); set(handles.plot_hold,'ForegroundColor',[0 0 1])
set(handles.plot_clc,'Value',0); set(handles.plot_clc,'ForegroundColor',[0 0 0])

function plot_clc_Callback(~, ~, handles)
set(handles.plot_new,'Value',0); set(handles.plot_new,'ForegroundColor',[0 0 0])
set(handles.plot_hold,'Value',0); set(handles.plot_hold,'ForegroundColor',[0 0 0])
set(handles.plot_clc,'Value',1); set(handles.plot_clc,'ForegroundColor',[0 0 1])

function disp_width_Callback(~, ~, ~)
function disp_width_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function disp_norm_Callback(~, ~, ~)
function disp_norm_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function disp_msize_Callback(~, ~, ~)
function disp_msize_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cek_normby_Callback(~, ~, handles)
set(handles.cek_norm,'Value',0); set(handles.cek_norm,'ForegroundColor',[0 0 0])
if get(handles.cek_normby,'Value') == 1
    set(handles.disp_norm,'BackgroundColor',[1 1 1])
    set(handles.disp_norm,'ForegroundColor',[0 0 0])
else set(handles.disp_norm,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.disp_norm,'ForegroundColor',[0.8 0.8 0.8])
end

function cek_norm_Callback(~, ~, handles)
set(handles.cek_normby,'Value',0)
set(handles.disp_norm,'BackgroundColor',[0.8 0.8 0.8])
set(handles.disp_norm,'ForegroundColor',[0.8 0.8 0.8])
set(handles.cek_normby,'ForegroundColor',[0 0 0])

function orb_s_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_px_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_py_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_pz_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_dxy_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_dyz_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_dxz_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_dz2_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_dx2_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_f_3_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_f_2_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_f_1_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_f0_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_f1_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_f2_Callback(~, ~, handles); set(handles.cek_all,'Value',0)
function orb_f3_Callback(~, ~, handles); set(handles.cek_all,'Value',0)

function cek_pband_Callback(~, ~, handles)
if get(handles.opt_dos,'Value') < 3 || get(handles.opt_dos,'Value') > 4
    set(handles.cek_pband,'Value',0)
end

function cek_all_Callback(~, ~, handles)
if get(handles.cek_all,'Value') == 1
    set(handles.orb_s,'Value',1); set(handles.orb_px,'Value',1)
    set(handles.orb_py,'Value',1); set(handles.orb_pz,'Value',1)
    set(handles.orb_dxy,'Value',1); set(handles.orb_dxz,'Value',1)
    set(handles.orb_dyz,'Value',1); set(handles.orb_dx2,'Value',1)
    set(handles.orb_dz2,'Value',1); set(handles.orb_f_3,'Value',1)
    set(handles.orb_f_2,'Value',1); set(handles.orb_f_1,'Value',1)
    set(handles.orb_f0,'Value',1); set(handles.orb_f1,'Value',1)
    set(handles.orb_f2,'Value',1); set(handles.orb_f3,'Value',1)
else set(handles.orb_s,'Value',0); set(handles.orb_px,'Value',0)
    set(handles.orb_py,'Value',0); set(handles.orb_pz,'Value',0)
    set(handles.orb_dxy,'Value',0); set(handles.orb_dxz,'Value',0)
    set(handles.orb_dyz,'Value',0); set(handles.orb_dx2,'Value',0)
    set(handles.orb_dz2,'Value',0); set(handles.orb_f_3,'Value',0)
    set(handles.orb_f_2,'Value',0); set(handles.orb_f_1,'Value',0)
    set(handles.orb_f0,'Value',0); set(handles.orb_f1,'Value',0)
    set(handles.orb_f2,'Value',0); set(handles.orb_f3,'Value',0)
end

function geometry = poscar(pathname)
fprintf('>> Reading POSCAR file \n')
fid = fopen([pathname 'POSCAR']);
fgetl(fid); scale = fscanf(fid, '%f',1);
geometry.lattice  = fscanf(fid, '%f %f %f', [3 3])'; 
geometry.lattice  = geometry.lattice*scale;
fgetl(fid); line  = fgetl(fid); cartesian = 0;
if sum(isstrprop(line, 'digit')) == 0;  line = fgetl(fid); end
geometry.atomcount = sscanf(line,'%d'); line = fgetl(fid);
if line(1) == 's' || line(1) == 'S';    line = fgetl(fid); end
if line(1) == 'C' || line(1) == 'c' ||  line(1) == 'K' || line(1) =='k'
	cartesian = 1;
end
for i = 1:sum(geometry.atomcount);
	line = fgetl(fid); geometry.coords(i,:) = sscanf(line, '%f %f %f');
end
if cartesian == 1
	geometry.coords = geometry.coords*scale;
    geometry.coords = geometry.coords/geometry.lattice;
end

function [energy, total_dos, efermi, pdos] = doscar(pathname)
fprintf('>> Reading DOSCAR file \n'); fid = fopen([pathname 'DOSCAR']);
buffer = fscanf(fid, '%f', 4)'; natoms = buffer(1); fgetl(fid); fgetl(fid); 
fgetl(fid); fgetl(fid); fgetl(fid); buffer = fscanf(fid, '%f', 5); 
energy = []; total_dos = []; efermi = []; pdos = []; 
if numel(buffer) ~=0
    nedos = buffer(3); efermi = buffer(4); fgetl(fid); position = ftell(fid);
    buffer = sscanf(fgetl(fid),'%f'); fseek(fid, position, 'bof');
    if max(size(buffer)) <= 5
        ispin = 1+(max(size(buffer))==5);
        total_dos = fscanf(fid, '%f', [(1+2*ispin) nedos])';
        energy = total_dos(:,1); total_dos = total_dos(:,2:(1+ispin));
        fgetl(fid); position = ftell(fid);
        if fgetl(fid)~=-1
            buffer = sscanf(fgetl(fid),'%f');
            sscanf(fgetl(fid),'%f'); sscanf(fgetl(fid),'%f');
            buffer2 = sscanf(fgetl(fid),'%f'); fseek(fid, position, 'bof');
            orb = max(size(buffer)); orb2 = max(size(buffer2)); 
            if orb == orb2; pdos = zeros([nedos orb-1 natoms]);
                for i = 1:natoms
                    if fgetl(fid)~=-1
                        buffer = fscanf(fid, '%f', [orb nedos])';
                        pdos(:,:,i) = buffer(:,2:end); fgetl(fid);
                    end
                end
            else pdos = zeros([nedos orb+orb2-1 natoms]);
                for i = 1:natoms
                    for j = 1:nedos
                        fgetl(fid); buffer = fscanf(fid, '%f', [orb 1])';
                        pdos(j,1:orb-1,i) = buffer(:,2:end);
                        buffer2 = fscanf(fid, '%f', [orb2 1])';
                        pdos(j,orb:end,i) = buffer2;
                    end
                end
            end
        end
    end
end
fclose(fid);

function [eigenvalues, kpoints, nelect] = eigenval(pathname)
fprintf('>> Reading EIGENVAL file \n')
fid = fopen([pathname 'EIGENVAL']); 
buffer = fscanf(fid, '%d', 4); ispin = buffer(4);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
nelect = fscanf(fid, '%d', 1); nkpoints = fscanf(fid, '%d', 1);
nbands = fscanf(fid, '%d', 1); kpoints = zeros(nkpoints,4);
eigenvalues = zeros(nkpoints, nbands, ispin); fgetl(fid);
for kpoint = 1:nkpoints
    fgetl(fid); kpoints(kpoint,:) = fscanf(fid, '%f', 4)'; fgetl(fid);
    for band = 1:nbands
        buffer = fscanf(fid,'%f',ispin+1); fgetl(fid);
        eigenvalues(kpoint,band,1:ispin) = buffer(2:(ispin+1));
    end
end
eigenvalues = sort(eigenvalues,2); fclose(fid);

function [pband1, pband2] = procar(pathname)
fprintf('>> Reading OUTCAR file \n'); fid = fopen([pathname 'OUTCAR']);
vasp = fscanf(fid, '%s', 1); vasp = str2double(vasp(6)); fclose(fid);
fid = fopen([pathname 'EIGENVAL']);
buffer = fscanf(fid, '%d', 4); ispin = buffer(4);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
buffer = fscanf(fid, '%d', 3); nkpoints = buffer(2); nbands = buffer(3);
fclose(fid); fid = fopen([pathname 'DOSCAR']);
buffer = fscanf(fid, '%f', 4)'; natoms = buffer(1); 
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); 
buffer = fscanf(fid, '%f', 5); fgetl(fid);
fscanf(fid, '%f', [(1+2*ispin) buffer(3)]); fgetl(fid); fgetl(fid);
buffer = sscanf(fgetl(fid),'%f'); orbital = (max(size(buffer))-1)/ispin;
fclose(fid); fprintf('>> Reading PROCAR file \n')
pband1 = zeros(nkpoints, nbands, natoms, orbital);
pband2 = zeros(nkpoints, nbands, natoms, orbital);
fid = fopen([pathname 'PROCAR']);
opt = fgetl(fid); opt = length(opt); fgetl(fid); fgetl(fid);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); 
for i = 1:nkpoints
    for j = 1:nbands
        buffer = fscanf(fid, '%f', [orbital+2 natoms])';
        pband1(i,j,:,:) = buffer(:,2:end-1);
        fgetl(fid); fgetl(fid); fgetl(fid); if vasp == 4; fgetl(fid); end
        if opt > 20; if vasp == 5 && orbital > 9; fgetl(fid); end;
            fscanf(fid, '%f', [2*natoms orbital+2]); if vasp == 4; fgetl(fid); end
        end
        fgetl(fid); fgetl(fid); 
        if vasp == 5; fgetl(fid); if orbital > 9 && opt < 28; fgetl(fid); end; end
    end
    fgetl(fid); 
    if vasp == 4 && isempty(fgetl(fid)); fgetl(fid); 
    elseif vasp == 5 && length(fgetl(fid)) == 1; fgetl(fid);
    end
end
if vasp == 4; fgetl(fid); end
if ispin == 2
    if length(fgetl(fid)) == 1; fgetl(fid); end
    for i = 1:nkpoints
        for j = 1:nbands
            buffer = fscanf(fid, '%f', [orbital+2 natoms])';
            pband2(i,j,:,:) = buffer(:,2:end-1);
            fgetl(fid); fgetl(fid); fgetl(fid); if vasp == 4; fgetl(fid); end
            if opt > 20; if vasp == 5 && orbital > 9; fgetl(fid); end;
                fscanf(fid, '%f', [2*natoms orbital+2]); if vasp == 4; fgetl(fid); end
            end 
            fgetl(fid); fgetl(fid);
            if vasp == 5; fgetl(fid); if orbital > 9 && opt < 28; fgetl(fid); end; end
        end
        fgetl(fid); 
        if vasp == 4 && isempty(fgetl(fid)); fgetl(fid); 
        elseif vasp == 5 && length(fgetl(fid)) == 1; fgetl(fid);
        end
    end
end
fclose(fid);

function pband = procar_lnon(pathname)
fprintf('>> Reading PROCAR file \n'); fid = fopen([pathname 'EIGENVAL']);
buffer = fscanf(fid, '%d', 4); ispin = buffer(4);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
buffer = fscanf(fid, '%d', 3); nkpoints = buffer(2); nbands = buffer(3);
fclose(fid); fid = fopen([pathname 'DOSCAR']);
buffer = fscanf(fid, '%f', 4)'; natoms = buffer(1);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
buffer = fscanf(fid, '%f', 5); fgetl(fid);
fscanf(fid, '%f', [(1+2*ispin) buffer(3)]); fgetl(fid); fgetl(fid);
buffer = sscanf(fgetl(fid),'%f'); sscanf(fgetl(fid),'%f'); 
sscanf(fgetl(fid),'%f'); buffer2 = sscanf(fgetl(fid),'%f');
orb = max(size(buffer)); orb2 = max(size(buffer2));
if orb == orb2; orbital = 9; else orbital = 16; end; fclose(fid);
pband = zeros(nkpoints, nbands, natoms, orbital); fid = fopen([pathname 'PROCAR']);
opt = fgetl(fid); opt = length(opt); fgetl(fid); fgetl(fid);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
for i = 1:nkpoints
    for j = 1:nbands
        buffer = fscanf(fid, '%f', [orbital+2 natoms])';
        pband(i,j,:,:) = buffer(:,2:end-1); fgetl(fid); fgetl(fid); 
        fscanf(fid, '%f', [orbital+2 natoms]); fgetl(fid); fgetl(fid); 
        fscanf(fid, '%f', [orbital+2 natoms]); fgetl(fid); fgetl(fid);
        fscanf(fid, '%f', [orbital+2 natoms]); fgetl(fid); fgetl(fid); 
        if opt > 20; fgetl(fid); for ion = 1:2*natoms; fgetl(fid); end; end
        fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); 
        if orbital > 9; fgetl(fid); end 
    end
    fgetl(fid); if length(fgetl(fid)) == 1; fgetl(fid); end
end
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
