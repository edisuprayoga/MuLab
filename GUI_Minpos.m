function varargout = GUI_Minpos(varargin)

% GUI_Minpos MATLAB code by Edi Suprayoga.
% feel free to contact me at suprayoga.edi@gmail.com
% last edit 3 Mar 2016

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Minpos_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Minpos_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout; [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else gui_mainfcn(gui_State, varargin{:})
end

function GUI_Minpos_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject; guidata(hObject, handles)
set(handles.figure1,'Color',[0.8 0.8 0.8])

function varargout = GUI_Minpos_OutputFcn(~, ~, handles) 
varargout{1} = handles.output; set(handles.autorized,'Value',1)

function Run_Callback(~, ~, handles); clc; tic
fprintf('   ============================================ \n')
fprintf('>> Reading input files \n')
Em = str2num(get(handles.Iso,'String'));
No = str2num(get(handles.no_atom,'String'));
filename = get(handles.disp_path,'String');
if get(handles.opt_loc,'Value') >= 3; i = 1; pot = []; pos = [];
    fprintf('\n>> Searching for minimum energy position \n')
    if get(handles.opt_loc,'Value') == 3 % folder
        fid = fopen([filename num2str(i) '/POSCAR']);
        while fid ~= -1
            fclose(fid); geo = poscar([filename num2str(i) '/POSCAR']);
            pot = [pot; outcar([filename num2str(i) '/OUTCAR'])];
            pos = [pos; geo.coords(end,:)]; i = i+1;
            fid = fopen([filename num2str(i) '/POSCAR']);
        end
        [en,j] = min(pot(:)); geoM = poscar([filename num2str(j) '/POSCAR']);
    elseif get(handles.opt_loc,'Value') == 4 % files
        fid = fopen([filename 'POSCAR/POSCAR_' num2str(i)]);
        while fid ~= -1
            fclose(fid); geo = poscar([filename 'POSCAR/POSCAR_' num2str(i)]);
            pot = [pot; outcar([filename 'OUTCAR/OUTCAR_' num2str(i)])];
            pos = [pos; geo.coords(end,:)];
            i = i+1; fid = fopen([filename 'POSCAR/POSCAR_' num2str(i)]);
        end
        [en,j] = min(pot(:)); geoM = poscar([filename 'POSCAR/POSCAR_' num2str(j)]);
    end
    pot = reshape(pot,Em); xm = reshape(pos(:,1),Em);
    ym = reshape(pos(:,2),Em); zm = reshape(pos(:,3),Em);
    posM = geoM.coords(end,:); latt = geo.lattice;
    pos0 = [xm(round(Em(1)/2),round(Em(2)/2),round(Em(3)/2)), ...
        ym(round(Em(1)/2),round(Em(2)/2),round(Em(3)/2)), ...
        zm(round(Em(1)/2),round(Em(2)/2),round(Em(3)/2))];
    if get(handles.corin,'Value') == 2;
        pos0 = pos0(1)*latt(1,:)+pos0(2)*latt(2,:)+pos0(3)*latt(3,:);
    end
    if get(handles.corout,'Value') == 2;
        posM = posM(1)*latt(1,:)+posM(2)*latt(2,:)+posM(3)*latt(3,:);
    end
    fprintf('     initial minimum position = %4.4f %4.4f %4.4f \n',pos0)
    fprintf('     minimum energy position  = %4.4f %4.4f %4.4f \n',posM)
    fprintf('     minimum energy in POSCAR_%1.0f \n',j)
    set(handles.Xinput,'String',num2str(pos0));
    set(handles.Xoutput,'String',num2str(posM));
    if get(handles.check_en,'Value') == 0
        set(handles.emin,'String',en-pot(round(Em(1)/2),round(Em(2)/2),round(Em(3)/2)))
    else set(handles.emin,'String',en)
    end
    if get(handles.check_dis,'Value') == 1
        set(handles.Emin,'String',pot(round(Em(1)/2),round(Em(2)/2),round(Em(3)/2)))
    end
    
elseif get(handles.opt_loc,'Value') <= 2
    fid = fopen(filename); geo = poscar(fid);
    latt = geo.lattice; latti = inv(latt);
    fgetl(fid); gridsize = fscanf(fid, '%d %d %d', [3 1])';
    locpot = fscanf(fid, '%f', [prod(gridsize,2) 1])';
    locpot = reshape(locpot,gridsize); fclose(fid);
    locpot(end+1,:,:) = locpot(1,:,:); locpot(:,end+1,:) = locpot(:,1,:); 
    locpot(:,:,end+1) = locpot(:,:,1);
    if get(handles.opt_loc,'Value') == 1
        fprintf('\n>> Searching for minimum potential position \n')
        loc = -locpot;
    elseif get(handles.opt_loc,'Value') == 2
        fprintf('\n>> Searching for Iso-fields position \n')
        loc = abs(locpot-Em);
    end
    Rd = str2double(get(handles.Rd,'String'));
    pos1 = str2num(get(handles.Xinput,'String'));
    x1 = pos1(1); x2 = pos1(2); x3 = pos1(3);
    fprintf('     predicted position = %4.4f %4.4f %4.4f \n',x1,x2,x3)
    if (get(handles.corin,'Value') == 2)
        pos = x1*latti(1,:)+x2*latti(2,:)+x3*latti(3,:);
        x1 = pos(1); x2 = pos(2); x3 = pos(3);
    end
    if get(handles.lock,'Value') == 1
        a1 = [x1 x2 x3];
        x = linspace(0,1,size(locpot,1)); a = abs(x-x1); [~,a] = min(a);
        y = linspace(0,1,size(locpot,2)); b = abs(y-x2); [~,b] = min(b);
        z = linspace(0,1,size(locpot,3)); c = abs(z-x3); [~,c] = min(c);
        en1 = loc(a,b,c);
    elseif x1 == 0 && x2 == 0 && x3 == 0
        [~,pos] = min(loc(:)); [i,j,k] = ind2sub(size(loc),pos);
        x = linspace(0,1,size(locpot,1));
        y = linspace(0,1,size(locpot,2)); 
        z = linspace(0,1,size(locpot,3));
        a1 = [x(i) y(j) z(k)]; en1 = loc(i,j,k);
    else [a1, en1] = minpos(geo,[x1 x2 x3],loc,Rd);
    end
    if get(handles.opt_loc,'Value') == 2
        x = linspace(0,1,size(locpot,1)); a = abs(x-a1(1)); [~,a] = min(a);
        y = linspace(0,1,size(locpot,2)); b = abs(y-a1(2)); [~,b] = min(b);
        z = linspace(0,1,size(locpot,3)); c = abs(z-a1(3)); [~,c] = min(c);
        en1 = locpot(a,b,c);
    end
    p2 = a1; p1 = a1(1)*latt(1,:)+a1(2)*latt(2,:)+a1(3)*latt(3,:);
    if (get(handles.corout,'Value') == 2); a1 = p1; end
    if get(handles.opt_loc,'Value') == 1
        fprintf('     minimum position   = %4.4f %4.4f %4.4f \n',a1)
        set(handles.emin,'String',en1)
        if get(handles.check_en,'Value') == 1; set(handles.emin,'String',en1)
        else set(handles.emin,'String',(en1-min(loc(:)))*1000)
        end
    elseif get(handles.opt_loc,'Value') == 2 && get(handles.check_en,'Value') == 0
        fprintf('     isofields position = %4.4f %4.4f %4.4f \n',a1)
        set(handles.emin,'String',en1)
    elseif get(handles.opt_loc,'Value') == 2 && get(handles.check_en,'Value') == 1
        fprintf('     isofields position = %4.4f %4.4f %4.4f \n',a1)
        fprintf('\n>> Calculate Field Distribution \n')
        x = linspace(0,1,size(locpot,1)); y = linspace(0,1,size(locpot,2)+1);
        z = linspace(0,1,size(locpot,3));
        pos = str2num(get(handles.dis_rad,'String'))*(latti(1,:)+latti(2,:)+latti(3,:));
        xa = pos(1); ya = pos(2); za = pos(3); if (get(handles.corout,'Value') == 2); a1 = p2; end
        pa = abs(x-a1(1)-xa); [~,pa] = min(pa); pb = abs(y-a1(2)-ya); [~,pb] = min(pb);
        pc = abs(z-a1(3)-za); [~,pc] = min(pc); qa = abs(x-a1(1)+xa); [~,qa] = min(qa);
        qb = abs(y-a1(2)+ya); [~,qb] = min(qb); qc = abs(z-a1(3)+za); [~,qc] = min(qc);
        en1 = locpot(qa:pa,qb:pb,qc:pc); 
        locT = en1(:); 
        if (get(handles.corout,'Value') == 2); a1 = p1; end
        Chg = linspace(str2double(get(handles.dis_min,'String')), ...
                str2double(get(handles.dis_max,'String')), ...
                str2double(get(handles.dis_N,'String')));
        Dist = zeros(1,length(locT));
        for i = 1:length(Chg)-1;
            n = (locT > Chg(i) & locT <= Chg(i+1));
            Dist(i+1) = sum(n(:));
        end
        Dist = Chg(Dist == max(Dist)); set(handles.emin,'String',mean(Dist))
    end
    if get(handles.check_dis,'Value') == 1
        if get(handles.opt_loc,'Value') == 1; set(handles.Emin,'String',min(loc(:)))
        elseif get(handles.opt_loc,'Value') == 2
            locT = locpot(:); Dist = zeros(1,length(locT));
            Chg = linspace(str2num(get(handles.dis_min,'String')), ...
                str2num(get(handles.dis_max,'String')),str2num(get(handles.dis_N,'String')));
            for i = 1:length(Chg)-1;
                n = (locT > Chg(i) & locT <= Chg(i+1));
                Dist(i+1) = sum(n(:));
            end
            Dist = Chg(Dist == max(Dist)); set(handles.Emin,'String',mean(Dist))
        end
    end
    set(handles.Xoutput,'String',num2str(a1))
end
fprintf('\n>> Drawing figure \n')
xb = str2num(get(handles.xb,'String'));
yb = str2num(get(handles.yb,'String')); zb = str2num(get(handles.zb,'String'));
if xb(1)<0; xb(1)=0; end; if xb(2)>1; xb(2)=1; end; if xb(2)<=xb(1); xa=xb; xb=[xa(2) xa(1)]; end
if yb(1)<0; yb(1)=0; end; if yb(2)>1; yb(2)=1; end; if yb(2)<=yb(1); ya=yb; yb=[ya(2) ya(1)]; end
if zb(1)<0; zb(1)=0; end; if zb(2)>1; zb(2)=1; end; if zb(2)<=zb(1); za=zb; xb=[za(2) za(1)]; end
if get(handles.cek_cls,'Value') == 1; figure(1), hold on; else figure; end
if get(handles.cek_cls,'Value') == 0; 
    col = get(handles.ion_col,'String'); rd = str2num(get(handles.ion_rad,'String'));
    if length(No) == 1; col_oct = col(1); else col_oct = col(2);end
    if get(handles.opt_bond,'Value') > 1; d0 = str2num(get(handles.bond_l,'String')); end
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
end
if get(handles.opt_loc,'Value') <= 2;
    xa = linspace(0,1,size(locpot,1)); 
    ya = linspace(0,1,size(locpot,2)); 
    za = linspace(0,1,size(locpot,3)); 
    [x,y,z] = ndgrid(xa,ya,za);
    x = x(:); y = y(:); z = z(:);
    xi = x*latt(1,1) + y*latt(2,1) + z*latt(3,1);
    yi = x*latt(1,2) + y*latt(2,2) + z*latt(3,2);
    zi = x*latt(1,3) + y*latt(2,3) + z*latt(3,3);
    x = reshape(xi,size(locpot)); y = reshape(yi,size(locpot)); z = reshape(zi,size(locpot));
    for i = 1:2
        a = abs(xa-xb(i)); [~,xb(i)] = min(a);
        a = abs(ya-yb(i)); [~,yb(i)] = min(a);
        a = abs(za-zb(i)); [~,zb(i)] = min(a);
    end
    x = x(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
    y = y(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
    z = z(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
    
    if get(handles.opt_loc,'Value') == 1
        loc = loc(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
        p = patch(isosurface(x,y,z,loc,min(loc(:))+Em*1E-3));
        set(p,'FaceColor','y','EdgeColor','none','facelighting','gouraud','FaceAlpha',0.5)
        draw_atom(p1,0.2,'r',1)
    elseif get(handles.opt_loc,'Value') == 2
        locpot = locpot(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
        if get(handles.points,'Value') == 1 % point
            x = x(:); y = y(:); z = z(:); locpot = locpot(:);
            cond = (locpot < Em-10); x(cond)=[]; y(cond)=[]; z(cond)=[]; locpot(cond) = [];
            cond = (locpot > Em+10); x(cond)=[]; y(cond)=[]; z(cond)=[];        
%             Cu = [0.5 0.5 0]; Cu = Cu(1)*latt(1,:)+Cu(2)*latt(2,:)+Cu(3)*latt(3,:);
%             R = sqrt((x-Cu(1)).^2+(y-Cu(2)).^2+(z-Cu(3)).^2);
%             cond = (R > 2.5); x(cond)=[]; y(cond)=[]; z(cond)=[];
            draw_atom([x y z],0.05,'r',1);
        else p = patch(isosurface(x,y,z,locpot,Em));
            set(p,'FaceColor','c','EdgeColor','none','facelighting','gouraud','FaceAlpha',0.5)
            draw_atom(p1,0.2,'r',1)
        end
    end
elseif get(handles.opt_loc,'Value') >= 3
    if get(handles.corin,'Value') == 1;
        pos0 = pos0(1)*latt(1,:)+pos0(2)*latt(2,:)+pos0(3)*latt(3,:);
    end
    if get(handles.corout,'Value') == 1;
        posM = posM(1)*latt(1,:)+posM(2)*latt(2,:)+posM(3)*latt(3,:);
    end
    draw_atom(pos0,0.3,'r',0.2); draw_atom(posM,0.3,'r',1)
    if get(handles.check_dis,'Value') == 1 && length(size(pot)) == 3;
        xm = xm(:); ym = ym(:); zm = zm(:);
        xi = xm*latt(1,1) + ym*latt(2,1) + zm*latt(3,1);
        yi = xm*latt(1,2) + ym*latt(2,2) + zm*latt(3,2);
        zi = xm*latt(1,3) + ym*latt(2,3) + zm*latt(3,3);
        xm = reshape(xi,size(pot)); ym = reshape(yi,size(pot));
        zm = reshape(zi,size(pot)); p = patch(isosurface(xm,ym,zm,pot));
        set(p,'FaceColor','y','EdgeColor','none','facelighting','gouraud','FaceAlpha',0.5)
    end
end
if get(handles.box,'Value') == 1; drbox(geo.lattice); end
daspect([1 1 1]); camlight right local; rotate3d on;
if get(handles.fix,'Value') == 1
    xlabel('a','FontSize',12); ylabel('b','FontSize',12); 
    zlabel('c','FontSize',12); axis tight
else axis off
end
t = toc; h = floor(t/3600); m = floor((t-h*3600)/60); t = t-h*3600 - m*60;
fprintf('   ============================================ \n');
if h >= 1; fprintf('   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n',h,m,t);
elseif m >= 1; fprintf('   Elapsed time is %1.0f min %1.4f sec.\n',m,t);
else fprintf('   '); toc
end

function Iso_Callback(~, ~, ~)
function Iso_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Xoutput_Callback(~, ~, ~)
function Xoutput_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Xinput_Callback(~, ~, ~)
function Xinput_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load_Callback(~, ~, handles)
global filename_loc pathname_loc
[filename_loc, pathname_loc] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname_loc)
    set(handles.disp_path,'String',[pathname_loc filename_loc])
end

function emin_Callback(~, ~, ~)
function emin_CreateFcn(hObject, ~, ~)
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
geometry.lattice = geometry.lattice*scale; cartesian = 0;
fgetl(fid); line = fgetl(fid); has_symbols_line = false;
if sum(isstrprop(line, 'digit')) == 0
    geometry.symbols = regexp(line, '([^ ]*)', 'match');
    line = fgetl(fid); has_symbols_line = true;
else geometry.symbols = {};
end
geometry.atomcount = sscanf(line,'%d'); natoms = sum(geometry.atomcount);
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

function result = outcar(filename)
fid = fopen(filename);
while ~feof(fid)
    buffer = ftell(fid); line = fgetl(fid);
    if numel(regexp(line,'energy\(sigma->0\)'))==1; pos = buffer; end
end; fseek(fid,pos,'bof'); line = fgetl(fid);
result = sscanf(line,'  energy  without entropy= %*f energy(sigma->0) = %f');
fclose(fid);

function [mpos, en] = minpos(geo,mpos,loc,Rd)
latt = geo.lattice;
mpos = mpos(1)*latt(1,:)+mpos(2)*latt(2,:)+mpos(3)*latt(3,:);
x = linspace(0,1,size(loc,1)); y = linspace(0,1,size(loc,2));
z = linspace(0,1,size(loc,3)); [x,y,z] = ndgrid(x,y,z);
x = x(:); y = y(:); z = z(:); loc = loc(:);
xr = x.*latt(1,1) + y*latt(2,1) + z*latt(3,1);
yr = x.*latt(1,2) + y*latt(2,2) + z*latt(3,2);
zr = x.*latt(1,3) + y*latt(2,3) + z*latt(3,3);
cond = (xr >= mpos(1)-Rd & xr <= mpos(1)+Rd);
xr = xr(cond); yr = yr(cond); zr = zr(cond); loc = loc(cond);
cond = (yr >= mpos(2)-Rd & yr <= mpos(2)+Rd);
xr = xr(cond); yr = yr(cond); zr = zr(cond); loc = loc(cond);
cond = (zr >= mpos(3)-Rd & zr <= mpos(3)+Rd);
xr = xr(cond); yr = yr(cond); zr = zr(cond); loc = loc(cond);
[~,pos] = min(loc(:));
mpos = [xr(pos) yr(pos) zr(pos)]; en = loc(pos);
latt = inv(latt);
mpos = mpos(1)*latt(1,:)+mpos(2)*latt(2,:)+mpos(3)*latt(3,:);

% function [mpos, en] = minpos(geo,mpos,loc,Rd)
% latt = geo.lattice;
% mpos = mpos(1)*latt(1,:)+mpos(2)*latt(2,:)+mpos(3)*latt(3,:);
% latt = latt(1,:)+latt(2,:)+latt(3,:); 
% x = linspace(0,latt(1),size(loc,1));
% y = linspace(0,latt(2),size(loc,2));
% z = linspace(0,latt(3),size(loc,3));
% xa = abs(x-mpos(1)+Rd); [~,xa] = min(xa); ya = abs(y-mpos(2)+Rd); [~,ya] = min(ya);
% za = abs(z-mpos(3)+Rd); [~,za] = min(za); xb = abs(x-mpos(1)-Rd); [~,xb] = min(xb);
% yb = abs(y-mpos(2)-Rd); [~,yb] = min(yb); zb = abs(z-mpos(3)-Rd); [~,zb] = min(zb);
% V = loc(xa:xb, ya:yb, za:zb); [~,pos] = min(V(:)); [i,j,k] = ind2sub(size(V),pos);
% gpos = [xa ya za]+[i j k]-1; 
% x = linspace(0,1,size(loc,1));
% y = linspace(0,1,size(loc,2));
% z = linspace(0,1,size(loc,3));
% mpos = [x(gpos(1)) y(gpos(2)) z(gpos(3))]; en = loc(gpos(1),gpos(2),gpos(3));

function drbox(latt)
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

function corin_Callback(~, ~, ~)
function corin_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function corout_Callback(~, ~, ~)
function corout_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function no_atom_Callback(~, ~, ~)
function no_atom_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function Emin_Callback(~, ~, ~)
function Emin_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function draw_atom(target,r,col,alp); hold on
[rx,ry,rz] = sphere(10); rx = rx*r; ry = ry*r; rz = rz*r;
for i = 1:size(target,1)
    h = surf((rx+target(i,1)),(ry+target(i,2)),(rz+target(i,3)));
    set(h,'edgecolor','none','facecolor',col,'facelighting','gouraud')
    alpha(h,alp)
end

function opt_loc_Callback(~, ~, handles)
set(handles.dis_rad,'ForegroundColor',[0.8 0.8 0.8])
set(handles.dis_rad,'BackgroundColor',[0.8 0.8 0.8])
set(handles.dis_min,'ForegroundColor',[0.8 0.8 0.8])
set(handles.dis_min,'BackgroundColor',[0.8 0.8 0.8])
set(handles.dis_max,'ForegroundColor',[0.8 0.8 0.8])
set(handles.dis_max,'BackgroundColor',[0.8 0.8 0.8])
set(handles.dis_N,'ForegroundColor',[0.8 0.8 0.8])
set(handles.dis_N,'BackgroundColor',[0.8 0.8 0.8])
if get(handles.opt_loc,'Value') <= 2
    set(handles.text5,'String','Input Pos: '); set(handles.text8,'String','Output Pos: ')
    set(handles.text5,'ForegroundColor',[0 0 0]); set(handles.Xinput,'ForegroundColor',[0 0 0]);
    set(handles.Rd,'ForegroundColor',[0 0 0]); set(handles.Rd,'BackgroundColor',[1 1 1])
end
if get(handles.opt_loc,'Value') == 1
    set(handles.uipanel4,'ForegroundColor',0.25*[1 1 1])
    if get(handles.check_en,'Value') == 1
        set(handles.check_en,'String','Energy (eV)')
    else set(handles.check_en,'String','dE (meV)')
    end
    set(handles.text3,'String','Isosurf (meV): ')
    set(handles.check_dis,'String','Min Energy')
    set(handles.load,'String','LOCPOT')
elseif get(handles.opt_loc,'Value') == 2
    set(handles.uipanel4,'ForegroundColor',[0 0 1])
    set(handles.text3,'String','Isosurf (G): ')
    if get(handles.check_en,'Value') == 1
        set(handles.check_en,'String','Dist (G)')
        set(handles.dis_rad,'ForegroundColor',[0 0 0])
        set(handles.dis_rad,'BackgroundColor',[1 1 1])
    else set(handles.check_en,'String','Fields (G)')
    end
    if get(handles.check_dis,'Value') == 1 || get(handles.check_en,'Value') == 1
        set(handles.dis_min,'ForegroundColor',[0 0 0])
        set(handles.dis_min,'BackgroundColor',[1 1 1])
        set(handles.dis_max,'ForegroundColor',[0 0 0])
        set(handles.dis_max,'BackgroundColor',[1 1 1])
        set(handles.dis_N,'ForegroundColor',[0 0 0])
        set(handles.dis_N,'BackgroundColor',[1 1 1])
    end
    set(handles.check_dis,'String','Dist Max'); set(handles.load,'String','DIPOLE')
else set(handles.uipanel4,'ForegroundColor',0.31*[1 1 1])
    set(handles.text3,'String','Npoints [x y z]: ')
    set(handles.check_dis,'String','E0 (eV)')
    if get(handles.check_en,'Value') == 1
        set(handles.check_en,'String','Min Energy')
    else set(handles.check_en,'String','dE (meV)')
    end
    if get(handles.opt_loc,'Value') == 3
        set(handles.load,'String','Folder_')
    else set(handles.load,'String','Filename_')
    end
    set(handles.lock,'ForegroundColor',[0 0 0]); set(handles.lock,'Value',0)
    set(handles.lock,'FontWeight','normal'); set(handles.Xinput,'ForegroundColor',[0 0 0])
    set(handles.text5,'String','Initial Pos: '); set(handles.text8,'String','Min Pos: ')
    set(handles.text5,'ForegroundColor',[0 0 1]); set(handles.Xinput,'ForegroundColor',[0 0 1]); 
    set(handles.Rd,'ForegroundColor',[0.8 0.8 0.8]); 
    set(handles.Rd,'BackgroundColor',[0.8 0.8 0.8])
end

function opt_loc_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function check_dis_Callback(~, ~, handles)
if get(handles.check_dis,'Value') == 1
    set(handles.Emin,'ForegroundColor',[0 0 0])
    set(handles.Emin,'BackgroundColor',[1 1 1])
else set(handles.Emin,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.Emin,'BackgroundColor',[0.8 0.8 0.8])
end
if get(handles.opt_loc,'Value') == 2
    if get(handles.check_dis,'Value') == 1 || get(handles.check_en,'Value') == 1
        set(handles.dis_min,'ForegroundColor',[0 0 0])
        set(handles.dis_min,'BackgroundColor',[1 1 1])
        set(handles.dis_max,'ForegroundColor',[0 0 0])
        set(handles.dis_max,'BackgroundColor',[1 1 1])
        set(handles.dis_N,'ForegroundColor',[0 0 0])
        set(handles.dis_N,'BackgroundColor',[1 1 1])
    else set(handles.dis_N,'BackgroundColor',[0.8 0.8 0.8])
        set(handles.dis_N,'ForegroundColor',[0.8 0.8 0.8])
        set(handles.dis_min,'ForegroundColor',[0.8 0.8 0.8])
        set(handles.dis_min,'BackgroundColor',[0.8 0.8 0.8])
        set(handles.dis_max,'ForegroundColor',[0.8 0.8 0.8])
        set(handles.dis_max,'BackgroundColor',[0.8 0.8 0.8])
    end
end

function disp_path_Callback(~, ~, ~)
function disp_path_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function dis_min_Callback(~, ~, ~)
function dis_min_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function dis_max_Callback(~, ~, ~)
function dis_max_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function dis_N_Callback(~, ~, ~)
function dis_N_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function dis_rad_Callback(~, ~, ~)
function dis_rad_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function check_en_Callback(~, ~, handles)
if get(handles.opt_loc,'Value') == 1
    if get(handles.check_en,'Value') == 1
        set(handles.check_en,'String','Energy (eV)')
    else set(handles.check_en,'String','dE (meV)')
    end
elseif get(handles.opt_loc,'Value') == 2
    if get(handles.check_en,'Value') == 1
        set(handles.check_en,'String','Dist (G)')
        set(handles.dis_rad,'ForegroundColor',[0 0 0])
        set(handles.dis_rad,'BackgroundColor',[1 1 1])
    else set(handles.check_en,'String','Fields (G)')
        set(handles.dis_rad,'ForegroundColor',[0.8 0.8 0.8])
        set(handles.dis_rad,'BackgroundColor',[0.8 0.8 0.8])
        set(handles.dis_min,'ForegroundColor',[0.8 0.8 0.8])
        set(handles.dis_min,'BackgroundColor',[0.8 0.8 0.8])
        set(handles.dis_max,'ForegroundColor',[0.8 0.8 0.8])
        set(handles.dis_max,'BackgroundColor',[0.8 0.8 0.8])
        set(handles.dis_N,'ForegroundColor',[0.8 0.8 0.8])
        set(handles.dis_N,'BackgroundColor',[0.8 0.8 0.8])
    end
    if get(handles.check_dis,'Value') == 1 || get(handles.check_en,'Value') == 1
        set(handles.dis_min,'ForegroundColor',[0 0 0])
        set(handles.dis_min,'BackgroundColor',[1 1 1])
        set(handles.dis_max,'ForegroundColor',[0 0 0])
        set(handles.dis_max,'BackgroundColor',[1 1 1])
        set(handles.dis_N,'ForegroundColor',[0 0 0])
        set(handles.dis_N,'BackgroundColor',[1 1 1])
    end
elseif get(handles.check_en,'Value') == 1
    set(handles.check_en,'String','Energy (eV)')
else set(handles.check_en,'String','dE (meV)')
end

function cell_Callback(~, ~, ~)
function cell_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function bond_l_Callback(~, ~, ~)
function bond_l_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
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

function opt_bond_Callback(~, ~, handles)
if get(handles.opt_bond,'Value') == 1
    set(handles.bond_l,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.bond_l,'BackgroundColor',[0.8 0.8 0.8])
else set(handles.bond_l,'ForegroundColor',[0 0 0])
    set(handles.bond_l,'BackgroundColor',[1 1 1])
end

function opt_bond_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function xb_Callback(~, ~, ~)
function xb_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function yb_Callback(~, ~, ~)
function yb_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function zb_Callback(~, ~, ~)
function zb_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function box_Callback(~, ~, ~)
function fix_Callback(~, ~, ~)
function cek_cls_Callback(~, ~, ~)

function lock_Callback(~, ~, handles)
if get(handles.opt_loc,'Value') <= 2
    if get(handles.lock,'Value') == 1; set(handles.lock,'ForegroundColor',[1 0 0])
        set(handles.lock,'FontWeight','bold'); set(handles.Xinput,'ForegroundColor',[1 0 0])
        set(handles.Rd,'ForegroundColor',[0.8 0.8 0.8]); 
        set(handles.Rd,'BackgroundColor',[0.8 0.8 0.8])
    else set(handles.lock,'ForegroundColor',[0 0 0])
        set(handles.lock,'FontWeight','normal'); set(handles.Xinput,'ForegroundColor',[0 0 0])
        set(handles.Rd,'ForegroundColor',[0 0 0]); set(handles.Rd,'BackgroundColor',[1 1 1])
    end
else set(handles.lock,'Value',0)
end

function ion_rad_Callback(~, ~, ~)
function ion_rad_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function ion_col_Callback(~, ~, ~)
function ion_col_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white')
end

function points_Callback(~, ~, handles)
if get(handles.points,'Value') == 1;
    set(handles.points,'ForegroundColor',[1 0 0])
    set(handles.points,'FontWeight','bold');
else set(handles.points,'ForegroundColor',[0 0 0])
    set(handles.points,'FontWeight','normal');
end

function Rd_Callback(~, ~, ~)
function Rd_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
