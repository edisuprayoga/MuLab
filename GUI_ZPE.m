function varargout = GUI_ZPE(varargin)

% GUI_ZPE MATLAB code by Edi Suprayoga
% feel free to contact me at suprayoga.edi@gmail.com
% last edit 25 July 2016

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ZPE_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ZPE_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout; [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else gui_mainfcn(gui_State, varargin{:});
end

function GUI_ZPE_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject; guidata(hObject, handles);
set(handles.figure1,'Color',[0.8 0.8 0.8])

function varargout = GUI_ZPE_OutputFcn(~, ~, handles) 
varargout{1} = handles.output; set(handles.autorized,'Value',1);

function Result_Callback(~, ~, ~)
function Result_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Run_Callback(~, ~, handles)
filename = get(handles.pathname,'String'); clc; tic
fprintf('   ============================================ \n')
fprintf('           Zero Point Energy Calculation \n')
fprintf('   ============================================ \n')
fprintf('>> Reading input files \n')
fid = fopen(filename); geo = poscar(fid);
latt = geo.lattice; fgetl(fid);
gridsize = fscanf(fid, '%d %d %d', [3 1])';
loc = fscanf(fid, '%f', [prod(gridsize,2) 1])';
loc = reshape(loc,gridsize); fclose(fid);
muon = str2num(get(handles.muon,'String')); n = str2double(get(handles.n,'String'));
Nx = str2double(get(handles.Nx,'String')); dx = str2double(get(handles.dx,'String'));
Ny = str2double(get(handles.Ny,'String')); dy = str2double(get(handles.dy,'String'));
Nz = str2double(get(handles.Nz,'String')); dz = str2double(get(handles.dz,'String'));
if get(handles.opt_grid,'Value') == 2;
    Nx = ceil(dx/Nx); Ny = ceil(dy/Ny); Nz = ceil(dz/Nz);
end
fprintf('     Hamiltonian matrix size: %1.0f x %1.0f \n',Nx*Ny*Nz,Nx*Ny*Nz)
fprintf('\n>> Writing Hamiltonian matrix \n')
if get(handles.opt_muon,'Value') == 1; m = 1.88353148e-28; % muon mass
else m = 9*1.88353148e-28; end % proton mass
if get(handles.opt_axis,'Value') == 1;
muon = muon(1)*latt(1,:) + muon(2)*latt(2,:) + muon(3)*latt(3,:); end
V = -loc; SX = []; NX = size(loc);
for i = -1:1; SX = [SX; loc]; end; loc = SX; SX = [];
for i = -1:1; SX = [SX loc]; end; loc = SX; SX = []; 
for i = -1:1
    if size(SX,3) == 1; SX(:,:,size(SX,3):size(SX,3)+NX(3)-1) = loc;
    else SX(:,:,size(SX,3)+1:size(SX,3)+NX(3)) = loc; end
end
loc = SX; clear SX; latt = inv(latt); 
loc(end+1,:,:) = loc(1,:,:); loc(:,end+1,:) = loc(:,1,:); loc(:,:,end+1) = loc(:,:,1);
muon = muon(1)*latt(1,:) + muon(2)*latt(2,:) + muon(3)*latt(3,:);
dxyz = (dx*latt(1,:) + dy*latt(2,:) + dz*latt(3,:))/2;
x = linspace(-1,2,size(loc,1)); 
y = linspace(-1,2,size(loc,2)); 
z = linspace(-1,2,size(loc,3));
xa = abs(x-muon(1)+dxyz(1)); [~,xa] = min(xa);
ya = abs(y-muon(2)+dxyz(2)); [~,ya] = min(ya);
za = abs(z-muon(3)+dxyz(3)); [~,za] = min(za);
xb = abs(x-muon(1)-dxyz(1)); [~,xb] = min(xb);
yb = abs(y-muon(2)-dxyz(2)); [~,yb] = min(yb);
zb = abs(z-muon(3)-dxyz(3)); [~,zb] = min(zb);
xi1 = linspace(x(xb),x(xa),Nx); yi1 = linspace(y(yb),y(ya),Ny);
zi1 = linspace(z(zb),z(za),Nz); [xi,yi,zi] = ndgrid(xi1,yi1,zi1); 
xi = xi(:); yi = yi(:); zi = zi(:); latt = geo.lattice;
xi1 = xi*latt(1,1) + yi*latt(2,1) + zi*latt(3,1);
yi1 = xi*latt(1,2) + yi*latt(2,2) + zi*latt(3,2);
zi1 = xi*latt(1,3) + yi*latt(2,3) + zi*latt(3,3);
xi = reshape(xi1,[Nx Ny Nz]); yi = reshape(yi1,[Nx Ny Nz]);
zi = reshape(zi1,[Nx Ny Nz]); [x,y,z] = ndgrid(x,y,z);
x = x(:); y = y(:); z = z(:);
xi1 = x*latt(1,1) + y*latt(2,1) + z*latt(3,1);
yi1 = x*latt(1,2) + y*latt(2,2) + z*latt(3,2);
zi1 = x*latt(1,3) + y*latt(2,3) + z*latt(3,3);
x = reshape(xi1,size(loc)); xi1 = [];  y = reshape(yi1,size(loc)); yi1 = [];
z = reshape(zi1,size(loc)); zi1 = []; xo = x(xa:xb, ya:yb, za:zb); x = [];
yo = y(xa:xb, ya:yb, za:zb);  y = []; zo = z(xa:xb, ya:yb, za:zb); z = [];
Vi = -loc(xa:xb, ya:yb, za:zb); Vi = interpn(xo,yo,zo,Vi,xi,yi,zi);

hn = figure; ver=version; if str2double(ver(1)) > 7; hn = hn.Number; end
if get(handles.plot_opt,'Value') == 1 % x-line
    Vo(1,:) = Vi(:,round(size(Vi,2)/2),round(size(Vi,3)/2));
    xi1(1,:) = xi(:,round(size(xi,2)/2),round(size(xi,3)/2));
    [En, phi1] = schrod(Vo,xi1,0,0,n,m); phi1 = phi1/max(phi1(:));
    fprintf('     Zero-point Energy   = %4.4f meV \n',(En+max(loc(:)))*1000)
    
    fprintf('\n>> Drawing figure \n');
    if get(handles.one,'Value') == 1; subplot(1,3,1); else figure(hn); end
    plot([min(xi1) max(xi1)],En*[1 1],'k--') % V(x)
    hold on; box on; plot(xi1,Vo, 'LineWidth',2)
    xlabel('r (Angs)','FontSize',18); ylabel('V (eV)','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,3,2); else figure(hn+1); end
    plot(xi1,real(phi1),'r', 'LineWidth',2); box on % Psi(x)
    xlabel('r (Angst)','FontSize',18); ylabel('|\psi(x)|^2','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,3,3); else figure(hn+2); end
    plot(xi1,Vo, 'b', xi1, real(phi1)+En, 'r', 'LineWidth',2) 
    box on; hold on; plot([min(xi1) max(xi1)],En*[1 1],'k--') % line plot 2 
    xlabel('r (Angst)','FontSize',18); ylabel('V (eV)','FontSize',18);
    legend('V(x)','|\psi(x)|^2')
    
elseif (get(handles.plot_opt,'Value') == 2) % y-line
    Vo(1,:) = Vi(round(size(Vi,1)/2),:,round(size(Vi,3)/2));
    yi1(1,:) = yi(round(size(yi,2)/2),:,round(size(yi,3)/2));
    [En, phi1] = schrod(Vo,yi1,0,0,n,m);  phi1 = phi1/max(phi1(:));
    fprintf('     Zero-point Energy   = %4.4f meV \n',(En+max(loc(:)))*1000)
    
    fprintf('\n>> Drawing figure \n');
    if get(handles.one,'Value') == 1; subplot(1,3,1); else figure(hn); end
    plot([min(yi1) max(yi1)],En*[1 1],'k--') % V(y)
    hold on; box on; plot(yi1,Vo, 'LineWidth',2)
    xlabel('r (Angs)','FontSize',18); ylabel('V (eV)','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,3,2); else figure(hn+1); end
    plot(yi1,real(phi1),'r', 'LineWidth',2); box on % Psi(y)
    xlabel('r (Angst)','FontSize',18); ylabel('|\psi(y)|^2','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,3,3); else figure(hn+2); end
    plot(yi1,Vo, 'b', yi1, real(phi1)+En, 'r', 'LineWidth',2)
    box on; hold on; plot([min(yi1) max(yi1)],En*[1 1],'k--'); % line plot 2
    xlabel('r (Angst)','FontSize',18); ylabel('V (eV)','FontSize',18);
    legend('V(y)','|\psi(y)|^2')
    
elseif (get(handles.plot_opt,'Value') == 3) % z-line
    Vo(1,:) = Vi(round(size(Vi,1)/2),round(size(Vi,2)/2),:);
    zi1(1,:) = zi(round(size(zi,2)/2),round(size(zi,3)/2),:);
    [En, phi1] = schrod(Vo,zi1,0,0,n,m); phi1 = phi1/max(phi1(:));
    fprintf('     Zero-point Energy   = %4.4f meV \n',(En+max(loc(:)))*1000)
    
    fprintf('\n>> Drawing figure \n')
    if get(handles.one,'Value') == 1; subplot(1,3,1); else figure(hn); end
    plot([min(zi1) max(zi1)],En*[1 1],'k--') % V(z)
    hold on; box on; plot(zi1,Vo, 'LineWidth',2)
    xlabel('r (Angs)','FontSize',18); ylabel('V (eV)','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,3,2); else figure(hn+1); end
    plot(zi1,real(phi1),'r', 'LineWidth',2); box on % Psi(z)
    xlabel('r (Angst)','FontSize',18); ylabel('|\psi(z)|^2','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,3,3); else figure(hn+2); end
    plot(zi1,Vo, 'b', zi1, real(phi1)+En, 'r', 'LineWidth',2)
    box on; hold on; plot([min(zi1) max(zi1)],En*[1 1],'k--') % line plot 2
    xlabel('r (Angst)','FontSize',18); ylabel('V (eV)','FontSize',18);
    legend('V(z)','|\psi(z)|^2')
    
elseif (get(handles.plot_opt,'Value') == 4) % XY-plane
    Vo(:,:) = Vi(:,:,round(size(Vi,3)/2)); x(:,:) = xi(:,:,round(size(xi,3)/2));
    y(:,:) = yi(:,:,round(size(yi,3)/2)); [En, phi1] = schrod(Vo,x,y,0,n,m);
    fprintf('     Zero-point Energy   = %4.4f meV \n',(En+max(loc(:)))*1000)
    phi1 = phi1/max(phi1(:));
    
    fprintf('\n>> Drawing figure \n')
    if get(handles.one,'Value') == 1; subplot(1,2,1); else figure(hn); end
    surf(x,y,Vo) % V(x,y)
    axis tight; xlabel('a','FontSize',18); rotate3d on
    ylabel('b','FontSize',18); zlabel('V(x,y)','FontSize',18)
    title('Potential Energy','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,2,2); else figure(hn+1); end
%     imagesc(xi1,-yi1,interp2(phi1',2)); % Psi(x,y)
    [C,h] = contour(x,y,phi1);
    clabel(C, 'FontSize',18)
    daspect([1 1 1]); % Psi(x,y)
    axis tight; xlabel('a'); ylabel('b');
%     title('Muon probability density','FontSize',18)
    
elseif (get(handles.plot_opt,'Value') == 5) % XZ-plane
    Vo(:,:) = Vi(:,round(size(Vi,2)/2),:); x(:,:) = xi(:,round(size(xi,2)/2),:);
    z(:,:) = zi(:,round(size(zi,2)/2),:); [En, phi1] = schrod(Vo,x,z,0,n,m);
    fprintf('     Zero-point Energy   = %4.4f meV \n',(En+max(loc(:)))*1000)
    phi1 = phi1/max(phi1(:));
    
    fprintf('\n>> Drawing figure \n')
    if get(handles.one,'Value') == 1; subplot(1,2,1); else figure(hn); end
    surf(x,z,Vo) % V(x,z)
    axis tight; xlabel('a','FontSize',18); rotate3d on
    ylabel('c','FontSize',18); zlabel('V(x,z)','FontSize',18)
    title('Potential Energy','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,2,2); else figure(hn+1); end
%     imagesc(xi1,-zi1,interp2(phi1',2)); daspect([1 1 1]); % Psi(x,z)
    [C,h] = contour(x,z,phi1);
    clabel(C,'FontSize',18)
    daspect([1 1 1]);
    axis tight; xlabel('a'); ylabel('c');
%     title('Muon probability density','FontSize',18)
    
elseif (get(handles.plot_opt,'Value') == 6) % YZ-plane
    Vo(:,:) = Vi(round(size(Vi,1)/2),:,:); y(:,:) = yi(round(size(yi,1)/2),:,:);
    z(:,:) = zi(round(size(zi,1)/2),:,:); [En, phi1] = schrod(Vo,y,z,0,n,m);
    fprintf('     Zero-point Energy   = %4.4f meV \n',(En+max(loc(:)))*1000)
    phi1 = phi1/max(phi1(:));
    
    fprintf('\n>> Drawing figure \n')
    if get(handles.one,'Value') == 1; subplot(1,2,1); else figure(hn); end
    surf(y,z,Vo) % V(y,z)
    axis tight; xlabel('b','FontSize',18); rotate3d on
    ylabel('c','FontSize',18); zlabel('V(y,z)','FontSize',18)
    title('Potential Energy','FontSize',18)
    
    if get(handles.one,'Value') == 1; subplot(1,2,2); else figure(hn+1); end
%     imagesc(yi1,-zi1,interp2(phi1',2)); daspect([1 1 1]); % Psi(y,z)
    [C,h] = contour(y,z,phi1);
    clabel(C,'FontSize',18)
    daspect([1 1 1]);
    axis tight; xlabel('b'); ylabel('c');
%     title('Muon probability density','FontSize',18)
    
elseif (get(handles.plot_opt,'Value') == 7) % XYZ
    [En, phi1] = schrod(Vi,xi,yi,zi,n,m); phi1 = phi1/max(phi1(:));
    fprintf('     Zero-point Energy   = %4.4f meV \n',(En+max(loc(:)))*1000)
    latr = [max(xi(:))-min(xi(:)) 0 0; 0 max(yi(:))-min(yi(:)) 0; 0 0 max(zi(:))-min(zi(:))];
    fprintf('\n>> Drawing figure \n')
    No = str2num(get(handles.atom,'String')); iso = str2double(get(handles.Isosurf,'String'));
    col = get(handles.ion_col,'String'); rd = str2num(get(handles.ion_rad,'String'));
    if get(handles.opt_bond,'Value') > 1; d0 = str2num(get(handles.bond_l,'String')); end
    xb = str2num(get(handles.xb,'String'));
    yb = str2num(get(handles.yb,'String')); zb = str2num(get(handles.zb,'String'));
    if xb(1)<0; xb(1)=0; end; if xb(2)>1; xb(2)=1; end; if xb(2)<=xb(1); xa=xb; xb=[xa(2) xa(1)]; end
    if yb(1)<0; yb(1)=0; end; if yb(2)>1; yb(2)=1; end; if yb(2)<=yb(1); ya=yb; yb=[ya(2) ya(1)]; end
    if zb(1)<0; zb(1)=0; end; if zb(2)>1; zb(2)=1; end; if zb(2)<=zb(1); za=zb; xb=[za(2) za(1)]; end
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
        for mn = 1:3
            if get(handles.one,'Value') == 1; subplot(1,3,mn); else figure(hn+mn-1); end
            draw_atom(Cu,rd(no),col(no),1)
            if get(handles.opt_bond,'Value') == 2 && no <= 2
                if no == 1; ion1 = Cu;
                elseif no == 2
                    if get(handles.one,'Value') == 1; subplot(1,3,mn); else figure(hn+mn-1); end
                    drbond(ion1,Cu,d0,col(1:2),2);
                end
            elseif get(handles.opt_bond,'Value') == 3 && no == 1
                if get(handles.one,'Value') == 1; subplot(1,3,mn); else figure(hn+mn-1); end
                drplane(Cu,d0,col_oct,0.4*[1 1 1],0.3)
            elseif get(handles.opt_bond,'Value') == 4 && no <= 2
                if no == 1; ion1 = Cu;
                    if get(handles.one,'Value') == 1; subplot(1,3,mn); else figure(hn+mn-1); end
                    drplane(Cu,d0,col_oct,0.4*[1 1 1],0.3)
                elseif no == 2
                    if get(handles.one,'Value') == 1; subplot(1,3,mn); else figure(hn+mn-1); end
                    drbond(ion1,Cu,d0,col(1:2),2);
                end
            end
        end
    end
    if get(handles.one,'Value') == 1; subplot(1,3,1); else figure(hn); end
    p = patch(isosurface(xi,yi,zi,Vi,En));
    set(p,'FaceColor','y','EdgeColor','none','facelighting','gouraud','FaceAlpha',0.5)
    camlight right local; daspect([1 1 1]); rotate3d on;
    if (get(handles.box,'Value') == 1); drbox(latt); end
    if (get(handles.box2,'Value') == 1)
        drbox(latr,[min(xi(:)) min(yi(:)) min(zi(:))],'k');
    end
    if get(handles.fix,'Value') == 1
        xlabel('a','FontSize',12); ylabel('b','FontSize',12); 
        zlabel('c','FontSize',12); axis tight
    else axis off
    end
    if get(handles.one,'Value') == 1; subplot(1,3,2); else figure(hn+1); end
    psi = patch(isosurface(xi,yi,zi,phi1,1-iso));
    set(psi,'FaceColor','m','EdgeColor','none','facelighting','gouraud','FaceAlpha',0.5)
    camlight right local; daspect([1 1 1]); rotate3d on;
    if (get(handles.box,'Value') == 1); drbox(latt); end
    if (get(handles.box2,'Value') == 1)
        drbox(latr,[min(xi(:)) min(yi(:)) min(zi(:))],'k');
    end
    if get(handles.fix,'Value') == 1
        xlabel('a','FontSize',12); ylabel('b','FontSize',12); 
        zlabel('c','FontSize',12); axis tight
    else axis off
    end
    if get(handles.one,'Value') == 1; subplot(1,3,3); else figure(hn+2); end
    psi = patch(isosurface(xi,yi,zi,phi1,1-iso));
    set(psi,'FaceColor','m','EdgeColor','none','facelighting','gouraud','FaceAlpha',0.6)
    if (get(handles.box,'Value') == 1); drbox(latt); end
    if (get(handles.box2,'Value') == 1)
        drbox(latr,[min(xi(:)) min(yi(:)) min(zi(:))],'k');
    end
    xa = linspace(0,1,size(V,1)); ya = linspace(0,1,size(V,2));
    za = linspace(0,1,size(V,3)); [xi,yi,zi] = ndgrid(xa,ya,za);
    xi = xi(:); yi = yi(:); zi = zi(:);
    x = xi*latt(1,1) + yi*latt(2,1) + zi*latt(3,1);
    y = xi*latt(1,2) + yi*latt(2,2) + zi*latt(3,2);
    z = xi*latt(1,3) + yi*latt(2,3) + zi*latt(3,3);
    x = reshape(x,size(V)); y = reshape(y,size(V));
    z = reshape(z,size(V)); 
    for i = 1:2
        a = abs(xa-xb(i)); [~,xb(i)] = min(a);
        a = abs(ya-yb(i)); [~,yb(i)] = min(a);
        a = abs(za-zb(i)); [~,zb(i)] = min(a);
    end
    x = x(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
    y = y(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
    z = z(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
    V = V(xb(1):xb(2), yb(1):yb(2), zb(1):zb(2));
    p = patch(isosurface(x,y,z,V,En));
    set(p,'FaceColor','y','EdgeColor','none','FaceAlpha',0.4)
    daspect([1 1 1]); rotate3d on; camlight right local;
    if get(handles.fix,'Value') == 1
        xlabel('a','FontSize',12); ylabel('b','FontSize',12); 
        zlabel('c','FontSize',12); axis tight
    else axis off
    end
end
set(handles.Result,'String',(En+max(loc(:)))*1000);
t = toc; h = floor(t/3600); m = floor((t-h*3600)/60); t = t-h*3600 - m*60;
fprintf('   ============================================ \n')
if h >= 1; fprintf('   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n',h,m,t);  
elseif m >= 1; fprintf('   Elapsed time is %1.0f min %1.4f sec.\n',m,t);
else fprintf('   '); toc
end

function load_Callback(~, ~, handles)
global filename pathname
[filename, pathname] = uigetfile({'*.*', 'All Files (*.*)'});
if ~isnumeric(pathname); set(handles.pathname,'String',[pathname filename]); end

function muon_Callback(~, ~, ~)
function muon_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_opt_Callback(~, ~, handles)
if get(handles.plot_opt,'Value') == 7
    set(handles.box2,'Value',1); set(handles.box,'Value',1); set(handles.fix,'Value',1)
    set(handles.uipanel2,'ForegroundColor',[0 0 1])
    set(handles.atom,'BackgroundColor',[1 1 1])
    set(handles.atom,'ForegroundColor',[0 0 0])
    set(handles.Isosurf,'BackgroundColor',[1 1 1])
    set(handles.Isosurf,'ForegroundColor',[0 0 0])
    set(handles.ion_rad,'BackgroundColor',[1 1 1])
    set(handles.ion_rad,'ForegroundColor',[0 0 0])
    set(handles.ion_col,'BackgroundColor',[1 1 1])
    set(handles.ion_col,'ForegroundColor',[0 0 0])
    set(handles.xb,'BackgroundColor',[1 1 1]); set(handles.xb,'ForegroundColor',[0 0 0])
    set(handles.yb,'BackgroundColor',[1 1 1]); set(handles.yb,'ForegroundColor',[0 0 0])
    set(handles.zb,'BackgroundColor',[1 1 1]); set(handles.zb,'ForegroundColor',[0 0 0])
    if get(handles.opt_bond,'Value') > 1
        set(handles.bond_l,'BackgroundColor',[1 1 1])
        set(handles.bond_l,'ForegroundColor',[0 0 0])
    end
else set(handles.box,'Value',0); set(handles.fix,'Value',0); set(handles.box2,'Value',0)
    set(handles.uipanel2,'ForegroundColor',0.31*[1 1 1])
    set(handles.atom,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.atom,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.Isosurf,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.Isosurf,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.ion_rad,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.ion_rad,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.ion_col,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.ion_col,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.xb,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.xb,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.yb,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.yb,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.zb,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.zb,'ForegroundColor',[0.8 0.8 0.8])
    set(handles.bond_l,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.bond_l,'ForegroundColor',[0.8 0.8 0.8])
end

function plot_opt_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nx_Callback(~, ~, ~)
function Nx_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dx_Callback(~, ~, ~)
function dx_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ny_Callback(~, ~, ~)
function Ny_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dy_Callback(~, ~, ~)
function dy_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nz_Callback(~, ~, ~)
function Nz_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dz_Callback(~, ~, ~)
function dz_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [En, phi] = schrod(V,x,y,z,n,mass)
h_bar = (6.62606957e-34)/(2*pi); joule_eV = 6.24150934e18;
k = h_bar^2/(2*mass)*joule_eV*1E20; % in eV
N = size(V); M = 1; L = 1; V = V(:); A = zeros(length(V));
if length(N) == 2 && min(N) == 1 % 1D
    dx = ((max(x(:))-min(x(:)))/length(x))^2; dx = k/dx;
    for i = 1:length(V)
        if i > 1; A(i,i-1) = A(i,i-1) + dx; end
        if i < length(V); A(i,i+1) = A(i,i+1) + dx; end
        A(i,i) = A(i,i) - (2*dx + V(i));
    end
elseif length(N) == 2 % 2D
    dx = ((max(x(:))-min(x(:)))/size(x,1))^2; dx = k/dx;
    dy = ((max(y(:))-min(y(:)))/size(y,2))^2; dy = k/dy;
    for i = 1:length(V)
        if i > 1+(M-1)*N(1); A(i,i-1) = A(i,i-1) + dx; end
        if i < M*N(1); A(i,i+1) = A(i,i+1) + dx; end
        if i > N(1); A(i,i-N(1)) = A(i,i-N(1)) + dy; end
        if i <= N(1)*(N(2)-1); A(i,i+N(1)) = A(i,i+N(1)) + dy; end
        A(i,i) = -(2*dx + 2*dy + V(i)); if i/N(1) == M; M = M+1; end
    end
elseif length(N) == 3 % 3D
    dx = ((max(x(:))-min(x(:)))/N(1))^2;
    dy = ((max(y(:))-min(y(:)))/N(2))^2;
    dz = ((max(z(:))-min(z(:)))/N(3))^2;
    for i = 1:length(V)
        if i > 1+(M-1)*N(1); A(i,i-1) = A(i,i-1) + k/dx; end
        if i < M*N(1); A(i,i+1) = A(i,i+1) + k/dx; end
        if i > N(1)*(N(2)*(L-1)+1); A(i,i-N(1)) = A(i,i-N(1)) + k/dy; end
        if i <= N(1)*(N(2)*L-1); A(i,i+N(1)) = A(i,i+N(1)) + k/dy; end
        if i > N(1)*N(2); A(i,i-N(1)*N(2)) = A(i,i-N(1)*N(2)) + k/dz; end
        if i <= N(1)*N(2)*(N(3)-1); A(i,i+N(1)*N(2))=A(i,i+N(1)*N(2))+k/dz; end
        A(i,i) = -(2*k/dx + 2*k/dy + 2*k/dz + V(i));
        if i/N(1) == M; M = M+1; end; if (M-1)/N(2) == L; L = L+1; end
    end
end
fprintf('\n>> Solving Schrodinger equation \n')
[psi,En] = eig(A); En = -En(end-n,end-n);
phi = reshape(psi(:,end-n).*conj(psi(:,end-n)),N);

function geometry = poscar(filename)
if ~isnumeric(filename); fid = fopen(filename); else fid = filename; end
fgetl(fid); scale = fscanf(fid, '%f',1);
geometry.lattice = fscanf(fid, '%f %f %f', [3 3])'; 
geometry.lattice = geometry.lattice*scale;
fgetl(fid); line = fgetl(fid); cartesian = 0; has_symbols_line = false;
if sum(isstrprop(line, 'digit')) == 0; has_symbols_line = true;
    geometry.symbols = regexp(line, '([^ ]*)', 'match'); line = fgetl(fid); 
else geometry.symbols = {}; 
end
geometry.atomcount = sscanf(line,'%d'); natoms = sum(geometry.atomcount);
line = fgetl(fid); geometry.selective = 0;
if line(1) == 's' || line(1) == 'S'; geometry.selective = 1; line = fgetl(fid); end
if line(1) == 'C' || line(1) == 'c' || line(1) == 'K' || line(1) =='k'; cartesian = 1; end
for i = 1:natoms
	line = fgetl(fid); geometry.coords(i,:) = sscanf(line, '%f %f %f');
    if ~has_symbols_line
        if numel(strfind(line,'!') == 1)
            str = regexp(line, '[^ ]*', 'match'); str = str{end};
            if numel(geometry.symbols) == 0; newelement = true;
            elseif strcmp(geometry.symbols{end},str) == 0; newelement = true; end
            if newelement; geometry.symbols{end+1} = str; end
        end
    end
end
if cartesian == 1; geometry.coords = geometry.coords*scale;
geometry.coords = geometry.coords/geometry.lattice; end 
if ~isnumeric(filename); fclose(fid); end

function n_Callback(~, ~, ~)
function n_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Isosurf_Callback(~, ~, ~)
function Isosurf_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function atom_Callback(~, ~, ~)
function atom_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_axis_Callback(~, ~, ~)
function opt_axis_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_muon_Callback(~, ~, handles)
set(handles.proton,'Value',0); set(handles.proton,'ForegroundColor',[0 0 0])
set(handles.opt_muon,'Value',1); set(handles.opt_muon,'ForegroundColor',[0 0 1])

function proton_Callback(~, ~, handles)
set(handles.proton,'Value',1); set(handles.proton,'ForegroundColor',[0 0 1])
set(handles.opt_muon,'Value',0); set(handles.opt_muon,'ForegroundColor',[0 0 0])

function pathname_Callback(~, ~, ~)
function pathname_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cell_Callback(~, ~, ~)
function cell_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_grid_Callback(~, ~, ~)
function opt_grid_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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

function draw_atom(target,r,col,alp); hold on
for i = 1:size(target,1)
    [rx,ry,rz] = sphere(10); rx = rx*r; ry = ry*r; rz = rz*r;
    h = surf((rx+target(i,1)),(ry+target(i,2)),(rz+target(i,3)));
    set(h,'edgecolor','none','facecolor',col,'facelighting','gouraud'); alpha(h,alp)
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

function fix_Callback(~, ~, handles)
if get(handles.plot_opt,'Value') < 7; set(handles.fix,'Value',0); end

function bond_l_Callback(~, ~, ~)
function bond_l_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_bond_Callback(~, ~, handles)
if get(handles.opt_bond,'Value') > 1 && get(handles.plot_opt,'Value') == 7
    set(handles.bond_l,'BackgroundColor',[1 1 1])
    set(handles.bond_l,'ForegroundColor',[0 0 0])
else set(handles.bond_l,'BackgroundColor',[0.8 0.8 0.8])
    set(handles.bond_l,'ForegroundColor',[0.8 0.8 0.8])
end

function opt_bond_CreateFcn(hObject, ~, ~)
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

function one_Callback(~, ~, handles)
set(handles.one,'Value',1); set(handles.one,'ForegroundColor',[0 0 1])
set(handles.sep,'Value',0); set(handles.sep,'ForegroundColor',[0 0 0])

function sep_Callback(~, ~, handles)
set(handles.one,'Value',0); set(handles.one,'ForegroundColor',[0 0 0])
set(handles.sep,'Value',1); set(handles.sep,'ForegroundColor',[0 0 1])

function box2_Callback(~, ~, handles)
if get(handles.plot_opt,'Value') < 7; set(handles.box2,'Value',0); end

function box_Callback(~, ~, handles)
if get(handles.plot_opt,'Value') < 7; set(handles.box,'Value',0); end
