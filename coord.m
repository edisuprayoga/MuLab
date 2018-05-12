function pos = coord(mode,pos,geo,grd)
if nargin == 2
    geo = [];
elseif nargin < 4
    grd = geo;
end
if isnumeric(geo)
    latt = geo;
elseif ischar(geo)
    geo = get_poscar(geo);
    latt = geo.lattice;
else
    latt = geo.lattice;
end
switch mode
    case 'd2r'
        mode = 'dir2real';
    case 'd2g'
        mode = 'dir2grid';
    case 'g2d'
        mode = 'grid2dir';
    case 'g2r'
        mode = 'grid2real';
    case 'r2d'
        mode = 'real2dir';
    case 'r2g'
        mode = 'real2grid';
end
switch mode
    case 'dir2real'
        for i = 1:size(pos,1)
            pos(i,1:3) = pos(i,1)*latt(1,:)+pos(i,2)*latt(2,:)+pos(i,3)*latt(3,:);
        end
    case 'dir2grid'
        x = linspace(0,1,grd(1));
        y = linspace(0,1,grd(2));
        z = linspace(0,1,grd(3));
        for i = 1:size(pos,1)
            a = abs(x-pos(i,1)); [~,a] = min(a);
            b = abs(y-pos(i,2)); [~,b] = min(b);
            c = abs(z-pos(i,3)); [~,c] = min(c);
            pos(i,:) = [a b c];
        end
    case 'dir2sup'
        if isnumeric(grd)
            geo2 = supercell(geo,grd);
        else
            geo2 = grd;
        end
        for i = 1:size(pos,1)
            pos(i,:) = pos(i,:)-0.5*[1 1 1];
        end
        pos = coord('dir2real',pos,geo);
        pos = coord('real2dir',pos,geo2);
        for i = 1:size(pos,1)
            pos(i,:) = pos(i,:)+0.5*[1 1 1];
        end
    case 'grid2dir'
        x = linspace(0,1,grd(1));
        y = linspace(0,1,grd(2));
        z = linspace(0,1,grd(3));
        for i = 1:size(pos,1)
            pos(i,:) = [x(pos(i,1)) y(pos(i,2)) z(pos(i,3))];
        end
    case 'grid2real'
        pos = coord('grid2dir',pos,grd);
        pos = coord('dir2real',pos,geo);
    case 'real2dir'
        latt = inv(latt);
        for i = 1:size(pos,1)
            pos(i,1:3) = pos(i,1)*latt(1,:)+pos(i,2)*latt(2,:)+pos(i,3)*latt(3,:);
        end
    case 'real2grid'
        pos = coord('real2dir',pos,geo);
        pos = coord('dir2grid',pos,grd);
    case 'real2sup'
        pos = coord('real2dir',pos,geo);
        if isnumeric(grd)
            geo2 = supercell(latt,grd);
        else
            geo2 = grd;
        end
        pos = coord('dir2sup',pos,geo,geo2);
        pos = coord('dir2real',pos,geo2);
    case 'reciproc'
        pos = 2*pi*inv(pos)';
end
end