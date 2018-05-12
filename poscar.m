function geometry = poscar(filename)
if ~isnumeric(filename); fid = fopen(filename); else fid = filename; end
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
            str = regexp(line, '[^ ]*', 'match'); str = str{end};
            newelement = false; % is this atom a new element?
            if numel(geometry.symbols) == 0; newelement = true;
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