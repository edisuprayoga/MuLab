%[chg, mx, my, mz, geo] =
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