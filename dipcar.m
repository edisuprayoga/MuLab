%[H, Hx, Hy, Hz, geo] =
function [H, Hx, Hy, Hz, geo] = dipcar(filename)
fid = fopen(filename); 
geo = poscar(fid); fgetl(fid); 
gridsize = fscanf(fid, '%d %d %d', [3 1])';
H = fscanf(fid, '%f', [prod(gridsize,2) 1])';
H = reshape(H,gridsize);
%if fgetl(fid) == ''; fgetl(fid); end
gridsize = fscanf(fid, '%d %d %d', [3 1])';
Hx = fscanf(fid, '%f', [prod(gridsize,2) 1])';
Hx = reshape(Hx,gridsize);
%if fgetl(fid) == ''; fgetl(fid); end
gridsize = fscanf(fid, '%d %d %d', [3 1])';
Hy = fscanf(fid, '%f', [prod(gridsize,2) 1])';
Hy = reshape(Hy,gridsize);
%if fgetl(fid) == ''; fgetl(fid); end
gridsize = fscanf(fid, '%d %d %d', [3 1])';
Hz = fscanf(fid, '%f', [prod(gridsize,2) 1])';
Hz = reshape(Hz,gridsize);
fclose(fid);