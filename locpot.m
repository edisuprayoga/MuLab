% [loc, geo] =
function [locpot, geometry] = locpot(filename)
fid = fopen(filename);
geometry = poscar(fid); fgetl(fid);
gridsize = fscanf(fid, '%d %d %d', [3 1])';
locpot = fscanf(fid, '%f', [prod(gridsize,2) 1])';
locpot = reshape(locpot,gridsize);
fclose(fid);
end