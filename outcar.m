function result = outcar(filename)
fid = fopen(filename);
while ~feof(fid); line = fgetl(fid);
    if numel(regexp(line,'magnetization \(x\)'))==1
        pos = ftell(fid);
    end
end
fseek(fid,pos,'bof'); fgetl(fid); line = fgetl(fid); fgetl(fid);
if length(line) > 40; sp = 6; else sp = 5; end
mx = fscanf(fid,'%f',[sp inf])'; mx = mx(:,end);
while ~feof(fid); line = fgetl(fid);
    if numel(regexp(line,'magnetization \(y\)'))==1
        pos = ftell(fid);
    end
end
fseek(fid,pos,'bof'); fgetl(fid); fgetl(fid); fgetl(fid);
my = fscanf(fid,'%f',[sp inf])'; my = my(:,end);
while ~feof(fid); line = fgetl(fid);
    if numel(regexp(line,'magnetization \(z\)'))==1
        pos = ftell(fid);
    end
end
fseek(fid,pos,'bof'); fgetl(fid); fgetl(fid); fgetl(fid);
mz = fscanf(fid,'%f',[sp inf])'; mz = mz(:,end);
fclose(fid); result = [mx my mz];