function status = main(opt_pos,opt_spn,opt_job,muon,coord,R_dip,file1,file2,ion_no, ...
    mag1,spin1,mag2,spin2,zpe_range,zpe_ngrid,R_min,N_cal,N_grid,lim,gpu)
clc; log = 1; opt = size(muon,1); if opt_job == 4; opt = 1; end
if isempty(file1); file1 = pwd; end; if isempty(file2); file2 = pwd; end
if opt > 1; tic; jam = clock; log = 0;
    fprintf('   ============================================ \n')
    fprintf('            Dipole Fields Calculation \n')
    fprintf('   ============================================ \n')
    fprintf(['>> Reading input files \n   ' file1 '\n'])
end; status = [];
for i = 1:opt;
    if opt > 1; fprintf('\n   run %1.0f/%1.0f ',i,opt); end
    if opt_pos == 1 % perfect
        if opt_spn == 1 % manual
            if opt_job == 1 % dipole
                Hdip = dip111(muon(i,:),coord,R_dip,file1,ion_no,mag1,spin1,gpu,log);
            elseif opt_job == 2 % zpe
                Hdip = dip112(muon(i,:),coord,R_dip,file1,ion_no,mag1,spin1,zpe_range,zpe_ngrid,gpu,log);
            elseif opt_job == 3 % test
                Hdip = dip113(muon(i,:),coord,R_dip,file1,ion_no,mag1,spin1,R_min,N_cal,gpu,log);
            elseif opt_job == 4 % makefile
                Hdip = dip114(R_dip,file1,ion_no,mag1,spin1,N_grid,gpu);
            end
        elseif opt_spn == 2 % vasp
            if opt_job == 1 % dipole
                Hdip = dip121(muon(i,:),coord,R_dip,file1,lim,gpu,log);
            elseif opt_job == 2 % zpe
                Hdip = dip122(muon(i,:),coord,R_dip,file1,zpe_range,zpe_ngrid,lim,gpu,log);
            elseif opt_job == 3 % test
                Hdip = dip123(muon(i,:),coord,R_dip,file1,R_min,N_cal,lim,gpu,log);
            elseif opt_job == 4 % makefile
                Hdip = dip124(R_dip,file1,N_grid,lim,gpu);
            end
        elseif opt_spn == 3 % chgcar
            if opt_job == 1 % dipole
                Hdip = dip131(muon(i,:),coord,R_dip,file1,gpu,log);
            elseif opt_job == 2 % zpe
                Hdip = dip132(muon(i,:),coord,R_dip,file1,zpe_range,zpe_ngrid,gpu,log);
            elseif opt_job == 3 % test
                Hdip = dip133(muon(i,:),coord,R_dip,file1,R_min,N_cal,gpu,log);
            elseif opt_job == 4 % makefile
                Hdip = dip134(R_dip,file1,N_grid,gpu);
            end
        end
    elseif opt_pos == 2 % defect
        if opt_spn == 1 % manual
            if opt_job == 1 % dipole
                Hdip = dip211(muon(i,:),coord,R_dip,file1,file2,ion_no,mag1,spin1,mag2,spin2,gpu,log);
            elseif opt_job == 2 % zpe
                Hdip = dip212(muon(i,:),coord,R_dip,file1,file2,ion_no,mag1,spin1,mag2,spin2,zpe_range,zpe_ngrid,gpu,log);
            elseif opt_job == 3 % test
                Hdip = dip213(muon(i,:),coord,R_dip,file1,file2,ion_no,mag1,spin1,mag2,spin2,R_min,N_cal,gpu,log);
            elseif opt_job == 4 % makefile
                Hdip = dip214(R_dip,file1,file2,ion_no,mag1,spin1,mag2,spin2,N_grid,gpu);
            end
        elseif opt_spn == 2 % vasp
            if opt_job == 1 % dipole
                Hdip = dip221(muon(i,:),coord,R_dip,file1,file2,lim,gpu,log);
            elseif opt_job == 2 % zpe
                Hdip = dip222(muon(i,:),coord,R_dip,file1,file2,zpe_range,zpe_ngrid,lim,gpu,log);
            elseif opt_job == 3 % test
                Hdip = dip223(muon(i,:),coord,R_dip,file1,file2,R_min,N_cal,lim,gpu,log);
            elseif opt_job == 4 % makefile
                Hdip = dip224(R_dip,file1,file2,N_grid,lim,gpu);
            end
        elseif opt_spn == 3 % chgcar
            if opt_job == 1 % dipole
                Hdip = dip231(muon(i,:),coord,R_dip,file1,file2,gpu,log);
            elseif opt_job == 2 % zpe
                Hdip = dip232(muon(i,:),coord,R_dip,file1,file2,zpe_range,zpe_ngrid,gpu,log);
            elseif opt_job == 3 % test
                Hdip = dip233(muon(i,:),coord,R_dip,file1,file2,R_min,N_cal,gpu,log);
            elseif opt_job == 4 % makefile
                Hdip = dip234(R_dip,file1,file2,N_grid,gpu);
            end
        end
    end
    if opt > 1; fprintf('--> %4.4f Gauss [%4.4f, %4.4f, %4.4f]',Hdip); end
    status = [status; Hdip];
end
if opt > 1; pathname = pwd; fprintf('\n\n>> Writing log file\n');
    fid = fopen([pathname '/log_main']);
    if fid ~= -1
        i = 1; fid = fopen([pathname '/log_main(1)']);
        while fid ~= -1; i = i+1; fclose(fid);
            fid = fopen([pathname '/log_main(' num2str(i) ')']);
        end
        fid = fopen([pathname '/log_main(' num2str(i) ')'],'w');
        fprintf(['     ' pathname '/log_main(' num2str(i) ') \n'])
    else fid = fopen([pathname '/log_main'],'w');
        fprintf(['     ' pathname '/log_main \n'])
    end
    fprintf(fid,'   ======================================================================== \n');
    fprintf(fid,'                        Dipole Fields Calculation \n');
    fprintf(fid,'   ======================================================================== \n');
    if jam(5) < 10; fprintf(fid,['   Running on: ' date '   %1.0f:0%1.0f \n'],jam(4:5));
    else fprintf(fid,['   Running on: ' date '   %1.0f:%1.0f \n\n'],jam(4:5));
    end
    if opt_pos == 1;     fprintf(fid,'   opt_pos = 1; %% Perfect crystal \n');
    elseif opt_pos == 2; fprintf(fid,'   opt_pos = 2; %% Defect crystal \n');
    end
    if opt_spn == 1;     fprintf(fid,'   opt_spn = 1; %% Manual input spin \n');
    elseif opt_spn == 2; fprintf(fid,'   opt_spn = 2; %% Input spin from VASP\n');
    elseif opt_spn == 3; fprintf(fid,'   opt_spn = 3; %% Spin distributed (from CHGCAR) \n');
    end
    if opt_job == 1;     fprintf(fid,'   opt_job = 1; %% Dipole fields \n');
    elseif opt_job == 2; fprintf(fid,'   opt_job = 2; %% Dipole fields with ZPE \n');
    elseif opt_job == 3; fprintf(fid,'   opt_job = 3; %% Converge test \n');
    elseif opt_job == 4; fprintf(fid,'   opt_job = 4; %% Makefile \n');
    end
    fprintf(fid,'\n   Input files : \n');
    fprintf(fid,['   ' file1 '\n']);
    if opt_pos == 2; fprintf(fid,['   ' file2 '\n']); end
    if opt_spn == 1; geo = poscar(file1);
        if isempty(geo.symbols);
            if length(ion_no) == 1; fprintf(fid,['\n   ion : '  num2str(ion_no) '\n']);
            else fprintf(fid,['\n   ion of : ['  num2str(ion_no) ']\n']); end;
        else ion = '';
            for i = 1:length(ion_no); ion = [ion cell2mat(geo.symbols(ion_no(i))) ' ']; end
            fprintf(fid,['\n   ion : ' ion '\n']);
        end
        if ischar(spin1); fid2 = fopen(spin1); spin1 = str2num(fgetl(fid2));
        while ~feof(fid2); spin1 = [spin1 str2num(fgetl(fid2))]; end; fclose(fid2);
        end
        if opt_pos == 2; fprintf(fid,['   Magmom1: '  num2str(mag1) ' muB \n']);
            fprintf(fid,'   Spins1: \n');
            fprintf(fid,'     %4.3f   %4.3f   %4.3f\n',spin1);
            if ischar(spin2); fid2 = fopen(spin2); spin2 = str2num(fgetl(fid2));
            while ~feof(fid2); spin2 = [spin2 str2num(fgetl(fid2))]; end; fclose(fid2);
            end
            fprintf(fid,['   Magmom2: '  num2str(mag2) ' muB \n']);
            fprintf(fid,'   Spins2: \n');
            fprintf(fid,'     %4.3f   %4.3f   %4.3f\n',spin2);
        else fprintf(fid,['   Magmom: '  num2str(mag1) ' muB \n']);
            fprintf(fid,'   Spins: \n');
            fprintf(fid,'     %4.3f   %4.3f   %4.3f\n',spin1);
        end
    elseif opt_spn == 2
        spin1 = outcar([file1 '/OUTCAR']); num = (1:size(spin1,1))';
        spin = sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2);
        if length(lim) == 1; cond = (spin(:,1) < lim);
            spin1(cond,:) = []; num(cond,:) = [];
        else cond = (spin(:,1) < lim(2)); spin1(cond,:) = []; num(cond,:) = [];
             cond = (spin(:,1) > lim(1)); spin1(cond,:) = []; num(cond,:) = [];
        end
        spinprt1 = [num spin1 sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2)]';
        fprintf(fid,'\n     No      Mx      My      Mz      Mtotal');
        fprintf(fid,'\n   --------------------------------------------\n');
        fprintf(fid,'   %4.0f    %6.3f  %6.3f  %6.3f    %6.3f \n',spinprt1(:));
        fprintf(fid,'   --------------------------------------------\n');
        fprintf(fid,'    sum    %6.3f  %6.3f  %6.3f    %6.3f \n',sum(spinprt1(2,:)), ...
            sum(spinprt1(3,:)),sum(spinprt1(4,:)),sum(spinprt1(5,:)));
        if opt_pos == 2 
            spin1 = outcar([file2 '/OUTCAR']); num = (1:size(spin1,1))';
            spin = sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2);
            if length(lim) == 1; cond = (spin(:,1) < lim);
                spin1(cond,:) = []; num(cond,:) = [];
            else cond = (spin(:,1) < lim(2)); spin1(cond,:) = []; num(cond,:) = [];
                 cond = (spin(:,1) > lim(1)); spin1(cond,:) = []; num(cond,:) = [];
            end
            spinprt2 = [num spin1 sqrt(spin1(:,1).^2+spin1(:,2).^2+spin1(:,3).^2)]';
            fprintf(fid,'\n     No2     Mx      My      Mz      Mtotal');
            fprintf(fid,'\n   --------------------------------------------\n');
            fprintf(fid,'   %4.0f    %6.3f  %6.3f  %6.3f    %6.3f \n',spinprt2(:));
            fprintf(fid,'   --------------------------------------------\n');
            fprintf(fid,'    sum    %6.3f  %6.3f  %6.3f    %6.3f \n',sum(spinprt2(2,:)), ...
                sum(spinprt2(3,:)),sum(spinprt2(4,:)),sum(spinprt2(5,:)));
        end
    elseif opt_spn == 3;
        [~,Sx2,Sy2,Sz2,geo2] = chgcar(file1); 
        latt2 = geo2.lattice; vol = abs(dot(latt2(1,:),cross(latt2(2,:),latt2(3,:))));
        sumx= [sum(Sx2(:)) sum(Sy2(:)) sum(Sz2(:))]/vol;
        if opt_pos == 1
            fprintf(fid,'\n   Magmom : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
        else fprintf(fid,'\n   Magmom1 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumx(1)^2+sumx(2)^2+sumx(3)^2),sumx);
            [~,Sx2,Sy2,Sz2,geo2] = chgcar(file2);
            latt2 = geo2.lattice; vol = abs(dot(latt2(1,:),cross(latt2(2,:),latt2(3,:))));
            sumy = [sum(Sx2(:)) sum(Sy2(:)) sum(Sz2(:))]/vol;
            fprintf(fid,'   Magmom2 : %4.4f muB [%4.4f %4.4f %4.4f]\n',sqrt(sumy(1)^2+sumy(2)^2+sumy(3)^2),sumy);
        end
    end
    fprintf(fid,'\n   Calculation range : %1.0f Angstrom \n',R_dip);
    if opt_job == 2
        if length(zpe_range) == 1
            fprintf(fid,['   ZPE Range (A) : ' num2str(zpe_range) ' \n']);
        else fprintf(fid,['   ZPE Range (A) : [' num2str(zpe_range) '] \n']);
        end
        if length(zpe_ngrid) == 1
            fprintf(fid,['   Num of points : ' num2str(zpe_ngrid) ' \n']);
        else fprintf(fid,['   Num of points : [' num2str(zpe_ngrid) '] \n']);
        end
    end
    out = [muon status]';
    fprintf(fid,'\n   muon position (x,y,z)       H(G)        [Hx, Hy, Hz] \n');
    fprintf(fid,'   %6.4f  %6.4f  %6.4f    %8.4f    [ %8.4f  %8.4f  %8.4f ] \n',out(:));
    fprintf(fid,'\n   ======================================================================== \n');
    fprintf(fid,'                          Calculations Completed! \n');
    fprintf(fid,'   ======================================================================== \n');
    t = toc; h = floor(t/3600); m = floor((t-h*3600)/60); t = t-h*3600 - m*60;
    if h >= 1; fprintf(fid,'   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n\n',h,m,t);
    elseif m >= 1; fprintf(fid,'   Elapsed time is %1.0f min %1.4f sec.\n\n',m,t);
    else fprintf(fid,'   Elapsed time is %1.4f sec.\n\n',t);
    end
    fclose(fid);
    fprintf('\n   ============================================ \n')
    fprintf('             Calculations Completed! \n')
    fprintf('   ============================================ \n')
    if h >= 1; fprintf('   Elapsed time is %1.0f hrs %1.0f min and %1.4f sec.\n',h,m,t);
    elseif m >= 1; fprintf('   Elapsed time is %1.0f min %1.4f sec.\n',m,t);
    else fprintf('   Elapsed time is %1.4f sec.\n',t);
    end
end
