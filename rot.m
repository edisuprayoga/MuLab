function r = rot(vek,ax,theta)
switch ax
    case 'a'; ax = 'x';
    case 'b'; ax = 'y';
    case 'c'; ax = 'z';
end
switch ax
    case 'x'
        R = [1 0 0; 0 cosd(theta) sind(theta); 0 -sind(theta) cosd(theta)];
    case 'y'
        R = [cosd(theta) 0 -sind(theta); 0 1 0; sind(theta) 0 cosd(theta)];
    case 'z'
        R = [cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; 0 0 1];
end
r = zeros(size(vek));
for i = 1:size(vek,1);
    r0 = R*vek(i,:)'; r(i,:) = r0';
end