% cell array to stock the coordinate of undivided lattice
l=100; % # of lattice
len=(l-1)*3.16;
MOS2=cell(l,l,3);  % Initial MoS2 coordiante
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Mo    %%%%%%%%%%%%%%%%%%%%
% reassign  coordinate
for i=1:l 
    for j=1:l
        if j<=(l-i+1)
            % relative coordinates
            XMo=(j-1)*3.16;
            YMo=(l-i)*3.16;
            % reassign to absolute coordinates
            XMo=XMo-YMo*sind(30);
            YMo=YMo*cosd(30);
            Mo=[XMo,YMo,0];
            MOS2{i,j,1}=Mo; 
        end
    end
end
%%%%%%%%%%%%%%%%%%    S1    %%%%%%%%%%%%%%%%%%%%%%%%%%%
% reassign  coordinate
for i=1:l 
    for j=1:l
        if j<(l-i+1)
            % reassign to absolute coordinates
            MOS2{i,j,2}=MOS2{i,j,1}+[1.5805,-0.9125,-1.5061];
        end
    end
end
%%%%%%%%%%%%%%%%%%    S2    %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:l 
    for j=1:l
        if j<(l-i+1)
            % reassign to absolute coordinates
            MOS2{i,j,3}=MOS2{i,j,1}+[1.5805,-0.9125,1.5061];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%    Monte Carlo Simulation   %%%%%%%%%%%%%%%
Simnum=160000;
temX=0;
temY=0;
temZ=0;
for i=1:Simnum
    tX=ceil(rand(1,1)*l);
    tY=ceil(rand(1,1)*l);
    tZ=ceil(rand(1,1)*3);
    while (all(ismember(temX,tX))==1 && all(ismember(temY,tY))==1 && all(ismember(temZ,tZ))==1) || tZ==1 && tY>l-tX+1 || tZ==2 && tY>=l-tX+1 || tZ==3 && tY>=l-tX+1    % not repeat & in the triangle
            tX=ceil(rand(1,1)*l);
            tY=ceil(rand(1,1)*l);
            tZ=ceil(rand(1,1)*3);
    end
    if tX==1 || tY==1 % atom on the first-row of right-angle side
        MOS2{tX,tY,tZ}=[];
    elseif isempty(MOS2{tX,tY-1,tZ}) || isempty(MOS2{tX-1,tY,tZ}) || isempty(MOS2{tX,tY+1,tZ}) %hypotenuse & side empty situation
        MOS2{tX,tY,tZ}=[];
    end
    temX=[temX;tX];
    temY=[temY;tY];
    temZ=[temZ;tZ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% StatRest# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num=0;
for i=1:l
    for j=1:l
        if isempty(MOS2{i,j,1})==0
            num=num+1;
        end
    end
end
% print S1 atom
for i=1:l
    for j=1:l
        if isempty(MOS2{i,j,2})==0
            num=num+1;
        end
    end
end
% print S2 atom
for i=1:l
    for j=1:l
        if isempty(MOS2{i,j,3})==0
            num=num+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%  Output  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print structure
fprintf('print structure')
MOS2
% print Mo atom
first=fopen('MOS2.txt','a');
fid=fopen('MOS2.txt','w');
fprintf(fid,'%d\n',num);
fprintf(fid,'%s\n','MoS2');
for i=1:l
    for j=1:l
        if isempty(MOS2{i,j,1})==0
            fprintf(fid,'%s\t','Mo');
            fprintf(fid,'%15.5g\t',MOS2{i,j,1}(1));
            fprintf(fid,'%15.5g\t',MOS2{i,j,1}(2));
            fprintf(fid,'%15.5g\n',MOS2{i,j,1}(3));
        end
    end
end
% print S1 atom
for i=1:l
    for j=1:l
        if isempty(MOS2{i,j,2})==0
            fprintf(fid,'%s\t','S');
            fprintf(fid,'%15.5g\t',MOS2{i,j,2}(1));
            fprintf(fid,'%15.5g\t',MOS2{i,j,2}(2));
            fprintf(fid,'%15.5g\n',MOS2{i,j,2}(3));
        end
    end
end
% print S2 atom
for i=1:l
    for j=1:l
        if isempty(MOS2{i,j,3})==0
            fprintf(fid,'%s\t','S');
            fprintf(fid,'%15.5g\t',MOS2{i,j,3}(1));
            fprintf(fid,'%15.5g\t',MOS2{i,j,3}(2));
            fprintf(fid,'%15.5g\n',MOS2{i,j,3}(3));
        end
    end
end
fclose(fid);