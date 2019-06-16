function A=read_snap(filenamestr,ngz,ngr)
%% read the grid data from snapshot file
% -inputs:
% filenamestr   string, the full path of the snapshot file
% ngz           double, number of grid in z direction
% ngr           double, number of grid in r direction
% -outputs:
% A             struct, the data passed from snapshot file

% Edited by Shuying SUN in 2019/06/12
% Contact: sunshuyingc@enn.cn, bangrsun@163.com
% ENN Sci. & Tech. Development Corporation, 2008-2019
% (c) All rights reserved.

fid=fopen(filenamestr,'r');
if(fid<=0)
    error(['Can not open file:', filenamestr]);
end
A={};
A.rgrid_zr=zeros(ngz,ngr);
A.zgrid_zr=zeros(ngz,ngr);
A.psi_zr=zeros(ngz,ngr);
A.Jzeta_zr=zeros(ngz,ngr);
for j=1:ngz
    A.rgrid_zr(j,:)=fscanf(fid,'%f\n',ngr);
end
for j=1:ngz
    A.zgrid_zr(j,:)=fscanf(fid,'%f\n',ngr);
end
for j=1:ngz
    A.psi_zr(j,:)=fscanf(fid,'%f\n',ngr);
end
for j=1:ngz
    A.psif_zr(j,:)=fscanf(fid,'%f\n',ngr);
end
for j=1:ngz
    A.psip_zr(j,:)=fscanf(fid,'%f\n',ngr);
end
for j=1:ngz
    A.pprim_zr(j,:)=fscanf(fid,'%f\n',ngr);
end

end