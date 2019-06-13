function [rgrid_zr, zgrid_zr]=read_grid(filenamestr,ngz,ngr)
%% read the grid data from grid file
% -inputs:
% filenamestr   string, the full path of the grid file
% ngz           double, number of grid in z direction
% ngr           double, number of grid in r direction
% -outputs:
% rgrid_zr      ngz*ngr array, the R coordinate of grid points
% zgrid_zr      ngz*ngr array, the Z coordinate of grid points

% Edited by Shuying SUN in 2019/06/12
% Contact: sunshuyingc@enn.cn, bangrsun@163.com
% ENN Sci. & Tech. Development Corporation, 2008-2019
% (c) All rights reserved.

fid=fopen(filenamestr,'r');
if(fid<=0)
    error(['Can not open file:', filenamestr]);
end
rgrid_zr=zeros(ngz,ngr);
zgrid_zr=zeros(ngz,ngr);
for j=1:ngz
    rgrid_zr(j,:)=fscanf(fid,'%f\n',ngr);
end
for j=1:ngz
    zgrid_zr(j,:)=fscanf(fid,'%f\n',ngr);
end

end