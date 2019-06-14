clear; clc; close all;
gridfile='../run/grid.dat';
snapfile='../run/snap.dat';
gfuncfile='../run/gfunc.dat';
ngz=64;
ngr=32;
nfcoil=13;
ngz1=ngz+1;
ngr1=ngr+1;
SNAP=read_snap(snapfile,ngz1,ngr1);

% gfunc_fzr=read_gfunc(gfuncfile,nfcoil,ngz,ngr);
figure('Unit','normalized',...
    'Position',[0.0,0.0,0.8,0.8],...
    'DefaultAxesFontSize',20,...
    'DefaultAxesFontWeight','normal',...
    'DefaultAxesLineWidth',3,...
    'DefaultAxesTickLength',[0.013,0.03]);
subplot(2,1,1);
surf(SNAP.zgrid_zr,SNAP.rgrid_zr,SNAP.psi_zr);
shading interp; view([0,0,1]); colorbar;
xlabel('$Z/Z_w$','Interpreter','latex');
ylabel('$R/R_w$','Interpreter','latex');
subplot(2,1,2);
surf(SNAP.zgrid_zr,SNAP.rgrid_zr,SNAP.Jzeta_zr);
shading interp; view([0,0,1]); colorbar;
xlabel('$Z/Z_w$','Interpreter','latex');
ylabel('$R/R_w$','Interpreter','latex');

