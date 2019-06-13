clear; clc; close all;
gridfile='../run/grid.dat';
snapfile='../run/snap.dat';
ngz=90;
ngr=65;
SNAP=read_snap(snapfile,ngz,ngr);

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

