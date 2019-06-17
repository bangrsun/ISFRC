clear; clc; close all;
gridfile='../run/grid.dat';
snapfile='../run/snap.dat';
gfuncfile='../run/gfunc.dat';
cn=[0.3,50,0.0];
ngz=64;
ngr=32;
nfcoil=13;
ngz1=ngz+1;
ngr1=ngr+1;
SNAP=read_snap(snapfile,ngz1,ngr1);
pres_zr=cn(1).*sech(cn(2).*(SNAP.psi_zr-cn(3)));

% gfunc_fzr=read_gfunc(gfuncfile,nfcoil,ngz,ngr);
figure('Unit','normalized',...
    'Position',[0.0,0.0,0.9,0.7],...
    'DefaultAxesFontSize',20,...
    'DefaultAxesFontWeight','normal',...
    'DefaultAxesLineWidth',3,...
    'DefaultAxesTickLength',[0.013,0.03]);
subplot(2,2,1);
surf(SNAP.zgrid_zr,SNAP.rgrid_zr,SNAP.psi_zr);
shading interp; view([0,0,1]); colorbar;
xlabel('$Z/Z_w$','Interpreter','latex');
ylabel('$R/R_w$','Interpreter','latex');
subplot(2,2,2);
surf(SNAP.zgrid_zr,SNAP.rgrid_zr,SNAP.psip_zr);
shading interp; view([0,0,1]); colorbar;
xlabel('$Z/Z_w$','Interpreter','latex');
ylabel('$R/R_w$','Interpreter','latex');
subplot(2,2,3);
surf(SNAP.zgrid_zr,SNAP.rgrid_zr,SNAP.pprim_zr);
shading interp; view([0,0,1]); colorbar;
xlabel('$Z/Z_w$','Interpreter','latex');
ylabel('$R/R_w$','Interpreter','latex');
subplot(2,2,4);
surf(SNAP.zgrid_zr,SNAP.rgrid_zr,pres_zr);
shading interp; view([0,0,1]); colorbar;
xlabel('$Z/Z_w$','Interpreter','latex');
ylabel('$R/R_w$','Interpreter','latex');
