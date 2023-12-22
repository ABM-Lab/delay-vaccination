
%%%%%%% Plot inci/hosp/death during the first 100 days %%%%%%%%%

clear all;
clc;

C0=[0.75 0.75 0.75];
C1=[0 0.45 0.74];
C2=[0.9 0.65 0.13];
C3=[0.49 0.18 0.56];
C4=[0.47 0.67 0.19];
C5=[0.64 0.08 0.18];
C6=[0.85 0.33 0.1];
C7=[200 0 100]/255;
C8=[0.3 0.5 0];
clr={C1,C2,C3};


TAUEND = 15;
TAU = 28:3.5:7*TAUEND; % waiting period for the 2nd dose: 1 weeks to 12 weeks
EPS1 = 0.3:0.025:0.6; % efficacy of the 1st dose against infection
EPS2 = 0.75:0.15:0.9; % 0.6:0.05:0.9; % efficacy of the 2nd dose against infection


m = length(TAU);
e1 = length(EPS1);
e2 = length(EPS2);
d1 = length(DLT1);
 % 
 load("R1_400.mat"); % waning parameter set I, R0=1.1, DSD period
 load("D1_400.mat");
 % 
 % load("RR100_R1400.mat"); % waning parameter set I, R0=1.1, first 100 days
 % load("RR100_D1400.mat");
 % 
 % load("R1_400_R18.mat"); % waning parameter set I, R0=1.8, DSD period
 % load("D1_400_R18.mat");
 % 
 % load("RR100_R1400_R18.mat"); % waning parameter set I, R0=1.8, first 100 days
 % load("RR100_D1400_R18.mat");
 %
%%

ee2 = 2; % fix epsilon2 = 0.75(1), 0.9(2)

% plot the relative reduction on inci/hosp/death in eps1 vs tau plane
RRhospP = zeros(e1, m);
RRdeathP = zeros(e1, m);
RRhospR= zeros(e1, m);
RRdeathR = zeros(e1, m);
RRinciP = zeros(e1, m);
RRinciR = zeros(e1, m);
RRinci = zeros(e1, m);
RRhosp = zeros(e1, m);
RRdeath = zeros(e1, m);

    for ee1 = 1:e1
        for mm = 1:m
            RRinci(ee1, mm) = (RinciP(ee1, ee2, mm)+RinciR(ee1, ee2, mm) - DinciP(ee1, ee2, mm)-DinciR(ee1, ee2, mm))./ (RinciP(ee1, ee2, mm)+RinciR(ee1, ee2, mm));
            RRhosp(ee1, mm) = (RhospP(ee1, ee2, mm)+RhospR(ee1, ee2, mm) - DhospP(ee1, ee2, mm)-DhospR(ee1, ee2, mm))./ (RhospP(ee1, ee2, mm)+RhospR(ee1, ee2, mm));
            RRdeath(ee1, mm) = (RdeathP(ee1, ee2, mm)+RdeathR(ee1, ee2, mm) - DdeathP(ee1, ee2, mm)-DdeathR(ee1, ee2, mm))./(RdeathP(ee1, ee2, mm)+RdeathR(ee1, ee2, mm));
            RRinciP(ee1, mm) = (RinciP(ee1, ee2, mm) - DinciP(ee1, ee2, mm))./ RinciP(ee1, ee2, mm);
            RRinciR(ee1, mm) = (RinciR(ee1, ee2, mm) - DinciR(ee1, ee2, mm))./ RinciR(ee1, ee2, mm);
            RRhospP(ee1, mm) = (RhospP(ee1, ee2, mm) - DhospP(ee1, ee2, mm))./ RhospP(ee1, ee2, mm);
            RRdeathP(ee1, mm) = (RdeathP(ee1, ee2, mm) - DdeathP(ee1, ee2, mm))./ RdeathP(ee1, ee2, mm);
            RRhospR(ee1, mm) = (RhospR(ee1, ee2, mm) - DhospR(ee1, ee2, mm))./ RhospR(ee1, ee2, mm);
            RRdeathR(ee1, mm) = (RdeathR(ee1, ee2, mm) - DdeathR(ee1, ee2, mm))./ RdeathR(ee1, ee2, mm);
        end
    end

%% PLOT of inciP / hospP / deathP

M1 = interpn(100*RRinciP, 6);
M2 = interpn(100*RRhospP, 6);
M3 = interpn(100*RRdeathP, 6);
% exportgraphics(gcf, 'R1400_d1075_e209_A9.pdf')

% figure()
hold on
subplot(1, 3, 1)
M = M1;
imagesc(M);
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-.','color',[0.5 0.5 0.5]); %%%% 'xcolor', 'k', 'ycolor', 'k');
n = 30;
% color map with red and blue on the edges
C = [1 0 0; 0 0.5 1];
% convert to HSV for interpolation
C_HSV = rgb2hsv(C);
% interpolate hue value
C_HSV_interp = interp1([0 n], C_HSV(:,1), 1:n);
% compose full HSV colormap
C_HSV = [C_HSV_interp(:), repmat(C_HSV(2:3), n, 1)];
% convert back to RGB
C = hsv2rgb(C_HSV);
% set colormap
colormap(C)
% colormap(flipud(C))
set(gca,'ydir','normal');
hh=colorbar;
hold on
[CC,h]=contour(M,[-6 -4 -2 -1  2 4],'linewidth', 1 ,'color','k');
CCL=clabel(CC,h,'FontSize',13);
hh.Ticks = [-6 -4 -2 0 2 4 ] ; %Create 8 ticks from zero to 1
hh.TickLabels = [{'-6'},{'-4'},{'-2'},{'0'}, {'2'}, {'4'}];
[CCh,Ch]=contour(M,[0 0],'linewidth', 4 ,'color',[204,76,2]/255);
CChL=clabel(CCh,Ch,'FontSize',15,'fontweight','bold','color',[204,76,2]/255);
clim([-6 4])
Ssize=size(M);
Vec1 = Ssize(2);
Vec2 = Ssize(1);
xlim([0 Vec1]);
ylim([0 Vec2]);

xTtick = (Vec1-1)/11;
xx=0:xTtick:Vec1; % when xlim = 769
xlab=1:1:12;
set(gca,'xtick',xx,'xticklabel',xlab);
xtickangle(0)
xlabel('Delay in second dose (weeks)')

yTtick=(Vec2-1)/6;
yy=0:yTtick:Vec2;
ylab=30:5:60;
set(gca,'ytick',yy,'yticklabel',ylab);
ylabel('Efficacy of first dose against infection (%)')

%clim([-14 6])
ylabel(hh, 'RR of incidence for primary infection (%)','FontSize',15)
set(gca,'fontsize',15);
%hh.Ticks = linspace(1, 6, 6) ; %Create 6 ticks 
%hh.TickLabels = num2cell([90 80 70 60 50 40]);
title('\rmD','position',[20, Vec2+10],'fontsize',20)
s =  [200 200 1500 400];
set(gcf,'Position',s)

hold on
subplot(1, 3, 2)
M = M2;
imagesc(M);
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-.','color',[0.5 0.5 0.5]); %%%% 'xcolor', 'k', 'ycolor', 'k');
n = 30;
% color map with red and blue on the edges
C = [1 0 0; 0 0.5 1];
% convert to HSV for interpolation
C_HSV = rgb2hsv(C);
% interpolate hue value
C_HSV_interp = interp1([0 n], C_HSV(:,1), 1:n);
% compose full HSV colormap
C_HSV = [C_HSV_interp(:), repmat(C_HSV(2:3), n, 1)];
% convert back to RGB
C = hsv2rgb(C_HSV);
% set colormap
colormap(C)
% colormap(flipud(C))
set(gca,'ydir','normal');
hh=colorbar;
hold on
[CC,h]=contour(M,[-0.5  2 4 8 12 13 14],'linewidth', 1 ,'color','k');
CCL=clabel(CC,h,'FontSize',13)
hh.Ticks = [-1 0 4 8 12 14] ; %Create 8 ticks from zero to 1
hh.TickLabels = [{'-1'},{'0'},{'4'},{'8'},{'12'},{'14'}];
[CCh,Ch]=contour(M,[0 0],'linewidth', 4 ,'color',[204,76,2]/255);
CChL=clabel(CCh,Ch,'FontSize',15,'fontweight','bold','color',[204,76,2]/255);
clim([-1 14])
Ssize=size(M);
Vec1 = Ssize(2);
Vec2 = Ssize(1);
xlim([0 Vec1]);
ylim([0 Vec2]);

xTtick = (Vec1-1)/11;
xx=0:xTtick:Vec1; % when xlim = 769
xlab=1:1:12;
set(gca,'xtick',xx,'xticklabel',xlab);
xtickangle(0)
xlabel('Delay in second dose (weeks)')

yTtick=(Vec2-1)/6;
yy=0:yTtick:Vec2;
ylab=30:5:60;
set(gca,'ytick',yy,'yticklabel',ylab);
ylabel('Efficacy of first dose against infection (%)')

%clim([-4 12])
ylabel(hh, 'RR of hospitalization for primary infection (%)','FontSize',15)
set(gca,'fontsize',15);
%hh.Ticks = linspace(1, 6, 6) ; %Create 6 ticks 
%hh.TickLabels = num2cell([90 80 70 60 50 40]);
title('\rmE','position',[20, Vec2+10],'fontsize',20)
s =  [200 200 1500 400];
set(gcf,'Position',s)


subplot(1, 3, 3)
M = M3;
imagesc(M);
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-.','color',[0.5 0.5 0.5]); %%%% 'xcolor', 'k', 'ycolor', 'k');
n = 30;
% color map with red and blue on the edges
C = [1 0 0; 0 0.5 1];
% convert to HSV for interpolation
C_HSV = rgb2hsv(C);
% interpolate hue value
C_HSV_interp = interp1([0 n], C_HSV(:,1), 1:n);
% compose full HSV colormap
C_HSV = [C_HSV_interp(:), repmat(C_HSV(2:3), n, 1)];
% convert back to RGB
C = hsv2rgb(C_HSV);
% set colormap
colormap(C)
% colormap(flipud(C))
set(gca,'ydir','normal');
hh=colorbar;
hold on
[CC,h]=contour(M,[-1 -0.5 2 4 6 8 10 12 14 16],'linewidth', 1 ,'color','k');
CCL=clabel(CC,h,'FontSize',13)
hh.Ticks = [-1 0 4 8 12 16] ; %Create 8 ticks from zero to 1
hh.TickLabels = [{'-1'},{'0'},{'4'},{'8'},{'12'},{'16'}];
[CCh,Ch]=contour(M,[0 0],'linewidth', 4 ,'color',[204,76,2]/255);
CChL=clabel(CCh,Ch,'FontSize',15,'fontweight','bold','color',[204,76,2]/255);
clim([-1 16])
Ssize=size(M);
Vec1 = Ssize(2);
Vec2 = Ssize(1);
xlim([0 Vec1]);
ylim([0 Vec2]);

xTtick = (Vec1-1)/11;
xx=0:xTtick:Vec1; % when xlim = 769
xlab=1:1:12;
set(gca,'xtick',xx,'xticklabel',xlab);
xtickangle(0)
xlabel('Delay in second dose (weeks)')

yTtick=(Vec2-1)/6;
yy=0:yTtick:Vec2;
ylab=30:5:60;
set(gca,'ytick',yy,'yticklabel',ylab);
ylabel('Efficacy of first dose against infection (%)')

%clim([-2 12])
ylabel(hh, 'RR of death for primary infection (%)','FontSize',15)
set(gca,'fontsize',15);
%hh.Ticks = linspace(1, 6, 6) ; %Create 6 ticks 
%hh.TickLabels = num2cell([90 80 70 60 50 40]);
title('\rmF','position',[20, Vec2+10],'fontsize',20)
s =  [200 200 1500 400];
set(gcf,'Position',s)

exportgraphics(gcf,'Sum100IncHospDed400_e90.pdf')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
ee2 = 1; % fix epsilon2 = 0.75(1), 0.9(2)

% plot the relative reduction on inci/hosp/death in eps1 vs tau plane
RRhospP = zeros(e1, m);
RRdeathP = zeros(e1, m);
RRhospR= zeros(e1, m);
RRdeathR = zeros(e1, m);
RRinciP = zeros(e1, m);
RRinciR = zeros(e1, m);
RRinci = zeros(e1, m);
RRhosp = zeros(e1, m);
RRdeath = zeros(e1, m);

    for ee1 = 1:e1
        for mm = 1:m
            RRinci(ee1, mm) = (RinciP(ee1, ee2, mm)+RinciR(ee1, ee2, mm) - DinciP(ee1, ee2, mm)-DinciR(ee1, ee2, mm))./ (RinciP(ee1, ee2, mm)+RinciR(ee1, ee2, mm));
            RRhosp(ee1, mm) = (RhospP(ee1, ee2, mm)+RhospR(ee1, ee2, mm) - DhospP(ee1, ee2, mm)-DhospR(ee1, ee2, mm))./ (RhospP(ee1, ee2, mm)+RhospR(ee1, ee2, mm));
            RRdeath(ee1, mm) = (RdeathP(ee1, ee2, mm)+RdeathR(ee1, ee2, mm) - DdeathP(ee1, ee2, mm)-DdeathR(ee1, ee2, mm))./(RdeathP(ee1, ee2, mm)+RdeathR(ee1, ee2, mm));
            RRinciP(ee1, mm) = (RinciP(ee1, ee2, mm) - DinciP(ee1, ee2, mm))./ RinciP(ee1, ee2, mm);
            RRinciR(ee1, mm) = (RinciR(ee1, ee2, mm) - DinciR(ee1, ee2, mm))./ RinciR(ee1, ee2, mm);
            RRhospP(ee1, mm) = (RhospP(ee1, ee2, mm) - DhospP(ee1, ee2, mm))./ RhospP(ee1, ee2, mm);
            RRdeathP(ee1, mm) = (RdeathP(ee1, ee2, mm) - DdeathP(ee1, ee2, mm))./ RdeathP(ee1, ee2, mm);
            RRhospR(ee1, mm) = (RhospR(ee1, ee2, mm) - DhospR(ee1, ee2, mm))./ RhospR(ee1, ee2, mm);
            RRdeathR(ee1, mm) = (RdeathR(ee1, ee2, mm) - DdeathR(ee1, ee2, mm))./ RdeathR(ee1, ee2, mm);
        end
    end

%% PLOT of inciP / hospP / deathP

M1 = interpn(100*RRinciP, 6);
M2 = interpn(100*RRhospP, 6);
M3 = interpn(100*RRdeathP, 6);
% exportgraphics(gcf, 'R1400_d1075_e209_A9.pdf')

% figure()
hold on
subplot(1, 3, 1)
M = M1;
imagesc(M);
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-.','color',[0.5 0.5 0.5]); %%%% 'xcolor', 'k', 'ycolor', 'k');
n = 30;
% color map with red and blue on the edges
C = [1 0 0; 0 0.5 1];
% convert to HSV for interpolation
C_HSV = rgb2hsv(C);
% interpolate hue value
C_HSV_interp = interp1([0 n], C_HSV(:,1), 1:n);
% compose full HSV colormap
C_HSV = [C_HSV_interp(:), repmat(C_HSV(2:3), n, 1)];
% convert back to RGB
C = hsv2rgb(C_HSV);
% set colormap
colormap(C)
% colormap(flipud(C))
set(gca,'ydir','normal');
hh=colorbar;
hold on
[CC,h]=contour(M,[-6 -4 -2 -1 2 4],'linewidth', 1 ,'color','k');
CCL=clabel(CC,h,'FontSize',13)
hh.Ticks = [-6 -4  -2 0  2 4] ; %Create 8 ticks from zero to 1
hh.TickLabels = [{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'}];
[CCh,Ch]=contour(M,[0 0],'linewidth', 4 ,'color',[204,76,2]/255);
CChL=clabel(CCh,Ch,'FontSize',15,'fontweight','bold','color',[204,76,2]/255);
clim([-6 4])
Ssize=size(M);
Vec1 = Ssize(2);
Vec2 = Ssize(1);
xlim([0 Vec1]);
ylim([0 Vec2]);

xTtick = (Vec1-1)/11;
xx=0:xTtick:Vec1; % when xlim = 769
xlab=1:1:12;
set(gca,'xtick',xx,'xticklabel',xlab);
xtickangle(0)
xlabel('Delay in second dose (weeks)')

yTtick=(Vec2-1)/6;
yy=0:yTtick:Vec2;
ylab=30:5:60;
set(gca,'ytick',yy,'yticklabel',ylab);
ylabel('Efficacy of first dose against infection (%)')

%clim([-10 12])
ylabel(hh, 'RR of incidence for primary infection (%)','FontSize',15)
set(gca,'fontsize',15);
%hh.Ticks = linspace(1, 6, 6) ; %Create 6 ticks 
%hh.TickLabels = num2cell([90 80 70 60 50 40]);
title('\rmA','position',[20, Vec2+10],'fontsize',20)
s =  [200 200 1500 400];
set(gcf,'Position',s)

hold on
subplot(1, 3, 2)
M = M2;
imagesc(M);
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-.','color',[0.5 0.5 0.5]); %%%% 'xcolor', 'k', 'ycolor', 'k');
n = 30;
% color map with red and blue on the edges
C = [1 0 0; 0 0.5 1];
% convert to HSV for interpolation
C_HSV = rgb2hsv(C);
% interpolate hue value
C_HSV_interp = interp1([0 n], C_HSV(:,1), 1:n);
% compose full HSV colormap
C_HSV = [C_HSV_interp(:), repmat(C_HSV(2:3), n, 1)];
% convert back to RGB
C = hsv2rgb(C_HSV);
% set colormap
colormap(C)
% colormap(flipud(C))
set(gca,'ydir','normal');
hh=colorbar;
hold on
[CC,h]=contour(M,[-0.5 2 4 8 12 14],'linewidth', 1 ,'color','k');
CCL=clabel(CC,h,'FontSize',13)
hh.Ticks = [-1 0 4 8 12 14] ; %Create 8 ticks from zero to 1
hh.TickLabels = [{'-1'},{'0'},{'4'},{'8'},{'12'},{'14'}];
[CCh,Ch]=contour(M,[0 0],'linewidth', 4 ,'color',[204,76,2]/255);
CChL=clabel(CCh,Ch,'FontSize',15,'fontweight','bold','color',[204,76,2]/255);
clim([-1 14])
Ssize=size(M);
Vec1 = Ssize(2);
Vec2 = Ssize(1);
xlim([0 Vec1]);
ylim([0 Vec2]);

xTtick = (Vec1-1)/11;
xx=0:xTtick:Vec1; % when xlim = 769
xlab=1:1:12;
set(gca,'xtick',xx,'xticklabel',xlab);
xtickangle(0)
xlabel('Delay in second dose (weeks)')

yTtick=(Vec2-1)/6;
yy=0:yTtick:Vec2;
ylab=30:5:60;
set(gca,'ytick',yy,'yticklabel',ylab);
ylabel('Efficacy of first dose against infection (%)')

%clim([-2 18])
ylabel(hh, 'RR of hospitalization for primary infection (%)','FontSize',15)
set(gca,'fontsize',15);
%hh.Ticks = linspace(1, 6, 6) ; %Create 6 ticks 
%hh.TickLabels = num2cell([90 80 70 60 50 40]);
title('\rmB','position',[20, Vec2+10],'fontsize',20)
s =  [200 200 1500 400];
set(gcf,'Position',s)


subplot(1, 3, 3)
M = M3;
imagesc(M);
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-.','color',[0.5 0.5 0.5]); %%%% 'xcolor', 'k', 'ycolor', 'k');
n = 30;
% color map with red and blue on the edges
C = [1 0 0; 0 0.5 1];
% convert to HSV for interpolation
C_HSV = rgb2hsv(C);
% interpolate hue value
C_HSV_interp = interp1([0 n], C_HSV(:,1), 1:n);
% compose full HSV colormap
C_HSV = [C_HSV_interp(:), repmat(C_HSV(2:3), n, 1)];
% convert back to RGB
C = hsv2rgb(C_HSV);
% set colormap
colormap(C)
% colormap(flipud(C))
set(gca,'ydir','normal');
hh=colorbar;
hold on
[CC,h]=contour(M,[-0.5 2 6 8 10 12 14 16],'linewidth', 1 ,'color','k');
CCL=clabel(CC,h,'FontSize',13)
hh.Ticks = [-1 0 2 4 8 12 16] ; %Create 8 ticks from zero to 1
hh.TickLabels = [{'-1'},{'0'},{'2'},{'4'},{'8'},{'12'},{'16'}];
[CCh,Ch]=contour(M,[0 0],'linewidth', 4 ,'color',[204,76,2]/255);
CChL=clabel(CCh,Ch,'FontSize',15,'fontweight','bold','color',[204,76,2]/255);
clim([-1 16])
Ssize=size(M);
Vec1 = Ssize(2);
Vec2 = Ssize(1);
xlim([0 Vec1]);
ylim([0 Vec2]);

xTtick = (Vec1-1)/11;
xx=0:xTtick:Vec1; % when xlim = 769
xlab=1:1:12;
set(gca,'xtick',xx,'xticklabel',xlab);
xtickangle(0)
xlabel('Delay in second dose (weeks)')

yTtick=(Vec2-1)/6;
yy=0:yTtick:Vec2;
ylab=30:5:60;
set(gca,'ytick',yy,'yticklabel',ylab);
ylabel('Efficacy of first dose against infection (%)')

%clim([-2 18])
ylabel(hh, 'RR of death for primary infection (%)','FontSize',15)
set(gca,'fontsize',15);
%hh.Ticks = linspace(1, 6, 6) ; %Create 6 ticks 
%hh.TickLabels = num2cell([90 80 70 60 50 40]);
title('\rmC','position',[20, Vec2+10],'fontsize',20)
s =  [200 200 1500 400];
set(gcf,'Position',s)

exportgraphics(gcf,'Sum100IncHospDed400_e75.pdf')