% Create daily plot from daily streamflow
%

tic

close all
clear
clc

%% setting
% -------------------------------------------------------------------------
% Starting and ending time period for River network routed flow
% -------------------------------------------------------------------------
byr_sim=1950; bmon_sim=1; bday_sim=1;
eyr_sim=1950; emon_sim=12; eday_sim=31;

% Create daily array
dnum_day_sim = day_stamp(byr_sim,bmon_sim,bday_sim,eyr_sim,emon_sim,eday_sim,1);
[yrs_day_sim, mon_day_sim, day_day_sim]=datevec(dnum_day_sim);

% -------------------------------------------------------------------------
% metadata for  evaluating gauge 
% -------------------------------------------------------------------------
gauge(1).huc='14';
gauge(1).vic='CAMEO';
gauge(1).usgs='009095500';
gauge(1).usgsname='Colorado River near Cameo';
gauge(1).rseg=14000645;
gauge(1).gseg=923;
% ------------------------------------------------------------------------
% Directory
% -------------------------------------------------------------------------
Qrnw_path = '/home/mizukami/mizuRoute/route.v1/output';
fig_path  = '/home/mizukami/mizuRoute/route.v1/verification';

%% Getting stremflow output from mizuRoute
% -------------------------------------------------------------------------
% Daily streamflow for all segments
% -------------------------------------------------------------------------
ncname = [Qrnw_path '/q_out.nc'];
Qkwt_day = netcdf2mat(ncname,'KWTroutedRunoff');
Qirf_day = netcdf2mat(ncname,'IRFroutedRunoff');
segID_rn = netcdf2mat(ncname,'reachID');
clear ncname

Qkwt_day_gauge = ones(length(dnum_day_sim),length(gauge))*NaN;
Qirf_day_gauge = ones(length(dnum_day_sim),length(gauge))*NaN;
for g = 1:length(gauge)
    i1=find(segID_rn==gauge(g).rseg);
    Qkwt_day_gauge(:,g)=Qkwt_day(:,i1);
    Qirf_day_gauge(:,g)=Qirf_day(:,i1);
    clear i1
end
%% -------------------------------------------------------------------------
% Make plot
% -------------------------------------------------------------------------
hFig=figure('Visible','on');
%figure size setting
X=15.0;
Y=8.0;
xMargin = 0.25;
yMargin = 0.25;
xSize = X - 2*xMargin; 
ySize = Y - 2*yMargin;
set(hFig, 'Units','centimeters','Position',[0 0 xSize ySize])
set(hFig, 'PaperUnits','centimeters')
set(hFig, 'PaperSize',[X Y])
set(hFig, 'PaperPosition',[xMargin yMargin xSize ySize])

d1 = find(dnum_day_sim >= datenum(1950,1,1) & dnum_day_sim <= datenum(1950,12,31) );
q_day_max = max([max(Qkwt_day_gauge(d1,1)) max(Qirf_day_gauge(d1,1))]);
q_day_lims = [0 q_day_max*1.05];
plot(dnum_day_sim(d1),Qkwt_day_gauge(d1,1),'b-','linewidth',1);hold on;
plot(dnum_day_sim(d1),Qirf_day_gauge(d1,1),'r-','linewidth',1);
set(gca, 'fontname','arial','fontsize',9);
set(gca,'xlim',[min(dnum_day_sim(d1)) max(dnum_day_sim(d1))],'ylim',q_day_lims);
datetick('x','yy/mm','keeplimits');
ylabel('DMQ [cms]');
legend('KWT','IRF-UH','Orientation','vertical','Location','best');
legend('boxoff','Color','white');
title(gauge(1).usgsname,'units','normalized');

svstr = sprintf('%s/daily_hydrograph_cameo',fig_path);
print('-dtiff','-r300',svstr);
