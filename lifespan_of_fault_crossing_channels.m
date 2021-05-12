% Find below the script used to produce figures 2 and 3 of "The lifespan of
% fault-crossing channels". The file DATA.txt must be in the same directory

% Author: Kelian Dascher-Cousineau
% Last edited: May 2021
% Written and tested with Matlab 2020a (academic liscence)


%% useful plotting functions:
ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);
yr2sec = @(yrs) yrs*60*60*24*365;

%% analysis

% Load data
dataTbl         = readtable('DATA.csv');

% Data removed in revision because the mechanism and interpretation of
% offset were not sufficiently straighforward.
dataTbl         = dataTbl(~isnan(dataTbl.z1),:);

Ndata           = height(dataTbl);

% Diffusivity calculated following Paola 1992 (see supplement section S3)
cf              = 0.01; % drag coefficient
Co              = 0.7;  % sediment concentration
s               = 2.7;  % specific density of sediment (unitless)
r               = 0.05; % m/yr rainfall Spotila et al., 2007
epsilon         = 0.4;  % unitless
A               = (epsilon/(1+epsilon))^(3/2); % channel geometry
q               = r*dataTbl.Reach;             % Basin River flux (m^2/s) (width normalized)

kappa           = 8*A*sqrt(cf)/(Co*(s-1)) * q / yr2sec(1); % diffusivity (m^2/s)
dataTbl.kappa   = kappa;

% Slip rate from Sieh & Jahn, 1984
dataTbl.vx      = 3.5*10^-2/yr2sec(1) ...       % m/s   slip rate
                    * ones(Ndata,1);

% Calculate the avulsion time scale (equation 3)
Tc = (4*(dataTbl.hc).^2)./(dataTbl.kappa.*(dataTbl.S0).^2);

% Calulate the ratio d_obs/dc (equation 4)
Tnorm = dataTbl.offset./dataTbl.vx./Tc;
                 

%% Figure 2

figure;  I = dataTbl.IsActive == 1;
sz = 20;
wt = 1;
al = 0.8;

r0 =[239423.969,3910925.726];   % This is the origin point to calculate distance along strike
rc =[dataTbl.xc,dataTbl.yc];
dr =sqrt(sum((ones(length(rc),1)*r0-rc).^2,2))./1000; % Distance along strike (in km)

% Offset, d_obs [m]
subplot(5,1,1); hold on;
scatter(dr(I), dataTbl.offset(I), sz,assign_TLC(dataTbl.stage(I)), 'filled','Markerfacealpha',al);
scatter(dr(~I),dataTbl.offset(~I),sz,assign_TLC(dataTbl.stage(~I)),'LineWidth',wt);
set(gca,'yscale','log')
xticks([])
ylabel({'Offset, d (m)'})

% Avulsion Threshold height, h_c [m]
subplot(5,1,2); hold on
scatter(dr(I), dataTbl.hc(I), sz,assign_TLC(dataTbl.stage(I)), 'filled','Markerfacealpha',al);
scatter(dr(~I),dataTbl.hc(~I),sz,assign_TLC(dataTbl.stage(~I)),'LineWidth',wt);
set(gca,'yscale','log')
xticks([])
ylabel({'Threshold','height, h_c (m)'})

% Initial slope, S_0 []
subplot(5,1,3); hold on
scatter(dr(I), dataTbl.S0(I), sz,assign_TLC(dataTbl.stage(I)), 'filled','Markerfacealpha',al);
scatter(dr(~I),dataTbl.S0(~I),sz,assign_TLC(dataTbl.stage(~I)),'LineWidth',wt);
xticks([])
ylabel({'Slope, S_o'})

% Catchement length [m]
subplot(5,1,4); hold on
scatter(dr(I), dataTbl.Reach(I), sz,assign_TLC(dataTbl.stage(I)), 'filled','Markerfacealpha',al);
scatter(dr(~I),dataTbl.Reach(~I),sz,assign_TLC(dataTbl.stage(~I)),'LineWidth',wt);
set(gca,'yscale','log')
xticks([])
ylabel({'Reach, L (m)'})


% d_obs/d_c []
subplot(5,1,5); hold on;
scatter(dr(I),Tnorm(I),sz, ...
    assign_TLC(dataTbl.stage(I)),'filled','Markerfacealpha',al);
scatter(dr(~I),Tnorm(~I),sz, ...
    assign_TLC(dataTbl.stage(~I)),'LineWidth',wt);
set(gca,'yscale','log','ytick',[1,100])
xlabel('Distance along strike (km)')
ylabel('^{d_{obs}}/_{d_c}')

set(gca,'YMinorTick','on')
set(findall(gcf,'-property','Fontsize'),'Fontsize',10)
setsize(gcf,3.5,6)


%% Figure 3

sz = 50;
I = Tnorm > 0;

% Fit logistic regression through the categorical data (Active/abandoned)
[B,sdev, stats]     = mnrfit(log10(Tnorm(I)),categorical(~dataTbl.IsActive(I)));
logisticReg         = @(x,b0,b1) 1./(1+exp(-(b0+b1*x)));

% 1 sigma Confidence intervals are computed directly from the logistic regression
prctRange = [50-66/2,50+66/2]/100;
interval = (1-prctRange)./prctRange;

TcNormFit = 10^(-B(1)/B(2));
TcConfInt = 10.^((log(interval)-B(1))/B(2));

disp([num2str(TcNormFit),' best separates active and abandoned channels'])
disp(['where p = ',num2str(stats.p(2)), ' is the probability that B(2) 0 (no division in data)'])

% measurements of Tc derived from incipitent or recent avulsions (yellow
% points)
Iyellow = Tnorm > 0 & strcmp(dataTbl.stage,'yellow');

figure; hold on

% Confidence interval
rh = rectangle('Position',[TcConfInt(1),0,diff(TcConfInt),1],...
          'FaceColor',[0.2 0.2 0.2 0.2], ...
          'EdgeColor','none');

% Data 
pth = scatter((Tnorm(I)),dataTbl.IsActive(I),sz,assign_TLC(dataTbl.stage(I)), ...
                'Filled', ...
                'MarkerFaceAlpha',0.7);

% Fit
tNormVecLog = linspace(min(log10(Tnorm(I))),max(log10(Tnorm(I))),100);
lh = plot(10.^tNormVecLog,logisticReg(tNormVecLog,B(1),B(2)),'--k','linewidth',2);
      
xlabel('d_{obs}/d_c')   
set(gca,'xscale','log',...
        'ytick',[0,1],...
        'yticklabel',{'Abandoned','Active'}, ...
        'xtick',[10^-2,10^0,10^2])

% Incipient/recent avulsions
yyaxis right
edges = logspace(log10(min(Tnorm(Iyellow))-eps(1)),log10(max(Tnorm(Iyellow))), 10);
hh = histogram(Tnorm(Iyellow),edges, ...
    'facecolor', [0.9290 0.6940 0.1250], ...
    'edgecolor','none', ...
    'facealpha',0.4);
YLIM = ylim; 
ylim([YLIM(1),YLIM(2)*2])

set(gca,'Xscale','log')

ylabel({'Recent/Incipient','Avulsion'})
grid
% legend(lh,['1/(1+e^{-(\beta_0+\beta_1 x)}', newline, 'p(\beta_1>0)=', num2str(stats.p(2),1)])

set(findall(gcf,'-property','Fontsize'),'Fontsize',10)
setsize(gcf,3.5,1.5)

%%
function c = assign_TLC(color)

for n = 1:length(color)
switch color{n}
    case 'green'
        c(n,:) = [0.4660 0.6740 0.1880];
        
    case 'yellow'
        c(n,:) = [0.9290 0.6940 0.1250];
        
    case 'red'
        c(n,:) = [0.6350 0.0780 0.1840];
end   
end

end
