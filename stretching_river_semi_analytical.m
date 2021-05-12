function tc = stretching_river_semi_analytical(S0,reach,r,v,w0,hc)

yr2sec  = @(yr) yr*(60*60*24*365);

if nargin == 0
    example = 8; % toggle through these examples (1-8) to get a walkthrough of the approach
    switch example
        case 1
            %% Example #1: Heaviside function

            % define the intial shape of the river: 
            x = 0:0.01:10;
            z = x<max(x)/2;
            N = 1000; % This defines the number of sinusoidal terms in the fourier expansion
            tarr = linspace(0,1,100);
            k = 1;

            figure; hold on
            ph0 = plot(x,z,'Color',[1 0 0]);
            for ti = tarr
                zi = diffuse(x,z,k,ti,N);
                plot(x,zi,'Color',[1-ti/max(tarr),0,ti/max(tarr)])    
            end
            legend(ph0,'Initial profile')
            xlabel('Along profile distance')
            ylabel('Elevation')

        case 2
            %% Example #2: Dam removal/normal fault scarp

            S0 = 0.1; % here we define an initial slope
            dx = 0.01; % spacing
            L  = 10; 
            Dam_height = 5;

            [x,z] = simple_profile(S0,dx,L);
            z = z + Dam_height*(x<max(x)/2);
            N = 500; % This defines the number of sinusoidal terms in the fourier expansion
            tarr = linspace(0,1,100);
            k = 1;

            figure; hold on
            ph0 = plot(x,z,'Color',[1 0 0]);
            for ti = tarr
                zi = diffuse(x,z,k,ti,N);
                plot(x,zi,'Color',[1-ti/max(tarr),0,ti/max(tarr)])    
            end
            legend(ph0,'Initial profile')
            xlabel('Along profile distance')
            ylabel('Elevation')


        case 3
            %% Example #3: Fault one offset
            S0 = 0.1; % here we define an initial slope
            dx = 0.01; % spacing
            L  = 10; 

            [x0,z0] = simple_profile(S0,dx,L);

            % here I introduce the strech function (see below), which stretches
            % an input profile by a factor (stretch), over a width (w) centered
            % on a point (xo) and resamples the stretched interval at a spacing
            % of dx
            offset = 10;
            fault_zone_width = 1;
            stretch_factor = sqrt(fault_zone_width^2 + offset^2);
            [x,z] = stretch(x0,z0,dx,fault_zone_width,max(x0)/2,stretch_factor);

            N = 1000; % This defines the number of sinusoidal terms in the fourier expansion
            tarr = linspace(0,10,5);
            k = 1;

            figure; hold on
            ph0 = plot(x,z,'Color',[1 0 0]);
            for ti = tarr
                zi = diffuse(x,z,k,ti,N);
                plot(x,zi,'Color',[1-ti/max(tarr),0,ti/max(tarr)])    
            end
            legend(ph0,'Initial profile')
            xlabel('Along profile distance')
            ylabel('Elevation')

        case 4
            %% Example #4: Fault multiple constant offsets with constant recurrence
            S0 = 0.1; % here we define an initial slope
            dx = 0.01; % spacing
            L  = 10; 

            [x,z] = simple_profile(S0,dx,L);

            % diffusive terms
            N = 1000; % This defines the number of sinusoidal terms in the fourier expansion
            k = 0.5;

            nEq = 10;
            recurranceInterval = 1*ones(1,nEq); % population of recurrence intervals
            offsetArray = 1*ones(1,nEq); % population of slips (charecteristic in this case)
            w0 = 1; % width of the fault zone
            hc = 0.1; % possible to make this related to water level
            avulseYN = false;

            tc = eq_cycle_channel_model(x,z,dx,w0, ...                         % Geometry
                                        k, N,      ...                         % Diffusion
                                        recurranceInterval, offsetArray, ...   % Earthquakes
                                        hc, avulseYN);                        % Channel

        case 5
            %% Example #5: Fault with mulitple offsets which can avulse
            S0 = 0.1; % here we define an initial slope
            dx = 0.01; % spacing
            L  = 10; 

            [x,z] = simple_profile(S0,dx,L);

            % diffusive terms
            N = 100; % This defines the number of sinusoidal terms in the fourier expansion
            k = 1;

            nEq = 100;
            v = 10;
            recurranceInterval = 1*ones(1,nEq);          % population of recurrence intervals
            offsetArray = v*(recurranceInterval);        % population of slips (charecteristic in this case)
            w0 = 1;   % width of the fault zone
            hc = 0.1; % possible to make this related to water level
            avulseYN = true;
            showLive = false;

            tc = eq_cycle_channel_model(x,z,dx,w0, ...                         % Geometry
                                        k, N,      ...                         % Diffusion
                                        recurranceInterval, offsetArray, ...   % Earthquakes
                                        hc, avulseYN, ...                    % Channel
                                        showLive);

        case 6 
            %% Example #6: Fault with mulitple random offsets and recurrance intervals can avulse
            S0 = 0.1; % here we define an initial slope
            dx = 0.01; % spacing
            L  = 10; 

            [x,z] = simple_profile(S0,dx,L);

            % diffusive terms
            N = 1000;                               % This defines the number of sinusoidal terms in the fourier expansion
            k = 0.05;

            nEq = 100;

            recurranceInterval = 2*10*rand(1,nEq);  % population of recurrence intervals
            offsetArray = 2*0.1*rand(1,nEq);        % population of slips (charecteristic in this case)
            w0 = 1;                                 % width of the fault zone
            hc = 0.1;                               % possible to make this related to water level
            avulseYN = true;
            showLive = true;

            tc = eq_cycle_channel_model(x,z,dx,w0, ...                         % Geometry
                                        k, N,      ...                         % Diffusion
                                        recurranceInterval, offsetArray, ...   % Earthquakes
                                        hc, avulseYN, ...                    % Channel
                                        showLive);

        case 7
            %% Example #7: Fault with "shutter ridges" veritical offset (hc(t))

            S0 = 0.1; % here we define an initial slope
            dx = 0.01; % spacing
            L  = 10; 

            [x,z] = simple_profile(S0,dx,L);

            % diffusive terms
            N = 1000; % This defines the number of sinusoidal terms in the fourier expansion
            k = 0.05*10^-7; 

            nEq = 500;
            v = 0.01;
            recurranceInterval = 10*ones(1,nEq);          % population of recurrence intervals
            offsetArray = v*(recurranceInterval);        % population of slips (charecteristic in this case)
            w0 = 1; % width of the fault zone
            hc = 0.1; % possible to make this related to water level
            avulseYN = true;

            shutterRidge_sp = 10;
            hcfh = @(t) 2*hc + hc*sin(2*pi*t/(shutterRidge_sp/v)); % e.g. shutter ridges
            showLive = true;

            tc = eq_cycle_channel_model(x,z,dx,w0, ...                         % Geometry
                                        k, N,      ...                         % Diffusion
                                        recurranceInterval, offsetArray, ...   % Earthquakes
                                        hcfh, avulseYN, ...                    % Channel
                                        showLive);
        case 8
            %% Example #8: Exploring realistic parameters

            S0      = 0.05; % slope (Rise over run)
            dx      = 1;    % point spacing (m) 
            L       = 3000; % domain size (m)

            r       = 0.05; % annual rainfall (m)

            [x,z]   = simple_profile(S0,dx,L);
            N       = 100;  % order of fourier expansion
            nEq     = 20; 
            vcmyr   = 3.3;  % slip velocity (cm/year)
            v       = vcmyr / 100 /yr2sec(1);   % m/s
            k       = 0.1*r*L/2/yr2sec(1);      % m^2/s

            recurranceInterval  = yr2sec(500)*ones(1,nEq);
            offsetArray         = v*(recurranceInterval);


            w0              = 3;        % width of the fault zone (m)
            hc              = 5;        % height to avulse (loosely based on the channel height at wallace creek
            avulseYN        = true;    % simulate avulsions
            shutterRidgeSp  = 100;      % shutter ridge spacing (m)
            shutterRidgeHeight = 0.00;  % add amplitude here to simulate shutter ridge
            upliftRate      = 0.000*v/yr2sec(1); % (m/s) add uplift here to simulate uplift
            hcfh            = @(t) ...     
                hc + shutterRidgeHeight*abs(sin(2*pi*v*t/(2*shutterRidgeSp))) + upliftRate*t; % hill-looking thing + uplift
            showLive        = true;
            tc = eq_cycle_channel_model(x,z,dx,w0, ...                         % Geometry
                                        k, N,      ...                         % Diffusion
                                        recurranceInterval, offsetArray, ...   % Earthquakes
                                        hcfh, avulseYN, ...                    % Channel
                                        showLive);

    end

elseif nargin == 6
    L = 2*reach;
    dx      = L/1000;
    N       = 100;
    [x,z]   = simple_profile(S0,dx,L);
    k       = 0.1*r*reach/yr2sec(1); % m^2/s
    
    charEQ  = 6; % characteristic slip (m)
    
    recurranceInterval  = charEQ/v;
    offsetArray         = v*(recurranceInterval);
    
    avulseYN = 'break';
    showLive = false;
    
    tc = eq_cycle_channel_model(x,z,dx,w0, ...                         % Geometry
                                k, N,      ...                         % Diffusion
                                recurranceInterval, offsetArray, ...   % Earthquakes
                                hc, avulseYN, ...                      % Channel
                                showLive);    
    
end

end

function ut = diffuse(x,fx,k,t,N)
    x = x(:); fx= fx(:);
    nx= length(x);

    L   = range(x);
    T1  = fx(1); 
    T2  = fx(end);
   
    u_E = T1 + (T2-T1)/L * x;
        
    %% determine coefficients:
    B = zeros(N,1);
    Narr = (1:N)';
    for n = Narr'
        B(n) = 2/L * trapz(x,(fx-u_E).* sin(n*pi*x/L));
    end
        
    %% solution at time t
    ut = u_E' +  B' * (sin(Narr*pi*x'/L) .* (exp(-k*(Narr*pi/L).^2*t)*ones(1,nx)));
    
end

function [x,z] = simple_profile(S0,dx,L)
        x = 0:dx:L;
        z = S0*L-S0*x;
end

function [nx,ny] = stretch(x,y,dx,w,xo,stretch_factor)

    % edge case to deal with: 
    % w<dx
    % fault right on/near the edge
    
    I  = x>xo-w/2 & x<xo+w/2;
    I0 = find(I,1,'first');
    If = find(I,1,'last'); 
    
    % initial stretched segment:
    yS_o = y(I);
    xS_o = x(I);
    
    % stretch
    % x -> x*stretch_factor
    xS = xS_o(1)+(xS_o-xS_o(1))*stretch_factor; 
    % densify
    xD = linspace(min(xS),max(xS),ceil(range(xS)/dx)); % NOTE! this results in uneven point spacing!
    
    % densify the y coord
    yD = interp1(xS,yS_o,xD);
    
    nx = [x(1:I0-1),xD,x(If+1:end)+(range(xD)-range(xS_o))];
    ny = [y(1:I0-1),yD,y(If+1:end)];

end
                                  
function tc = eq_cycle_channel_model(x,z,dx,w0, ...                         % Geometry
                                     k, N,      ...                         % Diffusion
                                     recurranceInterval, offsetArray, ...   % Earthquakes
                                     hc, avulseYN, ...                      % Channel 
                                     showLive)
        % turn hc into a function handle  if not already
        if ~isa(hc,'function_handle')
            hcfh = @(t) hc*ones(size(t));
        else 
            hcfh = hc;
        end
                                 
        % initial geometry
        [z0,zz] = deal(z);
        [x0,xx] = deal(x);
        nEq= length(recurranceInterval);
        
        figure; set(gcf,'color','w');
        
        subplot(2,2,1); 
        ph0 = plot(x0,z0,'Color',[1 0 0]); hold on
        
        wInProfile = w0;
        
        T = sum(recurranceInterval);
        
        fx0 = max(x)/2;
        xc = fx0 - w0/2;
        Ic = find(x0<xc,1,'last');
        plot(xc*ones(1,2),minmax(z0),'--');
        zxc = zeros(1,nEq);
        zxc0 = interp1(x0,z0,xc);
        tc = 0;
        NA = 1;
      
               
        for n = 1:nEq
            
            stretch_factor = sqrt(w0^2+(sum(offsetArray(NA:n))^2))/wInProfile;
            [x,z] = stretch(x,z,dx,wInProfile,fx0,stretch_factor);
            [xx,zz] = stretch(xx,zz,dx,wInProfile,fx0,stretch_factor);
            
            z = diffuse(x,z,k,recurranceInterval(n),N);
            
            % plot the profile after it has diffused (right before the next
            % earthquake)
            
            subplot(2,2,1); 
            lp = plot(x,z, ...
                'Color',[1-sum(recurranceInterval(1:n))/T,0,sum(recurranceInterval(1:n))/T]);
            I = x>xc & x<(xc + wInProfile);
            
            hold on
            
            lf = plot(x(I),z(I), ...
                'Color',[1-sum(recurranceInterval(1:n))/T,0,sum(recurranceInterval(1:n))/T], ...
                'LineWidth',1.5);
            
            lp.Color(4) = 0.5;
            lf.Color(4) = 0.5;
            
            if nargin == 11
                if showLive
                    plot(x0,z0)
                    plot(xx,zz,'k')
                    hold off
                    set(gca,'xlim',[min(x0),max(x0)*1.5])
                    set(gca,'ylim',minmax(z0))
                    drawnow
                end
            end
            
            subplot(2,2,3)
            
            curvature  = gradient(gradient(z, x),x);
            curvNorm   = curvature/max(abs(curvature));
            plot(x,curvNorm, ...
                'Color',[1-sum(recurranceInterval(1:n))/T,0,sum(recurranceInterval(1:n))/T]);
            I = x>xc & x<(xc + wInProfile);
            
            hold on
            
            lf = plot(x(I),curvNorm(I), ...
                'Color',[1-sum(recurranceInterval(1:n))/T,0,sum(recurranceInterval(1:n))/T], ...
                'LineWidth',1.5);
            
            lp.Color(4) = 0.5;
            lf.Color(4) = 0.5;
            
            if nargin == 11
                if showLive
                    plot(xlim,[0 0],'-','Color',[0.8 0.8 0.8])
                    hold off
                    set(gca,'Ylim',[-1,1],'Ytick',[-0.5,0.5],'yticklabel',{'Eroding', 'Aggrading'})
                    set(gca,'xlim',[min(x0),max(x0)*1.5])
                    drawnow                 
                end
            end
            
            ti = sum(recurranceInterval(1:n));
            zxc(n) = interp1(x,z,xc);
            
            if zxc(n)-zxc0 > hcfh(ti)
                tc = [tc,ti-sum(tc)];
                if avulseYN
                    % Resest the channel to the original form beyond xc
                    % this is not necessarily as straight forward as it may
                    % seem...
                    z = [z(1:Ic),z0((Ic+1):end)];                                   % reset to original channel after the critical point
                    % z = [z(1:Ic),z0((Ic+1):end) + (z(Ic)-z0(Ic+1))];              % reset to the original channel after the critical point, but shifted up by the aggradation at xc
                    % z = [z(1:Ic),z(Ic)-z(Ic)*(x0((Ic+1):end)-xc)/(max(x0)-xc)];   % reset to return to base level over the original domain size from xc
                    % z = [z(1:Ic), (z0(Ic+1)+hcfh(ti))-((z0(Ic+1)+hcfh(ti))-z0(end))/(max(x0)-xc)*(x0((Ic+1):end)-xc)];
                    x = [x(1:Ic),x0((Ic+1):end)];
                    
                    zz= z0;
                    xx= x0;
                    
                    wInProfile = w0;
                    NA = n;
                    
                elseif strcmp(avulseYN,'break')  
                    break        
                end
            end
            
            wInProfile = sqrt(w0^2+ sum(offsetArray(NA:n))^2);
            fx0        = xc + wInProfile/2;
            
        end
        
       
        
        tc = tc(2:end);
        yr2sec  = @(yr) yr*(60*60*24*365);
        t = cumsum(recurranceInterval)/yr2sec(1)/1000;
        ax1 = subplot(2,2,[2,4]); hold on
        plot(t,zxc-zxc0,'.-')
        ylabel('Aggradation [m]')
        xlabel('Time (kyr)')
        phhc = plot(t,hcfh(t),'--');
        legend(phhc,'h_c')
        set(gca,'Ylim',[0,max(zxc-zxc0)*1.5])

        t = title('c)'); 
        set(t,'Units','normalized','Position',[-0.22,0.95])
        
        subplot(2,2,1);
        try legend([ph0,lp],{'Initial profile','Final profile'}); catch; legend(lp,'Final profile'); end
        ylabel('Elevation [m]')
        t = title('a)'); 
        set(t,'Units','normalized','Position',[-0.22,0.95])
        
        subplot(2,2,3);
        xlabel('Along profile distance [m]')
        ylabel('Normalized Curvature')
        
        t = title('b)'); 
        set(t,'Units','normalized','Position',[-0.22,0.95])
        
        ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
        setsize = @(fh,dim1,dim2) set(fh,...
            'Units',        'Inches', ...
            'Position',     [0,0,dim1,dim2],...
            'PaperUnits',   'Inches',...
            'PaperSize',    [dim1,dim2]);
        ftsz(gcf,12)
        setsize(gcf,7,3.5)

end

% substitude for minmax
function OUT = minmax(A)

OUT = [min(A),max(A)];

end
