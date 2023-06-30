%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2023-- James O Thompson (Baylor University)
% Contributors: Kenneth Befus, Michael Manga , Chelsea Allison, Ben Black, 
% Ben Andrews, Anna Ruefer
% 1D diffusion of H2O secondary hydration through an Embayment
% Non-isothermal
% Michael Manga provided the conductive cooling model

% Requires: MATLAB 2019a or newer
%           Toolbox: curve_fitting_toolbox

% Last updated: 06/02/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    saveout = 1; %change to 0 for no, 1 for yes

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Only chnage the variables between the 0s and 1s
%000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
%%% DIFFUSION INPUTS

%%% FILENAME OF EMBAYMENT FILE
    filename = 'MesaFallsAllLines.xlsx'; %this is the data table
    outputname = 'MF35'; %name the sample
    
    position_pos = 271; %this is the column number of the position data in the table
    H2O_pos = 276; %this is the column number of the H2O data in the table
    
%%% H2O

    H2O_Concentration_in_Conduit = 2.02; %Measured far field concentration in wt%
    H2O_Concentration_in_Embayment = 0.83; %Concentration in embayment prior to eruption
    
%%% HYDRATION DEPTH
    Hydration_depth = 10; %meters
    
%%% ERUPTION TIMESCALE
    Hydration_Duration_yr = 10;      %years to run the model, choosing a close value to cooling time is critical for efficient program run 
    Hydration_Duration = Hydration_Duration_yr *(365*24); %Hours to run the model

%%% TEMPERATURE CHANGE (NON-ISOTHERMAL)
    % Temperature model speed
    RFactor = 0.5; %**** KENNY CHANGE THIS TO CHANGE COOLING SPEED *** LOWER NUMBER == FASTER

%%% DIFFUSION CHANGE (NON-STOICHI)
    % Diffusion change with distance
    DFactor = 1.5; %**** KENNY CHANGE THIS TO CHANGE CURVATURE *** HIGHER NUMBER == MORE CURVE

%000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

%111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
%%% CONDUCTION THERMAL MODEL INPUTS
    
    % Vertical spacing between grid points
    % in meters(n times dz is the vertical thickness of all units together)
    % WARNING: I would recommend leaving this alone
    dz = 0.05; 
    % Years to run model
    years_cooling = 30.0;
    % Porosity of substrate
    phi_sub = 0.2;
    % Porosity of fall
    phi_fall = 0.6; 
    % Porosity of ignimbrite 
    phi_ign = 0.6; 
    % k0 of substrate
    k0_sub = 1.6; 
    % k0 of fall
    k0_fall = 1.6; 
    % k0 of ignimbrite 
    k0_ign = 1.6; 
    % Specific heat of substrate in J/(kg K)
    c_sub = 1000.0;
    % Specific heat of fall in J/(kg K)
    c_fall = 1000.0; 
    % Specific heat of ignimbrite  in J/(kg K)
    c_ign = 1000.0; 
    % Temperature of substrate in degrees C
    t_sub = 5.0; 
    % Temperature of fall in degrees C
    t_fall = 5.0; 
    % Temperature of ignimbrite  in degrees C
    t_ign = 500.0; 
    % Thickness of substrate in m
    z_sub = 20.0; 
    % Thickness of fall in m
    z_fall = 5.0; 
    % Thickness of ignimbrite in m
    z_ign = 10.0; 
    % Depth of sample from top of fall in m
    z_sample_fall = 0.5; 

%111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
%%% END OF USER INPUTS
tic
%%% Rehydration Diffusion Model for Embayments %%%

%%% LOAD embayment file and extract axes, position on mount as columns of data into Matlab
    Validation_Data = readmatrix(filename);
    Validation_Data = Validation_Data(7:end,:);
    file_x = Validation_Data(:,position_pos); %second column x position
    file_x =  file_x(~isnan(file_x));
    vqH2O = (Validation_Data(:,H2O_pos));
    vqH2O =  vqH2O(~isnan(vqH2O));
    x_length_actual = max(file_x)-min(file_x); 
    x_length = length(file_x)/length(find(file_x==max(file_x)));
    dx = x_length_actual / (x_length-1);
    xq = file_x;
    minvqH2O = min(vqH2O);

%%% BUILD MASKS
    mask = zeros(x_length,1);
    conduit_mask = mask;
    conduit_mask(1,1) = 1;
    emb_mask = mask;
    emb_mask(2:end,1) = 1;
    emb_diff_mask = conduit_mask+emb_mask;
    conduit_mask(1,1) = 1;
   
%%% CONSTANTS
    MW_H2O = 18.01528;
    F_MDensity = 2380.0; % density in kg/m3
    
%%% FILL THE MESH AND DOMAIN (2D)
    Endy = x_length;
   %%H2O
    H2O_C_Con = conduit_mask .* H2O_Concentration_in_Conduit ; %Concentration_in_Conduit;
    H2O_C_Emb = emb_mask .* H2O_Concentration_in_Embayment ;
    H2O_C = H2O_C_Con + H2O_C_Emb ;
 
%%% TIME VARIABLES AND MANIPULATIONS FOR MODEL
    % calculate inital depth based on pressure
    GPa = (Hydration_depth * ( F_MDensity * 9.81 )) / 1e9;

    % Temperature change (non-isothermal)
    %old stuff
    %expmulti = 1000; %exponent multiplier
    %coolt = 30 * (365*24*60*60); %seconds
    %iTempK = (TempK_start-TempK_end)/expmulti *(10.^(log10(expmulti+1)/coolt*(coolt-0)) - 1)+ TempK_end; %in seconds
    %iTempK = (MFCoolingModel(0)) + 273.15;
    %%%
 
%%% TEMPERATURE MODEL (NON-ISOTHERMAL)
    % Run the Conductive_Cooling_Function to produce a cooling model
    [Sample_Conductive_Model,Temperature_Vertical_Profile] = ...
    Conductive_Cooling_Function(dz,years_cooling, ...
    phi_sub,phi_fall,phi_ign,k0_sub,k0_fall,k0_ign,c_sub,c_fall, ...
    c_ign,t_sub,t_fall,t_ign,z_sub,z_fall,z_ign,z_sample_fall);
    
    % Extract time (days) and fall sample cooling from the model results
    Time = Sample_Conductive_Model(:,1);
    FallSample = Sample_Conductive_Model(:,3);
    % Calculate cooling model fit from the Conductive_Cooling_Function
    MFCoolingModel = fit(Time,FallSample,'gauss3');
    %load MFFallCoolinModelGauss3.mat %model is in days
    
    % Temperature model fit constants
    a1 = MFCoolingModel.a1;
    b1 = MFCoolingModel.b1*RFactor; %move first peak left-right
    c1 = MFCoolingModel.c1*RFactor; %compress first peak 0.3
    a2 = MFCoolingModel.a2;
    b2 = MFCoolingModel.b2*RFactor; %move second peak left-right and up-down
    c2 = MFCoolingModel.c2*RFactor; %compress second peak 0.4
    a3 = MFCoolingModel.a3;
    b3 = MFCoolingModel.b3*RFactor; %move third peak left-right and up-down
    c3 = MFCoolingModel.c3*RFactor; %compress third peak
    
    % initial model results
    time0 = 0; %time in days
    iTempK = (a1.*exp(-((time0-b1)./c1).^2) + a2.*exp(-((time0-b2)./c2).^2) + ...
        a3.*exp(-((time0-b3)./c3).^2)) + 273.15;
    
%%% INTITIAL DIFFUSIVITY AND CONDITIONS
    H2O_Co = 1;  %C naught set to 1 wt % for the Ni and Zhang diffusion model
    time_step_diffusion = (exp((H2O_Concentration_in_Conduit/H2O_Co)^DFactor))* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)/iTempK)); 
    dt = 0.4*(dx^2/time_step_diffusion); %in seconds
    %dt=10000;
    Max_Time = 3600*Hydration_Duration;  % seconds %this is used as a max time, can be raised if needed
    dt_total=0; %time indicator

%%% MODEL HOW THE DISTURBANCE MIGRATES THRU TIME
    H2O_Cnew = H2O_C;     
    dr = dx;
    r_x = xq+(dx/2);      
    iter=1 ;
    z=1;
    gamma = 2; %Gamma - will be used to change between linear = 0; cylindrical = 1; spherical=2 %treatments; usually gamma will be set to 2 

%%% MAIN LOOP - Diffusion 1D %%% 
    figure()
    while dt_total < Max_Time
        
%          %H2Ot (<2wt%)
%          %D = C.*exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)./TempK)); %D in um2/s
%          %infinite D sink in conduit
%          H2O_D = ( exp(H2O_C.^DFactor).* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)./iTempK)) ).*emb_diff_mask; %D in um2/s in conduit - 10 fold increase in diffusion
         
         %%H2Ot
         %H2Ot concentration matrix
         H2Ot_le2 = double(H2O_C <= 2);
         H2Ot_gt2 = double(H2O_C > 2);
         %H2Ot (<2wt%)
         H2O_D2_le2 = ( exp(H2O_C.^DFactor).* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)./iTempK)) ).*emb_diff_mask; %D in um2/s in conduit - 10 fold increase in diffusion
         %Diffusion equation from Ni and Zhang (2008), eq.13 or in abstract, applies to H2Ototal less than 2wt.%
         %H2Ot (>2wt%)
         H2OSX = (H2O_C./MW_H2O)./ ( (H2O_C./MW_H2O) + ((100.0-H2O_C)./32.49) );
         H2O_D2_gt2 = ( H2OSX.* exp( 13.470 + (-49.996.*H2OSX) + (7.0827.*(H2OSX.^0.5)) + (1.8875.*GPa) - ...
               ( ( 9532.3 + (-91933.0.*H2OSX) + (13403.0.*(H2OSX.^0.5)) + (3625.6.*GPa) ) ./ iTempK ) ) ) .*emb_diff_mask;%
         %Diffusion equation from Ni and Zhang (2008), eq.12 or in abstract
         H2O_D = (H2O_D2_le2.*H2Ot_le2) + (H2O_D2_gt2.*H2Ot_gt2);

     % Simplifying variables in discretization
        
         up1x = x_length-1; low1x = 2; %starting and ending points in diffusion model x
         mult1x = dt./(xq(low1x:up1x).^gamma.*dx.^2);
         rgammax = r_x(low1x:up1x).^gamma; % same as "(r(i))^gamma"
        %%H2O
         H2O_Dmain = H2O_D(low1x:up1x);  %D(i)
         H2O_Cmain = H2O_C(low1x:up1x); %C(i)
          
     % Discretized Diffusion equation 
        %%H2O
         H2O_diffopperator = mult1x .* (rgammax.*((H2O_Dmain.*H2O_D(low1x+1:up1x+1))./(H2O_Dmain+H2O_D(low1x+1:up1x+1))).*(H2O_C(low1x+1:up1x+1)-H2O_Cmain)...
             -rgammax.*((H2O_Dmain.*H2O_D(low1x-1:up1x-1))./(H2O_Dmain+H2O_D(low1x-1:up1x-1)).*(H2O_Cmain-H2O_C(low1x-1:up1x-1))));
        
     % Update Concentration 
        %%H2O
         H2O_Cnew(low1x:up1x) = H2O_Cmain + H2O_diffopperator;
         H2O_Cnew(1) = H2O_Cnew(2);
         H2O_Cnew(end) = H2O_Cnew(end-1);
         H2O_C=H2O_Cnew;

      % Plotting Part of Model  
              
                H2O_textpositions = [ 0.05  0.1 0.15    0.2 0.25    0.30    0.3500    0.4    0.45    0.5    0.55    0.6   0.8 0.85    .9 ];
                H2O_ymax = H2O_Concentration_in_Embayment + 1; % ymin = 0;
                H2O_textloc = H2O_ymax.*H2O_textpositions;
         
                if iter==z  
                
                    
                plot(xq,H2O_C,'linewidth',3)
                xlabel('Position (\mum)')
                ylabel('dH_2O (wt%)')
                xlim([min(min(xq)),max(max(xq))])
                ylim([min(min(vqH2O))*0.7, max(max(vqH2O))*1.2])
                %title('Modeled')
                % Enlarge figure to full screen.
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
                hold on
                plot(xq,vqH2O,'linewidth',3)
                hold off
                title({append('Time Elapsed (yrs) = ',num2str(dt_total/(3600*24*365))),append('Diffusivity (um2/s) = ',num2str(H2O_D(2,1)))})

                drawnow
                             
                z=z+100; %This controls how often the plot updates, higher numbers is less frames

               end 

             iter=iter+1;


    
         % Update diss H2O in conduit based on above
         u_H2O_C_Con = conduit_mask .* H2O_Concentration_in_Conduit;
         u_H2O_C_Emb = emb_mask .* H2O_C;
         H2O_C = u_H2O_C_Con + u_H2O_C_Emb ;    
         %H2O_C(end) = H2O_Concentration_in_Embayment;  
         
         % Update temperature based on cooling 
         %iTempK = (MFCoolingModel((dt_total/(24*60*60)))) + 273.15;
         iTempK = (a1.*exp(-(((dt_total/(24*60*60))-b1)./c1).^2) + ...
             a2.*exp(-(((dt_total/(24*60*60))-b2)./c2).^2) + ...
             a3.*exp(-(((dt_total/(24*60*60))-b3)./c3).^2)) + 273.15;
         if iTempK < 273; iTempK = 273; end

         % Update time step for next round based on temperature
         dt = 0.4*(dx^2/H2O_D(1));
         if dt > 31536000 & Hydration_Duration_yr < 1e4; dt = 31536000; end
         dt_total = dt_total+dt;
         
         % Temporal variables to save
         H2O_D_Mouth_iter(iter) = H2O_D(end-10,1);
         dt_total_iter(iter) = dt_total;
         iTempK_iter(iter) = iTempK;
 
    end
    close()
    H2O_D_Mouth_iter(iter) = H2O_D(end-10,1);
    dt_total_iter(iter) = dt_total;
    iTempK_iter(iter) = iTempK;
    toc
%%
%%% PLOT FIGURES and SAVE
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
    plot(xq, H2O_C,'linewidth',2)
    xlabel('Position (\mum)')
    ylabel('dH_2O (wt%)')
    axis([0 max(xq) min(min(vqH2O))*0.7 max(max(vqH2O))*1.05])
    hold on
    plot(xq, vqH2O,'linewidth',2,'color','[1,0.563,0]')
    hold off
    leg1 = legend({'Model','FTIR'},'Location','southwest');
    leg1.FontSize = 16;
    title({append('Time Elapsed (yrs) = ',num2str(dt_total/(3600*24*365)), ...
        '    Cooling Rate: RFactor = ', num2str(RFactor)), ...
        append('H2O_{D}^{x} component: DFactor = ', num2str(DFactor), ...
        '    Diffusivity (\mum^{2}/s) = ',num2str(H2O_D(2,1))),[]})
    set(gca,'FontSize',16)
    saveas(gcf,append(outputname,'_Rehydration_Diffusion_1D_Profile.jpg'))
    
    drawnow

    %diffusivity over time
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
    plot(dt_total_iter(2:end)./31536000, H2O_D_Mouth_iter(2:end),'linewidth',2)
    xlabel('Time (Years)')
    ylabel('Diffusivity (um^{2}/s)')
    axis([0 max(dt_total_iter./31536000) min(min(H2O_D_Mouth_iter))*0.7 max(max(H2O_D_Mouth_iter))*1.5])
    set(gca, 'YScale', 'log')
    title({append('Time Elapsed (yrs) = ',num2str(dt_total/(3600*24*365)), ...
        '    Cooling Rate: RFactor = ', num2str(RFactor)), ...
        append('H2O_{D}^{x} component: DFactor = ', num2str(DFactor), ...
        '    Diffusivity (\mum^{2}/s) = ',num2str(H2O_D(2,1))),[]})
    set(gca,'FontSize',16)
    saveas(gcf,append(outputname,'_Rehydration_Diffusivity_Time.jpg'))

    %temperature over time
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
    plot(dt_total_iter(2:end)./31536000, iTempK_iter(2:end)-273.15,'linewidth',2)
    xlabel('Time (Years)')
    ylabel('Temperature (\circC)')
    axis([0 max(dt_total_iter./31536000) min(min(iTempK_iter(2:end)-273.15))*0.7 max(max(iTempK_iter(2:end)-273.15))*1.1])
    title({append('Time Elapsed (yrs) = ',num2str(dt_total/(3600*24*365)), ...
        '    Cooling Rate: RFactor = ', num2str(RFactor)), ...
        append('H2O_{D}^{x} component: DFactor = ', num2str(DFactor), ...
        '    Diffusivity (\mum^{2}/s) = ',num2str(H2O_D(2,1))),[]})
    set(gca,'FontSize',16)
    saveas(gcf,append(outputname,'_Rehydration_Temperature_Time.jpg'))

    %temperature vs diffusivity
    xtemp = 5:5:240;
    lowTyhat = 4E-6.*exp(-11.2.*(1000./(xtemp+273.15)))*1E12;
    highTyhat = 1E-6.*exp(-12.4.*(1000./(xtemp+273.15)))*1E12;
    
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
    plot(iTempK_iter(2:end)-273.15, H2O_D_Mouth_iter(2:end),'linewidth',2)
    hold on
    plot(xtemp, lowTyhat,'linewidth',2,'Color','k')
    plot(xtemp, highTyhat,'linewidth',2,'Color','r')
    xlabel('Temperature (\circC)')
    ylabel('Diffusivity (\mum^{2}/s)')
    axis([0 xtemp(end)+10 1E-14 1E-2])
    set(gca, 'YScale', 'log')
    set ( gca, 'xdir', 'reverse' )
    legend('This Study','Low-T Experiments','High-T Experiments')
    title({append('Time Elapsed (yrs) = ',num2str(dt_total/(3600*24*365)), ...
        '    Cooling Rate: RFactor = ', num2str(RFactor)), ...
        append('H2O_{D}^{x} component: DFactor = ', num2str(DFactor), ...
        '    Diffusivity (\mum^{2}/s) = ',num2str(H2O_D(2,1))),[]})
    set(gca,'FontSize',16)
    saveas(gcf,append(outputname,'_Rehydration_Diffusivity_Temperature.jpg'))

    if saveout == 1
%%% OUTPUT DATA
    
    header = {'Position (um)','H2O FTIR (wt.%)','H2O Model (wt.%)'};
    oned_data = [xq vqH2O H2O_C ];
    output_data_header = [header; num2cell(oned_data)];
    writecell(output_data_header,append(outputname,'_Rehydration_1D_Profile.xls'))

    header1 = {'Time Elapsed (years)','Temperature (\circC)','Diffusivity (um^{2}/s)'};
    oned_data1 = [(dt_total_iter(2:end)./31536000)' (iTempK_iter(2:end)-273.15)'  H2O_D_Mouth_iter(2:end)'];
    output_data_header1 = [header1; num2cell(oned_data1)];
    writecell(output_data_header1,append(outputname,'_Rehydration_Temporal_Parameters.xls'))
    end
    
    
%%%FUNCTIONS

% calculate the time-temperature history of particles in cooling deposits
function [Sample_Conductive_Model,Temperature_Vertical_Profile] = ...
    Conductive_Cooling_Function(dz,years_cooling, ...
    phi_sub,phi_fall,phi_ign,k0_sub,k0_fall,k0_ign,c_sub,c_fall, ...
    c_ign,t_sub,t_fall,t_ign,z_sub,z_fall,z_ign,z_sample_fall)

% outputs:
%           Sample_Conductive_Model (4xn array)
%               Column: 1 - time elapsed (days)
%                       2 - Temperature at fall base (C)
%                       3 - Temperature at fall sample (C)
%                       4 - Temperature at ignimbrite  center (C)
%           Temperature_Vertical_Profile (2xn array)
%               Row: 1 - Profile distance (meters)
%                    2 - Temperature through section (C)

    unit_thinkness = z_sub + z_fall + z_ign;
    imx = round(unit_thinkness/dz) ; %grid points x
    n = imx; % grid points n
    t = zeros(1,imx); %temperature in degrees C
    tn = zeros(1,imx); %temperature at the next time step
    rk = zeros(1,imx); %thermal conductivity in W/(m K)
    sat = zeros(1,imx); %saturation (I do nothing with this since I assume there was no liquid water)
    c = zeros(1,imx); %Specific heat in J/(kg K)
    z = zeros(1,imx); %Vertical position in m
    rho = zeros(1,imx); %density in kg/m^3
    seconds = years_cooling*365.0*24.0*3600.0; % total seconds to run model
    imx_sub = round(z_sub/dz); %grid points from base to top of sub
    imx_fall = round((z_sub/dz)+(z_fall/dz)); %grid points from base to top of fall
    imx_ign = round((z_sub/dz)+(z_fall/dz)+(z_ign/dz)); %grid points from base to top of ign
    imx_mid_ign = round((z_sub/dz)+(z_fall/dz)+(0.5*z_ign/dz)); %grid points from base to mid of ign
    imx_sample_fall = round(imx_fall-(z_sample_fall/dz)); %grid points from base to top of sample

    % define grid and reference values
    for i = 1:n
        c(i) = 1000.0;
        rk(i) = 2.0;
        rho(i) = 2000.0;
        sat(i) = 0.0;
        z(i) = dz*i;
    end

    % initial temperature
    % this is the substrate
    for i = 1:imx_sub
        t(i) = t_sub;
        c(i) = c_sub;
        rk(i) = k0_sub*(1.0-phi_sub)/(1.0+phi_sub);
        rho(i) = 2600.0*(1.0-phi_sub);
    end

    % fall deposit
    for i = imx_sub+1:imx_fall 
        t(i) = t_fall;
        c(i) = c_fall;
        rk(i) = k0_fall*(1.0-phi_fall)/(1.0+phi_fall);
        rho(i) = 2600.0*(1.0-phi_fall);
    end

    % ignimbrite
    for i = imx_fall+1:imx_ign 
        t(i) = t_ign;
        c(i) = c_ign;
        rk(i) = k0_ign*(1.0-phi_ign)/(1.0+phi_ign);
        rho(i) = 2600.0*(1.0-phi_ign);
    end
    t(n) = t_sub;

    % loop over time
    dt = 0.2*dz*dz*rho(1)*c(1)/rk(1);
    time = 0.0;
    it_time5000 = 0;
    n_time5000 = round((seconds/dt)/5000);
    Sample_Conductive_Model = zeros(n_time5000,4);

    for it = 1:round(seconds/dt)

        time = time + dt;

        for i = 2:n-1
            f = 1.0/(rho(i)*c(i)*dz*dz);
            a = 0.5*(rk(i+1)+rk(i));
            aa = 0.5*(rk(i-1)+rk(i));
            b = 0.5*(a+aa);
            tn(i) = t(i) + dt*f*(a*t(i+1)-2.0*b*t(i)+aa*t(i-1));
        end

        %update temperature profile
        t(2:end-1) = tn(2:end-1);

        %print*,time/(24.0*3600.0),t(1200),t(1450), t(1750)
        if  rem(it,5000) == 0  
            it_time5000 = it_time5000 + 1;
            Sample_Conductive_Model(it_time5000,:) = [time/(24.0*3600.0),t(imx_sub),t(imx_sample_fall),t(imx_mid_ign)];
        end

    end

    % compile temperature profile
    Temperature_Vertical_Profile = [z - z_sub ; t];

end

%%%References used in this code
    %Zhang, Y., & Ni, H. (2010). Diffusion of H, C, and O components in silicate melts. Reviews in Mineralogy and Geochemistry, 72(1), 171-225.
    %Ni, H., & Zhang, Y. (2008). H2O diffusion models in rhyolitic melt with new high pressure data. Chemical Geology, 250(1-4), 68-78.
    %Zhang, Y., Xu, Z., Zhu, M., & Wang, H. (2007). Silicate melt properties and volcanic eruptions. Reviews of Geophysics, 45(4).
    %Thermal conductivity of unsaturated tuffs from Sass et al. (1987) open file report 87-649
    %k0=1.6 and c=1000 from Lavallee et al. (Fronteirs 2015) and thermal conductivity model is Bagdassavov and Dingwell (1994)
       