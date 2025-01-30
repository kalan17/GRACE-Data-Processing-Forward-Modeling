%% P02. Forward Modeling
% This code is designed to improve the spatial resolution of GRACE data
% by using WATERGAP Total Water Storage (TWS) data.
% The forward modeling approach iteratively refines GRACE-derived mass changes.
%
% <Data Sources>
% Watergap_TWS_sh_2003.mat: WATERGAP TWS data from Jan 2003 to Dec 2003.
% Download and more info: https://doi.pangaea.de/10.1594/PANGAEA.918447

%% 1. Forward modeling
    clear,clc
    % 1.1 Load WATERGAP TWSmdata
    load ./Data/Watergap_TWS_sh_2003.mat Ini_tws Ini_tws_sh Ini_tws_sh_surfden lat lon
    % Ini_tws: WATERGAP TWS (grid format)
    % Ini_tws_sh: WATERGAP TWS (spherical harmonic coefficients)
    % Ini_tws_sh_surfden: WATERGAP TWS SH converted into surface density (kg/m^2 or mmH2O)

    % 1.2 Load processed GRACE GSM data (output from previous step P01)
    load ./output/p01_Slepian_GRACE_GSM_surfacemass_P4M8_DL60_cg_G400_mmH20.mat
    GSM_sl_D60_SL60=GSM_sl_D60_SL60(:,:,5:16); % from Jan 2003~ Dec 2003
    
    % 1.3 Variable settings
    N=60; % Spherical harmonic truncation degree
    load ./Data/mk_cg_05.mat; % Load land mask
    idx_cg=map_cg_05==1; % Identify land pixels

    diff_SH=SH_ord_eom(N); % Initialize SH difference array
    diff_RMS=nan(12,200); % Initialize RMS difference tracking array
    FM_iter=zeros(12,1);  % Store iteration counts
    FM_diff_grid=zeros(size(Ini_tws));
    FM_tws_sh_surfden=Ini_tws_sh_surfden; % Initialize forward modeled SH data (mmH2O,WATERGAP)
    FM_routing=Ini_tws; % Initialize forward modeled grid data (mmH2O,WATERGAP)
    
    %% 1.4 Start Forward Modeling
    % This process can take a long time.
    time_start=fix(clock)  % Record start time

    for t = 1:length(GSM_sl_D60_SL60(1,1,:))
        ii = 1;
        diff_grid_RMS = 200;
        while diff_grid_RMS > 0.1 && ii <= 200
            % Apply Slepian localization and Gaussian smoothing to SH coefficients
            FM_tws_sh_surfden(:,:,t)=Slepian_local(FM_tws_sh_surfden(:,:,t), G_60, sN_60);
            FM_tws_sh_surfden(:,:,t)=gaussian(FM_tws_sh_surfden(:,:,t), N, 400);

            % Compute difference between GRACE and WATERGAP SH data
            diff_SH(:,3:4)=GSM_sl_D60_SL60(:,3:4,t)-FM_tws_sh_surfden(:,3:4,t);

            % Convert SH difference to spatial grid
            diff_grid=plm2xyz(diff_SH, 0.5, [0.25 89.75 359.75 -89.75]);

            % Update forward-modeled grid with difference correction
            FM_routing(:,:,t)=FM_routing(:,:,t)+diff_grid;

            % Convert corrected grid back to SH domain
            FM_tws_sh_surfden(:,:,t)=xyz2plm(FM_routing(:,:,t), N, 'im');

            % Compute RMS of the difference
            diff_grid_RMS=rms(diff_grid(idx_cg));
            diff_RMS(t,ii)=diff_grid_RMS;
            FM_iter(t,1)=ii;

            % Convergence check: Stop if improvement is minimal
            if ii >= 2 && (diff_RMS(t,ii-1)-diff_RMS(t,ii)) < 0.01 && (diff_RMS(t,ii-1)-diff_RMS(t,ii)) > 0
                break;
            end
            ii=ii + 1;
        end
        FM_diff_grid(:,:,t)=diff_grid; % Store difference for visualization
        disp(t);
    end

    time_end=fix(clock); % Record end time
    %% save file
    time=time(5:16);
    save ./output/P02_FM_CSR_GSM_RL06_sp60_cg_P4M8_G400_deg1_c20_c30_slr_PE_ICE6GD_watergap.mat GSM_sl_D60_SL60 Ini_tws FM_tws_sh_surfden diff_SH diff_grid FM_routing FM_tws_sh_surfden  diff_grid_RMS diff_RMS FM_diff_grid  
%% 2. Visualization of Results
    Rm = 6.371e6; % Earth's mean radius (m)
    rho_w = 1000;  % Water density (kg/m^3)
    rho_e = 5500;  % Earth's mean density (kg/m^3)
    n = 1; % Select time step for visualization

    figure(1), set(gcf, 'color', 'w');

    subplot(2,2,1)
    imagesccont(plm2xyz(GSM_sl_D60_SL60(:,:,n), 0.5, [0.25 89.75 359.75 -89.75]))
    c = colorbar;
    c.Label.String = 'mmH2O';
    caxis([-100 100]);
    axis([5 40 70 120]);
    title('GRACE SH Water Mass Change');

    subplot(2,2,2)
    imagesccont(Ini_tws(:,:,n) .* map_cg_05);
    c.Label.String = 'mmH2O';
    caxis([-500 500]);
    axis([5 40 70 120]);
    colorbar;
    title('WATERGAP TWS');

    subplot(2,2,3)
    imagesccont((FM_routing(:,:,n) .* map_cg_05));
    caxis([-500 500]);
    c = colorbar;
    c.Label.String = 'mmH2O';
    axis([5 40 70 120]);
    title('Forward Modeling Result');

    subplot(2,2,4)
    imagesccont(FM_diff_grid(:,:,n));
    colorbar;
    axis([5 40 70 120]);
    title('Difference: GRACE SH vs. WATERGAP');