%% P01. GRACE Data Processing
% This code is for processing GRACE Level 2 (Spherical Harmonic Coefficients) data.
% The goal is to convert geoid/potential data into surface mass changes, apply necessary corrections,
% and filter the data for noise reduction and improved spatial localization.
% 
% <Data Sources>
% GRACE_200204_201706.mat: CSR Release-06 GRACE Level-2 Data from April 2002 to June 2017. 
% Download and more info:https://www2.csr.utexas.edu/grace/RL06.html

%% 1. Convert geoid/potential to surface mass (water)
% This step involves converting spherical harmonic coefficients to meaningful surface mass changes. 
% Necessary corrections such as degree 1 correction, C20/C30 replacement, and GIA corrections are applied.
    %% 1.1 Degree 1 correction
    clear,clc

    load ./Data/GRACE_200204_201706;
    FN1='./Data/deg1_coef.txt';
    % The degree-1 coefficients (geocenter corrections) are applied using TN-13a estimates.
    % For more details, refer to:
    % https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_CSR_RL0603.txt

    deg1=textread(FN1,'','headerlines',69);
    
    C10=deg1(1:2:end-1,4);
    C11=deg1(2:2:end,4); 
    S11=deg1(2:2:end,5);
    
    % Filling data gap for 2004.01
    % The solution for 2004.01 is based on less than 15 days of data and is not considered useful.
    % Therefore, the gap is filled by linearly interpolating between the December 2003 and February 2004 data.
    C10=[C10(1:16);(C10(16)+C10(17))/2;C10(17:end)];
    C11=[C11(1:16);(C11(16)+C11(17))/2;C11(17:end)];
    S11=[S11(1:16);(S11(16)+S11(17))/2;S11(17:end)];
    
    %% 1.2 Covert C20,C30
    load ./Data/c20_c30_slr.mat;
    % C20 & C30 Replacement:
    % C20 coefficients are replaced with SLR solutions from [Loomis et al., 2019, GRL, doi:10.1029/2019GL082929], 
    % ensuring consistency with GRACE/GRACE-FO data used in the Science community. 
    % C30 coefficients are replaced with SLR solutions from the same source, but only for GRACE-FO. 
    % Reference: TN-14 (https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt)

    % Extract valid time range and interpolate missing data
    GRACE_data_1=GRACE_data(:,:,3:end); % Start from August 2002
    C20=C20(3:end);
    C30_NaN=C30_NaN(3:end);
    
    for ii=1:length(C10)
        ii2=113:length(C10);
        % Apply Degree 1 corrections
        GRACE_data_1(2,3,ii)=C10(ii,1);
        GRACE_data_1(3,3,ii)=C11(ii,1);
        GRACE_data_1(3,4,ii)=S11(ii,1);

        % Replace C20 and C30
        GRACE_data_1(4,3,ii)=C20(ii,1);
        GRACE_data_1(7,3,ii2)=C30_NaN(ii2,1);
    end
    %% 1.3 Interpolate missing data (2002.09 to 2017.05)
    t0=datenum(2002,8,15,0,0,0);
    ddays=daysinmonth_eom(2002,8,2017,4);
    time=cumsum(ddays)+t0;

    I_GRACE_data=GRACE_data_1;
    thedates=time_GRACE(3:163,1);

    for d=1:size(I_GRACE_data,1)
        GSM_D60_C(d,:)=interp1(thedates,squeeze(I_GRACE_data(d,3,:)),time);
        GSM_D60_S(d,:)=interp1(thedates,squeeze(I_GRACE_data(d,4,:)),time);
    end

    GSM_data=repmat(GRACE_data(:,1:2),[1,1,size(time,1)]);
    GSM_data(:,3,:)=GSM_D60_C(:,:);
    GSM_data(:,4,:)=GSM_D60_S(:,:);

    %% 1.4  Apply GIA Correction
    load ('./Data/PE_ICE6GD_stokes_GRACE_monthly.mat')
    % GIA correction based on the ICE6G-D model (Peltier et al., 2018).
    PGR=PGR(1:1891,:);

    for ii=1:length(ddays)
        GSM_data_pgr(:,3,ii)=GSM_data(:,3,ii)-ii*PGR(:,3);
        GSM_data_pgr(:,4,ii)=GSM_data(:,4,ii)-ii*PGR(:,4);
    end
    
    GSM_data_pgr(:,1,:)=GSM_data(:,1,:);
    GSM_data_pgr(:,2,:)=GSM_data(:,2,:);
    %% 1.5 Convert geoid/potential to surface mass (water equivalent)
    Rm=6.371e6; % Earth's radius in meters
    rho_w=1000; % Density of water in kg/m^3
    rho_e=5500; % Average density of Earth in kg/m^3

    c=squeeze(GSM_data_pgr(:, 3, :));
    s=squeeze(GSM_data_pgr(:, 4, :));

    c(1, :)=0; % Remove degree 0 term

    N=60;
    [deg, ord]=deg_ord(N);

    c_mean=mean(c, 2);
    s_mean=mean(s, 2);

    c=c - c_mean; 
    s=s - s_mean;

    for i=1:size(c, 2)
        deg=squeeze(GSM_data_pgr(:, 1, i));
        k=lovenumbers(deg, 'k');

        c_sm_mmH2O(:, i)=c(:, i) * (Rm * rho_e) / 3 .* ((2 .* deg) + 1) ./ (1 + k);
        s_sm_mmH2O(:, i)=s(:, i) * (Rm * rho_e) / 3 .* ((2 .* deg) + 1) ./ (1 + k);
    end

    GRACE_sf_mass_mmH2O(:, 1, :)=GSM_data_pgr(:, 1, :);
    GRACE_sf_mass_mmH2O(:, 2, :)=GSM_data_pgr(:, 2, :);
    GRACE_sf_mass_mmH2O(:, 3, :)=c_sm_mmH2O;
    GRACE_sf_mass_mmH2O(:, 4, :)=s_sm_mmH2O;
%% 2. Filtering
% Filtering improves data quality by reducing noise and enhancing spatial localization. 
% This process includes:
% 1. De-correlation (PM Filtering)
% 2. Spatial localization (Slepian)
% 3. Gaussian smoothing

deg=GRACE_sf_mass_mmH2O(:,1,:);
ord=GRACE_sf_mass_mmH2O(:,2,:);
clm_e=GRACE_sf_mass_mmH2O(:,3,:);
slm_e=GRACE_sf_mass_mmH2O(:,4,:);
    %% 2.1 PM Filtering (De-correlation)
    % Reduces correlated errors due to aliasing and imperfect data processing.
    for ii=1:size(time,1)
        [cm(:,:,ii),sm(:,:,ii)]=sh_vec_mat(clm_e(:,ii),slm_e(:,ii),N);
        [clm_p(:,:,ii),slm_p(:,:,ii)]=decorr_cs(cm(:,:,ii),sm(:,:,ii),4,8,N,N);

        ccp(:,:,ii)=convt2low_plm(clm_p(:,:,ii));
        ssp(:,:,ii)=convt2low_plm(slm_p(:,:,ii));
    end
    % Convert filtered coefficients back to spherical harmonic format
      lmcosi_p(:,1,:)=deg;
      lmcosi_p(:,2,:)=ord;
      lmcosi_p(:,3,:)=ccp;
      lmcosi_p(:,4,:)=ssp;
    %% 2.2 Slepian Localization
    % Improves spatial resolution by focusing on the target region (e.g., Congo River Basin).
     load ./Data/mk_cg_05.mat;

     N=60;
    intv=0.5;
    maskgrid=map_cg_05;
    [G, sN, V, dllmm, N, intv_lat, intv_lon]=Slepian_basis(maskgrid,N,intv);

    for t=1:size(time,1)
        lmcosi_p_sp(:,:,t)=Slepian_local(lmcosi_p(:,:,t),G,sN);
    end
    %% 2.3 Gaussian Smoothing
    % Applies a low-pass filter to suppress high-frequency noise.
    for ii=1:size(time,1)
        [clm_g(:,ii),slm_g(:,ii)]=gauss(lmcosi_p_sp(:,3,ii),lmcosi_p_sp(:,4,ii),N,400e3);
    end

     lmcosi_g(:,1,:)=deg;
     lmcosi_g(:,2,:)=ord;
     lmcosi_g(:,3,:)=clm_g;
     lmcosi_g(:,4,:)=slm_g;
%% Save the file  
% Check if the 'output' folder exists, and create it if it does not.
if ~exist('./output', 'dir')
    mkdir('./output');
end

GSM_sl_D60_SL60=lmcosi_g;
G_60=G; 
sN_60=sN;

% Save the data
save('./output/P01_Slepian_GRACE_GSM_surfacemass_P4M8_DL60_cg_G400_mmH20.mat', 'GSM_sl_D60_SL60', 'time', 'G_60', 'sN_60');