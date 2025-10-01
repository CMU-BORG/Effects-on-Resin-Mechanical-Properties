%% Calculating and saving corrected stress data from Liao et al. 2025 rigid resin mechanical tests

clear all; close all; clc


% colors = {1.25*[0.2,0.8,0.6], [0.95,0.33,0.64], [0.09,0.45,0.6]};
colors = {[80,215,50]/255,[253,46,122]/255,[78,95,125]/255};
alpha = 0.3;
lines = {'-','--','-.'};
markers = {'o','v','diamond'};

resin_names = {{"3Dresyn Bioflex","A10 MB - IPA"},
               {"3Dresyn Bioflex","A10 MB - UNW2"},
               {"Asiga DentaGUIDE"},
               {"Asiga DentaGUM"},
               {"Formlabs Silicone","40A - IPA"},
               {"Formlabs Silicone","40A - IPA/BuOAc"},
               {"Liqcreate","Bio-Med Clear"},
               {"Phrozen AquaGray 8K"}};

sterilization = {'NS','Autoclave','EtOH'};


masterDataFile ="Liao et al. 2025 - Scientific Data Datasets\masterRawData_clean.xlsx";
geomOpts = detectImportOptions(masterDataFile,"Sheet","SpecificGeometries");
geomOpts.VariableTypes(5:6) = {'string','string'};
compOpts = detectImportOptions(masterDataFile,"Sheet","Compression_Rem");
compOpts.VariableTypes(4:5) = {'string','string'};
tensOpts = detectImportOptions(masterDataFile,"Sheet","Tensile_Rem");
tensOpts.VariableTypes(4:5) = {'string','string'};

geom_data = readtable(masterDataFile,geomOpts);
compress_data = readtable(masterDataFile,compOpts);
tensile_data = readtable(masterDataFile,tensOpts);

geom_data.Resin = lower(geom_data.Resin);
compress_data.Resin = lower(compress_data.Resin);
tensile_data.Resin = lower(tensile_data.Resin);

resins = unique(geom_data.Resin);
tests = unique(geom_data.TestType);



% stretch when accounting for tare load displacement
LAM = @(dL, a, L0) 1 + (dL + a)/L0;


%% Calculating stress and corrected Strain

RIGID_MATERIALS_DATA = cell(3,1);
resin_number = 0;
N_interp_points = 50;

peak_stress = zeros(8,1);
peak_stress(3) = 60;
peak_stress(7) = 60;
peak_stress(8) = 30;

peak_c_stress = zeros(8,1);
peak_c_stress(3) = 11;
peak_c_stress(7) = 11;
peak_c_stress(8) = 14;

for i=[3,7,8] % loop through only the rigid resins
    resin_number = resin_number + 1; % index into the SOFT_MATERIALS_MODEL cell array

    fprintf("\nResin: %s\n",resins{i})

    % storing the fits and statistics for the different sterilization
    % techniques
    stress_strain = cell(3,1);

    for j=1:length(sterilization)
        col = colors{j}; %rand(3,1);

        % Compression Data
        inds = ismember(compress_data.Resin,resins(i)) & ismember(compress_data.Sterilization,sterilization(j));
        subset = compress_data(inds,:);
        
        % loop through all three plates and both samples per plate
        YoungsMod = [];         % Young's modulus [MPa]
        CompMod = [];           % Compressive Modulus [MPa]
        strain_at_break = [];   % strain at failure [mm/mm]
        UTS = [];               % utilimate tensile strength [MPa]

        tensile_stress_ = [];
        tensile_strain_ = [];
        comp_stress_ = [];
        comp_strain_ = [];


        samp_count = 1;
        
        plates = unique(subset.PlateID);
        for k=1:length(plates)
            specimen = unique(subset{ismember(subset.PlateID,plates{k}),"SpecimenID"});
            for l=1:length(specimen)

                samp_ind = ismember(subset.PlateID,plates{k,1}) & ismember(subset.SpecimenID,specimen{l,1});
                samp_data = subset(samp_ind,:);

                % getting the rest geometry of the current sample
                geom_ind = ismember(geom_data.TestType,"Compression") & ...
                           ismember(geom_data.Resin,resins(i)) & ...
                           ismember(geom_data.Sterilization,sterilization(j)) & ...
                           ismember(geom_data.PlateID,plates{k}) & ...
                           ismember(geom_data.SpecimenID,specimen{l});

                use_sample = geom_data{geom_ind,"UseSample_NoForIfExperimentalError_"};
                sample_broke = geom_data{geom_ind,"SampleBroke_0Or1_"};

                if use_sample
                    L0 = geom_data{geom_ind,"ThicknessOrGageLength_mm_"};
                    A0 = geom_data{geom_ind,"Area_mm2_"};
                    
                    % pulling the data for the current sample
                    t = samp_data.Time_sec_;
                    dL = samp_data.Crosshead_mm_; 
                    F = samp_data.Load_N_;

                    % calculating the engineering stress
                    sig = F/A0;
                    dL = dL - dL(1);

                    sig = 1*sig; % negative because of compression
                    
                    % calculating the nominal stretch
                    lam = 1 + dL/L0;
                    
                    % interpolating the data to ensure equal weight between
                    % tensile and compression tests
                    lam_ = LAM(dL,0,L0);
                    lam_interp = linspace(min(lam_),max(lam_),N_interp_points);
                    sig_interp = interp1(lam_,sig,lam_interp);
    
                    e1 = diff(sig_interp)./diff(lam_interp);
                    e1 = smooth(e1,10,"moving");
                    
                    [~,ind_] = max(abs(e1));
    
                    if ind_+1 > length(lam_interp)
                        linear_region = length(lam_interp)-2:length(lam_interp);
                    elseif ind_<2
                        linear_region = 1:3;
                    else
                        linear_region = ind_-1:ind_+1;
                    end
    
                    eps_linear = lam_interp(linear_region)' - 1;
                    sig_linear = sig_interp(linear_region)';
    
                    p1 = polyfit(eps_linear,sig_linear,1);
    
                    E1 = p1(1);
                    delL = L0*(p1(2)/p1(1));
                    lam_ = LAM(dL,delL,L0);
                    lam_interp = linspace(min(lam_),max(lam_),N_interp_points);
                    
                    comp_strain_ = [comp_strain_; lam_interp-1];
                    comp_stress_ = [comp_stress_; sig_interp];
                    
    
                    samp_count = samp_count + 1;
                end

            end
        end

        % Tensile Data
        inds = ismember(tensile_data.Resin,resins(i)) & ismember(tensile_data.Sterilization,sterilization(j));
        subset = tensile_data(inds,:);
        plates = unique(subset.PlateID);
        for k=1:length(plates)
            
            specimen = unique(subset{ismember(subset.PlateID,plates{k}),"SpecimenID"});

            for l=1:length(specimen)
                
                % getting the data for the current sample
                samp_ind = ismember(subset.PlateID,plates{k,1}) & ismember(subset.SpecimenID,specimen{l,1});                
                samp_data = subset(samp_ind,:);
                
                % getting the rest geometry of the current sample
                geom_ind = ismember(geom_data.TestType,"Tensile") & ...
                           ismember(geom_data.Resin,resins(i)) & ...
                           ismember(geom_data.Sterilization,sterilization(j)) & ...
                           ismember(geom_data.PlateID,plates{k}) & ...
                           ismember(geom_data.SpecimenID,specimen{l});

                use_sample = geom_data{geom_ind,"UseSample_NoForIfExperimentalError_"};
                sample_broke = geom_data{geom_ind,"SampleBroke_0Or1_"};

                if use_sample
                    L0 = geom_data{geom_ind,"ThicknessOrGageLength_mm_"};
                    A0 = geom_data{geom_ind,"Area_mm2_"};
    
                    % pulling the data for the current sample
                    t = samp_data.Time_sec_;
                    dL = samp_data.Crosshead_mm_; 
                    F = samp_data.Load_N_;

                    % calculating the engineering stress
                    sig = F/A0;
                    dL = dL - dL(1);
                    
                    % calculating the nominal stretch
                    lam = 1 + dL/L0;
                    
                    % interpolating the data to ensure equal weight between
                    % tensile and compression tests
                    lam_ = LAM(dL,0,L0);
                    lam_interp = linspace(min(lam_),max(lam_),50);
                    sig_interp = interp1(lam_,sig,lam_interp);
                    
                    e1 = diff(sig_interp)./diff(lam_interp);
                    e1 = smooth(e1,5,"moving");
                    
                    [~,ind_] = max(e1);
                    if ind_ < 2
                        % highest slope is right at the beginning
                        linear_region = 1:3;
                    else
                        linear_region = ind_-1:ind_+1;
                    end
                    eps_linear = lam_interp(linear_region)' - 1;
                    sig_linear = sig_interp(linear_region)';
    
                    p1 = polyfit(eps_linear,sig_linear,1);
    
                    E1 = p1(1);
                    delL = L0*(p1(2)/p1(1));
                    lam_ = LAM(dL,delL,L0);
                    lam_interp = linspace(min(lam_),max(lam_),50);
    
                    tensile_strain_ = [tensile_strain_; lam_interp-1];
                    tensile_stress_ = [tensile_stress_; sig_interp];
                    


                end
            end
        end


               
        stress_strain{j} = {comp_strain_,comp_stress_,tensile_strain_,tensile_stress_};
               
    end
    
    RIGID_MATERIALS_DATA{resin_number} = [resins(i);stress_strain];

end

save("Liao et al. 2025 - Scientific Data Datasets\RigidResinStressStrainData.mat","RIGID_MATERIALS_DATA")