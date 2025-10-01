%% Plotting Model fits using data from Scientific Reports

clear all; close all; clc


%% Loading the soft resin data

soft_resin_data = load("Liao et al. 2025 - Scientific Data Datasets\SoftMaterialModels.mat").SOFT_MATERIALS_MODELS;
rigid_resin_data = load("Liao et al. 2025 - Scientific Data Datasets\RigidMaterialModels.mat").RIGID_MATERIALS_MODELS;
rigid_resin_stress_strain = load("Liao et al. 2025 - Scientific Data Datasets\RigidResinStressStrainData.mat").RIGID_MATERIALS_DATA;



%{
                        Elastomer Data
    The data from each resin is stored in a cell array of structs. Asiga
    DentaGUM is in position 3, FormLabs Silicone (IPA) is position 4, and
    FormLabs Silicone (Mix) is position 5. Within each of those cell arrays
    (in position 2) is a cell array of structs containing the data.
    Position 1 is the control data, Position 2 is the Autoclave data, and
    Position 3 is the Ethanol data. 

                        Rigid Resin Data
    This data is also stored cell arrays. Position 1 is Asiga DentaGUIDE,
    Position 2 is Liqcreate BioMed Clear, and Position 3 is Phrozen Aqua
    Gray. In these cell arrays, the first position is a tag with the resin
    name, and the follow three positions are cell arrays for the control,
    autoclave, and ethanol groups, respectively. In each
    sterilization-specific array, the data is organized as {YoungsMod,
    100*strain_at_break,UTS,CompMod}.
%}

sterilization = {'NS','Autoclave','EtOH'};

% Breaking out the data by resin
asiga_dentaGUM = soft_resin_data{3};
formlabs_IPA = soft_resin_data{4};
formlabs_mix = soft_resin_data{5};

asiga_dentaGUIDE = rigid_resin_data{1};
liqcreate_biomed = rigid_resin_data{2};
phrozen_aqua_gray = rigid_resin_data{3};

asiga_dentaGUIDE_data = rigid_resin_stress_strain{1};
liqcreate_biomed_data = rigid_resin_stress_strain{2};
phrozen_aqua_gray_data = rigid_resin_stress_strain{3};

%% Model code for Yeoh model


% uniaxial deformation ( lam1 = lam, lam2 = lam3 = 1/sqrt(lam) )
I1 = @(lam) lam.^2 + 2./lam;

% Yeoh model Engineering stress
sig_Yeoh = @(lam,C1,C2,C3) 2*(lam - lam.^-2).*( C1 + 2*C2*(I1(lam) - 3) + 3*C3*(I1(lam)-3).^2);

% colors = {1.25*[0.2,0.8,0.6], [0.95,0.33,0.64], [0.09,0.45,0.6]};
colors = {[80,215,50]/255,[253,46,122]/255,[78,95,125]/255};
alpha = 0.3;
lines = {'-','--','-.'};



%% Plotting the mean model fits for the soft resin under each condition

% asiga dentaGUM
PlotElastomerModelFits(asiga_dentaGUM,sterilization,lines,alpha,colors)
saveas(gcf,"Output Figures\AsigaDentaGum_ModelFits.svg")

% formlabs silicone IPA
PlotElastomerModelFits(formlabs_IPA,sterilization,lines,alpha,colors)
saveas(gcf,"Output Figures\FormlabsSiliconeIPA_ModelFits.svg")

% formlabs silicone mix
PlotElastomerModelFits(formlabs_mix,sterilization,lines,alpha,colors)
saveas(gcf,"Output Figures\FormlabsSiliconeMIX_ModelFits.svg")


%% PLotting the mean model fits for the rigid resins under each condition
PlotRigidTensileModelFits(asiga_dentaGUIDE,asiga_dentaGUIDE_data,sterilization,lines,alpha,colors)
saveas(gcf,"Output Figures\AsigaDentaGUIDE_ModelFits.svg")




function PlotElastomerModelFits(dataset,sterilization,lines,alpha,colors)

    figure("Position",[100,100,600,400],"Color","w"); hold all
    stats = dataset{3};

    for k=1:length(sterilization)
        fill([stats{k}.lam_all,fliplr(stats{k}.lam_all)],[stats{k}.sig_mean + stats{k}.sig_std, fliplr(stats{k}.sig_mean - stats{k}.sig_std)],colors{k},'FaceAlpha',alpha,'EdgeColor','none',"HandleVisibility","off");
        plot(stats{k}.lam_all,stats{k}.sig_mean,lines{k},'color',colors{k},"DisplayName",sterilization{k} + "  (R^2="+num2str(stats{k}.R2,3)+")","LineWidth",1.5)
    end
    xlabel("Stretch [mm/mm]")
    ylabel("Engineering Stress [MPa]")
    legend("Location","northwest")
    %     xlim([0.5,3])
    %     ylim([-1,5])
    set(gca,"FontSize",15)
end

function PlotRigidTensileModelFits(dataset,raw_data,sterilization,lines,alpha,colors)

    figure("Position",[100,100,600,400],"Color","w"); hold all
    
    strain = linspace(0,0.2,50);
    
    % Tensile data
    for k=1:length(sterilization)
        E = dataset{1+k}{1};
        stress = E'*strain;
    
        stress_mean = mean(stress,1);
        stress_std = std(stress,[],1);

        subplot(3,4,[1,2,5,6,9,10]); hold all
        fill(100*[strain,fliplr(strain)],[stress_mean + stress_std, fliplr(stress_mean - stress_std)],colors{k},'FaceAlpha',alpha,'EdgeColor','none',"HandleVisibility","off");
        plot(100*strain,stress_mean,lines{k},'color',colors{k},"DisplayName",sterilization{k},"LineWidth",1.5)
        
        ylim([0,60])
        xlim([0,15])
        xlabel("Strain [%]")
        ylabel("Engineering Stress [MPa]")
        legend("Location","southeast")
        set(gca,"FontSize",15)

        E_mean = mean(E);
        [~,closest_to_mean] = min(abs(E - E_mean)); % find the slope that is closest to the mean to plot as representative sample
        strain_data = raw_data{1+k}{3}(closest_to_mean,:);
        stress_data = raw_data{1+k}{4}(closest_to_mean,:);

        subplot(3,4,[4*k-1,4*k]); hold all
        fill(100*[strain,fliplr(strain)],[stress_mean + stress_std, fliplr(stress_mean - stress_std)],colors{k},'FaceAlpha',alpha,'EdgeColor','none',"HandleVisibility","off");
        plot(100*strain,stress_mean,lines{k},'color',colors{k},"DisplayName",sterilization{k},"LineWidth",1.5)
        
        % plot(100*strain_data,stress_data,'o','color',colors{k},"LineWidth",1.5,"HandleVisibility","off")
        for m=1:length(E)
            closest_to_mean = m;
            strain_data = raw_data{1+k}{3}(closest_to_mean,:);
            stress_data = raw_data{1+k}{4}(closest_to_mean,:);

            plot(100*strain_data,stress_data,'.','color',colors{k},"HandleVisibility","off","LineWidth",1.5)
        end
        ylim([0,60])
        xlim([0,15])
        if k==3
            xlabel("Strain [%]")
        end
        ylabel({"Engineering","Stress [MPa]"})
        % legend("Location","northwest")
        set(gca,"FontSize",11)
        
    end

    





end