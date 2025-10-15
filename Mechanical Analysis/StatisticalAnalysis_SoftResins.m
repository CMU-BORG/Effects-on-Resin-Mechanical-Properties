%% Statistical comparison of the bootstrapped Elastomer model parameters

clear all; close all; clc


%% Loading the soft resin data

soft_resin_data = load("Liao et al. 2025 - Scientific Data Datasets\SoftMaterialModels.mat").SOFT_MATERIALS_MODELS;

sterilization = {'NS','Autoclave','EtOH'};

% Breaking out the data by resin
asiga_dentaGUM = soft_resin_data{3};
formlabs_IPA = soft_resin_data{4};
formlabs_mix = soft_resin_data{5};

figure("color","w","Position",[100,100,680,300]); hold all
xline(0,'--k','LineWidth',2);   % null difference
comparison_ID = 0;

y_ticks = [];
y_tick_labels = {};


%% Formlabs Silicone 40A IPA
disp("Formlabs Silicone 40A IPA Comparisons")
formlabs_IPA_comparison = struct();
% non-sterile to autoclave
[interval,samples] = BootstrapDifference(formlabs_IPA{2}{1}.E,formlabs_IPA{2}{2}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Formlabs IPA NS-Autoclave"};
comparison_ID = comparison_ID - 1;
fprintf("\tNS-Autoclave: [%.3f,%.3f] MPa\n",interval(1),interval(2))

% non-sterile to EtOH
[interval,samples] = BootstrapDifference(formlabs_IPA{2}{1}.E,formlabs_IPA{2}{3}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Formlabs IPA NS-EtOH"};
comparison_ID = comparison_ID - 1;
fprintf("\tNS-EtOH: [%.3f,%.3f] MPa\n",interval(1),interval(2))

% autoclave to EtOH
[interval,samples] = BootstrapDifference(formlabs_IPA{2}{2}.E,formlabs_IPA{2}{3}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Formlabs IPA Autoclave-EtOH"};
comparison_ID = comparison_ID - 2;
fprintf("\tAutoclave-EtOH: [%.3f,%.3f] MPa\n",interval(1),interval(2))


%% Formlabs Silicone 40A Mix
disp("Formlabs Silicone 40A Mix Comparisons")

formlabs_mix_comparison = struct();
% non-sterile to autoclave
[interval,samples] = BootstrapDifference(formlabs_mix{2}{1}.E,formlabs_mix{2}{2}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Formlabs IPA NS-Autoclave"};
comparison_ID = comparison_ID - 1;
fprintf("\tNS-Autoclave: [%.3f,%.3f] MPa\n",interval(1),interval(2))

% non-sterile to EtOH
[interval,samples] = BootstrapDifference(formlabs_mix{2}{1}.E,formlabs_mix{2}{3}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Formlabs IPA NS-EtOH"};
comparison_ID = comparison_ID - 1;
fprintf("\tNS-EtOH: [%.3f,%.3f] MPa\n",interval(1),interval(2))

% autoclave to EtOH
[interval,samples] = BootstrapDifference(formlabs_mix{2}{2}.E,formlabs_mix{2}{3}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Formlabs IPA Autoclave-EtOH"};
comparison_ID = comparison_ID - 2;
fprintf("\tAutoclave-EtOH: [%.3f,%.3f] MPa\n",interval(1),interval(2))



%% Asiga DentaGUM
disp("Asiga DentaGUM Comparisons")

asiga_dentaGUM_comparison = struct();
% non-sterile to autoclave
[interval,samples] = BootstrapDifference(asiga_dentaGUM{2}{1}.E,asiga_dentaGUM{2}{2}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Asiga DentaGUM NS-Autoclave"};
comparison_ID = comparison_ID - 1;
fprintf("\tNS-Autoclave: [%.3f,%.3f] MPa\n",interval(1),interval(2))

% non-sterile to EtOH
[interval,samples] = BootstrapDifference(asiga_dentaGUM{2}{1}.E,asiga_dentaGUM{2}{3}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Asiga DentaGUM NS-EtOH"};
comparison_ID = comparison_ID - 1;
fprintf("\tNS-EtOH: [%.3f,%.3f] MPa\n",interval(1),interval(2))

% autoclave to EtOH
[interval,samples] = BootstrapDifference(asiga_dentaGUM{2}{2}.E,asiga_dentaGUM{2}{3}.E,[0.025,0.975]);
PlotTestResult(interval,samples,comparison_ID)
y_ticks = [y_ticks, comparison_ID];
y_tick_labels = {y_tick_labels{:}, "Asiga DentaGUM Autoclave-EtOH"};
comparison_ID = comparison_ID - 2;
fprintf("\tAutoclave-EtOH: [%.3f,%.3f] MPa\n",interval(1),interval(2))


yticks(flip(y_ticks));
yticklabels(flip(y_tick_labels))

xlabel("Young's Modulus Difference [MPa]")



function [interval, samples] = BootstrapDifference(data1,data2,interval_range)
%{
    Calculate the bootstrapped difference between sampling distributions
    for model parameters. NOTE: data1 and data2 should be estimates of
    model parameters and not raw data. For raw data comparisons, an
    additional step is required where a model parameter (eg, a mean) value
    is calculated
%}


samples = data1 - data2'; 
samples = samples(:);     % full bootstrapp distribution of difference (ie. every pairwise difference in parameter values). Not efficient for large datasets

interval = quantile(samples,interval_range);

end

function PlotTestResult(interval,samples,ID)
    
    mu = mean(samples);
    sig = std(samples);

    plot(interval,ID*[1,1],'-k','LineWidth',2)
    plot(mu,ID,'sk',"MarkerFaceColor",'k','MarkerSize',10)
end