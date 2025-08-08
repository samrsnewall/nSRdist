%% Create table to show chi2Stats for each pairing
numPDFs = 4;
numDatas = 4;
chi2Table = table('Size', [numPDFs, numDatas], 'VariableTypes', ["double", "double", "double", "double"],'VariableNames',["BMhist", "LinCoresLinMethBMode", "NewCoresLinMethBMode", "NewCoresNewMethBMode"], 'RowNames',["BMpdf", "LCLMBModeMLN", "NCLMBModeMLN", "NCNMBModeMLN"]);
chi2Table.BMhist(1) = chiStatBMvBM.chi2stat;
chi2Table.LinCoresLinMethBMode(1) = LinCoreLinMeth.chiStatBM.chi2stat;
chi2Table.LinCoresLinMethBMode(2) = LinCoreLinMeth.chiStat.chi2stat;
chi2Table.NewCoresLinMethBMode(1) = NewCoreLinMeth.chiStatBM.chi2stat;
chi2Table.NewCoresLinMethBMode(3) = NewCoreLinMeth.chiStat.chi2stat;
chi2Table.NewCoresNewMethBMode(1) = NewCoreNewMeth.chiStatBM.chi2stat;
chi2Table.NewCoresNewMethBMode(4) = NewCoreNewMeth.chiStat.chi2stat;

figure;
hold on
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(x, LinCoreLinMeth.mixLogBMode(:,2), '-b', 'DisplayName', "LinCores,Marine09", 'LineWidth', 1)
plot(x, NewCoreLinMeth.mixLogBMode(:,2), '-g', 'DisplayName', "NewCores,Marine09", 'LineWidth', 1)
plot(x, NewCoreNewMeth.mixLogBMode(:,2), '-r', 'DisplayName', "NewCores,Marine20", 'LineWidth', 1)
legend()
xlim([0 6])
xlabel("nSR")
ylabel("pdf")
title("Bchronology Mode Results")

%% Plot mix log normal from individual random sampling runs (note, random sampling of Bchronology at given depths - i.e. Lin Meth

commonYLim = [0 max([LinCoreLinMeth.MLN1R.c95up, NewCoreLinMeth.MLN1R.c95up, NewCoreNewMeth.MLN1R.c95up])];

figure;
hold on
plot(x, LinCoreLinMeth.MLN1R.pdfs, 'Color', [0,0,0,0.1], 'HandleVisibility', 'off')
plot(x, LinCoreLinMeth.mixLogBMode(:,2), '-r', 'DisplayName', "LinCores,Marine09", 'LineWidth', 1)
plot(NaN, NaN, 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(LinCoreLinMeth.MLN1R.chiStat1RunT.h1R).*100, 3) + "%")
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(NaN, NaN, 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(LinCoreLinMeth.MLN1R.chiStat1RunT_BM.h1R).*100, 3) + "%")
xlim([0 6])
ylim(commonYLim)
legend()
title("LinCores; Marine09; R = 0±0")

figure;
hold on
plot(x, NewCoreLinMeth.MLN1R.pdfs, 'Color', [0,0,0,0.1], 'HandleVisibility', 'off')
plot(x, NewCoreLinMeth.mixLogBMode(:,2), '-r', 'DisplayName', "NewCores,Marine09", 'LineWidth', 1)
plot(NaN, NaN, 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NewCoreLinMeth.MLN1R.chiStat1RunT.h1R).*100, 3) + "%")
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BIGMACS", 'LineWidth', 1)
plot(NaN, NaN, 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NewCoreLinMeth.MLN1R.chiStat1RunT_BM.h1R).*100, 3) + "%")
xlim([0 6])
xlim([0 6])
ylim(commonYLim)
legend()
title("NewCores; Marine09; R = 0±0")

figure;
hold on
plot(x, NewCoreNewMeth.MLN1R.pdfs, 'Color', [0,0,0,0.1], 'HandleVisibility', 'off')
plot(x, NewCoreNewMeth.mixLogBMode(:,2), '-r', 'DisplayName', "BchronMode: NewCores,Marine20", 'LineWidth', 1)
plot(NaN, NaN,'LineStyle', "none", 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NewCoreNewMeth.MLN1R.chiStat1RunT.h1R).*100, 3) + "%")
plot(x, MLN_BIGMACS(:,2), '-k', 'DisplayName', "BchronMode: BIGMACS", 'LineWidth', 1)
plot(NaN, NaN, 'LineStyle', "none", 'DisplayName', "Random Runs that reject H0 = " + num2str(mean(NewCoreNewMeth.MLN1R.chiStat1RunT_BM.h1R).*100, 3) + "%")
xlim([0 6])
%ylim(commonYLim)
legend()
xlabel("nSR")
ylabel("PDF")
title("NewCores; Marine20; R = 0±200")