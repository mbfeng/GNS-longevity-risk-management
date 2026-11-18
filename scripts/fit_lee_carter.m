function lc = fit_lee_carter(export)
% FIT_LEE_CARTER  Fig. 1: Lee-Carter fit to England & Wales female data (Section 4).
%   Fits the Lee-Carter model (ages 60-89, years 1962-2011, Human Mortality Database)
%   by Poisson maximum likelihood, plots a_x, b_x and k_t, and prints the parameters
%   of the random walk for kappa that the paper reports.
%
%   fit_lee_carter()       fits, plots and saves figures/Fig_abk.png
%   fit_lee_carter(false)  does not save the figure

if nargin < 1, export = true; end
lc = kappa_model();

fig = figure;
panels = {lc.a, '\alpha_x', 'Age (x)'; lc.b, '\beta_x', 'Age (x)'; lc.k, '\kappa_t', 'Year (t)'};
for j = 1:3
    subplot(1, 3, j)
    hold on
    box on
    plot(panels{j, 1}, 'LineWidth', 2)
    title(panels{j, 2})
    ylabel('Value')
    xlabel(panels{j, 3})
    if j < 3
        xticks(1:5:30)
        xticklabels(lc.Age(1):5:lc.Age(end))
    else
        xticks(1:10:50)
        xticklabels(lc.Year(1):10:lc.Year(end))
        ylim([-15 12])
    end
end
if export, export_figure(fig, 'Fig_abk'); end

fprintf('Random walk for kappa (paper: drift -0.4971, variance 0.7467, kappa_0 -14.2328):\n');
fprintf('  drift %.4f, variance %.4f, kappa_0 %.4f\n', lc.kmu, lc.kVar, lc.k0);
end
