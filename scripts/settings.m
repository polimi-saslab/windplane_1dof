%% Figures settings
set(groot,'defaultAxesFontSize',20)
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
plt.LineStyle = '-';
plt.LineStyle2 = '--';
plt.dim = [0 3 20 11.5];
plt.b = [0, 0.4470, 0.7410];
plt.r = [0.8500, 0.3250, 0.0980];
plt.y = [0.9290, 0.6940, 0.1250];
plt.v = [0.4940 0.1840 0.5560];

%% Import casadi
addpath('C:\Users\trevi\OneDrive - Politecnico di Milano\Desktop\PhD\Working\casadi-3.7.1-windows64-matlab2018b')
import casadi.*
