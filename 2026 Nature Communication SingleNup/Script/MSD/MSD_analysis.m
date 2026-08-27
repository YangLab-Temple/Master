SPACE_UNITS = 'nm';
TIME_UNITS = 's';

N_PARTICLES = 10;
N_TIME_STEPS = 100;
N_DIM = 2; % 2D

% Typical values taken from studies of proteins diffusing in membranes:
% Diffusion coefficient
D  = 1e-3; % µm^2/s
% Time step between acquisition; fast acquisition!
dT = 0.05; % s,

% Area size, just used to disperse particles in 2D. Has no impact on
% analysis.
SIZE = 2; % µm

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

ma = ma.addAll(tracks);

ma = ma.computeMSD;
imsd = ma.msd
%% 
figure
hmsd = ma.plotMeanMSD(gca, true);
mmsd = ma.getMeanMSD;

%% 
ma = ma.fitLogLogMSD(0.5);
ma.loglogfit
alphas = ma.loglogfit.alpha;
gammas = ma.loglogfit.gamma;
%% 
newarr = [imsd{:}];

writematrix(newarr,'TPRC_inditrack.csv');
writematrix(mmsd,'TPRC_sumtrack.csv');
writematrix(alphas,'TPRC_alpha.csv');
writematrix(gammas,'TPRC_gamma.csv');