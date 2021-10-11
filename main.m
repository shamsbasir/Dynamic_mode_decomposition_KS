% Shamsulhaq Basi 
clear all; close all; clc;
FS = 'FontSize';LW = 'LineWidth';

%% Data generated using KS.m
t = 0.0:.001:4.0;
dt = t(2)-t(1);
data = importdata('Data.mat');

%% Create DMD data matrices
X1 = data(:,1:1999);
X2 = data(:,2:2000);

%% SVD and rank -r truncation
r = 21;
[U S V] = svd(X1,'econ');
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);

%% Build Atilde and DMD Modes
Atilde = Ur'*X2*Vr/Sr; 
[W, eigs] = eig(Atilde);   
Phi = X2*Vr/Sr*W; % DMD Modes 

%% DMD Spectra
lambda = diag(eigs);                       % discrete-time DMD eigenvalues
omega = log(lambda)/dt;                    % the continuous-time DMD eigenvalues
[omega,I]=sort(omega,'descend','ComparisonMethod','real'); % sort the omegas
Phi= Phi(:,I(1:r));                        % sort the DMD Modes accordingly

%% 5 DMD Modes
figure
plot(real(Phi(:,1:5)),LW,1.6) % can set this to r
title('First 5 DMD Modes')
legend('mode 1','mode 2','mode 3','mode 4','mode 5')
set(gca,FS,10)

%% EigenValues 
figure
scatter(real(lambda),imag(lambda),'ob');
title('21 Eigenvalue, Real{\Lambda} Vs. Imag{\Lambda}')
xlabel("Real{\Lambda}");
ylabel("Imag{\Lambda}");
set(gca,FS,10)

%% Compute DMD Solution
x1 = X1(:,1);
b = Phi\x1;
time_dynamics = zeros(r,length(t));
for iter = 1:length(time_dynamics)
    time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Phi * time_dynamics;

%% POD vs EXACT 
figure
subplot(1,2,1)
surfl(real(X_dmd));
shading interp; colormap(gray);
title("DMD solution");

subplot(1,2,2)
surfl(real(data));
shading interp;colormap(gray);
title("Exact solution");

%% DMD forcast 
t = [2.5 3.0 3.5 4.0];
figure
for i = 1:4
    subplot(2,2,i)
    hold on
    plot(real(data(:,t(i)/dt+1)),'-k',LW,1.6)
    plot(real(X_dmd(:,t(i)/dt+1)),'-.r',LW,1.6)
    hold off
    title(['Time = ',num2str(t(i)),' (s)'])
    legend('PDE','DMD')
end



