%% Progetto Cibernetica Fisiologica -- Febbre Gialla

%% Definizione parametri

beta_1_fg = 0.6; 
beta_2_fg = 0.5;
gamma_h_fg = 0.31;
theta_fg = 0.15;
%b_h_fg = 4.94*10^-5;
b_h_fg = 0;
b_v_fg = 0.051;
a_fg = 3;
delta_fg = 0.143;
alpha_fg = 3.5*10^-1;
epsilon_fg = 0.01;
rho_fg = 0.5;
sigma_fg = 10^-6;
mu_h_fg = 4.94*10^-5;
mu_v_fg = 0.051;
sigma_v = 1/3.5;
p_fg = 0.5;
%% Condizioni iniziali

Sh_fg_0 = 0.62;
Eh_fg_0 = 3.4*10^-2;
Th_fg_0 = 6*10^-3;
Ih_fg_0 = 4*10^-2;
Rh_fg_0 = 0.3;
Sv_fg_0 = 0.95;
Ev_fg_0 = 0;
Iv_fg_0 = 0.05;
Nh_fg_0 = Sh_fg_0 + Eh_fg_0 + Th_fg_0 + Ih_fg_0 + Rh_fg_0;
Nv_fg_0 = Sv_fg_0 + Iv_fg_0;

%% Vettore da passare a Simulink

Par = [beta_1_fg beta_2_fg gamma_h_fg theta_fg b_h_fg b_v_fg a_fg delta_fg... 
       alpha_fg epsilon_fg rho_fg sigma_fg mu_h_fg mu_v_fg sigma_v p_fg];

%% Ross Number

Ross_fg = (gamma_h_fg*a_fg^2*beta_1_fg*beta_2_fg*theta_fg*(mu_h_fg-rho_fg ...
    *sigma_fg))/(mu_v_fg*(mu_h_fg+epsilon_fg)*(mu_h_fg+alpha_fg+delta_fg) ...
    *(gamma_h_fg+mu_h_fg));

