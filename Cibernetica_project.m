%% Progetto di Cibernetica Fisiologica

%% Angelo Massara, Mirco Calzaretta, Antonio Nakula Brigida
clear all
%% Dinamica dell'epidemia

syms alpha rho_b rho_w p r q_a q_i q_r beta_w sigma_w gamma_w mu_w pi_b ...
    beta_b sigma_b gamma_b mu_b pi_v beta_v sigma_v mu_v 
syms S_b E_b A_b I_b I_bm R_b S_w E_w A_w I_w I_wm R_w S_v E_v I_v
syms b_v eta N_b N_v N_w 

% disease force of infection rates

lambda_w = (beta_w*b_v*I_v)/N_w;
lambda_v = beta_v*b_v*((I_w + rho_w*A_w + eta*(I_b + rho_b*A_b))/(N_w + eta*N_b));
lambda_b = (eta*beta_b*b_v*I_v)/N_b;

S_b_dot = pi_b - q_a*pi_b*A_w - q_i*pi_b*I_w - q_r*pi_b*R_w - ...
    lambda_b*S_b - (alpha + mu_b)*S_b;
E_b_dot = lambda_b*S_b - (alpha + sigma_b + mu_b)*E_b;
A_b_dot = q_a*pi_b*A_w + (1 - p)*sigma_b*E_b - (alpha + gamma_b + mu_b)*A_b;
I_b_dot = q_i*pi_b*I_w + p*sigma_b*E_b - (alpha +gamma_b + mu_b)*I_b;
I_bm_dot = r*q_r*pi_b*R_w - (alpha + mu_b)*I_bm;
R_b_dot = (1-r)*q_r*pi_b*R_w + gamma_b*A_b + gamma_b*I_b - (alpha + mu_b)*R_b;
S_w_dot = alpha*S_b - lambda_w*S_w - mu_w*S_w;
E_w_dot = lambda_w*S_w - (sigma_w + mu_v)*E_w;
A_w_dot = (1-p)*sigma_w*E_w - (gamma_w + mu_w)*A_w;
I_w_dot = p*sigma_w*E_w - (gamma_w + mu_w)*I_w;
I_wm_dot = alpha*I_bm - mu_w*I_wm;
R_w_dot = alpha*R_b + gamma_w*A_w + gamma_w*I_w - mu_w*R_w;
S_v_dot = pi_v - lambda_v*S_v - mu_v*S_v;
E_v_dot = lambda_v*S_v - (mu_v + sigma_v)*E_v;
I_v_dot = sigma_v*E_v - mu_v*I_v;

csi_dot = [S_b_dot E_b_dot A_b_dot I_b_dot I_bm_dot R_b_dot S_w_dot E_w_dot ...
    A_w_dot I_w_dot I_wm_dot R_w_dot S_v_dot E_v_dot I_v_dot];
csi_eq = [pi_b/(alpha + mu_b) 0 0 0 0 0 alpha*pi_b/(mu_w*(alpha + mu_b)) ...
    0 0 0 0 0 pi_v/mu_v 0 0];


%% Valori parametri

% N_b = S_b + E_b + A_b +I_b + I_bm + R_b;
% N_w = S_w + E_w + A_w + I_w + I_wm + R_w;
% N_v = S_v + E_v + I_v;

%% Calcolo matrici Jacobiane 

csi_dot = [ E_b_dot A_b_dot I_b_dot I_bm_dot  E_w_dot ...
    A_w_dot I_w_dot I_wm_dot  E_v_dot I_v_dot];


F = jacobian(csi_dot,[E_b A_b I_b I_bm E_w ...
    A_w I_w I_wm  E_v I_v]);

F = subs(F, [S_b E_b A_b I_b I_bm R_b S_w E_w ...
   A_w I_w I_wm R_w S_v E_v I_v], csi_eq);

V = jacobian([S_b_dot E_b_dot A_b_dot R_b_dot S_w_dot E_w_dot ...
    A_w_dot R_w_dot S_v_dot E_v_dot], [S_b E_b A_b R_b S_w E_w ...
    A_w R_w S_v E_v]);


%% Definizione parametri

mos_red = input(['Mosquito Reduction: Select 1 for no control, 2 for low ' ...
    'control, 3 for high control: ']); 
pers_prot = input(['Personal Protection: Select 1 for no control, 2 for ' ...
    'low control, 3 for high control: ']);
del_preg = input(['Delayed Pregnancy: Select 1 for no control, 2 for ' ...
    'low control, 3 for high control: ']);

alpha = 1/(16*365);     % Tasso di maturazione
rho_b = 0.5;            % Modula l'infettività dei bambini
rho_w = 0.5;            % Modula l'infettività degli adulti 
beta_w = 0.33;          % Probabilità di trasmissione degli adulti
sigma_w = 1/7.5;        % Tasso di progressione degli adulti esposti
gamma_w = 1/8.5;        % Tasso di guarigione adulti
mu_w = 1/(70*365);      % Tasso di morte naturale degli adulti 
pi_b = 1/(15*365);      % Tasso di natalità dei bambini/
beta_b = 0.33;          % Probabilità di trasmissione dei bambini
sigma_b = 1/7.5;        % Tasso di progressione dei bambini esposti 
gamma_b = 1/8.5;        % Tasso di guarigione bambini
mu_b = 1/(18.6*365);    % Tasso di morte naturale dei bambini 
pi_v = 250;             % Tasso di nascita zanzare
beta_v = 0.33;          % Probabilità di infezione tramite zanare
b_v = 0.5;              % Tasso di infettività per morsi di zanzare
sigma_v = 1/3.5;        % Tasso di progressione nelle zanzare
mu_v = 1/14;            % Tasso di morte delle zanzare
p = 0.5;                % Frazione di adulti e bambini nati asintomatici        
eta = 0.8;              % Modification parameter
q_a = 0.5;              % Frazione di bambini nati asintomatici 
q_r = 0.5;              % Frazione di bambini nati rimossi
q_i = 0.5;              % Frazione di bambini nati infetti
r = 0.5;                % Frazione di bambini nati senza microcefalia

if(mos_red == 1)
    mu_v = 1/21;
    pi_v = 500;
elseif(mos_red == 2)
    mu_v = 1/14;
    pi_v = 250;
elseif(mos_red == 3)
    mu_v = 1/8;
    pi_v = 125;
end

if(pers_prot == 1)
    b_v = 0.5;
elseif(pers_prot == 2)
    b_v = 0.25;
elseif(pers_prot == 3)
    b_v = 0.125;
end

if(del_preg == 1)
    pi_b = 1/(15*365); 
elseif(del_preg == 2)
    pi_b = (1/(15*365))/2; 
elseif(del_preg == 3)
    pi_b = (1/(15*365))/10; 
end
csi_eq = [pi_b/(alpha + mu_b) 0 0 0 0 0 alpha*pi_b/(mu_w*(alpha + mu_b)) ...
    0 0 0 0 0 pi_v/mu_v 0 0];

%% Condizioni iniziali 

S_b_0 = 8000;
E_b_0 = 0;
A_b_0 = 0;
I_b_0 = 1;
I_bm_0 = 0;
R_b_0 = 0;
S_w_0 = 8000;
E_w_0 = 0;
A_w_0 = 0;
I_w_0 = 1;
I_wm_0 = 0;
R_w_0 = 0;
S_v_0 = 0.225*10^6;
E_v_0 = 0;
I_v_0 = 1;
N_b_0 = S_b_0 + E_b_0 + A_b_0 + I_b_0;
N_w_0 = S_w_0 + A_w_0 + I_w_0;
N_v_0 = S_v_0 + I_v_0;

N_b_eq = csi_eq(1);
N_w_eq = csi_eq(7);
N_v_eq = csi_eq(13);

% Vettori dei parametri

P = [alpha, rho_b, rho_w, beta_w, sigma_w, gamma_w, mu_w, pi_b, beta_b, ...
    sigma_b, gamma_b, mu_b, pi_v, beta_v, b_v, sigma_v, mu_v, p, eta, ...
    q_a, q_r, q_i, r];


%% Calcolo numero di ROSS

k1 = alpha + mu_b;
k2 = alpha + sigma_b + mu_b;
k3 = alpha + gamma_b + mu_b;
k4 = k3;
k8 = sigma_w + mu_w;
k9 = gamma_w + mu_w;
k10 = k9;
k14 = mu_v + sigma_v;
k5 = alpha + mu_b;
k12 = mu_w;

F_ross = [0 0 0 0 0 0 0 0 0 eta*beta_b*b_v
          0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 beta_w*b_v
          0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 
          0 0 0 0 0 0 0 0 0 0 
          0 eta*rho_b*beta_v*b_v*csi_eq(13)/(csi_eq(7)+eta*csi_eq(1))...
          eta*beta_v*b_v*csi_eq(13)/(csi_eq(7)+eta*csi_eq(1)) 0 0 ...
          rho_w*beta_v*b_v*csi_eq(13)/(csi_eq(7)+eta*csi_eq(1)) ...
          beta_v*b_v*csi_eq(13)/(csi_eq(7)+eta*csi_eq(1)) 0 0 0
          0 0 0 0 0 0 0 0 0 0];

V_ross = [k2 0 0 0 0 0 0 0 0 0
       -(1-p)*sigma_b k3 0 0 0 -q_a*pi_b 0 0 0 0
       -p*sigma_b 0 k4 0 0 0 -q_i*pi_b 0 0 0
       0 0 0 k5 0 0 0 0 0 0 
       0 0 0 0 k8 0 0 0 0 0 
       0 0 0 0 -(1-p)*sigma_w k9 0 0 0 0
       0 0 0 0 -p*sigma_w 0 k10 0 0 0
       0 0 0 -alpha 0 0 0 k12 0 0
       0 0 0 0 0 0 0 0 k14 0 
       0 0 0 0 0 0 0 0 -sigma_v mu_v];

Ross_0 = max(abs(eig(F_ross*V_ross^-1)))


ross_v = (csi_eq(13)*beta_v*b_v*sigma_v)/(k14*mu_v);
ross_b = (eta^2 * beta_b*b_v*sigma_b*(rho_b*(1-p)*k4 + p*k3))/(k2*k3*k4*(csi_eq(7) + eta*csi_eq(1)));
ross_w = (beta_w*b_v*sigma_w)/(k8*(csi_eq(7)+eta*csi_eq(1)))*...
    ((1-p)*(k3*rho_w + eta*rho_b*q_a*pi_b)/(k3*k9) + p*(k4 + eta*q_i*pi_b)/(k4*k10));

ross_0 = sqrt(ross_v*(ross_w + ross_b));



% % % Definizione parametri
% 
% alpha =  1/(18*365) + (1/(15*365) - 1/(18*365))*rand();               % Tasso di maturazione
% rho_b =  0.1 + (1 - 0.1)*rand();                                      % Modula l'infettività dei bambini
% rho_w =  0.1 + (1 - 0.1)*rand();                                      % Modula l'infettività degli adulti 
% beta_w = 0.1 + (0.75 - 0.1)*rand();                                   % Probabilità di trasmissione degli adulti
% sigma_w = 1/12 + (1/3 - 1/12)*rand();                                 % Tasso di progressione degli adulti esposti
% gamma_w = 1/14 + (1/3 - 1/14)*rand();                                 % Tasso di guarigione adulti
% mu_w = 1/(76*365) + (1/(68*365) - (1/(76*365)))*rand();               % Tasso di morte naturale degli adulti 
% pi_b = 1/(18*365) + (1/(12*365) - 1/(18*365))*rand();                 % Tasso di natalità dei bambini
% beta_b = 0.001 + (0.54 - 0.001)*rand();                               % Probabilità di trasmissione dei bambini
% sigma_b = 1/12 + (1/3 - 1/12)*rand();                                 % Tasso di progressione dei bambini esposti
% gamma_b = 1/14 + (1/3 - 1/14)*rand();                                 % Tasso di guarigione bambini
% mu_b = 1/(22.32*365) + (1/(14.88*365) - (1/(22.32*365)))*rand();      % Tasso di morte naturale dei bambini 
% pi_v = 50 + (5000 - 50)*rand();                                       % Tasso di nascita zanzare
% beta_v = 0.1 + (0.75 - 0.1)*rand();                                   % Probabilità di infezione tramite zanare
% b_v = 0.33 + (1 - 0.33)*rand();                                       % Tasso di infettività per morsi di zanzare
% sigma_v = 1/6 + (1/2 - 1/6)*rand();
% mu_v = 1/42 + (1/8 - 1/42)*rand();                                    % Tasso di morte delle zanzare
% p = 0.2;
% eta = 0.2;
% q_a = 0.5;
% q_r = 0.5;
% q_i = 0.5;
% r = 0.5;
