%% Grafici progetto Cib Fis
% 
% cum_low_i_b = load("Data no control.mat","cum_i_b").cum_i_b;
% cum_mod_i_b = load("Data moderate control.mat","cum_i_b").cum_i_b;
% cum_max_i_b = load("Data max control.mat","cum_i_b").cum_i_b;
% 
% figure
% plot(cum_mod_i_b);
% hold on;
% plot(cum_low_i_b);
% plot(cum_max_i_b);
% xlim([0 275]);
cum_i_b = out.New_inf_b;                     % Bambini infetti cumulativi
cum_i_bm = out.New_inf_bm;                   % Bambini con microcefalia
cum_i_w = out.New_inf_w;                     % Adulti infetti cumulativi

%% Plot
if mu_b ~= 0

    figure
    plot(cum_i_b, LineWidth=1.5);
    title('Bambini Infetti Cumulativi')
    xlabel('Giorni')
    ylabel('Numero Casi')
    xlim([0 275]);
    grid on
    figure
    plot(cum_i_w, LineWidth=1.5);
    title('Adulti Infetti Cumulativi')
    xlabel('Giorni')
    ylabel('Numero Casi')
    xlim([0 275]);
    grid on
    figure
    plot(cum_i_bm, LineWidth=1.5);
    title('Bambini con Microcefalia Cumulativi')
    xlabel('Giorni')
    ylabel('Numero Casi')
    xlim([0 275]);
    grid on
    
    if mos_red == 1
        save("Data no control.mat","cum_i_b");
    elseif mos_red == 2
        save("Data low control.mat","cum_i_b");
    elseif mos_red == 3
        save("Data max control.mat","cum_i_b");
    end
else
    figure
    plot(cum_i_b);
    title('Bambini Infetti Cumulativi senza parametri vitali')
    xlim([0 275]);
    grid on
    figure
    plot(cum_i_w);
    title('Adulti Infetti Cumulativi senza parametri vitali')
    xlim([0 275]);
    grid on
    figure
    plot(cum_i_bm);
    title('Bambini con Microcefalia Cumulativi senza parametri vitali')
    xlim([0 275]);
    grid on
    
    if mos_red == 1
        save("Data no control.mat","cum_i_b");
    elseif mos_red == 2
        save("Data low control.mat","cum_i_b");
    elseif mos_red == 3
        save("Data max control.mat","cum_i_b");
    end

end