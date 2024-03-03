% S_b_plot = out.S_b_plot;
% E_b_plot = out.E_b_plot;
% A_b_plot = out.A_b_plot;
% I_b_plot = out.I_b_plot;
% I_bm_plot = out.I_bm_plot;
% R_b_plot = out.R_b_plot;
% 
% S_w_plot = out.S_w_plot;
% E_w_plot = out.E_w_plot;
% A_w_plot = out.A_w_plot;
% I_w_plot = out.I_w_plot;
% I_wm_plot = out.I_wm_plot;
% R_w_plot = out.R_w_plot;
% 
% S_v_plot = out.S_v_plot;
% E_v_plot = out.E_v_plot;
% I_v_plot = out.I_v_plot;
% 
% t = 275;
% dt = 1;
% 
% %% Andamento Suscettibili
% figure
% plot(0:dt:t, S_b_plot, LineWidth=1.4)
% grid on
% hold on
% title ("Andamento Suscettibili")
% xlabel("Giorni")
% ylabel("Numero Suscettibili")
% xlim([0 275]);
% plot(0:dt:t, S_w_plot, LineWidth=1.4)
% legend ('Bambini', 'Adulti')
% hold off
% %% Andamento Esposti
% figure
% plot(0:dt:t, E_b_plot, LineWidth=1.4)
% grid on
% hold on
% title ("Andamento Esposti")
% xlabel("Giorni")
% ylabel("Numero Esposti")
% xlim([0 275]);
% plot(0:dt:t, E_w_plot, LineWidth=1.4)
% legend ('Bambini', 'Adulti')
% hold off
% 
% %% Andamento Asintomatici 
% 
% figure
% plot(0:dt:t, A_b_plot, LineWidth=1.4)
% grid on
% hold on
% title ("Andamento Asintomatici")
% xlabel("Giorni")
% ylabel("Numero Asintomatici")
% xlim([0 275]);
% plot(0:dt:t, A_w_plot, LineWidth=1.4)
% legend ('Bambini', 'Adulti')
% hold off
% %% Andamento Infetti  
% 
% figure
% plot(0:dt:t, I_b_plot, LineWidth=1.4)
% grid on
% hold on
% title ("Andamento Infetti")
% xlabel("Giorni")
% ylabel("Numero Infetti")
% xlim([0 275]);
% plot(0:dt:t, I_w_plot, LineWidth=1.4)
% legend ('Bambini', 'Adulti')
% hold off
% %% Andamento Infetti Microcefalia 
% 
% figure
% plot(0:dt:t, I_bm_plot, LineWidth=1.4)
% grid on
% hold on
% title ("Andamento Infetti con Microcefalia")
% xlabel("Giorni")
% ylabel("Numero Infetti con microcefalia")
% xlim([0 275]);
% plot(0:dt:t, I_wm_plot, LineWidth=1.4)
% legend ('Bambini', 'Adulti')
% hold off
% %% Andamento Rimossi 
% 
% figure
% plot(0:dt:t, R_b_plot, LineWidth=1.4)
% grid on
% hold on
% title ("Andamento Rimossi")
% xlabel("Giorni")
% ylabel("Numero Rimossi")
% xlim([0 275]);
% plot(0:dt:t, R_w_plot, LineWidth=1.4)
% legend ('Bambini', 'Adulti')
% hold off
% 
% %% Andamento Zanzare
% 
% figure
% grid on
% title('Andamento Zanzare Suscettibili')
% xlabel('Giorni')
% ylabel('Numero Zanzare')
% xlim([0 275])
% plot(0:dt:t, S_v_plot, LineWidth=1.4)
% figure
% grid on
% hold on 
% title('Andamento Zanzare Esposte/Infette')
% xlim([0 275])
% plot(0:dt:t, E_v_plot, LineWidth=1.4)
% plot(0:dt:t, I_v_plot, LineWidth=1.4)
% legend('Esposte', 'Infette')
% 

%% Confronto Zika-Febbre Gialla

cum_i_h_fg = out.cum_inf_hum_fg;             % Umani infetti cumulativi 
                                             % febbre gialla

figure
plot(cum_i_w, LineWidth=1.5);
title('Umani Infetti Cumulativi')
xlim([0 275]);
ylabel('Numero Casi')
xlabel('Giorni')
grid on
hold on 
plot(cum_i_h_fg, LineWidth=1.5);
legend('Zika', 'Febbre Gialla');

dead_zika = out.Dead_Zika;
dead_fg = out.dead_fg;

% figure
% plot(dead_zika, LineWidth=1.5);
% title('Andamento morti')
% xlim([0 275]);
% ylabel('Numero Casi')
% xlabel('Giorni')
% grid on
figure
plot(dead_fg, LineWidth=1.5);
title('Andamento morti')
xlim([0 275]);
ylabel('Numero Casi')
xlabel('Giorni')
grid on
