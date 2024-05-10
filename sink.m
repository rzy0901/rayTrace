power_transmit = 0;
power_receive = 0;
power_relay_receive = 0;
power_relay_transmit = 0;
for i = 1:length(time)
    power_transmit = power_transmit + abs(tx_signal(i)).*abs(tx_signal(i));
    power_receive = power_receive + abs(rx_signal_receive(i)).*abs(rx_signal_receive(i));
    power_relay_receive = power_relay_receive + abs(re_signal_receive(i)).*abs(re_signal_receive(i));
    power_relay_transmit = power_relay_transmit + abs(re_signal_transmit(i)).*abs(re_signal_transmit(i));
end
power_transmit = power_transmit/length(time);
power_receive = power_receive/length(time);
power_relay_receive = power_relay_receive/length(time);
power_relay_transmit = power_relay_transmit/length(time);

SNR = power_receive/(noise_var*noise_var);
SNR_dB = 10*log10(SNR);

path_loss_tx_re = power_transmit/power_relay_receive;
path_loss_tx_re_dB = 10*log10(path_loss_tx_re);
path_loss_re_rx = power_relay_transmit/power_receive;
path_loss_re_rx_dB = 10*log10(path_loss_re_rx);

delay_mean = sum(ray_matrix(:,1).*ray_matrix(:,2))/sum(ray_matrix(:,2));
delay_rms = sqrt(sum(ray_matrix(:,2).*(ray_matrix(:,1)-delay_mean).^2)/sum(ray_matrix(:,2)));

fprintf('Signal-to-noise ratio:\t%.2f\tdB\n',SNR_dB)
fprintf('Path loss (Tx to R):\t%.2f\tdB\n',path_loss_tx_re_dB)
fprintf('Path loss (R to Rx):\t%.2f\tdB\n',path_loss_re_rx_dB)
fprintf('Delay spread:\t\t\t%.3f\tns\n',1e9*delay_rms)

% for i = 1:9
%     figure(i)
% end