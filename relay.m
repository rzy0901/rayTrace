tx_re_attenuation = abs(sum(fade_coefficients));
tx_re_phase_shift = angle(sum(fade_coefficients));

re_signal_receive = signal_receive*exp(-1j*tx_re_phase_shift);

relay_real = relay_receive_gain*real(re_signal_receive);
relay_imag = relay_receive_gain*imag(re_signal_receive);

% figure(6)
% subplot(2,1,1)
% plot(1000*time,relay_real,'Color',blue)
% title('Real relay received baseband signal')
% xlabel('Time (ms)')
% ylabel('Amplitude')
% subplot(2,1,2)
% plot(1000*time,relay_imag,'Color',red)
% title('Imaginary relay received baseband signal')
% xlabel('Time (ms)')
% ylabel('Amplitude')

signal_transmit = relay_amplitude*relay_transmit_gain*(relay_real+1j.*relay_imag);
re_signal_transmit = signal_transmit;

message_receive_real = zeros(1,modulated_length);
message_receive_imag = zeros(1,modulated_length);

for i = 1:modulated_length
    message_receive_real(i) = relay_real(i*frequency_sample/bandwidth);
    message_receive_imag(i) = relay_imag(i*frequency_sample/bandwidth);
end

re_message_receive = message_receive_real+1j*message_receive_imag;

% figure(8)
% hold on
% plot(real(message_modulated),imag(message_modulated),'x','Color',blue,'MarkerSize',10,'LineWidth',2)
% plot(real(re_message_receive),imag(re_message_receive),'x','Color',red,'MarkerSize',10)
% title('Signal Constellations')
% legend('Transmitted constellation','Received constellation','Location','SouthEast')
% axis([min(real(message_modulated))-2,max(real(message_modulated))+2,min(imag(message_modulated))-2,max(imag(message_modulated))+2])
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';

origin = relay_pos;
termination = receive_pos;

figure_counter = figure_counter + 1;