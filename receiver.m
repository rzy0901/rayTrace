re_rx_attenuation = abs(sum(fade_coefficients));
re_rx_phase_shift = angle(sum(fade_coefficients));

rx_signal_receive = signal_receive*exp(-1j*re_rx_phase_shift);

rx_receive_real = receive_gain*real(rx_signal_receive);
rx_receive_imag = receive_gain*imag(rx_signal_receive);
rx_signal = rx_receive_real+1j*rx_receive_imag;

% figure(7)
% subplot(2,1,1)
% plot(1000*time,rx_receive_real,'Color',blue)
% title('Real received baseband signal')
% xlabel('Time (ms)')
% ylabel('Amplitude')
% subplot(2,1,2)
% plot(1000*time,rx_receive_imag,'Color',red)
% title('Imaginary received baseband signal')
% xlabel('Time (ms)')
% ylabel('Amplitude')

message_receive_real = zeros(1,modulated_length);
message_receive_imag = zeros(1,modulated_length);

for i = 1:modulated_length
    message_receive_real(i) = rx_receive_real(i*frequency_sample/bandwidth);
    message_receive_imag(i) = rx_receive_imag(i*frequency_sample/bandwidth);
end

rx_message = message_receive_real+1j*message_receive_imag;

% figure(9)
% hold on
% plot(real(message_modulated),imag(message_modulated),'x','Color',blue,'MarkerSize',10,'LineWidth',2)
% plot(real(rx_message),imag(rx_message),'x','Color',red,'MarkerSize',10)
% title('Signal Constellations')
% legend('Transmitted constellation','Received constellation','Location','SouthEast')
% axis([min(real(message_modulated))-2,max(real(message_modulated))+2,min(imag(message_modulated))-2,max(imag(message_modulated))+2])
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';