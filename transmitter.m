frequency_sample = sample_factor*bandwidth;
time = 0:1/frequency_sample:(modulated_length+1)/bandwidth;

message_real = real(message_modulated);
message_imag = imag(message_modulated);

tx_transmit_real = zeros(1,length(time));
tx_transmit_imag = zeros(1,length(time));

for i = 1:modulated_length
    tx_transmit_real = tx_transmit_real + message_real(i)*sinc(bandwidth*time-i);
    tx_transmit_imag = tx_transmit_imag + message_imag(i)*sinc(bandwidth*time-i);
end

signal_transmit = transmit_amplitude*transmit_gain*(tx_transmit_real+1j.*tx_transmit_imag);
tx_signal = signal_transmit;

% figure(5)
% subplot(2,1,1)
% hold on
% title('Real transmitted baseband signal')
% plot(1000*time,tx_transmit_real,'Color',blue)
% plot(1000*(1:modulated_length)/bandwidth,message_real,'.','Color',blue,'MarkerSize',10)
% xlabel('Time (ms)')
% ylabel('Amplitude')
% subplot(2,1,2)
% hold on
% title('Imaginary transmitted baseband signal')
% plot(1000*time,tx_transmit_imag,'Color',red)
% plot(1000*(1:modulated_length)/bandwidth,message_imag,'.','Color',red,'MarkerSize',10)
% xlabel('Time (ms)')
% ylabel('Amplitude')

origin = transmit_pos;
termination = relay_pos;