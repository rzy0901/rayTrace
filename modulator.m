% Convert message to passband signal using specified modulation scheme.

switch modulation
    case 1
        modulated_length = message_length;
        message_modulated = zeros(1,modulated_length);
        for i = 1:modulated_length
            switch message(i)
                case 1
                    message_modulated(i) = 1+0j;
                case 0
                    message_modulated(i) = -1+0j;
                otherwise
                    message_modulated(i) = 0+0j;
            end
        end
    case 2
        modulated_length = message_length/2;
        message_modulated = zeros(1,modulated_length);
        for i = 1:modulated_length
            switch int2str(message(2*i-1:2*i))
                case '0  0'
                    message_modulated(i) = 1/sqrt(2)+1j/sqrt(2);
                case '0  1'
                    message_modulated(i) = -1/sqrt(2)+1j/sqrt(2);
                case '1  1'
                    message_modulated(i) = -1/sqrt(2)-1j/sqrt(2);
                case '1  0'
                    message_modulated(i) = 1/sqrt(2)-1j/sqrt(2);
                otherwise
                    message_modulated(i) = 0+0j;
            end
        end
    case 3
        modulated_length = message_length/4;
        message_modulated = zeros(1,modulated_length);
        for i = 1:modulated_length
            switch int2str(message(4*i-3:4*i))
                case '0  0  0  0'
                    message_modulated(i) = 1/sqrt(11/2);
                case '0  0  0  1'
                    message_modulated(i) = 1j/sqrt(11/2);
                case '0  0  1  0'
                    message_modulated(i) = -1/sqrt(11/2);
                case '0  0  1  1'
                    message_modulated(i) = -1j/sqrt(11/2);
                case '0  1  0  0'
                    message_modulated(i) = 3/(2*sqrt(11/2))+1j*3/(2*sqrt(11/2));
                case '0  1  0  1'
                    message_modulated(i) = -3/(2*sqrt(11/2))+1j*3/(2*sqrt(11/2));
                case '0  1  1  0'
                    message_modulated(i) = -3/(2*sqrt(11/2))-1j*3/(2*sqrt(11/2));
                case '0  1  1  1'
                    message_modulated(i) = 3/(2*sqrt(11/2))-1j*3/(2*sqrt(11/2));
                case '1  0  0  0'
                    message_modulated(i) = 2/sqrt(11/2);
                case '1  0  0  1'
                    message_modulated(i) = (1j*2)/sqrt(11/2);
                case '1  0  1  0'
                    message_modulated(i) = -2/sqrt(11/2);
                case '1  0  1  1'
                    message_modulated(i) = (-1j*2)/sqrt(11/2);
                case '1  1  0  0'
                    message_modulated(i) = 5/(2*sqrt(11/2))+1j*5/(2*sqrt(11/2));
                case '1  1  0  1'
                    message_modulated(i) = -5/(2*sqrt(11/2))+1j*5/(2*sqrt(11/2));
                case '1  1  1  0'
                    message_modulated(i) = -5/(2*sqrt(11/2))-1j*5/(2*sqrt(11/2));
                case '1  1  1  1'
                    message_modulated(i) = 5/(2*sqrt(11/2))-1j*5/(2*sqrt(11/2));
                otherwise
                    message_modulated(i) = 0+0j;
            end
        end
    otherwise
        fprintf('Specify correct modulation scheme.')
end