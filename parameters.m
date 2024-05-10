% Set user-defined simulation parameters.

% Distance between points along ray:
ray_resolution = 0.1;
% Number of points to draw objects:
object_resolution = 40;

% Maximum number of ray collisions:
ray_collisions = 3;

% Modulation scheme:
% 1 = BPSK, 2 = 4-QPSK, 3 = 16-QAM.
modulation = 2;
% Message length:
message_length = 16;           % 256 bits

% Bandwidth:
bandwidth = 1*10^6;             % 1 MHz
% Sample factor:
sample_factor = 10000;          % 10 kS

% Carrier frequency:
frequency_carrier = 2.4*10^9;   % 2.4 GHz

% Transmit signal amplitude:
transmit_amplitude = 1;         % 1 V
% Receive signal amplitude:
relay_amplitude = 1;            % 1 V

% Transmitter gain:
transmit_gain = 10;             % 10 linear == 10 dB
% Receiver gain:
receive_gain = 0;              % 10 linear == 10 dB
% Relay receiver gain:
relay_receive_gain = 10;        % 10 linear == 10 dB
% Relay transmitter gain:
relay_transmit_gain = 10;       % 10 linear == 10 dB

% Noise variance:
noise_var = 0.01;

% Set environment geometry and material properties.

% Transmitter and receiver positions:
transmit_pos = [3,5.5,1];
receive_pos = [3,0.5,1];
relay_pos = receive_pos;

% Geometric properties of objects:
% Normal, point and bounds.

object_geometry(1,1:12) = [[0,0,1],[3,3,0],[0,6],[0,6],[-99,99]];
%object_geometry(2,1:12) = [[0,1,0],[3,0,1],[0,6],[-99,99],[0,2]];
object_geometry(3,1:12) = [[0,1,0],[3,6,1],[0,6],[-99,99],[0,2.5]];
object_geometry(4,1:12) = [[1,0,0],[0,3,1],[-99,99],[0,6],[0,2.5]];
object_geometry(5,1:12) = [[1,0,0],[6,3,1],[-99,99],[0,6],[0,2.5]];
%object_geometry(6,1:12) = [[1,0,0],[3,3,1],[-99,99],[2.5,3.5],[0.8,1.2]];

% object_geometry(1,1:12) = [[1,0,0],[4,0,0],[-99,99],[0,4],[0,4]];
% object_geometry(2,1:12) = [[0,1,0],[0,4,0],[0,4],[-99,99],[0,4]];
% object_geometry(3,1:12) = [[0,0,1],[0,0,0],[0,4],[0,4],[-99,99]];

% Matrial properties of objects;
% Reflectance.
object_material(1,1) = 0.4;
object_material(2,1) = 0.4;
object_material(3,1) = 0.4;
object_material(4,1) = 0.4;
object_material(5,1) = 0.4;
object_material(6,1) = 0.4;

[object_number,~] = size(object_geometry);

% indicident points
start_points = [];
end_points = [];

% Generate a human.
ped = backscatterPedestrian('Height',1.7,'WalkingSpeed',1,'InitialPosition',[1;3;0],'OperatingFrequency',frequency_carrier);
% ped2 = backscatterPedestrian('Height',1.7,'WalkingSpeed',1,'InitialPosition',[5;5;0],'OperatingFrequency',frequency_carrier);
dt = 0.1;