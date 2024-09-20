%% calculate the range in meters of four targets with respective measured beat frequencies [0 MHz, 1.1 MHz, 13 MHz, 24 MHz]
%% Give: The radar maximum range = 300m and range resolution = 1m

%% Range resolution solely dependent on the bandwidth of the chirp B_sweep:
% dist_res = c / (2 * B_sweep)


% TODO : Find the Bsweep of chirp for 1 m resolution

% Speed of light
c = 3*10^8;
% Range resolution 
dist_res = 1;

B_sweep = c /(2*dist_res);


% TODO : Calculate the chirp time based on the Radar's Max Range
R_max = 300;
Ts = 5.5 * 2 * R_max / c;


% TODO : define the frequency shifts
beat_freq = [0, 1.1e6, 13e6, 24e6];
calculated_range = c * Ts * beat_freq / (2 * B_sweep);


% Display the calculated range
disp(calculated_range);