%% Task: calculate the velocity in m/s of four targets with the following doppler frequency shifts: [3 KHz, -4.5 KHz, 11 KHz, -3 KHz]

% Doppler Velocity Calculation
c = 3*10^8;         %speed of light
frequency = 77e9;   %frequency in Hz

% TODO: Calculate the wavelength
% lamda = c / f

lambda = c / frequency;


% TODO: Define the doppler shifts in Hz using the information from above 
fd = [3e3, -4.5e3, 11e3, -3e3];

% TODO: Calculate the velocity of the targets  fd = 2*vr/lambda
vr = (fd * lambda) / 2;


% TODO: Display results
disp(vr);