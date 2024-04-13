mystartdefaults

tic

recipunit= 1.0E+10;
ekinscale= ((hbar*recipunit)^2/(2*elm))/qel;

%Ex:1.1:

% defining the step size:

step_size=0.5; %Angstroms

% descritizing the first interval bound of U [0, 15]:

bound1_start=0;
bound1_end=15;
bound_range=bound1_start:step_size:bound1_end;

%discretizing the second interval bound of U: [65,80]:
start_bound2=65;
end_bound2=80;
discrete_bound2= start_bound2:step_size:end_bound2;

% conversion
convert_to_j=1.60218e-19; % conversion factor from ev to joules
% convert grid step to meters for calculations:

grid_step_m=step_size* 1e-10; %1 angstrom

% % finding maximum momentum to avoid aliasing:
% % correspond to a de Broglie wavelength of twice the grid step size
% lambda_min=2 * grid_step_m; %minimum wavelength
% p_max=hbar*(2*pi) /lambda_min; % convert energy to ev
% 
% % maximum energy (upper bound Eb):
% Eb_j=p_max^2 / (2*elm); % energy in joules
% Eb_ev=Eb_j/ convert_to_j; % convert to ev
% 
% % lower bound Ea:
% Ea_ev=0;
% 
% fprintf('Lower energy bound (Ea) = %f eV\n', Ea_ev);
% fprintf('Lower energy bound (Ea) = %f eV\n', Eb_ev);


%% 
% Ex: 1.2:

% Energy range and step size
E_min = 0; % in eV
E_max = 0.3; % in eV
deltaE = 0.0005; % in eV
step_size = 0.5;
x1_min = -20;
x1_max = 100;
tau = 1e-9; % 1 ns in seconds
gam=(hbar*2*pi/tau)/qel; % damping factor
recipunit=1.0E+10;
ekinscale=(hbar*recipunit)^2/(2*elm)/qel;
n = (x1_max - x1_min)/step_size;
%% 
% Creating potential U

U = zeros(1,n); % Potential barrier
x1 = zeros(1,n);
for i = 1:n
    x1(i) = x1_min + step_size/2 + (i - 1) * step_size;
    if x1(i) <= 15 && x1(i)>=0
        U(i) = 0.2;
    elseif x1(i) >= 65 && x1(i) <= 80
        U(i) = 0.2;
    else
        U(i) = 0;
    end
end
plot(x1,U,'LineWidth',2)
xlabel('X(Angsroms)');
ylabel('E(eV)');
title('Potential Barrier')
ylim([0,0.22])
hold on;