clear all;
length = 1000000; % Number of points in wave



%% cantilever properties 

k = 0.5 % N/m
f_0 = 25000     ; % Resonant frequency
Q = 3; % Qualiy factor
temperature = i_temp; % Kelvin
m = k/(2*pi*f_0)^2;
b=sqrt(k*m)/Q;
% timestep (dt) depends on the resonance frequency of the cantilever:
sample_pts =500; %sample_pts is the number of points in each period 
dt = 1/(sample_pts*f_0);
dt=500e-9;
time_length=dt*length;

%% equilibrium_position
equilibrium_position_begin = 50e-10; %starting position of the simulation
equilibrium_position_end = 3e-10; % ending position of the simulation
equilibrium_position=equilibrium_position_begin;

%% force field 
x=linspace(equilibrium_position_end,equilibrium_position_begin,5000);
for i=1:numel(x)
simulation_force(i)=100*(1e-11/x(i))^6-.0006*(1e-11/x(i))^3; 
%simulation_energy(i)=energytrace(x(i))/(1.38e-23*300);
end
%forcetrace()=simulation_force()
%%  driving force (set as 0 in this case)
driving_frequency = f_0;
driving_amplitude = 0; %  driving amplitude if at resonance (5e-9)
drive_force=driving_amplitude*k*(driving_frequency/f_0)*sqrt((f_0/driving_frequency-driving_frequency/f_0)^2+(1/Q)^2);

%% initial position and velocity

% position=zeros(length,1);
% impact=zeros(length/sample_pts,1); % velocity
vel = 0; 
acc = 0;
position_temp = equilibrium_position;

%% Loop
vel_turn=0; % vel_meas=0
counter=0;
oversampling=10;
oversamploop=0;
counter2=0;

while (counter2<(length-1))
    vel_1=vel;
    %%xindex=find(x,position_temp);
    [~,mindex]=min(abs(x-position_temp));
    forcetrace_temp=simulation_force(mindex);
    %%forcetrace_temp=forcetrace(position_temp);
    equilibrium_position=equilibrium_position_begin+(equilibrium_position_end-equilibrium_position_begin)*counter/(length*oversampling);
    F_noise=normrnd(0,(sqrt(temperature*1.38e-23*2*k/(pi*Q*2*dt*f_0))));
    Fdrive=drive_force*sin(2*pi*counter*dt*driving_frequency); % a sinewave at the driving frequency with amplitude
    acc=1/m *(-k*(position_temp-equilibrium_position) - b*vel + Fdrive+forcetrace_temp+F_noise);
    vel=vel+acc*dt;
    position_temp=position_temp+vel*dt;
    counter=counter+1;
    if(vel*vel_1<0 & vel>0)
        impact(vel_turn+1)=forcetrace_temp;
        vel_turn=vel_turn+1;
    end
    if(oversamploop==oversampling)
        position(counter2+1)=position_temp;
        eqb_pos(counter2+1)=equilibrium_position;
        force_profile(counter2+1)=forcetrace_temp;
        counter2=counter2+1;
        oversamploop=0;
    end
    oversamploop=oversamploop+1;
end

for p=1:numel(position)
    scaled_position(p)=position(p)-((equilibrium_position_end-equilibrium_position_begin)*p/length + equilibrium_position_begin);
end

%save('Simulation_single_well_variabless','scaled_position','eqb_pos','x','simulation_force');
scaled_position=scaled_position*10^9;
eqb_pos=eqb_pos*10^9;

csvwrite(strcat('pos_k0.5_T',int2str(i_temp),'_DFS',int2str(i_rep),'.csv'),eqb_pos.');
csvwrite(strcat('defl_k0.5_T',int2str(i_temp),'_DFS',int2str(i_rep),'.csv'),scaled_position.');

%figure(1);
%plot(eqb_pos,scaled_position,'red'); % plots deflection versus equilibrium z-position
%set(gca,'FontSize',36,'FontWeight','Bold');
%xlabel('Position (nm)','FontSize',36,'FontWeight','Bold')
%ylabel('Deflection (nm)','FontSize',36,'FontWeight','Bold')