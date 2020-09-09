% Script to generate corrugated horn profiles for simulation with openEMS.
%
% Formulations and example from "Design of Corrugated Horns: A Primer", IEEE 
% Antennas and Propagation Magazine,Vol.47,No 2,April 2006, by Christophe Granet
% Author: Paul Klasmann
% Date: March 2019
%
% 1) Enter flags for RUN_SIMULATION and RUN_CSXCAD, 1 to run and 0 not to run.
%    It's advised not to run the simulation while setting up the horn's geometry.
% 2) Enter values in the user editable parameters.  
% 3) If using the ring loaded slot mode converter, enter a value for delta2, 
%    otherwise ignore this variable.
% 4) If using the variable pitch to slot width mode converter, enter a value for
%    delta_min, otherwise ignore this variable.
% 5) Select an integer number from 1 to 8 for the horn_profile.
% 6) Select an integer number from 1 to 3 for the mode converter type.
% 7) Length of horn profile is chosen in the length = xx*p line. Enter a value
%    that is multiplied by the pitch "p".  This will give the length of the
%    profiled length of the horn.  Line 81.
% 
% Note 1: If you wish to only generate the complete horn profile you may comment
%       out all lines from 408 if you do not have openEMS installed.  
% Note 2: Set circular waveguide length wgl to something sensible for a given diameter.

close all
clear
clc

RUN_SIMULATION = 0;   % Set to 1 to run openEMS, 0 not to run
RUN_CSXCAD = 1;       % Set to 1 to open CSXCAD to show 3D geometry, 0 not to open CSXCAD

% USER EDITABLE PARAMETERS
fmin = 10.7                % Minimum frequency in GHz
fmax = 14.5                % Maximum frequency in GHz
pitch_fraction = 8  % Choose a fraction between 10 to 5 (lambda_c / pitch_fraction)
delta = 0.8               % Pitch to width ratio 0.7 to 0.9
sigma = 0.42               % Percentage factor for first slot depth, 0.4 to 0.5 
NMC = 5                    % Number of corrugations in mode converter
wgl = 30;                  % Length of circular feeding waveguide
num_of_corrugations = 60;
% delta2 is only used for case 2, ring loaded slot mode converter
delta2 = 0.1;

% delta_min is only used for case 3,variable pitch to width slot mode converter
delta_min = 0.125;           % Should be greater than 0.125 and less than delta

% Choose profile, 1=LINEAR, 2=SINUSOID, 3=ASYMMETRIC SINE-SQUARED, 4=TANGENTIAL,
% 5=x.rho, 6=EXPONENTIAL, 7=HYPERBOLIC, 8=POLYNOMIAL
horn_profile = 7;          % Must be an integer from 1 to 8

% Choose the type of mode converter, 1=VARIABLE SLOT DEPTH MC,
% 2=RING LOADED SLOTS MC, 3=VARIABLE PITCH TO WIDTH SLOT MC
mode_converter_type = 1;   % Must be an integer from 1 to 3

% END OF USER EDITABLE PARAMETERS

% Calculate center frequency fc based on narrow or wide bandwidth.
fratio = fmax/fmin             % ratio of fmax/fmin
if (fratio >= 2.4);            % check fmax/fmin is less than 2.4
  disp('Error, fmax/fmin is greater than 2.4!');
  fc = 0
  elseif (fratio <= 1.4);      % Use Narrowband formula if fmax <= 1.4fmin
    fc = sqrt(fmin*fmax)
     elseif (fmax >= 1.4*fmin && fmax <= 2.4*fmin); % Use wideband formula for 1.4fmin<=fmax<=2.4fmin
       fc = 1.2*fmin
endif

if (fratio <= 1.4);
  fo = 1.02*fc    % For narrow band choose fc <= fo <= 1.05fc
    elseif (fmax >= 1.4*fmin && fmax <= 2.4*fmin);
      fo = 1.10*fc    % For wideband choose 1.05fc <= fo <= 1.15fc
endif

unit = 1e-3;                    % Units in mm
lambda_c = 300/fc               % Center frequency wavelength
lambda_o = 300/fo               % Output frequency
depth_nominal = lambda_c/4      % Nominal slot depth at center frequency
ai = (3 * lambda_c)/(2*pi)      % Radius of input waveguide in mm
ao = 1.95*lambda_c              % Radius of output waveguide in mm
p = lambda_c/pitch_fraction;    % Pitch in mm, lambda_c/10 to lambda_c/5
length = num_of_corrugations*p  % Length of horn profile
N = length/p                    % Total number of corrugations
kc = (2*pi)/lambda_c            % Wave number at center frequency
ko = (2*pi)/lambda_o            % Wave number at output frequency
z = 0:p:length;                 % z index distance array from 0 to length of horn
r_app = ao*1e-3;                % Aperture radius
A_app = pi*(r_app)^2;           % Aperture area for gain calculation

% One of the following horn profiles will be selected depending on value of
% horn_profile chosen above.
switch(horn_profile)
   case 1
   %%% Linear profile %%%
   a = ai+(ao-ai)*z/length;
   
   plot(z, a);    
   set(gca, "linewidth",2, "fontsize", 14 )
   axis equal;
   xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
   ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
   title( 'Linear Horn Profile', 'FontSize', 16 );

   case 2
   %%% Sinusoid profile %%%
   A = 1;     % Amplitude factor 'A' should be between 0 and 1
   rho = 2;     % rho should be between 0.5 and 5, default is 2
   a = ai+(ao-ai)*((1-A)*(z/length)+A*power(sin((pi*z)/(2*length)),rho));
   
   plot(z, a);    
   set(gca, "linewidth",2, "fontsize", 14 )
   axis equal;
   xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
   ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
   title( 'Sinusoidal Horn Profile', 'FontSize', 16 );

   case 3
   %%% Asymmetric Sine-Squared profile %%%
   L1 = length/3;    % Choose a value for L1, must be less than the horns length
   L2 = length-L1;   % L2 is the length between L1 and the end of the horn
   gamma = L2/L1;
   idx = find(z <= L1)     % Find the index of z corresponding to L1
   zelements = size(z,2)   % Total number of points in z axis of the horn
   
   za = z(1: max(idx))
   aa = ai+((2*(ao-ai))/(1+gamma))*sin((pi*za)/(4*L1)).^2;
    
   zb = z(max(idx)+1 : zelements)
   ab = ai+((2*(ao-ai))/(1+gamma))*(gamma*sin(((pi*(zb+L2-L1))/(4*L2))).^2+((1-gamma)/2));
   
   a = [aa,ab];
   z = [za,zb];
   
   plot(z, a);    
   set(gca, "linewidth",2, "fontsize", 14 )
   axis equal;
   xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
   ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
   title( 'Asymmetric Sine Squared Horn Profile', 'FontSize', 16 );
   
   case 4
   %%% Tangential profile %%%
   A = 1;
   rho = 2;
   a = ai+(ao-ai)*((1-A)*(z/length)+A*power(tan((pi*z)/(4*length)),rho));
   
   plot(z, a);    
   set(gca, "linewidth",2, "fontsize", 14 )
   axis equal;
   xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
   ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
   title( 'Tangential Horn Profile', 'FontSize', 16 );

   case 5
   %%% x.rho profile %%%
   A = 1;
   rho = 2;
   a = ai+(ao-ai)*((1-A)*(z/length)+A*power(z/length,rho));
   
   plot(z, a);    
   set(gca, "linewidth",2, "fontsize", 14 )
   axis equal;
   xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
   ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
   title( 'xp Horn Profile', 'FontSize', 16 );

   case 6
   %%% Exponential profile %%%
   a=ai*exp(log(ao/ai)*(z/length));
   
   plot(z, a);    
   set(gca, "linewidth",2, "fontsize", 14 )
   axis equal;
   xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
   ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
   title( 'Exponential Horn Profile', 'FontSize', 16 );

   case 7
   %%% Hyperbolic profile %%%
   a = sqrt(ai^2 + (power(z,2) * (ao^2-ai^2) / length^2));
   
   plot(z, a);    
   set(gca, "linewidth",2, "fontsize", 14 )
   axis equal;
   xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
   ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
   title( 'Hyperbolic Horn Profile', 'FontSize', 16 );

   case 8
   %%% POLYNOMIAL Profile %%%
   rho = 3;
   a=ai+(rho+1)*(ao-ai)*(1-((rho*z)/((rho+1)*length))).*power(z/length,rho);
   
   plot(z, a);    
   set(gca, "linewidth",2, "fontsize", 14 )
   axis equal;
   xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
   ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
   title( 'Polynomial Horn Profile', 'FontSize', 16 );
   
endswitch

switch(mode_converter_type)
  
  case 1                            % case 1=VARIABLE SLOT DEPTH MODE CONVERTER
    % Mode Converter depths for element j
    ajmc = a(1:NMC);                     % Index range for mode converter
    idx = 1:NMC;
    djmc = (sigma-((idx-1)./NMC).*(sigma-(0.25.*exp(1./(2.114.*(kc*ajmc).^1.134)))))*lambda_c;
    % Depth of remaining corrugations
    aj = a(NMC+1:end);
    idx = NMC+1:N+1;
    dj = ((lambda_c/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-NMC-1)/(N-NMC-1)).*((lambda_c/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lambda_o/4).*exp(1./(2.114.*(ko*ao).^1.134)));
    d = [djmc, dj];       % Combining the mode converter and horn depth values

    % Generate z,y coordinates as len and rad vector
    n = 0;
    len(1) = 0;
    len(2) = 0;
      for i = 1:N;
       rad(i+n) = a(i);
       rad(i+n+1) = a(i)+d(i);
       rad(i+n+2) = a(i)+d(i);
       rad(i+n+3) = a(i+1);
       rad(i+n+4) = a(i+1);
       len(i+n+2) = len(i+n)+delta*p;
       len(i+n+3) = len(i+n+2);
       len(i+n+4) = len(i+n+3)+(1-delta)*p;
       len(i+n+5) = len(i+n+4);
       n = n+3;
      endfor
    
    z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
    len = len(1:z_number); % Truncate z axis data points to equal rad vector length
    
    figure
    plot(len,rad);
    set(gca, "linewidth",2, "fontsize", 14 )
    xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
    ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
    title( 'Mode Converter and Corrugation Profile', 'FontSize', 16 );
    axis equal;   % Scale axis equally for aspect ratio 1:1

    case 2                             % case 2=RING LOADED SLOT MODE CONVERTER
      % Mode Converter depths for element j
      ajmc = a(1:NMC);                 % Index range for mode converter
      idx = 1:NMC;
      djmc = (lambda_c/4).*exp(1./(2.114.*(kc*ajmc).^1.134));
      % Width of bjth slot for mode converter
      bj = (0.1+(idx-1).*((delta2-0.1)./NMC)).*p;
      % Height of hjth slot for mode converter
      hj = (2/3).*djmc;
    
      % Depth of remaining corrugations
      aj = a(NMC+1:end);
      idx = NMC+1:N+1;
      dj = ((lambda_c/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-NMC-1)/(N-NMC-1)).*((lambda_c/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lambda_o/4).*exp(1./(2.114.*(ko*ao).^1.134)));
      d = [djmc, dj];       % Combining the mode converter and horn depth values
    
      % Generate z,y coordinates as len,rad vector
      n = 5;
      len = [0, 0, (-delta*p)+bj(1), (-delta*p)+bj(1), bj(1), bj(1)];
      rad = [a(1), a(1)+d(1)-hj(1), a(1)+d(1)-hj(1), a(1)+d(1), a(1)+d(1), a(2)];
        for i = 2:NMC;
          rad(i+n) = a(i);
          rad(i+n+1) = a(i)+d(i)-hj(i);
          rad(i+n+2) = a(i)+d(i)-hj(i);
          rad(i+n+3) = a(i)+d(i);
          rad(i+n+4) = a(i)+d(i);
          rad(i+n+5) = a(i+1);
          len(i+n) = i*p-p;
          len(i+n+1) = i*p-p;
          len(i+n+2) = len(i+n+1)-(delta*p)+bj(i);
          len(i+n+3) = len(i+n+1)-(delta*p)+bj(i);
          len(i+n+4) = len(i+n+3)+(delta*p)+bj(i);
          len(i+n+5) = len(i+n+3)+(delta*p)+bj(i);
          n = n+5;
        endfor
      % Add extra coordinate points before remaining corrugations
      len(NMC*(NMC+1)+1) = len(NMC*(NMC+1))+(1-delta)*p; 
      rad(NMC*(NMC+1)+1) = a(NMC+1);  

      figure
      plot(len,rad);
      set(gca, "linewidth",2, "fontsize", 14 )
      xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
      ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
      title( 'Ring Loaded Mode Converter Profile', 'FontSize', 16 );
      axis equal;   % Scale axis equally for aspect ratio 1:1

      n = n+NMC+1
      for i = NMC+1:N;
        rad(n) = a(i);
        rad(n+1) = a(i)+d(i);
        rad(n+2) = a(i)+d(i);
        rad(n+3) = a(i+1);
        rad(n+4) = a(i+1);
        len(n+1) = len(n);
        len(n+2) = len(n+1)+delta*p;
        len(n+3) = len(n+2);
        len(n+4) = len(n+3)+(1-delta)*p;
        n = n+4;
      endfor

      z_number = (NMC*2)+(N*4)+1; % Number of coordinate points for corrugated length of horn
   
      figure
      plot(len,rad);
      set(gca, "linewidth",2, "fontsize", 14 )
      xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
      ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
      title( 'Ring Loaded Mode Converter Horn Profile', 'FontSize', 16 );
      axis equal;   % Scale axis equally for aspect ratio 1:1

    case 3                 % case 3=VARIABLE PITCH TO WIDTH SLOT MODE CONVERTER
      % Mode Converter depths for element j
      ajmc = a(1:NMC);                     % First indexes for mode converter
      idx = 1:NMC;
      djmc = (sigma*(lambda_c/1.15)+((idx-1)./(NMC-1)).*(lambda_c/4-(sigma*lambda_c/1.15))).*exp(1./(2.114.*(kc*ajmc).^1.134))
      % Depth of remaining corrugations
      aj = a(NMC+1:end);
      idx = NMC+1:N+1;
      dj = ((lambda_c/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-NMC-1)/(N-NMC-1)).*((lambda_c/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lambda_o/4).*exp(1./(2.114.*(ko*ao).^1.134)))
      d = [djmc, dj];       % Combining the mode converter and horn depth values

      % Generate z,y coordinates as len,rad vector
      n = 0;
      len(1) = 0;
      len(2) = 0;
        for i = 1:NMC;
          rad(i+n) = a(i);
          rad(i+n+1) = a(i)+d(i);
          rad(i+n+2) = a(i)+d(i);
          rad(i+n+3) = a(i+1);
          rad(i+n+4) = a(i+1);
          len(i+n+2) = len(i+n)+(delta_min+(((i-1)./(NMC-1))*(delta-delta_min)))*p;
          len(i+n+3) = len(i+n+2);
          len(i+n+4) = len(i+n+3)+(1-(delta_min+(((i-1)./(NMC-1))*(delta-delta_min))))*p;
          len(i+n+5) = len(i+n+4);
          n = n+3;
        endfor

      figure
      plot(len(1:end-1),rad);
      set(gca, "linewidth",2, "fontsize", 14 )
      xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
      ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
      title( 'Variable Pitch To Width Mode Converter Profile', 'FontSize', 16 );
      axis equal;   % Scale axis equally for aspect ratio 1:1 
 
        for i = NMC+1:N;
          rad(i+n) = a(i);
          rad(i+n+1) = a(i)+d(i);
          rad(i+n+2) = a(i)+d(i);
          rad(i+n+3) = a(i+1);
          rad(i+n+4) = a(i+1);
          len(i+n+2) = len(i+n)+delta*p;
          len(i+n+3) = len(i+n+2);
          len(i+n+4) = len(i+n+3)+(1-delta)*p;
          len(i+n+5) = len(i+n+4);
          n = n+3;
        endfor
         
      z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
      len = len(1:z_number); % Truncate z axis data points to equal rad vector length

      figure
      plot(len,rad);
      set(gca, "linewidth",2, "fontsize", 14 )
      xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
      ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
      title( 'Variable Pitch To Width Mode Converter Horn Profile', 'FontSize', 16 );
      axis equal;   % Scale axis equally for aspect ratio 1:1
      
endswitch

% Add the rest of the geometry to create a closed path
% a_offset is the inner horn profile shifted up to give the horn a thickness
a_offset = a.+(lambda_c/2+2);
%figure;                    % Uncomment these three lines for debugging
%plot(z, a_offset);
%axis equal;

% Add vertical surface at horn aperture
len = [len, len(z_number)];
rad = [rad, a_offset(N)];
radmsh=rad;                 % radmesh to fix mesh lines to corrugations
%figure;                    % Uncomment these three lines for debugging
%plot(len, rad);
%axis equal;

% Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
outer_surface = fliplr(a_offset);
z_flip = fliplr(z);
extent = len(end);  % Fudge to make horn aperture planar for ring loaded slot MC
z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
% Add outer profile and circular waveguide to horn
len = [len, z_flip, -wgl, -wgl, 0];
rad = [rad, outer_surface, ai+(lambda_c/2+2), ai, ai];

figure
plot(len,rad);
set(gca, "linewidth",2, "fontsize", 14 )
xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
title( 'Complete Corrugated Horn Profile', 'FontSize', 16 );
axis equal;   % Scale axis equally for aspect ratio 1:1

% openEMS setup begins here
% EM related physical constants
physical_constants;

% frequency range of interest
f_start =  fmin*1e9;
f_stop  =  fmax*1e9;

% frequency to calculate fields
f0 = 12.46*1e9;

%% setup FDTD parameter & excitation function
FDTD = InitFDTD( 'NrTS', 50000, 'EndCriteria', 0.5e-3 );
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond(FDTD, BC);

%% setup CSXCAD geometry & mesh
max_res = c0/(f_stop)/unit/20; % cell size: lambda/20
CSX = InitCSX();               % Initialise CSX structure

% Calculate lambda/4 at lowest frequency to use as distance to nf2ff surfaces
lambda_max = c0/f_start/unit/4;

% Create fixed lines for the simulation box, structure and port
mesh.x = [(-a_offset(end)-(9*max_res)-lambda_max) -radmsh(1:4:end) 0 radmsh(1:4:end) (a_offset(end)+(9*max_res)+lambda_max)];
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.5); % Create a smooth mesh between specified fixed mesh lines
mesh.y = mesh.x;                                 % Same as x mesh
% Create fixed lines for the simulation box,port and given number of lines inside the horn
mesh.z = [-wgl-lambda_max-(9*max_res) -wgl-1 -wgl -wgl+10 0 len(1:2:z_number) length+2*lambda_max+(9*max_res)];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% create horn
% horn + waveguide, defined by a rotational polygon
CSX = AddMetal(CSX, 'Corrugated_Horn');
coords = [rad; len];
CSX = AddRotPoly(CSX,'Corrugated_Horn',10,'x',coords,'z');

% End cap to prevent the radiation coming out of the back of the horn
CSX = AddMetal(CSX, 'Cap');
CSX = AddCylinder(CSX,'Cap',10,[0 0 -wgl],[0 0 (-wgl-1)],a_offset(1));
% End of model geometry 

%% Apply the excitation %%
start=[-ai -ai -wgl];
stop =[ai ai -wgl+10];
[CSX, port] = AddCircWaveGuidePort( CSX, 0, 1, start, stop, ai*unit, 'TE11', 0, 1);

% Dump box for Electric field at Phi=0 (vertical cut)
CSX = AddDump(CSX,'Et_V_dump', 'SubSampling', '4,4,4');
start=[0 (-a_offset(end)-lambda_max) (-wgl-lambda_max)];
stop =[0 (a_offset(end)+lambda_max) (length+2*lambda_max)];
CSX = AddBox(CSX,'Et_V_dump',0,start,stop);

% Dump box for Electric field at Phi=90 (horizontal cut)
CSX = AddDump(CSX,'Et_H_dump', 'SubSampling', '4,4,4');
start=[(-a_offset(end)-lambda_max) 0 (-wgl-lambda_max)];
stop =[(a_offset(end)+lambda_max) 0 (length+2*lambda_max)];
CSX = AddBox(CSX,'Et_H_dump',0,start,stop);

% nf2ff calc
start = [mesh.x(9) mesh.y(9) mesh.z(9)];
stop  = [mesh.x(end-8) mesh.y(end-8) mesh.z(end-8)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions', [1 1 1 1 1 1], 'OptResolution', max_res*4);

% Prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'Corrugated_Horn.xml';
[status, message, messageid] = rmdir( Sim_Path, 's'); % Clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % Create empty simulation folder

% Write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

% Show the structure
if(RUN_CSXCAD == 1)
% CSXGeomPlot([Sim_Path '/' Sim_CSX], ['--export-polydata-vtk=tmp']);
CSXGeomPlot([Sim_Path '/' Sim_CSX], ['--export-STL=tmp']);
end

% Run openEMS
if(RUN_SIMULATION == 1)
%openEMS_opts = '--debug-PEC --no-simulation';   % Uncomment to visualise mesh in Paraview
%RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);
RunOpenEMS(Sim_Path, Sim_CSX, '--numThreads=4');
end

% Postprocessing & do the plots
freq = linspace(f_start,f_stop,201);
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;

% Plot reflection coefficient S11
figure
plot(freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2);
xlim([fmin fmax]);
ylim([-40 0]);
set(gca, "linewidth",2, "fontsize", 14)
grid on
title('Reflection Coefficient S_{11}', 'FontSize', 16);
xlabel('Frequency (GHz)','FontSize', 14);
ylabel('Reflection Coefficient |S_{11}| (dB)','FontSize', 14);
drawnow

% NFFF plots

% Calculate the far field at phi=0, 45 and at phi=90 degrees
thetaRange = (0:0.2:359) - 180;
disp('calculating far field at phi=[0 45 90] deg...');
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, [0 45 90]*pi/180);

Dlog=10*log10(nf2ff.Dmax);      % Calculate maximum Directivity in dB
G_a = 4*pi*A_app/(c0/f0)^2;     % Calculate theoretical gain for given aperture
e_a = nf2ff.Dmax/G_a;           % Calculate Efficiency

% Display some antenna parameters from above calculations
disp(['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp(['directivity: Dmax = ' num2str(Dlog) ' dBi']);
disp(['aperture efficiency: e_a = ' num2str(e_a*100) '%']);

% Directivity
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 2 3]);
ylim([-30 25]);
xlim([-180 180]);
grid on
set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40)
title('Farfield Directivity @ 12.46GHz','FontSize', 16);
xlabel('Theta (degrees)','FontSize', 14);
ylabel('Directivity (dBi)','FontSize', 14);
drawnow

% Plot Ludwig3 cross polar
plotFFcocx(nf2ff,'xaxis','theta','param',[2]);
ylim([-30 25]);
xlim([-180 180]);
grid on
set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40)
title('Farfield Directivity with Ludwig3 XPOL @ 12.46GHz','FontSize', 16);
xlabel('Theta (degrees)','FontSize', 14);
ylabel('Directivity (dBi)','FontSize', 14);
drawnow

% Polar plot
figure
leg=[];   %legend
polarFF(nf2ff,'xaxis','theta','param',[1 2 3],'logscale',[-30 35], 'xtics', 12);
title('Farfield Directivity @ 12.46GHz','FontSize', 16);
xlabel('Theta (degrees)','FontSize', 14);
ylabel('Directivity (dBi)','FontSize', 14);
drawnow

%% Calculate 3D pattern
%phiRange = sort(unique([-180:5:-100 -100:2.5:-50 -50:1:50 50:2.5:100 100:5:180]));
%thetaRange = sort(unique([0:1:50 50:2:100 100:5:180]));
phiRange = sort(unique([-180:1:-100 -100:1:-50 -50:1:50 50:1:100 100:1:180]));
thetaRange = sort(unique([0:1:50 50:1:100 100:1:180]));

disp('calculating 3D far field...');
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, phiRange*pi/180, 'Verbose',2,'Outfile','nf2ff_3D.h5');

figure
colormap jet;
plotFF3D(nf2ff, 'logscale', -40);        % plot 3D far field in dB

% Save far field in VTK to plot in ParaView
E_far_normalized = nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:));
DumpFF2VTK([Sim_Path '/Farfield.vtk'],E_far_normalized,thetaRange,phiRange,'scale', 0.008, 'logscale', -30, 'maxgain', Dlog);

%%% END OF SCRIPT %%%
