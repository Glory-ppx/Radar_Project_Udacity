clear all
clc
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
Frequence_of_operation= 77e9;
% Max Range = 200m
MaxRange = 200;
% Range Resolution = 1 m
RangeResolution = 1;
% Max Velocity = 100 m/s
MaxVelocity = 100;
% Velocity Resolution = 1m/s
VelResolution = 1;
%speed of light = 3e8
c = 3e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
Bandwidth = c/(2*VelResolution);
Tchirp = 5.5*2*MaxRange/c;
Slope = Bandwidth/Tchirp;

%Operating carrier frequency of Radar 
             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
N_doppler=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
N_range=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
time=linspace(0,N_doppler*Tchirp,N_range*N_doppler); %total amount of time samples

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
target_initial_position = 100;
target_velocity = 50;
 
%% Generating the signals
% *%TODO* :
%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Transmitted_signal=cos(2*pi*(Frequence_of_operation*time+Slope*time.^2/2)); %transmitted signal
Range_t = target_initial_position + target_velocity*time; %target position
TimeDelay = 2*Range_t/c; % time delay
Received_signal=cos(2*pi*(Frequence_of_operation*(time-TimeDelay)+Slope*(time-TimeDelay).^2/2)); %received signal
Mixed_signal = Transmitted_signal.*Received_signal; %beat/mixed signal



%% RANGE MEASUREMENT

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mixed_signal = reshape(Mixed_signal,[N_range,N_doppler]) ;

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
fft_beat = fft(Mixed_signal, N_range,1 )/N_range;

 % *%TODO* :
% Take the absolute value of FFT output
fft_beat_abs = abs(fft_beat);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
fft_beat_abs_first_half = fft_beat_abs(1:N_range/2,:);


%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
 plot(fft_beat_abs_first_half)

 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.


% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mixed_signal,N_range,N_doppler);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:N_range/2,1:N_doppler);
sig_fft2 = fftshift (sig_fft2);
RadarDopplerMap = abs(sig_fft2);
RadarDopplerMap = 10*log10(RadarDopplerMap) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,N_doppler);
range_axis = linspace(-200,200,N_range/2)*((N_range/2)/400);
figure,surf(doppler_axis,range_axis,RadarDopplerMap);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
N_Training_cell_Range = 8; % Number of training cells range
N_Training_cell_Doppler = 8; % Number of training cells doppler

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
N_Guard_cell_Range = 2; % Number of guard cells range
N_Guard_cell_Doppler = 2; % Number of guard cells doppler
% *%TODO* :
% offset the threshold by SNR value in dB
Offset=1.7;



% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR

   for i = 1:(length(RadarDopplerMap)-(2*N_Guard_cell_Range+2*N_Training_cell_Range+1))
    for j = 1:(N_doppler-(2*N_Guard_cell_Doppler+2*N_Training_cell_Doppler+1))
        % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
        % CFAR
        Noise_Level = get_training_mean(RadarDopplerMap(i:i-1+2*N_Guard_cell_Range+2*N_Training_cell_Range+1,j:j-1+2*N_Guard_cell_Doppler+2*N_Training_cell_Doppler+1),N_Training_cell_Range,N_Training_cell_Doppler,N_Guard_cell_Range,N_Guard_cell_Doppler);
        threshold = Noise_Level*Offset;
        if (RadarDopplerMap(i+N_Training_cell_Range+N_Guard_cell_Range, j+N_Training_cell_Doppler+N_Guard_cell_Doppler) < threshold)
            RadarDopplerMap(i+N_Training_cell_Range+N_Guard_cell_Range, j+N_Training_cell_Doppler+N_Guard_cell_Doppler) = 0;
        else
            RadarDopplerMap(i+N_Training_cell_Range+N_Guard_cell_Range, j+N_Training_cell_Doppler+N_Guard_cell_Doppler) = 1;
        end    
    end
   end

% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
RadarDopplerMap(1:N_Training_cell_Range+N_Guard_cell_Range,:) = 0;
RadarDopplerMap(end-N_Training_cell_Range-N_Guard_cell_Range:end,:)= 0;
RadarDopplerMap(N_Training_cell_Range+N_Guard_cell_Range+1:end-(N_Training_cell_Range+N_Guard_cell_Range),1:N_Training_cell_Doppler+N_Guard_cell_Doppler) = 0;
RadarDopplerMap(N_Training_cell_Range+N_Guard_cell_Range+1:end-(N_Training_cell_Range+N_Guard_cell_Range),end-N_Training_cell_Doppler-N_Guard_cell_Doppler:end)= 0;



% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RadarDopplerMap);


 
 