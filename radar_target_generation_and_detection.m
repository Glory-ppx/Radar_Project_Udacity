clear all
clc
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
MaxRange = 200;
% Range Resolution = 1 m
rangeResolution = 1;
% Max Velocity = 100 m/s
MaxVelocity = 100;
% Velocity Resolution = 1m/s
velResolution = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3e8;
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
target_initial_position = 100;
target_velocity = 50;
 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
Bandwidth = c/(2*velResolution);
Tchirp = 5.5*2*MaxRange/c;
Slope = Bandwidth/Tchirp;

%Operating carrier frequency of Radar 
carrier_frequence= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Ndoppler=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nrange=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Ndoppler*Tchirp,Nrange*Ndoppler); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=cos(2*pi*(carrier_frequence*t+Slope*t.^2/2)); %transmitted signal
r_t = target_initial_position + target_velocity*t; %target position
td = 2*r_t/c; % time delay
Rx=cos(2*pi*(carrier_frequence*(t-td)+Slope*(t-td).^2/2)); %received signal
Mix = Tx.*Rx; %beat/mixed signal


%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix,[Nrange,Ndoppler]) ;

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
fft_beat = fft(Mix, Nrange,1 )/Nrange;

 % *%TODO* :
% Take the absolute value of FFT output
fft_beat_abs = abs(fft_beat);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
fft_beat_abs_first_half = fft_beat_abs(1:Nrange/2);


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

Mix=reshape(Mix,[Nrange,Ndoppler]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nrange,Ndoppler);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nrange/2,1:Ndoppler);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Ndoppler);
range_axis = linspace(-200,200,Nrange/2)*((Nrange/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 8; % Number of training cells range
Td = 8; % Number of training cells doppler

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 2; % Number of guard cells range
Gd = 2; % Number of guard cells doppler
% *%TODO* :
% offset the threshold by SNR value in dB
offset=1.7;
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(length(RDM)-(2*Gr+2*Tr+1),Ndoppler-(2*Gd+2*Td+1));


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

   for i = 1:(length(RDM)-(2*Gr+2*Tr+1))
    for j = 1:(Ndoppler-(2*Gd+2*Td+1))
        % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
        % CFAR
        noise_level(i,j) = get_training_mean(RDM(i:i-1+2*Gr+2*Tr+1,j:j-1+2*Gd+2*Td+1),Tr,Td,Gr,Gd);
        threshold = noise_level(i,j)*offset;
        if (RDM(i+Tr+Gr, j+Td+Gd) < threshold)
            RDM(i+Tr+Gr, j+Td+Gd) = 0;
        else
            RDM(i+Tr+Gr, j+Td+Gd) = 1;
        end    
    end
   end

% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
RDM(1:Tr+Gr,:) = 0;
RDM(end-Tr-Gr:end,:)= 0;
RDM(Tr+Gr+1:end-(Tr+Gr),1:Td+Gd) = 0;
RDM(Tr+Gr+1:end-(Tr+Gr),end-Td-Gd:end)= 0;



% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RDM);
colorbar;


 
 