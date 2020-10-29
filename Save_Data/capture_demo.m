%  NOTE:
%      (C) Copyright 2016 Texas Instruments, Inc.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions
%  are met:
%
%    Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
%
%    Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the
%    distribution.
%
%    Neither the name of Texas Instruments Incorporated nor the names of
%    its contributors may be used to endorse or promote products derived
%    from this software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
%  capture_demo.m
%
%  Description:
%  This file is a sample MATLAB script to process the data captured from
%  Capture Demo.
%  It reads data from a .dat file which is saved from CCS or 
%  a .csv file which is saved from LVDS + TSW setup. 
%  It then generates first FFT graph of the captured data.
%  Please refer to the capture demo script for the parameters used in the
%  demo. Such as number of samples, number of chirps etc.
%  
%  Param[in]: fileName   
%             The name of input file to be processed. 
%             It is either a .dat file (CCS memory save data file) or 
%             a .csv file (LVDS captured data file). 
%  Param[in]: dataFormat 
%             Apply when input file is .csv file. Select one of 3 data formats 
%             1: CBUFF_DataFmt_ADC_DATA  (Only ADC Data)
%             2: CBUFF_DataFmt_CP_ADC    (Chirp Parameters + ADC Data)
%             3: CBUFF_DataFmt_CP_ADC_CQ (Chirp Parameters + ADC Data + Chirp Quality)
%             Don't care when input file is .dat file.
% 

function capture_demo(fileName, dataFormat)
close all;

% Global variables
% Please align with the actual demo for these variables.

% Number of ADC samples per chirp
K = 256;
% Number of chirps
Nchirp = 128;
% Slope const (MHz/usec)
freqSlope = 40 * 1e12;
% Velosity of light
C = 3e8;
% ADC sampling rate
sampleRate = 8000 * 1e3;
% ADC start time (usec)
adcStartTime = 5 * 1e-6;
% Start frequency (GHz)
startFreq = 77 * 1e9;
% ramp End time
rampEndTime = 40 * 1e-6;

CPLen = 1;
magicWordsLen = 4;
CQLen = 167;

% Create a array to hold pure ADC samples of one frame
ADCSamplesFrame = zeros(K*Nchirp, 2);

datFile = strfind(fileName,'.dat');
csvFile = strfind(fileName,'.csv');

if (numel(datFile) > 0) && (numel(csvFile) == 0)
    fprintf('Input is a .dat file!\n');
    
    % Read in CCS memory file and convert to real number
    fid = fopen(fileName, 'r');
    rawData = textscan(fid, '%s', 'Delimiter', '', 'headerLines',1);
    data = cell2mat(cat(2, rawData{:}));
    % Convert data to decimal, and reshape in two columns for I and Q
    dataDec = hex2dec(data);
    dataDec = dataDec - ( dataDec >=2.^15).* 2.^16;
    rad_in = reshape(dataDec, [2, length(dataDec)/2]).';
    fclose(fid);
    
    % pure ADC samples of one frame 
    ADCSamplesFrame(1:K*Nchirp,:) = rad_in(1:K*Nchirp,:);

elseif (numel(datFile) == 0) && (numel(csvFile) > 0)
    fprintf('Input is a .csv file!\n');
    
    % Read in LVDS csv file
    data = csvread(fileName);
    % Reshape in two columns for I and Q
    [row col] = size(data);
    rad_in = reshape(data.', [2, row*col/2]).';

    if dataFormat == 1 %CBUFF_DataFmt_ADC_DATA (Only ADC Data is to sent out)
        capturedChirpLen = magicWordsLen + K;
        adcSamplesStart  = magicWordsLen + 1;
        adcSamplesEnd    = magicWordsLen + K;
    elseif dataFormat == 2 %CBUFF_DataFmt_CP_ADC (Chirp Parameters + ADC Data)
        capturedChirpLen = CPLen + magicWordsLen + K;
        adcSamplesStart  = CPLen + magicWordsLen + 1;
        adcSamplesEnd    = CPLen + magicWordsLen + K;
    elseif dataFormat == 3 %CBUFF_DataFmt_CP_ADC_CQ (Chirp Parameters + ADC Data + Chirp Quality)
        capturedChirpLen = CPLen + magicWordsLen + K + CQLen;
        adcSamplesStart  = CPLen + magicWordsLen + 1;
        adcSamplesEnd    = CPLen + magicWordsLen + K;
    else
        fprintf('Wrong data format!\n');
        fprintf('Support 3 data formats\n');
        fprintf('1: CBUFF_DataFmt_ADC_DATA  (Only ADC Data)\n');
        fprintf('2: CBUFF_DataFmt_CP_ADC    (Chirp Parameters + ADC Data)\n');
        fprintf('3: CBUFF_DataFmt_CP_ADC_CQ (Chirp Parameters + ADC Data + Chirp Quality)\n');
        return;
    end

    % pure ADC samples of one frame 
    for chirpIdx = 1:Nchirp
        ADCSamplesFrame(K*(chirpIdx-1)+1:K*chirpIdx,:) = rad_in(capturedChirpLen*(chirpIdx-1)+adcSamplesStart:capturedChirpLen*(chirpIdx-1)+adcSamplesEnd,:);  
    end
else
    fprintf('Wrong file!\n');
    fprintf('Support .dat file (CCS memory data file) or .csv file (LVDS data file)!\n');
    return;
end

% Plot I and Q
figure(1);
plot(ADCSamplesFrame(:,1));grid on;
hold on;
plot(ADCSamplesFrame(:,2),':r');grid on;
xlabel('ADC sample index');
ylabel('Sample value');
title('ADC raw data');

% Change to complex data and Calculate average
% To handle static object only, use average as an example filter.
frame = ADCSamplesFrame(1:K*Nchirp,1) + i*ADCSamplesFrame(1:K*Nchirp,2);
chAve = zeros(K,1);
avgChirp = Nchirp;
for nn = 1:avgChirp
    chAve = chAve + frame((nn-1)*K+1:(nn-1)*K+K);
end
chAve = chAve/avgChirp;

% Calculate range
rangeResolutionsInMeters = C * sampleRate / (2 * freqSlope * K);
rangeIdxToMeters = C * sampleRate / (2 * freqSlope * K);
rangeBin = linspace(0, K * rangeIdxToMeters, K);

% Do FFT on chirp average and plot results
fftResult = fft(chAve);

% Do FFT on chirp average and plot results
figure(2);
plot(rangeBin, abs(fftResult), '-m'); grid on;
hold on;
xlabel('range (m)');
title('Range FFT');



