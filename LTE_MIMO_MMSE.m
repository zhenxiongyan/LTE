
clc;
close;
clear all;
warning off;

%%  system parameter

FFTsize = 2048;
CPsize = 144;
Nt = 2;% # of tx antennas
Nr = 2;% # of rx antennas
modOrderVec = [2 4 6];% M-QAM
snrVec = [0:2:30];
nframes = 100;% number of ofdm frames
RENum = 1200;
SampleRate = 30.72e6;
GuardBandSize = (FFTsize - RENum)/2;

Len = length(snrVec);
SimLen = length(modOrderVec);
BER = zeros(SimLen,Len);


%% init the modules

% Create an OFDM modulator and demodulator
ofdmMod = comm.OFDMModulator('FFTLength',FFTsize,...
    'PilotInputPort',true,...
    'NumGuardBandCarriers',[GuardBandSize;GuardBandSize-1],...
    'InsertDCNull',true,...
    'PilotInputPort',false, ...
    'NumTransmitAntennas',Nt,...
    'CyclicPrefixLength',CPsize);
ofdmDemod = comm.OFDMDemodulator(ofdmMod);
ofdmDemod.NumReceiveAntennas = 2;

%%
lteChannel = comm.LTEMIMOChannel(...
    'Profile',              'ETU 70HZ',...
    'AntennaConfiguration', '2x2',...
    'CorrelationLevel',     'Low',...
    'AntennaSelection',     'Off',...
    'RandomStream',         'mt19937ar with seed',...
    'Seed',                 1000,...
    'PathGainsOutputPort',  true,...
    'SampleRate', SampleRate);

%%
% Determine the dimensions of the OFDM modulator
ofdmModDim = info(ofdmMod);
numData = ofdmModDim.DataInputSize(1);   % Number of data subcarriers
numSym = ofdmModDim.DataInputSize(2);    % Number of OFDM symbols
numTxAnt = ofdmModDim.DataInputSize(3);  % Number of transmit antennas


 for index = 1:SimLen
    modOrder = modOrderVec(index);
    
    %% Generate data symbols to fill 100 OFDM frames.
    data = randi([0 1],nframes*numData*modOrder*2,1);
    
    %% M QAM
    hmod = comm.RectangularQAMModulator('ModulationOrder',2^modOrder,'BitInput',true,'NormalizationMethod','Average power');
    Moddata = step(hmod,data(:));
    %scatterplot(Moddata)
    dataMod = reshape(Moddata,nframes*numData,numSym,numTxAnt);
    %
    %% OFDM Modulation
    data_all = zeros(nframes*(numData-400),numSym,numTxAnt);
    
    for sim_index = 1:Len
        snr = snrVec(sim_index);
        
        hAWGN = comm.AWGNChannel('NoiseMethod',...
            'Signal to noise ratio (SNR)','SNR',snr);
        
        for k = 1:nframes
            
            % Find row indices for kth OFDM frame
            indData = (k-1)*numData+1:k*numData;
            indData_pilotRem = (k-1)*(numData-400)+1:k*(numData-400);
            
            % Generate random OFDM pilot index
            pilotIdx = [(2:6:RENum).' (5:6:RENum).'];
            
            % pilot and data 
            dataSym = dataMod(indData,:,:);
            % the first stream
            dataSym(pilotIdx(:,1),:,1) = 1;
            dataSym(pilotIdx(:,2),:,1) = 0;
            % the second stream
            dataSym(pilotIdx(:,1),:,2) = 0;
            dataSym(pilotIdx(:,2),:,2) = 1;
            
            % remove the pilot
            pilotVec = reshape(pilotIdx,[],1);
            DataTx = dataSym;
            DataTx(pilotVec,:,:) = [];
            data_tx(indData_pilotRem,:,:) = DataTx;
            
            % Modulate  symbols using OFDM
            dataOFDM = step(ofdmMod,dataSym)*sqrt(FFTsize);
            
            %% LTE channel
            [LTEChanOut,LTEPathGain] = step(lteChannel,dataOFDM);
            %         LTEChanOut = dataOFDM;
            %% AWGN Channel
            AWGNOut = step(hAWGN,LTEChanOut);
            
            %%
            % Demodulate OFDM data
            
            [receivedOFDMData] = step(ofdmDemod,AWGNOut)/sqrt(FFTsize);
            y = squeeze(receivedOFDMData);
            %%
            %%%   channel estimate
            H_est = zeros(numData,2,2);
            for tx = 1:2
                for rx = 1:2
                    % LS 信道估计
                    H_LS = squeeze(receivedOFDMData(pilotIdx(:,tx),1,rx));
                    % 插值
                    p_idx = pilotIdx(:,tx);
                    inter_p_idx = 1:RENum;
                    % 内插
                    H_est(:,tx,rx) = interp1(p_idx,H_LS,inter_p_idx,'linear','extrap');
                end
            end
            ChanG = H_est;
            
            %%
            %MMSE detection
            numData = size(receivedOFDMData, 1);
            Variance =  1/10^(hAWGN.SNR);
            r=zeros(numData,2);
            for i = 1:numData
                h = squeeze(ChanG(i,:,:)).';
                w = pinv(h'*h+Variance*eye(2))*h';
                r(i,:) = w*y(i,:).';
            end
            % remove the pilot
            r_rx = r;
            r_rx(pilotVec,:,:) = [];
            data_all(indData_pilotRem,:,:)= r_rx;
            
        end
        %% 16 demodulation
        rxSig = data_all(:);
         %%scatterplot(rxSig)
        hDemod = comm.RectangularQAMDemodulator('ModulationOrder',2^modOrder,'BitOutput',true,'NormalizationMethod','Average power');
        dataout = step(hDemod,rxSig);
        
        
        txSig = data_tx(:);
        % scatterplot(rxSig)
        hDemod = comm.RectangularQAMDemodulator('ModulationOrder',2^modOrder,'BitOutput',true,'NormalizationMethod','Average power');
        dataSource = step(hDemod,txSig);
        %
        %%
        
        [number,ber]= biterr(dataSource,dataout);
        
        BER(index,sim_index) = ber;
        
    end;
    
end;

%%

figure(1);
semilogy(snrVec,BER(1,:),'k-s','linewidth',1);
hold on;
semilogy(snrVec,BER(2,:),'k-o','linewidth',1);
hold on
semilogy(snrVec,BER(3,:),'k-o','linewidth',1);
xlabel('SNR(dB)');
ylabel('BER');
grid on;
legend('QPSK','16QAM','64QAM');