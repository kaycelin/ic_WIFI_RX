%% WiFi Receiver
% History
%% 
% * 2024-10-28, Draft, scripts from Example_WIFI_TRx_240912_v4_draft.mlx
% * 
% Set Receiver Line-up

close all
clear all
clc
clear rxAntSig lnaSig % clear signal

Rohm    = 1 % set resistor in ohm
RbwKHz  = 10 % set resolution bw for aclr plot

format shortG

simulation_trx = 'rx'
switch simulation_trx
    case {'rx','trx'} % set simulation case as 'rx'
        is_sim_AdjacentChannelRejection     = 1 % simulate adjacent channel rejection
        is_set_rxAntWantedPower             = 2; % set by 1: Sensitiviy, 2: AdjacentChannelRejection, 3: FreeSpaceDistance, 4: MaximumInputLevel, otherwise: Sensitiviy+isAntWantedPowerSet
        is_sweep_AdjChIntferPower           = 0 % sweep power of adj. ch. interferer
        isAdjChIntferSource_NR              = 0 % 0. Based on WiFI ACR specifications, 1. -21dBm of NR max. input level, otherwise: adjInterferer_dBm_set

        isTRSwitch      = 1; % apply TRSwitch
        isLNA_Bypass    = 0; % set LNA mode
        isReceiverNoise = 1; % apply receiver noise 

        is_sweep_lnaNF  = 0 % sweep LNA NF
        is_sweep_lnaTOI = 1 % sweep LNA TOI
        isLNA_Flatness  = 1 % apply LNA fltaness
        isPropagationChannelModel = 0 % apply propoagation channel model
end

% Load Transmitter Signal

if 1 % import - tx signal
    load('WIFI_VHT80_MCS9_23.80dBm_NumPkts4_320MHz.mat') % load waveform
    % load('WIFI_EHT320_MCS13_15.80dBm_NumPkts4_1280MHz')
    wvParams;
    aclrVal_wv = wvParams.PlotACLR; % plot
    pltComm = wvParams.PlotComm; % plot
    dlSlots = wvParams.usefulWvfm;
    fsMHz = wvParams.SampleRateMHz;
    signalFormat = wvParams.SignalFormat;
    MCS_index = wvParams.MCS;
    wv = wvParams.finalSignal;
    configParamsWG = wvParams.configWG;
    cfg = wvParams.cfgFormat;
    cbw = wvParams.currentBw;
else
    wv = wvParams.finalSignal;
end

wv_org = wv; % store wv
fs = fsMHz*1e6;
transmit_dBm = powerDbm(wv,'rms',[],'dBm',dlSlots,Rohm); % check power

% get sensitiviy spec.
getSensLvlVsCBW = "Sensitivity_"+configParamsWG.bw+"_dBm";
[SENS_tab,SENS_Requ_dBm] = wlan_spec('sensitivity',signalFormat,{'MCS',getSensLvlVsCBW},{MCS_index,[]});

% plot
pltComm.sa(pltComm, wv(dlSlots), {091801,'wv(dlSlots)'});
pltComm.sa(pltComm, wv(:), {091801,'wv(:)'});
% Receiver - <https://www.mathworks.com/help/wlan/propagation-channel-models.html?s_tid=CRUX_lftnav *Propagation Channel Models*>
% This example uses a TGax non-line-of-sight (NLOS) indoor channel model with 
% delay profile Model-B. Model-B is considered NLOS when the distance between 
% transmitter and receiver is greater than or equal to 5 meters. For more information 
% about the TGax channel model, see <https://www.mathworks.com/help/wlan/ref/wlantgaxchannel-system-object.html 
% |wlanTGaxChannel|>.        

switch simulation_trx
    case {'rx','receiver','both','trx'}

        TransmitReceiveDistance = 10 % set init. tx-to-rx distance
        PassDecod_ch = 0;
        PassDecod_rx = nan;
        SIMULATING_1 = 1;

        distance_vec = [];
        power_vec = [];

        while SIMULATING_1 % STARTING...
            pause(1e-3)
            wv = wv_org;
            % isDecodedPass_ch = 0;

            if isPropagationChannelModel

                TransmitReceiveDistance;
                if logical(TransmitReceiveDistance < 2)
                    error("TransmitReceiveDistance < 2")
                end

                if isnan(PassDecod_rx) % Configure TGax channel
                    clear tgaxChannel

                    % Erceg, V., L. Schumacher, P. Kyritsi, et al. TGn Channel Models. Version 4. IEEE 802.11-03/940r4, May 2004.
                    % Model-A: For a office environmnet, NLOS conditions.
                    % Model-B: For a large open space and office environment, NLOS conditions.
                    % Model-C: For a large open space (indoor and outdoor), NLOS conditions.
                    % Model-D: Same as Model-C, LOS conditions.
                    % Model-E: For a large open space (indoor and outdoor).
                    % Doppler spread fd = vo/lamda
                    % vo is the environmental speed and is proposed equal to 1.2km.hr
                    % and lamda is the wavelenght defined by lamda=c/fc.
                    tgaxChannel = wlanTGaxChannel; % set propagation channel model
                    tgaxChannel.DelayProfile            = 'Model-F';
                    tgaxChannel.TransmitReceiveDistance = TransmitReceiveDistance;
                    tgaxChannel.NumTransmitAntennas     = configParamsWG.NumTransmitAntennas;
                    tgaxChannel.NumReceiveAntennas      = configParamsWG.NumTransmitAntennas;
                    tgaxChannel.ChannelBandwidth        = cfg.ChannelBandwidth;
                    tgaxChannel.LargeScaleFadingEffect  = 'Pathloss and shadowing';
                    tgaxChannel.SampleRate              = fsMHz * 1e6;

                    if 1 % set CarrierFrequency
                        switch configParamsWG.bw
                            case {"20MHz","40MHz"}
                                tgaxChannel.CarrierFrequency = 2.484e9 - wvParams.currentBw;
                            case {"80MHz","160MHz"}
                                tgaxChannel.CarrierFrequency = 5.25e9 - wvParams.currentBw;
                            otherwise
                                tgaxChannel.CarrierFrequency = 7.125e9 - wvParams.currentBw;
                        end
                    else
                        tgaxChannel.CarrierFrequency = 7.125e9 - wvParams.currentBw;
                    end
                else
                    TransmitReceiveDistance = TransmitReceiveDistance - 1; % decrease distance
                    release(tgaxChannel)
                    tgaxChannel.TransmitReceiveDistance = TransmitReceiveDistance;
                end

                if 0
                    % add trailing zeros to allow for channel delay
                    wv_zeroPad = [wv; zeros(50,cfg.NumTransmitAntennas)];
                else
                    wv_zeroPad = wv;
                end

                reset(tgaxChannel); % reset channel for different realization
                wv_propChanModl = tgaxChannel(wv_zeroPad); % pass through fading indoor TGax channel
                wv_propChanModl_dBm = powerDbm(wv_propChanModl,'rms',[],'dBm',dlSlots,Rohm); % check power
                if isnan(PassDecod_rx)
                    sprintf('if simulation_trx is rx, continue to tune the TransmitReceiveDistance to meet the Sensitivity power requeirement');
                    if SENS_Requ_dBm > wv_propChanModl_dBm + 5  % margin 5dB
                        TransmitReceiveDistance = TransmitReceiveDistance - 1; % decrease distance
                        continue
                    elseif SENS_Requ_dBm < wv_propChanModl_dBm -1 % margin 1dB
                        TransmitReceiveDistance = TransmitReceiveDistance + 1; % increase distance
                        continue
                    else
                        TransmitReceiveDistance;
                    end
                end

                try % demod
                    wvParams.finalSignal = wv_propChanModl;
                    demod_propChanModl = wlanDemodulation(wv_propChanModl, wvParams, 'rx')
                    PassDecod_ch = mean(demod_propChanModl.rxPER_percent, 'all') < 10;
                catch
                    PassDecod_ch = 0;
                end

                if ~PassDecod_ch
                    TransmitReceiveDistance = TransmitReceiveDistance - 1*isnan(PassDecod_rx); % decrease distance
                    close all
                    continue
                end

                if 0 % debug - power
                    wv_dBm = powerDbm(wv,'rms',[],'dBm',[],Rohm);
                    wvPad_dBm = powerDbm(wv_zeroPad,'rms',[],'dBm',[],Rohm);
                    wv_propChanModl_dBm;
                end

                % export
                wv = wv_propChanModl;
            end

% Receiver - Assign Sensitivity Power

            % set test case: adjacent channel rejection
            is_sim_AdjacentChannelRejection;

            % set receiver input signal
            if isPropagationChannelModel
                rxAntSig = wv_propChanModl;
            else
                % set ant. input power
                switch is_set_rxAntWantedPower
                    case {1,'Sensitivity'}
                        rxAntWantedPower_dBm_set = SENS_Requ_dBm;
                    case {2,'AdjacentChannelRejection'}
                        rxAntWantedPower_dBm_set = SENS_Requ_dBm + 3*is_sim_AdjacentChannelRejection + -0;
                    case {3,'FreqSpaceDistance'}
                        trxDistanceMeter = 0.1
                        freqSpaceLossDb = 20*log10(trxDistanceMeter)+20*log10(2.4e9)+20*log10(4*pi/3e8)
                        rxAntWantedPower_dBm_set = transmit_dBm - freqSpaceLossDb;
                    case {4,'MaximumInputLevel'}
                        rxAntWantedPower_dBm_set = -30;
                    otherwise
                        sprintf('Sensitivity + isAntWantedPowerSet')
                        rxAntWantedPower_dBm_set = SENS_Requ_dBm + is_set_rxAntWantedPower;

                end
                rxAntSig = powerDbm(wv,'set',rxAntWantedPower_dBm_set,'dBm',dlSlots);
            end

            if 0 % demod
                wvParams.finalSignal = rxAntSig;
                demod_rxant = wlanDemodulation(rxAntSig, wvParams, 'rx');
            end

            % export
            rxAntSig_org = rxAntSig; % store
            rxAnt_dBm = powerDbm(rxAntSig,'rms',[],'dBm',dlSlots);

% Receiver - Adjacent Channel Interferer

            if 1 % Initialized
                IterationNo = 1
                SIMULATING_2 = 1; % simulating...
                lna_iip3Step_dB = 1;
                lna_ip1dBStep_dB = 1;
                lna_nfStep_dB = 0.5;
                RESETPARMS = 1; % 2024-05-15, Reset parameter simulation

                PER_iter = [];
                NFdB_Cal_iter = [];
                NFdB_Set_iter = [];
                IIP3_dBm_iter = [];
                IIP3_L_iter = [];
                IIP3_H_iter = [];
                ACR_dB_iter = [];
                AdjInterferer_dBm_iter = [];
                EVM_RMS_iter = [];
                EVM_Subcarrier_iter = [];
            end

            if is_sim_AdjacentChannelRejection && ~exist('wv_adjChInf', 'var')
                ACR_tab = wlan_spec("acr", signalFormat);
                [ACR_tab, ACR_Requ_dB] = wlan_spec("acr", signalFormat, {'MCS','AdjacentChannelRejection_dB'}, {MCS_index,[]})
                if isnan(ACR_Requ_dB)
                    error('check the input format')
                end

                dispTile_adjCh = 'ADJACENT CHANNEL REJECTION';
                dispLegd_adjChSource = 'ADJ. Interferer';
                % get wv parameters - interferer wv
                configParamsWG_acr = configParamsWG;
                configParamsWG_acr.rng_seed = 'shuffle';
                configParamsWG_acr.bbFilter_en = 1;
                configParamsWG_acr.bbFilterSbAtten = 150;
                configParamsWG_acr.isACR_interferer = 1;
                wvParams_adjChInf = waveformGenWLAN.runWG(configParamsWG_acr); % create adj. channel interferer waveform
                wv_adjChInf = wvParams_adjChInf.finalSignal;
                wv_adjChInf = dfe_nco(wv_adjChInf,cbw,fs); % apply frequency shift to adjacent channel

                % set markers
                ftone1 = [];
                ftone2 = [];

                if 0
                    cond.markers = [10e6 20e6];
                    plot_comm(wv_adjChInf(:), fs, 'spectra', pltParms, {041502, [dispLegd_adjChSource]}, [ftone1 ftone2]); hold on
                    pltComm.sa(pltComm, wv_adjChInf(:), {091801, dispLegd_adjChSource},[],[],cond);

                    wv_2tones_dBm = powerDbm(wv_2tones,'rms',[],'dBm',[],Rohm)
                    plot_comm_sa(wv_adjChInf,fs,10e3)

                    [wv_inf_2, noise0] = rfSim_noise(wv_adjChInf, 'NoisePowerDbmHz', -150, fs);
                    plot_comm_sa(wv_inf_2,fs,10e3)
                end
            end

            isLNA_IIP3_debug = 0;
            if isLNA_IIP3_debug
                SIMULATING_2 = 0;
                rxAntWantedPower_dBm_set = -20
            end

            is_sweep_params = any([is_sweep_lnaNF,is_sweep_lnaTOI,is_sweep_AdjChIntferPower]);
            while SIMULATING_2

                if ~PassDecod_ch*isPropagationChannelModel
                    break % SIMULATING_2
                else
                    rxAntSig = rxAntSig_org;

                    if is_sim_AdjacentChannelRejection
                        % set ant. input power - adj. ch. interferer
                        switch isAdjChIntferSource_NR
                            case {0,'WiFi Interferer','wifi'}
                                adjChInf_add_dB = 32 + (IterationNo-1)*is_sweep_AdjChIntferPower;
                                adjChInf_dBm_set = rxAnt_dBm + ACR_Requ_dB + adjChInf_add_dB;
                            case {1,'NR Interferer','nr'}
                                nrMRInterferer_dBm = 24;
                                nrMRMinCouplingLoss_dB = 45;
                                adjChInf_dBm_set = nrMRInterferer_dBm - nrMRMinCouplingLoss_dB + (-5)*(0) + 1.5*0
                            otherwise
                                adjChInf_dBm_set = isAdjChIntferSource_NR;
                        end
                        adjChInfSig = powerDbm(wv_adjChInf,'set',adjChInf_dBm_set,'dBm',dlSlots);
                        adjChInf_dBm = powerDbm(adjChInfSig,'rms',[],'dBm',dlSlots);

                        % plus rxAntSignal
                        rxAntSig_adjChInf = rxAntSig*(~isLNA_IIP3_debug) + adjChInfSig;
                        rxAntSig_adjChInfer_dBm = powerDbm(rxAntSig_adjChInf,'rms',[],'dBm',dlSlots)

                        % export - rxAntSignal
                        rxAntSig = rxAntSig_adjChInf;
                    else
                        adjChInf_dBm = [];
                    end

                    if 0 % demod
                        demod_antIn = wlanDemodulation(rxAntSig, wvParams, 'rx')
                        pltComm.sa(pltComm, rxAntSig, {102901,'rxAntSignal'});
                    end
% Receiver - Apply Noise Floor to LNA Input Signal

                    % *** get signal combined with noise floor for evm's reference ***
                    if 0
                        t0_deg_set = 25^1; % set init. noise temperature
                        [rxAntSig_n0, noise0] = rfSim_noise(rxAntSig, 'NoiseTemperature', t0_deg_set, fs);
                    else
                        n0_dBmHz_set = -174 - 10*log10(size(rxAntSig,2)) % set init. noise psd in dBm/Hz
                        [rxAntSig_n0, noise0] = rfSim_noise(rxAntSig, 'NoisePowerDbmHz', n0_dBmHz_set, fs);
                    end
                    noise0_dBmHz = powerDbm(sum(noise0,2),'rms',[],fs,dlSlots);

                    if n0_dBmHz_set >= -174
                        lnaInSig = rxAntSig_n0;
                    else
                        % *** get input LNA signal lnaInSignal from rxAntSignal_n0 with adding thermal noise ***
                        if 0
                            [lnaInSig, noiseIn] = rfSim_lna(rxAntSig_n0, {fs,0,0});
                        elseif 1*1
                            t0_deg_set_2 = 0
                            [lnaInSig, noiseIn] = rfSim_noise(rxAntSig_n0, 'NoisePowerDbmHz', 10*log10(1000*1.380649e-23*(273+t0_deg_set_2)), fs);
                        else
                            lnaInSig = rxAntSig_n0;
                        end
                    end

                    try % demod
                        wvParams.finalSignal = lnaInSig;
                        demod_lnaIn = wlanDemodulation(lnaInSig, wvParams, 'rx')
                        PER_percent = demod_lnaIn.rxPER_percent;
                        PassDecod_rx = mean(PER_percent, 'all') < 10;
                    catch
                        PassDecod_rx = 0;
                    end

                    if isPropagationChannelModel * ~PassDecod_rx
                        continue
                    end

                    if 0
                        pltComm.sa(pltComm, rxAntSig, {0418, ['wanted']});
                        pltComm.sa(pltComm, rxAntSig_adjChInf, {0418, ['wanted plus w/ adj.ch.interferer']});
                        pltComm.sa(pltComm, rxAntSig_n0, {0418, ['ant w/ thermal noise ', num2str(noise0_dBmHz,5), 'dBm/Hz']});
                        pltComm.sa(pltComm, lnaInSig, {0418, ['LNA input w/ thermal noise ', num2str(noise0_dBmHz,5), 'dBm/Hz']});
                        powerDbm(rxAntSig(dlSlots))
                    end
                end
% Receiver - Characterize LNA

                if 1 % demod
                    demod_antIn = wlanDemodulation(lnaInSig, wvParams, 'rx')
                end
                if isTRSwitch % apply TRSwitch model
                    trswModel.LinearGain    = -0.5;
                    trswModel.IIP3          = 20 + trswModel.LinearGain;
                    trswModel.NFdB          = - trswModel.LinearGain;
                    trswModel.SampleRate    = fs;
                    [lnaInSig, ~] = rfSim_lna(lnaInSig, trswModel);
                    if 0 % demod
                        demod_trsw = wlanDemodulation(lnaInSig, wvParams, 'rx')
                    end
                end

                switch isLNA_Bypass % set LNA model
                    case {0,'HG','High Gain Mode'}
                        if IterationNo==1
                            % clear lnaModel % removed clear
                            lna_NF_requ_dB = 2+2;
                            lna_IIP3_requ_dBm = 7 - 0.5;
                            lna_IP1dB_requ_dBm = 0;
                            lnaModel.LinearGain = 16;
                            if 1
                                lnaModel.TOISpecification = 'IIP3';
                            else
                                % lnaModel.TOISpecification = 'IP1dB';
                            end

                            if 0 % add margin
                                lna_toi_dB_margin = - 1;
                                lna_nf_dB_margin = 1;
                            else
                                lna_toi_dB_margin = 0;
                                lna_nf_dB_margin = 0;
                            end
                        end

                        switch lnaModel.TOISpecification
                            case "IIP3"
                                lnaModel.IIP3 = lna_IIP3_requ_dBm + (0 + lna_iip3Step_dB*(IterationNo-1))*(-is_sweep_lnaTOI) + lna_toi_dB_margin;
                            case "IP1dB"
                                lnaModel.IP1dB = lna_IP1dB_requ_dBm + (0 + lna_ip1dBStep_dB*(IterationNo-1))*(-is_sweep_lnaTOI) + lna_toi_dB_margin;
                        end
                        % if RESETPARMS==1
                        lnaModel.NFdB = (lna_NF_requ_dB+lna_nfStep_dB*(1*IterationNo-1)*is_sweep_lnaNF) + lna_nf_dB_margin;
                        % else
                        %     lnaModel.NFdB = lnaModel.NFdB - 0.5;
                        % end
                        if lnaModel.NFdB<=0
                            break % SIMULATING_2
                        end
                        lnaModel.SampleRate = fs;
                        [lnaSig, lnaNoise] = rfSim_lna(lnaInSig, lnaModel); % create LNA signal
                    case {1,'BP','ByPass Mode'}
                        lnaModel.LinearGain = -3;
                        lnaModel.NFdB = - lnaModel.LinearGain;
                        lnaModel.SampleRate = fs;
                        lnaModel.IIP3 = 1000;
                        [lnaSig, lnaNoise] = rfSim_lna(lnaInSig, lnaModel); % create LNA signal

                end

                if 0 % demod
                    demod_lna = wlanDemodulation(lnaSig, wvParams, 'rx')
                    pltComm.sa(pltComm, lnaSig, {0418, ['LNA output']});
                    powerDbm(lnaNoise(dlSlots)) - 10*log10(fs)
                    powerDbm(lnaSig(dlSlots)) - powerDbm(lnaNoise(dlSlots))
                    evm(lnaInSig(1:12e-6*fs), lnaSig(1:12e-6*fs))

                    [SNR_lnaOut_dB, NFcal_lna_dB] = cal_snr(rxAntSig_n0,lnaSig,lnaInSig,cbw*[-0.5 0.5],fs)
                    [SNR_lnaOut_dB, NFcal_lna_dB] = cal_snr(rxAntSig_n0,lnaSig,lnaInSig)
                end

                % export
                wvParams.finalSignal = lnaSig;

% Receiver - Apply Flatness

                if isLNA_Flatness
                    if IterationNo==1 % set flatness response
                        lna_flatnessMax_dB = -1 - 0.5*1;
                        if configParamsWG.osr==1
                            lna_flatMagn_dB = [0 -0.5]
                            lna_flatFreq_Hz = [0 cbw/2]
                        elseif configParamsWG.osr>2
                            lna_flatMagn_dB = [0 -0.5 lna_flatnessMax_dB lna_flatnessMax_dB]
                            lna_flatFreq_Hz = [0 cbw/2 cbw/4+fs/4 fs/2]
                        else
                            error('check the oversampling ratio')
                        end
                        % end
                    end

                    [lnaSig_flat, b_flat] = rfSim_fltaness(lnaSig,fs,lna_flatMagn_dB,lna_flatFreq_Hz,10e3); % create flatness signal

                    % export
                    wvParams.finalSignal = lnaSig_flat;

                    if 0 % demod
                        demod_lnaFlat = wlanDemodulation(lnaSig_flat,wvParams,'rx')
                        powerDbm(lnaSig_flat(dlSlots))
                        powerDbm(lnaSig(dlSlots))
                        pltComm.aclr(pltComm, lnaSig_flat, {0415,'lna signal flatness','ACLR'})
                    end
                end
% Receiver - Apply wlanDemodulation and Measurements


                if isLNA_IIP3_debug
                    PER_percent = 0
                end

                if 1 % calculate SNR of LNA
                    [EVM_lnaIn, SNR_lnaIn_dB] = evm(rxAntSig_n0, lnaInSig);
                    [EVM_lnaOut, SNR_lnaOut_dB] = evm(rxAntSig_n0, lnaSig);
                    NFcal_lna_dB = SNR_lnaIn_dB - SNR_lnaOut_dB
                end

                % import signal
                receiverInSig = wvParams.finalSignal;

                if isReceiverNoise
                    % modulate signal with noise
                    if 0
                        % set receiver noise in dBm/Hz
                        receiverNoiseDbmHz = -160;
                        [receiverOutSig, receiverNoise] = rfSim_noise(receiverInSig, 'NoisePowerDbmHz', receiverNoiseDbmHz, fs);
                    else
                        receiverModel.LinearGain = 0;
                        receiverModel.NFdB = 12;
                        receiverModel.SampleRate = fs;
                        [receiverOutSig, receiverNoise] = rfSim_lna(receiverInSig, receiverModel);
                    end

                    receiverNoiseDbmHz = powerDbm(receiverNoise) - 10*log10(fs); % check noise psd
                else
                    receiverOutSig = receiverInSig;
                end

                % export signal
                wvParams.finalSignal = receiverOutSig;

                if 0
                    pltComm.sa(pltComm, receiverOutSig, {0418, ['Receiver output']});
                end

                try % demod
                    demod_final = wlanDemodulation(receiverOutSig, wvParams, 'rx')
                    PER_percent = demod_final.rxPER_percent;
                    PassDecod_rx = mean(PER_percent, 'all') < 10;
                catch
                    PassDecod_rx = 0;
                end

                PER_iter(IterationNo) = PER_percent;
                EVM_RMS_iter(IterationNo) = max(demod_final.evmRMS,[],'all');
                EVM_Subcarrier_iter(IterationNo) = max(cell2mat(demod_final.evmPerSc),[],'all');
                ACR_dB_iter(IterationNo) = adjChInf_dBm - rxAnt_dBm;
                AdjInterferer_dBm_iter(IterationNo) = adjChInf_dBm;
                IIP3_dBm_iter(IterationNo) = lnaModel.IIP3;
                NFdB_Cal_iter(IterationNo) = NFcal_lna_dB;
                NFdB_Set_iter(IterationNo) = lnaModel.NFdB;

                if ~is_sweep_params
                    SIMULATING_2 = 0;
                elseif ~PassDecod_rx
                    SIMULATING_2 = 0;
                else
                    IterationNo = IterationNo + 1;
                    close all
                    continue % SIMULATING_2
                end

            end % SIMULATING_2 end

            if isPropagationChannelModel && ~PassDecod_ch
                continue % SIMULATING_1
            else
                SIMULATING_1 = 0;
            end

        end % SIMULATING_1 end

% Receiver - Summary

        if 0
            mappingPass = containers.Map({1,0},{"Pass","Failed"}); % mapping
            summary.Format = dispLegd_format;
            summary.LnaGainDb = lnaModel.LinearGain;
            summary.LnaIIP3Dbm = lnaModel.IIP3;
            summary.LnaNFDb = lnaModel.NFdB;
            summary.FlatnessResponseCell = {b_flat};
            summary.RxSensitivityResult = string(arrayfun(@(x)mappingPass(x), demod_final.isRxSensitivityPass, 'UniformOutput', false))';
            summary.TotalPERPerct = demod_final.rxPER_percent
            summary.EVM_RMS = mean(demod_final.evmRMS');
        end
        summary_rx.Format = configParamsWG.wlanPHY;
        summary_rx.MCS = MCS_index;
        summary_rx.CBW = string(cfg.ChannelBandwidth);
        summary_rx.TransmitterDbm = transmit_dBm;
        summary_rx.ReceiverDbm = rxAnt_dBm;
        summary_rx.SensitivityDbm_Req = SENS_Requ_dBm;
        summary_rx.WantedSignalDbm = rxAnt_dBm;
        if isPropagationChannelModel
            summary_rx.CH_DelayProfile = tgaxChannel.DelayProfile;
            summary_rx.CH_TransmitReceiveDistanceMeter = tgaxChannel.TransmitReceiveDistance;
            summary_rx.CH_NumTransmitAntennas = tgaxChannel.NumTransmitAntennas;
            summary_rx.CH_NumReceiveAntennas = tgaxChannel.NumReceiveAntennas;
            summary_rx.CH_CarrierFrequencyGHz = tgaxChannel.CarrierFrequency/1e9;
            summary_rx.CH_SNRDb = -mean(demod_propChanModl.evmRMS,'all');
        end
        if is_sim_AdjacentChannelRejection
            summary_rx.AdjChannelIntefererDbm = AdjInterferer_dBm_iter(max(1, IterationNo-1));
            summary_rx.ACRDb_Req = ACR_Requ_dB;
        else
            summary_rx.AdjChannelIntefererDbm = nan;
            summary_rx.ACRDb_Req = nan;
        end
        if 1
            summary_rx.ANTin_SNRDb = - max(demod_antIn.evmRMS);
        end
        if isTRSwitch
            summary_rx.TRSW_lossDb = trswModel.LinearGain;
        else
            summary_rx.TRSW_lossDb = 0;
        end
        if isLNA_Bypass
            summary_rx.LNA_mode = 'BP';
        else
            summary_rx.LNA_mode = 'HG';
        end
        summary_rx.LNA_GainDb = lnaModel.LinearGain;
        summary_rx.LNA_IIP3 = IIP3_dBm_iter(max(1, IterationNo-1));
        summary_rx.LNA_NFDb = NFdB_Set_iter(max(1, IterationNo-1));

        if isLNA_Flatness
            summary_rx.LNA_FltanessMaxDb = lna_flatnessMax_dB;
        else
            summary_rx.LNA_FltanessMaxDb = nan;
        end
        if isReceiverNoise
            % summary_rx.RxNoiseDbmHz = receiverNoiseDbmHz;
            summary_rx.ReceiverGainDb = receiverModel.LinearGain;
            summary_rx.ReceiverNFDb = receiverModel.NFdB;
        else
            summary_rx.ReceiverGainDb = nan;
            summary_rx.ReceiverNFDb = nan;
        end

        summary_rx.Decoded_SNRDb = -1*round(EVM_RMS_iter(max(1, IterationNo-1)),2);
        summary_rx.Decoded_PER_percent = PER_iter(max(1, IterationNo-1));
        summary_rx.Decoded_SENS_Pass = (summary_rx.Decoded_PER_percent<10);
        summary_tab_tmp = struct2table(summary_rx); % convert struct to table
        if exist("summary_tab_rx",'var') % vertcat
            summary_tab_rx = vertcat(summary_tab_rx,summary_tab_tmp)
        else
            summary_tab_rx = summary_tab_tmp
        end
end % switch is_simulation_case


%%
function [SNROut_dB, NFcal_dB] = cal_snr(ref,mea_output,mea_input,bwInband,fs)
isFFTFreqSelect = 1;
if ~exist('bwInband','var')||isempty(bwInband)
    isFFTFreqSelect = 0;
elseif isscalar(bwInband)
    freq_1 = -bwInband/2;
    freq_end = bwInband/2;
else
    freq_1 = bwInband(1);
    freq_end = bwInband(2);
end

if isFFTFreqSelect
    freq_vec = fs/numel(ref) * [-numel(ref)/2+1:numel(ref)/2];
    [~,idx_cbw_1] = min(abs(freq_vec-freq_1));
    [~,idx_cbw_end] = min(abs(freq_vec-freq_end));
    idx_cbw = [idx_cbw_1:idx_cbw_end]';
    REF = fftshift(fft(ref));
    MEA_OUT = fftshift(fft(mea_output));
    MEA_IN = fftshift(fft(mea_input));
    ref = ifft(fftshift(REF(idx_cbw)));
    mea_output = ifft(fftshift(MEA_OUT(idx_cbw)));
    mea_input = ifft(fftshift(MEA_IN(idx_cbw)));
end

[~, SNRIn_dB] = evm(ref, mea_input);
[~, SNROut_dB] = evm(ref, mea_output);
NFcal_dB = SNRIn_dB - SNROut_dB;
end