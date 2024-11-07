%% 2024-03-26, draft
% 2024-05-02, compensate type of complex signal
% 2024-07-17, update n0DbmHz for case where gnDb<0
% 2024-07-17, update nAddGnDbmHz for case where gnDb<0

function [y, noiseAdd] = rfSim_lna(x,parms_struct)
    if 1 
        % init.
        pm = parms_struct;
        Rohm = 1;
        bw = 1;
        nfDb = 0;
        gnDb = 0;
        isTOI = 0;
        isNF = 0;
        noInDbmHz = nan;
        if isreal(x)
            isComplex = 0;
        else
            isComplex = 1;
        end

        if ~isstruct(pm)
            try
                bw = pm{1}; % cell
            catch
                bw = pm(1); % numeric
            end
            isNF = 1;
            try
                nfDb = pm{2};
            end
            try
                gnDb = pm{3};
            end
        elseif isstruct(pm)
            if isfield(pm,'SampleRate')&&~isempty(pm.SampleRate)
                bw = pm.SampleRate;
            end
            if isfield(pm,'NoiseInDbmHz')&&~isempty(pm.NoiseInDbmHz)&bw~=1
                noInDbmHz = pm.NoiseInDbmHz;
            elseif isfield(pm,'NoiseInDbm')&&~isempty(pm.NoiseInDbm)&&bw==1
                noInDbmHz = pm.NoiseInDbm;
            end
            if isfield(pm,'NFdB')&&~isempty(pm.NFdB)
                isNF = 1;
                nfDb = pm.NFdB;
            end
            if isfield(pm,'LinearGainDb')&&~isempty(pm.LinearGainDb)
                gnDb = pm.LinearGainDb;
            end
            if isfield(pm,'TOISpecification')&&~isempty(pm.TOISpecification)
                isTOI = 1;
            end
        end

        % get
        gn = 10^(gnDb/10);
        a1 = sqrt(gn);
        nSamps = size(x,1);
        nAnts = size(x,2);
    end

    if isTOI % p1dB and IP3 modeling
        switch pm.TOISpecification
            case 'IIP3'
                iip3 = 10^((pm.IIP3Dbm-30)/10);
                a3 = -2*a1^1 / (3*Rohm) / iip3;
            case 'OIP3'
                oip3 = 10^((pm.OIP3Dbm-30)/10);
                a3 = -2*a1^3 / (3*Rohm) / oip3;
            case 'IP1dB'
                ip1dB = 10^((pm.IP1dBDbm-30)/10);
                a3 = - 2*(10^-0.1)*(1-10^-0.05) / (3*Rohm) *a1^1 / ip1dB;
            case 'OP1dB'
                op1dB = 10^((pm.OP1dBDbm-30)/10);
                a3 = - 2*(10^-0.1)*(1-10^-0.05) / (3*Rohm) *a1^3 / op1dB;
        end
        c5 = 0;
        if isComplex % 2024-05-02, compensate type of complex signal
            a3 = a3*sqrt(2);
        end
        if 1
            yAmp = (a1)*x + a3*x.^3 + c5*x.^5;
        else
            yAmp = (a1)*x + a3*abs(x).^2.*x + c5*abs(x).^4.*x;
        end
        if 1
            a2 = 2
            yAmp = yAmp + a2*x.^2;
        end
        y = yAmp;
    else
        y = (a1)*x;
    end

    if isNF
        nf = 10^(nfDb/10);
        if nf == 0 
            nf = 1+eps;
        end
        kT0 = 1.380649e-23 * 290;
        if gnDb > 0
            n0DbmHz = 10*log10(1000*(nf-1)*kT0*gn);
        else
            % 2024-07-17, update n0DbmHz for case where gnDb<0
            n0DbmHz = 10*log10(1000*(nf)*kT0*1);
        end
        if gnDb > 0
            nAddGnDbmHz = 10*log10(1000*(nf-1)*gn*kT0);
        else
            % 2024-07-17, update nAddGnDbmHz for case where gnDb<0
            nAddGnDbmHz = 10*log10(1000*(nf)*kT0*1);
        end
        if ~isnan(noInDbmHz)
            nAddGnDbmHz = 10*log10(10^((noInDbmHz + gnDb)/10) + 10^(n0DbmHz/10));
        end
        [~,noiseAdd] = rfSim_noise([], 'NoisePowerDbmHz', nAddGnDbmHz, {bw, nSamps, isComplex, nAnts});
        powerDbm(noiseAdd) - 10*log10(bw);
        powerDbm(y) - powerDbm(x);
        y = y + noiseAdd;
    else
        noiseAdd = [];
    end

end