function LORLINE = Lorenz(T2,T1, H0, Hcentr)
    % ------------------------------
    % numpt is the number of points in the spectrum
    numpt = 2^10;
    % sweep is the spectral range in G
    sweep = 100;
    % H0 is the center of the line - corresponds to the Larmor frequency
    % H0=3320
    % Hcentr is the center of the scan
    step = sweep/(numpt-1);
    H = Hcentr-sweep/2:step:Hcentr+sweep/2;
    % T2 relaxation time (s)
    % T2=1E-7
    % T1 relaxation time (s)
    % T1=1E-6
    gammaE = 1.76E7;
    % h1 is the amplitude of the microwave field (G)
    h1 = 0.1;
    % DELTA is the difference the RF frequency 
    % and the Larmor resonance frequency
    DELTA = gammaE*(H-H0);
    R2SATSQ = (1/T2^2)*(1+T1*T2*gammaE^2*h1^2);
    % SP=(gammaE*h1/T2)./(DELTA.^2+R2SATSQ)
    SPDER = -2*(DELTA*gammaE*h1/T2)./(DELTA.^2+R2SATSQ).^2;
    LORLINE = [H;SPDER]';
end
