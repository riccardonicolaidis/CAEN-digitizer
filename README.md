# Data Acquisition Silicon Detector Sr90_300um_500um Logbook

## Instrumentation

- Power Supply:
    - Rode & Schwarz NGA100 for the bias voltage
    - Rigol DP831A for feeding the preamplifiers

- Silicon Detector:
    - Mirion PIPS 300 um (preamplifier manufactured by AGE scientific)
    - Mirion PIPS 500 um (preamplified embedded)

- Data Acquisition System:
    - Hardware:
        - CAEN DT1081B Four Fold Programmable Logic Unit
        - CAEN DT5730B digitizer
        - CEN Fani In - Fan Out module (N625)
        - CAEN 16 Channels Fast Amplifier (N979)
        - CAEN 4 channels Variable Gain Amplifier (N470)
    - Firmware:
        - Waveform Recording Firmware
    - Software
        - Wavedump
        - Web Browser for configuring the DT1081B Programmable Logic Unit


## Experimental Setup

- High Voltages:
    - Voltage: 70 V
    - Current Limit: 0.04 A

- Preamplifiers:
    - AGE Preamplifier: +5V, GND
    - Embedded Preamplifier: +12V, -12V, GND


### Channel configuration - Cabling

- PIPS 300 um:
    - Connected to the AGE Preamplifier
    - Preamplifier is connected to an AC decoupling capacitor
    - Connection to CH1 Prgogrammable gain amplifier G = 2
    - To CH0 of the digitizer

- PIPS 500 um:
    - Embedded preamplifier
    - Connected to CH0 of the fast amplifier (N979)
    - Connected to CH0 of the programmable gain amplifier G = 4 (total gain of x40)
    - To Fan In - Fan Out module (N625) CH0
        - One output to the digitizer CH1
        - One output to the programmable logic unit (DT1081B) CHA: The output is connected to the TRIG-IN port of the Digitizer (DT5730B)

### Configuration of the Programmable Logic Unit
- Wire mode
- Discriminator mode
- Input 50 Ohm
- Output monostabled TTL 1000 ns
- Threshold 108 mV

### Configuration of the Digitizer
See the configuration file:
- Frequency: 500 MHz
- Total number of samples 5000
- CH0: PIPS 300 um - AGE preamplifier (negative polarity)
- CH1: PIPS 500 um - Embedded preamplifier (positive polarity)



# Fit results - Calibration with Am 241

## Americium source near the PIPS 300 um
### PIPS 300 um
The model is a double gaussian with a linear background.
```
1  p0           4.34705e+02   5.13893e+00   1.99631e-02  -8.24896e-05
2  p1           3.81086e+01   1.50040e+00   5.67214e-03   6.01980e-04
3  p2           1.15323e+02   1.46675e+00   3.63943e-03  -4.19998e-04
4  p3           2.49819e+02   4.04548e+00   1.46587e-02   5.83921e-05
5  p4           4.75022e+02   2.31214e+00   8.38836e-03   3.78241e-06
6  p5           1.26240e+02   2.11223e+00   5.85899e-03   5.78542e-04
7  p6          -2.79570e+00   1.82051e+00   4.11245e-03  -3.67764e-04
```

1. `p0` = 434.705 +/- 5.13893 | Amplitude of the noise part (pedestal)
2. `p1` = 38.1086 +/- 1.50040 | Mean value of the noise part (pedestal)
3. `p2` = 115.323 +/- 1.46675 | Sigma of the noise part (pedestal)
4. `p3` = 249.819 +/- 4.04548 | Amplitude of the peak (Am 241 - 59.5 keV)
5. `p4` = 475.022 +/- 2.31214 | Mean value of the peak (Am 241 - 59.5 keV)
6. `p5` = 126.240 +/- 2.11223 | Sigma of the peak (Am 241 - 59.5 keV)
7. `p6` = -2.79570 +/- 1.82051 | Offset


### PIPS 500 um
The model is a single gaussian.
```
1  Constant     2.92032e+02   2.61387e+00   1.33932e-02  -1.45910e-06
2  Mean        -5.80919e-01   9.08120e-01   5.71560e-03   3.04754e-06
3  Sigma        1.24237e+02   6.49743e-01   8.88216e-06   2.63713e-03
```

1. `Constant` = 292.032 +/- 2.61387 | Amplitude of the pedestal
2. `Mean` = -0.580919 +/- 0.908120 | Mean value of the pedestal
3. `Sigma` = 124.237 +/- 0.649743 | Sigma of the pedestal


## Americium source near the PIPS 500 um
### PIPS 500 um
The model is a double gaussian with a linear background.
```
1  p0           2.48144e+02   6.53293e+00   1.96450e-02  -5.43737e-07
2  p1           9.22064e+01   2.01427e+00   8.11574e-03  -2.24308e-07
3  p2           8.31250e+01   2.58134e+00   6.15505e-03  -2.17902e-06
4  p3           6.65089e+02   8.25039e+00   2.60406e-02  -2.82972e-07
5  p4           7.04201e+02   1.53419e+00   5.37403e-03   1.00073e-06
6  p5           1.07718e+02   1.76933e+00   4.07410e-03  -2.08205e-06
7  p6           8.21911e+01   3.56876e+00   6.67527e-03  -1.49814e-06
```

1. `p0` = 248.144 +/- 6.53293 | Amplitude of the noise part (pedestal)
2. `p1` = 92.2064 +/- 2.01427 | Mean value of the noise part (pedestal)
3. `p2` = 83.1250 +/- 2.58134 | Sigma of the noise part (pedestal)
4. `p3` = 665.089 +/- 8.25039 | Amplitude of the peak (Am 241 - 59.5 keV)
5. `p4` = 704.201 +/- 1.53419 | Mean value of the peak (Am 241 - 59.5 keV)
6. `p5` = 107.718 +/- 1.76933 | Sigma of the peak (Am 241 - 59.5 keV)
7. `p6` = 82.1911 +/- 3.56876 | Offset 

### PIPS 300 um
The model is a single gaussian.
```
1  Constant     5.54824e+02   4.44118e+00   1.69917e-02   8.52710e-06
2  Mean        -7.74966e-01   6.47802e-01   3.11520e-03  -1.83426e-05
3  Sigma        9.90323e+01   5.08731e-01   6.63102e-06   7.32235e-02
```

1. `Constant` = 554.824 +/- 4.44118 | Amplitude of the pedestal
2. `Mean` = -0.774966 +/- 0.647802 | Mean value of the pedestal
3. `Sigma` = 99.0323 +/- 0.508731 | Sigma of the pedestal