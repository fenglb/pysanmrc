# Stereochemical Assignment of Diastereoisomers Using NMR Shift Calculation

###Introduction

DFT calculations for prediction of NMR spectra can decide which of the predictions match the experimental NMR data best. This technique has scince played key roles in the stereostructure assignment, many different ways to quantitatively evaluate goodness of fix between calculational and exprimental data have been developed, including mean absolute error, correlation coefficient, CP3[smith2009], DP4[smith2010] and DP4+[grimblat2015].

###Project
Build a web site for the stereochemical assignment using the quantum chemical calcaulations of nmr shifts and experimental data input by the user. The ways quantitatively evaluated goodness of fix include **mean absolute error**, **correlation coefficient**, **CP3**, **DP4**, and **DP4+**. 

The statical parameters referenced the supporting information of these papers[smith2009][smith2010][grimblat2015].

Table 1: Expectation values and standard deviations in normal distribution of CP3 for both a correct assignment and an incorrect assignment.[smith2009]

|       |C13     |H1      |C13+H1   |
|--------|-----|-------|------|
|Correct | 0.547(0.253) | 0.478(0.305) | 0.512(0.209) |
|Incorrect | -0.487(0.533) | -0.786(0.835) | -0.637(0.499) |

Table 2: Values of std and degree for the fitted t-distribution in a scaled calculated shift for DP4[smith2010].

|       |C13     |H1    |
|--------|-----|-------|
|`\mu`  |0      |0      |
|`\sigma`  |2.306      |0.185     |
|`\zeta`  |11.38    |14.18    |

Table 3: Values of the [`\mu`, `\sigma`, `\zeta`] statistical parameters of DP4+.

| Function          |     | B3LYP    | B3LYP    | B3LYP    | B3LYP     | B3LYP    | B3LYP     | B3LYP     | B3LYP    | B3LYP      | B3LYP    | B3LYP    | B3LYP    | mPW1PW91 | mPW1PW91  | mPW1PW91 | mPW1PW91  | mPW1PW91  | mPW1PW91 | mPW1PW91   | mPW1PW91 | mPW1PW91 | mPW1PW91 | mPW1PW91 | mPW1PW91  |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|-------------------|-----|----------|----------|----------|-----------|----------|-----------|-----------|----------|------------|----------|----------|----------|----------|-----------|----------|-----------|-----------|----------|------------|----------|----------|----------|----------|-----------|-----|-----------|-----------|-----|------------|-----|----------|----------|-----|-----------|-----|-----------|-----------|-----|------------|-----| 
| Solvent           |     | GasPhase | GasPhase | GasPhase | GasPhase  | GasPhase | GasPhase  | PCM       | PCM      | PCM        | PCM      | PCM      | PCM      | GasPhase | GasPhase  | GasPhase | GasPhase  | GasPhase  | GasPhase | PCM        | PCM      | PCM      | PCM      | PCM      | PCM       |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
| Method |     | 6-31G(d) | 6-31G(d , p)      | 6-31+G(d , p)      | 6-311G(d) | 6-311G(d , p)      | 6-311+G(d , p)   | 6-31G(d) | 6-31G(d , p)      | 6-31+G(d , p)      | 6-311G(d) | 6-311G(d , p)      | 6-311+G(d , p)      | 6-31G(d) | 6-31G(d , p)      | 6-31+G(d , p) | 6-311G(d) | 6-311G(d , p) | 6-311+G(d , p) | 6-31G(d) | 6-31G(d , p) | 6-31+G(d , p) | 6-311G(d) | 6-311G(d , p) | 6-311+G(d , p) | 
| TMS               | C   | 189.7149 | 191.55   | 192.5898 | 183.6135  | 183.8904 | 183.4253  | 190.0717  | 191.8936 | 193.0714   | 184.1054 | 184.3504 | 183.8888 | 193.8603 | 195.5355  | 196.0955 | 188.8451  | 189.1589  | 188.6787 | 194.2673   | 195.9285 | 196.6095 | 188.2374 | 188.4876 | 188.0185  |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | H   | 32.1826  | 31.7142  | 31.6377  | 32.2012   | 31.9154  | 31.9017   | 32.1769   | 31.7082  | 31.6354    | 32.1985  | 31.9135  | 31.9016  | 32.1159  | 31.6431   | 31.5661  | 32.3737   | 32.0829   | 32.0661  | 32.1073    | 31.6342  | 31.56    | 32.0903  | 31.795   | 31.7794   |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
| nonscaledC13(sp2) | μ   | -6.16    | -4.62    | -2.68    | 5.92      | 6.54     | 7.37      | -5.23     | -3.7     | -1.42      | 7.28     | 7.81     | 8.92     | -4.79    | -3.35     | -2.18    | 6.08      | 6.78      | 7.52     | -3.79      | -2.37    | -0.92    | 6.32     | 6.88     | 7.91      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | σ   | 2.49     | 2.46     | 2.29     | 2.35      | 2.24     | 2.25      | 1.84      | 1.84     | 1.99       | 2.65     | 2.41     | 3.06     | 2.41     | 2.3       | 2.1      | 2.04      | 2.05      | 1.97     | 1.72       | 1.67     | 1.75     | 2.31     | 2.14     | 2.7       |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | ν   | 6.53     | 6.71     | 8.05     | 5.86      | 5.03     | 6.19      | 4.62      | 4.86     | 6.56       | 14.74    | 8.78     | 77.02    | 7.03     | 6.53      | 6.55     | 4.28      | 4.12      | 4.84     | 4.48       | 4.67     | 5.36     | 9.6      | 6.66     | 46.56     |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | Σ-n | 2.96     | 2.91     | 2.64     | 2.88      | 2.85     | 2.71      | 2.41      | 2.37     | 2.38       | 2.85     | 2.74     | 3.11     | 2.83     | 2.74      | 2.51     | 2.78      | 2.82      | 2.51     | 2.27       | 2.19     | 2.19     | 2.59     | 2.55     | 2.76      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
| nonscaledC13(sp3) | μ   | 1.3      | 2.05     | 3.58     | 4.73      | 5.26     | 5.16      | 1.6       | 2.34     | 4.04       | 5.16     | 5.68     | 5.6      | 0.75     | 1.41      | 2.43     | 4.28      | 4.85      | 4.75     | 1.09       | 1.74     | 2.91     | 3.64     | 4.15     | 4.07      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | σ   | 1.65     | 1.76     | 1.9      | 2.28      | 2.21     | 2.29      | 1.59      | 1.75     | 1.9        | 2.35     | 2.27     | 2.4      | 1.48     | 1.52      | 1.61     | 1.74      | 1.69      | 1.81     | 1.43       | 1.51     | 1.6      | 1.86     | 1.79     | 1.94      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | ν   | 6.29     | 7.32     | 7.18     | 15.73     | 14.15    | 18.58     | 6.06      | 7.41     | 7.36       | 20.7     | 17.41    | 23.24    | 5.69     | 6.07      | 6.04     | 6.3       | 6.18      | 8.74     | 5.56       | 6.22     | 6.27     | 9.4      | 8.59     | 11.16     |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | Σ-n | 1.99     | 2.07     | 2.25     | 2.45      | 2.39     | 2.43      | 1.94      | 2.05     | 2.23       | 2.48     | 2.42     | 2.51     | 1.83     | 1.85      | 1.97     | 2.12      | 2.07      | 2.07     | 1.79       | 1.83     | 1.96     | 2.1      | 2.05     | 2.15      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
| scaledC13         | μ   | 0        | 0        | 0        | 0         | 0.02     | 0.01      | 0         | 0        | 0          | 0        | 0.01     | 0        | 0        | 0         | 0        | -0.02     | -0.01     | -0.01    | 0          | 0        | 0        | -0.01    | 0        | -0.01     |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | σ   | 1.58     | 1.85     | 2.09     | 1.92      | 1.84     | 1.83      | 1.53      | 1.78     | 2.03       | 1.84     | 1.74     | 1.8      | 1.27     | 1.45      | 1.6      | 1.6       | 1.55      | 1.58     | 1.23       | 1.39     | 1.56     | 1.56     | 1.49     | 1.57      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | ν   | 5.97     | 8.74     | 10.39    | 7.16      | 6.34     | 7.39      | 6.09      | 8.78     | 10.17      | 7.32     | 6.37     | 7.94     | 4.19     | 5.33      | 6.11     | 4.9       | 4.48      | 5.76     | 4.43       | 5.58     | 6.23     | 5.8      | 5.23     | 6.67      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | Σ-n | 1.94     | 2.12     | 2.33     | 2.27      | 2.22     | 2.15      | 1.87      | 2.04     | 2.28       | 2.17     | 2.1      | 2.09     | 1.73     | 1.83      | 1.96     | 2.08      | 2.06      | 1.96     | 1.65       | 1.75     | 1.9      | 1.93     | 1.9      | 1.9       |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
| nonscaledH1(sp2)  | μ   | -0.17    | 0.02     | 0.14     | -0.04     | 0.11     | 0.16      | -0.09     | 0.1      | 0.24       | 0.06     | 0.21     | 0.27     | -0.04    | 0.14      | 0.24     | 0.31      | 0.48      | 0.52     | 0.05       | 0.23     | 0.35     | 0.13     | 0.3      | 0.34      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | σ   | 0.19     | 0.18     | 0.16     | 0.17      | 0.15     | 0.14      | 0.17      | 0.15     | 0.13       | 0.15     | 0.12     | 0.12     | 0.18     | 0.17      | 0.15     | 0.17      | 0.14      | 0.14     | 0.16       | 0.14     | 0.12     | 0.15     | 0.12     | 0.11      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | ν   | 13.81    | 9.2      | 6.84     | 10.93     | 7.27     | 5.78      | 10.34     | 7.27     | 5.35       | 8.85     | 5.99     | 5.2      | 13.04    | 8.6       | 5.84     | 10.21     | 6.25      | 5.39     | 10.58      | 6.95     | 4.91     | 8.76     | 5.32     | 4.88      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | σ   | 0.21     | 0.2      | 0.2      | 0.19      | 0.18     | 0.18      | 0.19      | 0.18     | 0.17       | 0.17     | 0.16     | 0.16     | 0.2      | 0.19      | 0.18     | 0.19      | 0.17      | 0.17     | 0.18       | 0.17     | 0.16     | 0.17     | 0.15     | 0.15      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
| nonscaledH1(sp3)  | μ   | -0.06    | -0.07    | -0.06    | -0.03     | -0.07    | -0.05     | -0.01     | -0.03    | 0          | 0.02     | -0.02    | 0.01     | -0.06    | -0.08     | -0.07    | 0.24      | 0.19      | 0.21     | -0.02      | -0.04    | -0.02    | 0.01     | -0.04    | -0.01     |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | σ   | 0.17     | 0.15     | 0.13     | 0.15      | 0.13     | 0.14      | 0.15      | 0.13     | 0.12       | 0.13     | 0.1      | 0.11     | 0.16     | 0.14      | 0.13     | 0.15      | 0.12      | 0.13     | 0.14       | 0.12     | 0.11     | 0.13     | 0.1      | 0.11      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | ν   | 6.95     | 4.83     | 4.05     | 3.92      | 3.34     | 3.53      | 7.44      | 4.72     | 3.97       | 3.47     | 2.76     | 3.03     | 6.13     | 4.22      | 3.68     | 3.88      | 3.12      | 3.39     | 6.26       | 4.03     | 3.65     | 3.54     | 2.6      | 2.89      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | Σ-n | 0.2      | 0.19     | 0.18     | 0.21      | 0.19     | 0.2       | 0.18      | 0.17     | 0.16       | 0.18     | 0.17     | 0.18     | 0.19     | 0.18      | 0.18     | 0.2       | 0.19      | 0.19     | 0.16       | 0.16     | 0.16     | 0.18     | 0.17     | 0.18      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
| scaledH1          | μ   | 0        | 0        | 0        | 0         | 0        | 0         | 0         | 0        | 0          | 0        | 0        | 0        | 0        | 0         | 0        | 0         | 0         | 0        | 0          | 0        | 0        | 0        | 0        | 0         |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | σ   | 0.14     | 0.13     | 0.13     | 0.12      | 0.11     | 0.12      | 0.13      | 0.12     | 0.11       | 0.1      | 0.1      | 0.1      | 0.13     | 0.12      | 0.12     | 0.11      | 0.11      | 0.12     | 0.12       | 0.11     | 0.1      | 0.1      | 0.09     | 0.1       |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | ν   | 5.94     | 4.81     | 4.7      | 3.29      | 3.2      | 3.41      | 6.41      | 4.9      | 4.27       | 3.23     | 2.85     | 2.93     | 5.03     | 4.3       | 4.32     | 3.16      | 3.06      | 3.39     | 5.54       | 4.32     | 3.89     | 3.16     | 2.78     | 2.89      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 
|                   | Σ-n | 0.17     | 0.17     | 0.17     | 0.18      | 0.18     | 0.18      | 0.16      | 0.15     | 0.15       | 0.16     | 0.16     | 0.16     | 0.16     | 0.16      | 0.16     | 0.17      | 0.17      | 0.18     | 0.15       | 0.14     | 0.14     | 0.16     | 0.16     | 0.16      |     |           |           |     |            |     |          |          |     |           |     |           |           |     |            |     | 



###Aims

This web site will be provided more compatible and more conventient stereochemical assignment using NMR shift.

[The Goodman group](http://www-jmg.ch.cam.ac.uk/tools/nmr/) also provide the web site for stereostructure assignment using calculation of CP3 and DP4. The site is not compatible, because they uses the Java program run in the browser that can not work in popular browser like as firefox, chrome and IE defaultly.

the Grimblat use the Excel file to compute the DP4+ probability. But the `TDIST(x, deg_freedom, tails)` return the value of the t-distribution in Excel, the **Deg_freedom** is integer, that make some different with the Scipy/Python. You can download the Excel from [DP4+](http://sarotti-nmr.weebly.com/).

[smith2009]: Assigning the Stereochemistry of Pairs of Diastereoisomers Using GIAO NMR Shift Calculation, Steven G. Smith and Jonathan M. Goodman, J. Org. Chem. 2009, 74, 4597-4607
[smith2010]: Assigning Stereochemistry ot Single Diastereoisomers by GIAO NMR Calculation: The DP4 Probability, Steven G. Smith and Jonathan M. Goodman, J. Am. Chem. Soc. 2010, 132, 12946-12959
[grimblat2015]: Beyond DP4: an improved probability for the stereochemical assignment of isomeric compounds using quantum chemical calculations of NMR shifts, Grimblat, Nicolás Zanardi, María Marta Sarotti, Ariel Marcelo, J. Org. Chem. 2015, 80, 12526-12534
