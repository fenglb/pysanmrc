# Stereochemical Assignment of Diastereoisomers Using NMR Shift Calculation

###Introduction

DFT calculations for prediction of NMR spectra can decide which of the predictions match the experimental NMR data best. This technique has scince played key roles in the stereostructure assignment, many different ways to quantitatively evaluate goodness of fix between calculational and exprimental data have been developed, including mean absolute error, correlation coefficient, [CP3][smith2009], [DP4][smith2010] and [DP4+][grimblat2015].

###Project
Build a web site for the stereochemical assignment using the quantum chemical calcaulations of nmr shifts and experimental data input by the user. The ways quantitatively evaluated goodness of fix include **mean absolute error**, **correlation coefficient**, **CP3**, **DP4**, and **DP4+**. 

The statical parameters referenced the supporting information of [these papers][smith2009,smith2010,grimblat2015].

Table 1: Expectation values and standard deviations in normal distribution of CP3 for both a correct assignment and an incorrect assignment.

|       |C13     |H1      |C13+H1   |
|--------|-----|-------|------|
|Correct | 0.547(0.253) | 0.478(0.305) | 0.512(0.209) |
|Incorrect | -0.487(0.533) | -0.786(0.835) | -0.637(0.499) |

Table 2: Values of std and degree for the fitted t-distribution in a scaled calculated shift for DP4

|       |C13     |H1    |
|--------|-----|-------|
|`\mu`  |0      |0      |
|`\sigma`  |2.306      |0.185     |
|`\zeta`  |11.38    |14.18    |

Values of the [`\mu`, `\sigma`, `\zeta`] statistical parameters of DP4+ in `data/dp4_plus.csv` file.

###Aims

This web site will be provided more compatible and more conventient stereochemical assignment using NMR shift.

[The Goodman group](http://www-jmg.ch.cam.ac.uk/tools/nmr/) also provide the web site for stereostructure assignment using calculation of CP3 and DP4. The site is not compatible, because they uses the Java program run in the browser that can not work in popular browser like as firefox, chrome and IE defaultly.

the Grimblat use the Excel file to compute the DP4+ probability. But the `TDIST(x, deg_freedom, tails)` return the value of the t-distribution in Excel, the **Deg_freedom** is integer, that make some different with the Scipy/Python. You can download the Excel from [DP4+](http://sarotti-nmr.weebly.com/).

[smith2009]: Assigning the Stereochemistry of Pairs of Diastereoisomers Using GIAO NMR Shift Calculation, Steven G. Smith and Jonathan M. Goodman, J. Org. Chem. 2009, 74, 4597-4607

[smith2010]: Assigning Stereochemistry ot Single Diastereoisomers by GIAO NMR Calculation: The DP4 Probability, Steven G. Smith and Jonathan M. Goodman, J. Am. Chem. Soc. 2010, 132, 12946-12959

[grimblat2015]: Beyond DP4: an improved probability for the stereochemical assignment of isomeric compounds using quantum chemical calculations of NMR shifts, Grimblat, Nicolás Zanardi, María Marta Sarotti, Ariel Marcelo, J. Org. Chem. 2015, 80, 12526-12534
