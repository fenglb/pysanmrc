from sanmrc.comparemethods import *
from sanmrc.predata import *
# DP4 parameters
meanC = 0.0
meanH = 0.0
stdevC = 2.306
stdevH = 0.18731058105269952
degreeC = 11.38
degreeH = 14.18
 
calc_a_c13 = ([74.25543650739220000000, 48.44569810240660000000, 176.85339877593200000000, 11.67264309668400000000, 138.21012855822000000000, 121.89729430807700000000, 122.50376880179900000000, 121.05681827872700000000, 58.83075033676600000000, 15.62731956011590000000])
calc_b_c13 = ([76.54606326018590000000, 48.31514100521790000000, 174.56828551337200000000, 18.42221331063060000000, 139.58787236106400000000, 121.50248362594900000000, 122.56072415709800000000, 121.01379975085700000000, 58.34041422828900000000, 15.42125704559360000000])
expValue = ([76.1218 ,47.0828, 175.6654, 14.2679 , 141.5826, 126.5307, 128.2366, 127.7682, 60.5621, 13.9824])
isomers = [ calc_a_c13, calc_b_c13 ]

cdp4_s = []
output = {}
i = 0
for isomer in isomers:
    scaled_value = scaledValue( isomer, expValue )
    #cdp4_s.append( calculateCDP4( scaled_value, expValue, meanC, stdevC ) )
    cdp4_s.append( calculateTDP4( scaled_value, expValue, meanC, stdevC, degreeC ) )
    output[i] = []
    output[i].append("\nAssigned " + "C" + " shifts: (label, calc, corrected, exp, error)" )
    output[i] += [format(isomer[j], "6.2f") + ' ' +
                    format(scaled_value[j], "6.2f") + ' ' +
                    format(expValue[j], "6.2f") + ' ' +
                    format(expValue[j]-scaled_value[j], "6.2f") for j in range(0, len(isomer) ) ]
    i = i+1
results = list( 100*x/sum(cdp4_s) for x in cdp4_s )
for i, res in enumerate(results):
    output[i].append( "\nResults of DP4 using t distrubution:" )
    output[i].append( "Isomer :" + str(i+1) + ": " + format(res, "4.1f") + "%")
    for line in output[i]:
        print( line )
