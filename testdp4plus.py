from sanmrc.comparemethods import *
from sanmrc.predata import *
# DP4+ parameters
mean_s = 0.0
stdev_s = 1.58055
degree_s = 5.9698
mean_u_sp = -6.15715
stdev_u_sp = 2.4881
degree_u_sp = 6.53169
mean_u_sp3 = 1.30169
stdev_u_sp3 = 1.64843
degree_u_sp3 = 6.28871

 
sp = [0, 0, 1, 0, 1, 1, 1, 1, 0, 0]
calc_a_c13 = ([74.25543650739220000000, 48.44569810240660000000, 176.85339877593200000000, 11.67264309668400000000, 138.21012855822000000000, 121.89729430807700000000, 122.50376880179900000000, 121.05681827872700000000, 58.83075033676600000000, 15.62731956011590000000])
calc_b_c13 = ([76.54606326018590000000, 48.31514100521790000000, 174.56828551337200000000, 18.42221331063060000000, 139.58787236106400000000, 121.50248362594900000000, 122.56072415709800000000, 121.01379975085700000000, 58.34041422828900000000, 15.42125704559360000000])
isomers = [ calc_a_c13, calc_b_c13 ]
expValue = ([76.1218 ,47.0828, 175.6654, 14.2679 , 141.5826, 126.5307, 128.2366, 127.7682, 60.5621, 13.9824])
cdp4_s = []
cdp4_u = []
for isomer in isomers:
    scaled_value = scaledValue( isomer, expValue )
    cdp4_s.append( calculateTDP4( scaled_value, expValue, mean_s, stdev_s, degree_s ) )

    usv_sp = [mean_u_sp, stdev_u_sp, degree_u_sp]
    usv_sp3 = [mean_u_sp3, stdev_u_sp3, degree_u_sp3]
    cdp4_u.append( calculateDP4plus( isomer, expValue, sp, usv_sp, usv_sp3) )

print( list(100*x/sum(cdp4_s) for x in  cdp4_s ) )
print( list(100*x/sum(cdp4_u) for x in cdp4_u ) )
cdp4 = list( x*y for x, y in zip(cdp4_s, cdp4_u) )
print( list(100*x/sum(cdp4) for x in cdp4 ) )
