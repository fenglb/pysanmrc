def GenDSInchis(inchi):

    ilist = list(inchi)
    #Inchis of all diastereomers, including the parent structure
    resinchis = []

    #get the number of potential diastereomers
    layers = inchi.split('/')

    sterelayer = filter( lambda item: item.startswith("t"), layers )[0]

    numds = 2**( len( sterelayer.translate(None, 't,1234567890')) - 1 )

    print "Number of diastereomers to be generated: " + str(numds)

    #find configuration sites (+ and -)
    bs = ilist.index('t')
    es = ilist[bs:].index('/')
    spos = []
    for s in range(bs, bs+es):
        if ilist[s] == '+' or ilist[s] == '-':
            spos.append(s)

    temps = []
    #Generate inversion patterns - essentially just binary strings
    for i in range(0, numds):
        template = bin(i)[2:].zfill(len(spos)-1)
        temps.append(template)

    #For each 1 in template, invert the corresponding stereocentre
    #and add the resulting diastereomer to the list
    invert = {'+': '-', '-': '+'}

    for ds in range(0, numds):
        t = list(ilist)
        for stereocentre in range(1, len(spos)):
            if temps[ds][stereocentre-1] == '1':
                t[spos[stereocentre]] = invert[t[spos[stereocentre]]]
        resinchis.append(''.join(t))

    return resinchis
   
if __name__ == "__main__": 
    print GenDSInchis("InChI=1/C6H8O6/c7-1-2(8)5-3(9)4(10)6(11)12-5/h2,5,7-10H,1H2/t2-,5+/m0/s1")
