from itertools import chain, product, permutations, combinations
from functools import partial
from nmrdata import NMRData, Molecular
import calculatemethods as cm
from predata import scaledValue
#import pandas
class StatisticalNMR:
    def __init__(self):
        self.nmrdata = []
        self.molecular = Molecular()
        self.statpara = cm.StatPara("../data/statparas.db")
        self.calculate_results = {}

        self.DFT_software = "Gaussian09"
        self.function     = "B3LYP"
        self.method       = "6-31G(d,p)"
        self.solvent      = "GasPhase"

    def setFunction(self, function):
        self.function = function
    def setMethod(self, method):
        self.method = method
    def setSolvent(self, solvent):
        self.solvent = solvent

    def getStatPara(self, name, dtype, dist='t'):
        usv = self.statpara.getStatParaData(name=name, 
            soft=self.DFT_software, 
            function=self.function,
            method=self.method,
            solvent=self.solvent, 
            dtype=dtype, dist=dist)
        if not usv:
            raise TypeError("Your Settings donot match for %s. No Statistical paraments found!" % name)
        return usv

    def readDataFromFile(self, filename):
        raw_calc_data = {}
        raw_exp_data = {}
        with open(filename, "r") as f:
            while True:
                line = f.readline()
                if not line: break
                if line.startswith("#"): continue
                if line.startswith("NAME"):
                    self.molecular.setName( line.lstrip("NAME").strip() )
                    continue
                if line.startswith("DIAS"):
                    items = line.strip().split()
                    self.molecular.setNameofDias(items[1], items[2])
                    continue
                if line.startswith("INDX"):
                    Labels = line.strip().split()[1:]
                    calcItems = filter( lambda line: line.islower(), Labels )
                    expItems  = filter( lambda line: line.isupper(), Labels )
                    continue
                if line.startswith("%BEGIN"):
                    _exchange_list = []
                    _label = []
                    _hybs  = []
                    for item in calcItems:
                        raw_calc_data[item] = []
                    for item in expItems:
                        raw_exp_data[item] = []
                    if not Labels: 
                        raise SyntaxError("the format of input file has error!")
                    if( len(line.split()) != 2 ): 
                        raise SyntaxError("the format of input file has error!")
                    dtype = line.split()[1]
                    while True:
                        line = f.readline()
                        if line.startswith("%END"): break
                        if line.startswith("EXCH"):
                            items = line.strip().split()[1:]
                            _exchange_list.append(items)
                            continue
                        items = line.strip().split()
                        _label.append(items[0])
                        for i, item in enumerate(Labels):
                            if item.isupper(): raw_exp_data[item].append(items[i+1])
                            elif item.islower(): raw_calc_data[item].append(items[i+1])
                            else: _hybs.append(items[i+1])
                    _nmrdata = NMRData(dtype)
                    _nmrdata.setLabel( _label )
                    self.molecular.setLabel(dtype, _label)
                    self.molecular.setHybs(dtype, _hybs)
                    avg = lambda x: 1.0*sum(x)/len(x)
                    str2float = lambda x: avg(map(float, x.split("-"))) if "-" in x else float(x)
                    for key in raw_calc_data:
                        _nmrdata.setCalcData(key, map(float, raw_calc_data[key]))
                    for key in raw_exp_data:
                        _nmrdata.setExpData(key, map(str2float, raw_exp_data[key]))
                    _nmrdata.setExchangeLabel(_exchange_list)
                    self.nmrdata.append(_nmrdata)
                    continue

    def calculateMae(self):
        _calculate_method = "Meam Average Error"
        for index, dtype, label, exp_data, calc_data in self.productRawData():
            mae_results = {}
            for exp in exp_data:
                for calc in calc_data:
                    scale_calc = scaledValue(calc_data[calc], exp_data[exp])
                    mae_results[exp+calc] = cm.calculateMae(
                          scale_calc, exp_data[exp])
            self._restoreResults(dtype, index, _calculate_method, mae_results)

    def calculateCC(self):
        _calculate_method = "correlation"
        for index, dtype, label, exp_data, calc_data in self.productRawData():
            cc_results = {}
            for exp in exp_data:
                for calc in calc_data:
                    scale_calc = scaledValue(calc_data[calc], exp_data[exp])
                    cc_results[exp+calc] = cm.calculateCC(
                        scale_calc, exp_data[exp])
            self._restoreResults(dtype, index, _calculate_method, cc_results)

    def _getStatParaofCP3(self, dtype):
        self.setMethod("6-31G")
        if dtype == "13C":
            uv_c = self.getStatPara("CP3", "correctC13", "n")
            uv_i = self.getStatPara("CP3", "incorrectC13", "n")
        if dtype == "1H":
            uv_c = self.getStatPara("CP3", "correctH1", "n")
            uv_i = self.getStatPara("CP3", "incorrectH1", "n")
        if dtype == "13C1H":
            uv_c = self.getStatPara("CP3", "correctC13H1", "n")
            uv_i = self.getStatPara("CP3", "incorrectC13H1", "n")
        return uv_c, uv_i

    def _restoreResults(self, dtype, index, method, value):
        if not self.calculate_results.has_key(dtype):
            self.calculate_results[dtype] = {}
        if not self.calculate_results[dtype].has_key(index):
            self.calculate_results[dtype][index] = {}
        self.calculate_results[dtype][index][method] = value
        
    def calculateCP3(self):
        _calculate_method = "CP3"
        for index, dtype, label, exp_data, calc_data in self.productRawData():
            cp3_results = {}
            uv_correct, uv_incorrect = self._getStatParaofCP3(dtype)
            for exp_label, exp_data in self.productPairData(exp_data,2):
                for calc_label, calc_data in self.productPairData(calc_data,2):
                    calc_data = map(scaledValue, calc_data, exp_data)
                    AaBb, AbBa = cm.calculateCP3(exp_data, calc_data, uv_correct, uv_incorrect)
                    cp3_results[exp_label+calc_label] = AaBb*100.
                    calc_label = calc_label[::-1]
                    cp3_results[exp_label+calc_label] = AbBa*100.

            self._restoreResults(dtype, index, _calculate_method, cp3_results)

    def _getStatParaofDP4(self, dtype):
        if dtype == "13C":
            usv = self.getStatPara("DP4", "C13")
        if dtype == "1H":
            usv = self.getStatPara("DP4", "H1")
        return usv

    def calculateDP4(self):
        _calculate_method = "DP4"
        dp4_results = {}
        for index, dtype, label, exp_data, calc_data in self.productRawData():
            usv = self._getStatParaofDP4(dtype)
            for exp in exp_data:
                results = []; keys = []
                for calc in calc_data:
                    scale_calc = scaledValue(calc_data[calc], exp_data[exp])
                    keys.append(exp+calc)
                    results.append(cm.calculateDP4(scale_calc,
                                 exp_data[exp], usv))
                results = [100.*x/sum(results) for x in results]
                dp4_results = dict(zip(keys, results))
            self._restoreResults(dtype, index, _calculate_method, dp4_results)

    def _getStatParaofDP4pwithSP2(self, dtype):
        if dtype == "13C":
            return self.getStatPara("DP4+", "nonscaledC13(sp2)")
        if dtype == "1H":
            return self.getStatPara("DP4+", "nonscaledH1(sp2)")
    def _getStatParaofDP4pwithSP3(self, dtype):
        if dtype == "13C":
            return self.getStatPara("DP4+", "nonscaledC13(sp3)")
        if dtype == "1H":
            return self.getStatPara("DP4+", "nonscaledH1(sp3)")
    def _getStatParaofsDP4(self, dtype):
        if dtype == "13C":
            return self.getStatPara("DP4+", "scaledC13")
        if dtype == "1H":
            return self.getStatPara("DP4+", "scaledH1")

    def calculateDP4p(self):
        _calculate_method = "DP4+"
        dp4p_results = {}
        udp4_results = {}

        for index, dtype, label, exp_data, calc_data in self.productRawData():
            if not self.molecular.hybs[dtype]: 
                raise TypeError("No hybs in %s, %s" % (dtype, index))
            spX = self.molecular.hybs[dtype]
            usv = self._getStatParaofsDP4(dtype)
            usv_sp = self._getStatParaofDP4pwithSP2(dtype)
            usv_sp3 = self._getStatParaofDP4pwithSP3(dtype)
            for exp in exp_data:
                dresults = []; uresults = []; keys = []
                for calc in calc_data:
                    keys.append(exp+calc)
                    udp4 = cm.calculateuTDP4(
                        calc_data[calc], exp_data[exp], spX, usv_sp, usv_sp3)
                    scale_calc = scaledValue(calc_data[calc], exp_data[exp])
                    sdp4 = cm.calculateDP4(scale_calc, exp_data[exp], usv)
                    uresults.append(udp4); dresults.append(udp4*sdp4)

                uresults = [100.*x/sum(uresults) for x in uresults]
                dresults = [100.*x/sum(dresults) for x in dresults]
                dp4p_results = dict(zip(keys, dresults))
                udp4_results = dict(zip(keys, uresults))
            self._restoreResults(dtype, index, _calculate_method, dp4p_results)
            self._restoreResults(dtype, index, "uDP4", udp4_results)

    def printNMR(self, label, calc, exp, dtype):
        calc_items = sorted(calc.keys())
        exp_items = sorted(exp.keys())
        data_list = calc_items + exp_items

        label_str = "\t".join([format(item, "<4") for item in ['',] + data_list])
        print("\t"+label_str)

        for i in range(len(label)):
            items = [label[i],] + [format(calc[item][i], ".2f") for item in calc_items ] + \
                [format(exp[item][i], ".2f") for item in exp_items]
            label_str = "\t".join([format(item, "<4") for item in items])
            print("\t"+label_str)

    def checkExchangeItem(self, raw_item ,item):
        exchange_items = []
        for x, y in zip(raw_item, item):
            if x == y: continue
            else: exchange_items.append("%s => %s" % ( x, y ))
        if not exchange_items: return "Raw data"
        return ", ".join( exchange_items )

    def productPairData(self, data, num):
        for idata in combinations(data, num):
            yield "".join(idata), [data[key] for key in idata]
       
    def productRawData(self):
        for _nmrdata in self.nmrdata:
            for i, item in enumerate(_nmrdata.all_exchange_list):
                temp_list = lambda y: [_nmrdata.label[x] for x in chain.from_iterable(y)]
                data_label = self.checkExchangeItem(temp_list(_nmrdata.all_exchange_list[0]), temp_list(item))
                temp_exp_data = _nmrdata.exchangeNMR(item)
                yield data_label, _nmrdata.dtype, _nmrdata.label, temp_exp_data, _nmrdata.calc_data

    def report(self):
        for index, dtype, label, exp_data, calc_data in self.productRawData():
            self.printNMR(label, calc_data, exp_data, dtype)

            results = self.calc(exp_data, calc_data, 
                    ['calculateDP4','calculateCP3', 'calculateCC', 'calculateMae', 'calculateDP4plus'],
                    dtype)
            new_res = {'calculateDP4plus':{}, 'calculateMae':{}, 'calculateCC':{}, 'calculateDP4':{}, 'calculateCP3':{}}

            for key in new_res:
                for item in results[key]:
                    cdp4_s = results[key][item]
                    cdp4_s = map(lambda x: (item+x, cdp4_s[x]),  cdp4_s)
                    new_res[key].update( dict(cdp4_s) )
            #print(pandas.DataFrame(new_res))
            print new_res

if __name__ == "__main__":
    import sys
    s = StatisticalNMR()
    filename = sys.argv[1]
    s.readDataFromFile(filename)
    s.calculateCC()
    s.calculateMae()
    s.calculateDP4()
    s.calculateDP4p()
    s.calculateCP3()
    print s.calculate_results
    #s.report()
