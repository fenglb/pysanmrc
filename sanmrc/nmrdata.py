from itertools import permutations, product, chain
class NMRData:
    def __init__(self, dtype=None):
        self.dtype = dtype
        self.calc_data = {}
        self.exp_data = {}
        self.label = []
        self.all_exchange_list = []
    def setDataType(self, dtype):
        self.dtype = dtype
    def setLabel(self, label):
        self.label = label
    def setExpData(self, index, value):
        self.exp_data[index] = value
    def setCalcData(self, index, value):
        self.calc_data[index] = value
    def setExchangeLabel(self, exchange_list ):
        if not self.label: raise ValueError("don't exist labe!")
        # if exchange_list is none, then all_exchange_list = [()]
        self.all_exchange_list = list(product(*(
                                permutations(self.label.index(x) for x in item)
                                for item in exchange_list)))
    def exchangeNMR(self, exchange_item):
        if not exchange_item: return self.exp_data
        temp_exp_data = {}
        for key in self.exp_data:
            temp_exp_data[key] = self.exp_data[key][:]
            for x, y in zip(chain.from_iterable(self.all_exchange_list[0]), chain.from_iterable(exchange_item)):
                temp_exp_data[key][x] = self.exp_data[key][y]
        return temp_exp_data

class Molecular:
    def __init__(self, name=None):
        self.isomers = []
        self.name = name
        self.hybs = {}
        self.label = {}
    def setName(self, name):
        self.name = name
    def setHybs(self, key, value):
        self.hybs[key] = value
    def setLabel(self, key, value):
        self.label[key] = value
    def setNameofDias(self, index, name):
        self.isomers.append( {name: index} )
