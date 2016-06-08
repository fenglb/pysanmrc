class StatisticalNMR:
    def setLabel(self, label):
        self.label = label
    def setExpData(self, value):
        self.exp_data = value
    def setCalcData(self, value):
        self.calc_data = value
    def readDataFromFile(self, filename):
        with open(filename, "r") as f:
            lines = f.readlines()
    def report(self):
