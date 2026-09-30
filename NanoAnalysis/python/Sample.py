class Sample:
    def  __init__(self, dbblock, origin, year):

        # fixme check on origin
        
        self.year    = year
        self.origin  = origin
        self.name    = dbblock.name
        self.process = dbblock.process
        self.setLuminosity()
        #fixme: to be expandend
        

    def isMC(self):
        return self.origin == 'MC'

    def __str__(self):
        return f"origin={self.origin}, year={self.year}, name={self.name}, process={self.process}"

    def __repr__(self):
        return self.__str__()

    
    def path(self):
        return f"samples/{self.year}/{self.name}.root"


    def setLuminosity(self):
        
        if   self.year == 2026:
            self.luminosity = 28.06
        elif self.year == 2025:
            self.luminosity = 110.63
        elif self.year == 2024:
            self.luminosity = 109.95
        elif self.year == 2023:
            self.luminosity = 28.28
        elif self.year == 2022:
            self.luminosity = 34.76
        elif self.year == 2018:
            self.luminosity = 59.56
        elif self.year == 2017:
            self.luminosity = 42.07
        elif self.year == 2016:
            self.luminosity = 36.31
        else:
            self.luminosity = -1.

    def explainYourself(self):
        return f"sample: {self.name}, process tag: {self.process}, year: {self.year}, corresponding luminosity: {self.luminosity}/fb"

    def print(self):
        print(self.explainYourself)
            
            
            
            
            
            
