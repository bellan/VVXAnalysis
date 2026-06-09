class Sample:
    def  __init__(self, dbblock, origin, year):

        # fixme check on origin
        
        self.year    = year
        self.origin  = origin
        self.name    = dbblock.name
        self.process = dbblock.process
        #fixme: to be expandend
        

    def isMC(self):
        return self.origin == 'MC'

    def __str__(self):
        return f"origin={self.origin}, year={self.year}, name={self.name}, process={self.process}"

    def __repr__(self):
        return self.__str__()

    
    def path(self):
        return f"samples/{self.year}/{self.name}.root"
