class EpiSolver:
    def sensitivity(self, A, C):
        AC = A + C
        return (A / AC) 
    def specificity(self, D, B):
        BD = B + D
        return (D / BD)
    def PVP(self, A, B):
        AB = A + B
        return (A / AB)

test = EpiSolver()

# print (test.sensitivity(380, 20))


test1 = 34
test2 = 325

testlist = [test1, test2]

for i in testlist:
    format