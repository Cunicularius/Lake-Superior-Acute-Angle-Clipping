from InitProb import *
from collections import defaultdict

PNFile = "lsup"

OVERRIDE = 1

vertices = [[1, 1], [0, 3], [-1, 1], [-1, -1], [1, -1], [0.5, 0], [0, 0.75], [-0.5, 0], [0, 1]] 
edgelists = [[0,1,2,3,4], [5,6,7,8]]
trpl = [[1, 1], [0, 3], [-1, 1]]


c = .1
box =  [np.array([  0.08*c, 0.08*c]),\
        np.array([ -0.08*c, 0.08*c]),\
        np.array([ -0.08*c,-0.08*c]),\
        np.array([  0.08*c,-0.08*c])]

#def clip(PND, aL):
    
    

#Instead of modifying the polychain data, just reproduce it based on PND.PCC

if __name__ == "__main__" or OVERRIDE:
    OFF = PolyNodeData()
    OFF.loadOFF("lsup.off")
    acuteAngles = [] #This will be a list of lists that we fill with vertex idx for each polychain
    
    PND = PolyNodeData()
    PND.loadPND(PNFile + ".poly", PNFile + ".node")

    aL = identifyAcute(PND)

    
#    L = aL[0]
#    A = aL[-3][0]
#    out = []
#
#    for aI in L:
#        for b in box:
#            out.append((PND.nodes[aI]+b).tolist())
#
#    for b in box:
#        out.append((PND.nodes[A]+b).tolist())
#        
#
#    OFFOut = OFFData()
#
#    OFFOut.vertices = out
#    OFFOut.elements = [[i*4 + j + 518 for j in range(4)] for i in range(len(OFFOut.vertices)//4)]
#    #OFFOut.elements = [[i*4 + j for j in range(4)] for i in range(len(OFFOut.vertices)//4)]
#    OFFOut.NV = len(out)
#    OFFOut.NE = len(OFFOut.elements)
#
#    OFFOut.export("marks.off")

    

#   exportToOFF(vertices, edgelists, "House.off") 



    
#    #VVV All meant to determine boundary from a tri/quadrangulation
#    
#    D = defaultdict(lambda: 0)
#    OFF = importOFF("superior.off")
#    
#    #Determine Boundary and Holes
#
#    for ele in OFF.elements:
#        for idx in range(len(ele)):
#            D[str({ele[idx], ele[(idx+1) % len(ele)]})] += 1
#            
#    unsortedBoundary = set()
#
#    for edge in D:
#        if D[edge] == 1:
#            unsortedBoundary.add(edge)

