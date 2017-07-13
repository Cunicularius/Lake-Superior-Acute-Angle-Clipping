from __future__ import division
from copy import deepcopy
import random as r
import math as m
import numpy as np
import matplotlib.pyplot as plt
import copy
import sys


########################################
# General
########################################


class PolyNodeData:
    def __init__(self):
        self.PCC = []
        self.polychains = [[]]
        self.nodes = []
        self.NC = None
        

    def loadPND(self, polyname, nodename):
        with open(polyname) as p:
            polylines = p.readlines()
        

        ###Main Polychain###
        self.PCC.append(int(polylines[1].split()[0])) #PCC : PolyChainCount : A list of integers
        #Integers describe number of edges in polychains, with first being main and subsequent
        #being holes

        self.polychains = [[]]
        idx = 0
        for j in range(2, 2+self.PCC[0]):
            self.polychains[0].append(int(polylines[j].split()[1]))


        if len(polylines) == self.PCC[0] + 2:
            return 


        ###Hole Polychains###      
        zer0 = self.PCC[0]+3
        self.PCC += [None for i in range(int(polylines[zer0-1].split()[0]))] #Number of Holes
        for j in range(0, len(self.PCC)-1):

            p0s = sum(self.PCC[1:j+1]) + j
            if zer0+p0s == len(polylines):
                break
            self.PCC[j+1] = int(polylines[zer0 + p0s].split()[0])

            self.polychains.append([])
            for k in range(zer0 + p0s + 1, zer0 + p0s + 1 + self.PCC[j+1]):
                self.polychains[j+1].append(int(polylines[k].split()[1]))

        
        ###Node Data###
        with open(nodename) as p:
            nodelines = p.readlines()


        for i in range(1, len(nodelines)):
            self.nodes.append([float(j) for j in nodelines[i].split()])
            
        self.NC = len(self.nodes)

        
    def exportPND(self, fileprefix):
        with open(fileprefix + ".node", "w") as f:
            f.write("{} 2 0 1\n".format(self.NC))
            for i in self.nodes:
                f.write("{} {}\n".format(*[float(j) for j in i]))

        with open(fileprefix + ".poly", "w") as f:
            f.write("0 2 0 1\n")
            f.write("{} 1\n".format(int(self.PCC[0])))

            CC = self.polychains[0] #CurrentChain
            for i in range(1, len(self.polychains[0])):
                f.write("{1} {0} {1}\n".format(CC[i-1], CC[i]))
            f.write("TEST\n")
            f.write("{1} {0} {2}\n".format(CC[i], CC[i]+1, CC[0]))

            f.write("{}\n".format(int(len(self.PCC)-1)))
            
            for j in range(1, len(self.PCC)):
                f.write("{} 1\n".format(self.PCC[j]))
                CC = self.polychains[j] #CurrentChain
                for i in range(1, len(CC)):
                    f.write("{1} {0} {1}\n".format(CC[i-1], CC[i]))
                f.write("TEST\n")
                f.write("{1} {0} {2}\n".format(CC[i], CC[i]+1, CC[0]))
                
        return 0

    def loadOFF(self, filename):
        with open(filename) as f:
            lines = f.readlines()
        
        header = [int(i) for i in lines[1].split()]
        self.NC = header[0]
        NE = header[1]

        self.nodes = []
        for j in range(2, 2+self.NC):
            self.nodes.append([float(i) for i in lines[j].split()[0:2]])
        
        self.polychains = []
        for j in range(2+self.NC, 2+self.NC+NE):
            line = lines[j].split()
            self.polychains.append([int(i) for i in line[1:]])
            self.PCC.append(int(line[0]))

    def exportOFF(self, filename):
            exportToOFF(self.nodes, self.polychains, filename)
                    
    def rechain(self):
        self.polychains = []
        for idx, item in enumerate(self.PCC):
            self.polychains.append([sum(self.PCC[0:idx]) + j for j in range(item)])

    def manual(self, nodeslist): #nodeslist being list of lists, lists being list of nodes repr chains
        for i in nodeslist:
            self.PCC.append(len(i))
            self.nodes += i
        self.NC = len(nodeslist)
        self.rechain()
    
    def nearestNode(self, tIdx): #target [Node] Index
        sD = float("inf") #shortest Distance
        sIdx = None #shortest [Node by] Index
        for cIdx, cNode in enumerate(self.nodes): #current Index current Node
            if cIdx == tIdx:
                continue
            D = dist(self.nodes[cIdx], self.nodes[tIdx])
            if D < sD:
                sD = D
                sIdx = cIdx
        return([sD, sIdx])

    def clipAcute(self, C=0.5):
        aL = identifyAcute(self)
        rO = 0 #rolling Offset
        for idx, angles in enumerate(aL):
            for A in angles:
                D = self.nearestNode(A)[0]
                trpl = [self.nodes[(A-i+rO)%sum(self.PCC[0:idx+1])] for i in range(-1,2)]
                quad = clipTrpl(trpl, D, C)
                print("clipTrpl({}, {}, {}) = \n\t{}\n".format(trpl, D, C, quad))
                self.nodes[A+rO] = quad[1]
                self.nodes.insert(A+rO, quad[2])
                self.PCC[idx]+=1
                rO+=1
        self.rechain()
            

def clipTrpl(trpl, D, C=0.5):
    trpl = [np.array(i) for i in deepcopy(trpl)]
    v = [np.array([i[1],i[3]]) for i in vectorFormat(trpl)]
    trpl1A = trpl[1] + C*D/dist(trpl[0], trpl[1])*v[1]
    trpl1B = trpl[1] + C*D/dist(trpl[1], trpl[2])*v[0]
    quad = [trpl[0], trpl1A, trpl1B, trpl[2]]
    return [i.tolist() for i in quad]
    
    
            

def printAngles(PND):
    angleList = [[] for i in range(len(PND.PCC))]
    
    for i, chain in enumerate(PND.polychains):
        for j in range(len(chain)):
            trpl = [PND.nodes[j] for j in (chain[j:]+chain[:j])[0:3]]
            angleList[i].append(angleDebug(trpl))
    
    with open("angles_lsup.txt", "w") as f:
        for chain in angleList:
            for a in chain:
                f.write(str(180 - (a/(np.pi) * 180)) + "\n")
            f.write("\n")
    return angleList

def printAngleIndices(PND):
    acuteList = identifyAcute(PND)
    
    with open("anglesIndices_lsup.txt", "w") as f:
        for chain in acuteList:
            for a in chain:
                f.write(str(a) + "\n")
            f.write("\n")
    return acuteList

def identifyAcute(PND):
    acuteList = [[] for i in range(len(PND.PCC))]

    for i, chain in enumerate(PND.polychains):
        for j in range(len(chain)):
            trpl = [PND.nodes[j] for j in (chain[j:]+chain[:j])[0:3]]
            if isAcute(trpl):
                S = sum([PND.PCC[i] for i in range(i)])
                acuteList[i].append((j+1+S)%PND.NC)
            #acuteList[i].append(angleDebug(trpl))

    return acuteList

def averageDist(PND):
    S = 0 #Shortest
    for idx, item in enumerate(PND.nodes):
        for jdx, jtem in enumerate(PND.nodes):
            if idx == jdx:
                continue
            D = dist(item, jtem)
            S += D
    return S/PND.NC
            
def averageEdgeDist(PND):
    S = [] #Shortest
    for polychain in PND.polychains:
        for edgeIdx in polychain:
            D = dist(PND.nodes[edgeIdx[0]],PND.nodes[edgeIdx[1]])
            S.append(D)
    return sum(S)/len(S)

def longestDist(PND):
    S = 0 #Shortest
    for idx, item in enumerate(PND.nodes):
        for jdx, jtem in enumerate(PND.nodes):
            D = dist(item, jtem)
            if D > S:
                S = D
    return S
            
def dist(P1, P2):
    return ( (P2[0] - P1[0])**2 + (P2[1] - P1[1])**2 )**.5
    
        

def angleDebug(trpl):
    #LAST FOCUS
    """Takes list len=3 of nodes and returns True if acute angle, false otherwise"""
    vd = [np.array(i) for i in vectorFormat(trpl)]
    return angle_between(*vd)
           
def isAcute(trpl):
    """Takes list len=3 of nodes and returns True if acute angle, false otherwise"""
    vd = vectorFormat(trpl)
    if angle_between(*vd) < np.pi/2:
        return True
    else:
        return False

def vectorFormat(trpl):
    """trpl: [[p0x, p0y], [p1x, p1y], [p2x, p2y]]"""
    v = [[trpl[1][0], trpl[2][0], trpl[1][1], trpl[2][1]], [trpl[1][0], trpl[0][0], trpl[1][1], trpl[0][1]]]
    vd = [dispOrigin(vi) for vi in v]
    return vd
     



def robinsonAspectRatio(QVL):
    """ QuadVertexList : [[X0,Y0],...,[X3,Y3]] """
    bisectors = [ [.5*(QVL[i][0]+QVL[(i+1)%4][0]),.5*(QVL[i][1]+QVL[(i+1)%4][1])] for i in range(4) ]
    
    centroid = [0,0]
    for i in range(4):
        centroid[0] += .25*QVL[i][0]
        centroid[1] += .25*QVL[i][1]

    r1_h1 = ((centroid[0] - bisectors[0][0])**2 + (centroid[1] - bisectors[0][1])**2)**.5
    r2_h1 = ((centroid[0] - bisectors[1][0])**2 + (centroid[1] - bisectors[1][1])**2)**.5

    n1 = ((bisectors[0][0] - centroid[0]) / r1_h1, (bisectors[0][1] - centroid[1]) / r1_h1)
    n2 = ((bisectors[1][0] - centroid[0]) / r2_h1, (bisectors[1][1] - centroid[1]) / r2_h1)

    theta = 180/m.pi * m.acos(n1[0] * n2[0] + n1[1] * n2[1])

    sin = m.sin(theta * m.pi / 180)
    
    r1_h2 = sin * r2_h1
    r2_h2 = sin * r1_h1    

    return max(max(r1_h1,r1_h2)/min(r1_h1,r1_h2), max(r2_h1,r2_h2)/min(r2_h1,r2_h2))


def exportToOFFDebug( vertices, filename ):
    try:
        f = open(filename, "w")
        f.write("OFF\n")
        f.write("{} {} 0\n".format(str(len(vertices)), str( len(edgelists) )))
        for i, vertex in enumerate(vertices):
            f.write("{} {} 0.0\n".format(float(vertex[0]), float(vertex[1])))
        for i, vertex in enumerate(vertices):
            f.write("{}".format(len(edges)))
            f.write((" {}"*len(edges)).format(*edges))
            f.write("\n")
    except:
        print("exportQuadsToOFF(): FileError\n")
        return 1
    finally:
        f.close()
        return 0

def exportToOFF( vertices, edgelists, filename ):
    try:
        f = open(filename, "w")
        f.write("OFF\n")
        f.write("{} {} 0\n".format(str(len(vertices)), str( len(edgelists) )))
        for i, vertex in enumerate(vertices):
            f.write("{} {} 0.0\n".format(float(vertex[0]), float(vertex[1])))
        for edges in edgelists:
            f.write("{}".format(len(edges)))
            f.write((" {}"*len(edges)).format(*edges))
            f.write("\n")
    except:
        print("exportQuadsToOFF(): FileError\n")
        return 1
    finally:
        f.close()
        return 0


def exportQuadsToOFF( quads, filename ):
    """ quads is a list of QVL's such as those given by unitQuad and assessed by robinsonAspectRatio : quads = [ [[X0,Y0],...,[X3,Y3]], ... , [[U0,V0],...,[U3,V3]] ] """
    try:
        f = open(filename, "w")
        f.write("OFF\n")
        f.write("{} {} 0\n".format(str(4*len(quads)), str(len(quads))))
        for i, quad in enumerate(quads):
            for vertex in quad:
                f.write("{} {} 0.0\n".format(float(vertex[0] + 3*i), float(vertex[1])))
        for i in range(len(quads)):
            f.write("4 {} {} {} {}\n".format(*[4*i + j for j in range(4)]))
    except:
        print("exportQuadsToOFF(): FileError\n")
        return 1
    finally:
        f.close()
        return 0



def convertOFFtoELENODE( offname ):
    """Converts an off in the same directory to ele node files"""
    with open(offname, "r") as OFF:
        OFFLines = OFF.readlines()

    OFFData = []
    for line in OFFLines:
        OFFData.append(line.split())
        
    numVertices = int(OFFData[1][0])
    numFaces = int(OFFData[1][1])
    numPerFace = int(OFFData[2+numVertices+1][0])

    outname = offname.split(".")[0] #To name the output files

    with open( outname + ".ele", "w") as ELE:
        ELE.write( "{}\t{}\t0\n".format(numFaces, numPerFace)) #Placing the number of elements, and taking the number of vertices in an element from the first element that appears in the off
        
        for i in range(2 + numVertices, 2 + numVertices + numFaces):
            temp = []
            for j in range( 1, 1+numPerFace):
                temp.append( int(OFFData[i][j]) + 1 )

            template = "{}\t" + "{}\t"*numPerFace + "\n"
            ELE.write( template.format( i-numVertices-1, *temp))

    with open( outname + ".node", "w") as NODE:
        NODE.write( "{}\t2\t0\t0\n".format(numVertices)) #Placing the number of elements, and taking the number of vertices in an element from the first element that appears in the off
        
        for i in range(2, 2 + numVertices):

            template = "{}\t{}\t{}\n"
            NODE.write( template.format( i-1, *OFFData[i]))
    
    return



def listFilenameFormat( L ):
    return str(L).strip("[]").replace(",","_").replace(" ","").replace("\'", "")
 


#############################
#Perimarea
#############################



def PerimareaRatio(QVL):
    scale(QVL, perimeter(QVL)**-1)
    return 16 * quadArea(QVL) / perimeter(QVL)

def perimeter(PVL):
    P = 0
    n = len(PVL)
    for i in range(n):
        P += ( (PVL[(i)%n][0] - PVL[(i+1)%n][0])**2 + (PVL[(i)%n][1] - PVL[(i+1)%n][1])**2 )**.5
    return P

def quadArea(QVL):
    area = 0
    for i in [0,2]:
        TVL = []
        for j in range(3): #size of Triangle
            TVL.append(QVL[(i+j)%4])
        area += abs(.5 * (( (TVL[1][0] - TVL[0][0]) * (TVL[2][1] - TVL[0][1]) ) - ( (TVL[2][0] - TVL[0][0]) * (TVL[1][1] - TVL[0][1]) )) )
    return area


def scale(PVL, S):
    for i, Point in enumerate(PVL): #len(PVL)
        for j, Coord in enumerate(Point):
            PVL[i][j] = S*Coord
    return PVL


########################################
# CompC
########################################

def Area2( a, b, c ):
    X, Y = 0, 1
    return \
            (b[X] - a[X]) * (c[Y] - a[Y]) - \
            (c[X] - a[X]) * (b[Y] - a[Y])

def Left( a, b, c ):
    return Area2( a, b, c ) > 0

########################################
# Angle Problem
########################################



def randomConvexQuad( Min = 1, Max = 179 ):
    """ Random Convex Quad
        Produces the angles for a valid, convex quad by 'splitting' 2 180's
        takes the components of the splits and there it is 
    """
    ang = []
    for i in range( 2 ):
        Slice = r.randint( Min, Max )
        ang.extend( [Slice, 180-Slice] )
    ang[1], ang[2] = ang[2], ang[1]
    return ang



def SPAlt( ang ): 
    """
    Returns list ang cycled such that the ang[2] is the smallest angle
    """
    indexMin = 0    
    itemMin = 360
    for i, item in enumerate( ang ):
        if (item < itemMin):
            indexMin = i
            itemMin = item
    return ( ang[(indexMin-2)%4:] + ang[:(indexMin-2)%4] )



def intersection(v1, v2):
    """
    Returns the intersection of v1 and v2. Note however that these are not 'bezier lines', x1, y1, x3, and y3
    are all *changes* in x, not describing a point. So, rather than x0 + (x1 - x0)*t, its just x0 + x1*t.
    It just made the algebra slightly easier.
    list v1 = [x0, x1, y0, y1], list v2 = [x2, x3, y2, y3]
    """
    x = v1[0:2] + v2[0:2]
    y = v1[2:4] + v2[2:4]
    if( x[3] == 0 ): #To avoid a divide by zero, if x[3] is 0 then we just solve for where lineA equals x[2]
        t1 =    (x[2] - x[0])/\
                (x[1])
        return [ v1[0] + v1[1]*t1, v1[2] + v1[3]*t1 ]

    else: 
        t1 =    ( y[0] - y[2] + (y[3]/x[3])*(x[2] - x[0]) )/\
                ( (y[3]*x[1])/x[3] - y[1] )
        return [ v1[0] + v1[1]*t1, v1[2] + v1[3]*t1 ]


def unitQuad(angOriginal, offset = 1):
    """ Returns list of points that describe a unit quad based on ang, which contains angles
        describing a valid convex quad
        For example, [90,90,90,90] returns [Point(0,0),Point(1,0),Point(1,1),Point(0,1)]
    """
    ang = angOriginal.copy()

    for i, item in enumerate(ang):
        ang[i] = m.radians(item)

    points = []
    points.append( [0,0] )#P0
    points.append( [offset,0] )#P1
    points.append( None )#P2 
    points.append( [m.cos(ang[0]), m.sin(ang[0])] )#P3

    # !!! listA = [x0, x1, y0, y1], listB = [x2, x3, y2, y3]
    lineA = [  points[3][0], m.cos( ang[3] - (m.pi - ang[0]) ), points[3][1], m.sin( ang[3] - (m.pi - ang[0]) )  ]
    lineB = [  points[1][0], m.cos( m.pi - ang[1] ), points[1][1], m.sin( m.pi - ang[1] )  ]
    points[2] = (  intersection( lineA, lineB )  )

    return points



########################################
# Edge Problem
########################################



def randomEdgeLengths(): 
    """
    Uses a convoluted method of generating a quad with the Angle code and then taking the edge lengths from that. 
    """

    angles = randomConvexQuad()                     #Random angles
    quad = unitQuad(angles, r.uniform(0.05, 0.95))  #For a random quad, with varied base edge length

    edges = []
    n = len(quad)
    for i in range(n):  #For each edge in the quad, record the edge's length
        edges.append(( (quad[(i)%n][0] - quad[(i+1)%n][0])**2 + (quad[(i)%n][1] - quad[(i+1)%n][1])**2 )**.5) 
   
    maxItem = max(edges)
    maxIndex = edges.index(maxItem) 

    edges = edges[maxIndex:] + edges[:maxIndex] #Reorder the edge lengths so that the greatest one comes first
    for i, item in enumerate(edges):            #And normalize them so that the greatest edge has length 1
        edges[i] = item / maxItem
 
    return edges



def circleIntersection( p0, r0, p1, r1 ):
    circle0 = (*p0, r0) # Reformat the data for our 2 circles to
    circle1 = (*p1, r1) # (xi, yi, ri)

    #Send it out to CCIPP.py's circle intersection function
    intersections = Geometry.circle_intersection(None, circle0, circle1) 
    
    if len(intersections) != 2: #There's 1 or no intersections, in either case the excetpion is handled later
        raise ValueError

    elif intersections[0][1] > intersections[1][1]: #Return the point of intersection with greater y value
        return intersections[0]
    else:
        return intersections[1]
    

    
class Geometry(object):
    def circle_intersection(self, circle1, circle2):
        '''
        Source: https://gist.github.com/xaedes/974535e71009fa8f090e
        @summary: calculates intersection points of two circles
        @param circle1: tuple(x,y,radius)
        @param circle2: tuple(x,y,radius)
        @result: tuple of intersection points (which are (x,y) tuple)
        '''
        # return self.circle_intersection_sympy(circle1,circle2)
        d2r = m.pi/180
        x1,y1,r1 = circle1
        x2,y2,r2 = circle2
        # http://stackoverflow.com/a/3349134/798588
        dx,dy = x2-x1,y2-y1
        d = m.sqrt(dx*dx+dy*dy)
        if d > r1+r2:
            print("#1")
            return None # no solutions, the circles are separate
        if d < abs(r1-r2):
            print("#2")
            return None # no solutions because one circle is contained within the other
        if d == 0 and r1 == r2:
            print("#3")
            return None # circles are coincident and there are an infinite number of solutions

        a = (r1*r1-r2*r2+d*d)/(2*d)
        h = m.sqrt(r1*r1-a*a)
        xm = x1 + a*dx/d
        ym = y1 + a*dy/d
        xs1 = xm + h*dy/d
        xs2 = xm - h*dy/d
        ys1 = ym - h*dx/d
        ys2 = ym + h*dx/d

        return (xs1,ys1),(xs2,ys2)



def unitCircPt(theta):
    return np.array([np.cos(theta),np.sin(theta)])



def angle_between(v2, v1):
    """
    Note however that these are not 'bezier lines', x1, y1, x3, and y3
    are all *changes* in x, not describing a point. So, rather than x0 + (x1 - x0)*t, its just x0 + x1*t.
    It just made the algebra slightly easier.
    list v1 = [x0, x1, y0, y1], list v2 = [x2, x3, y2, y3]
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    result = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    if np.isnan(result):
        if abs(v1_u + v2_u) < .5 * (abs(v1_u) + abs(v2_u)):
            return np.pi
        else:
            return 0.0
    if Left( [v2[1],v2[3]], [0,0], [v1[1],v1[3]] ):
        return 2*np.pi - result
    return result



def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def dispOrigin(vector):
    vector[1] -= vector[0]
    vector[3] -= vector[2]
    vector[0] = vector[2] = 0
    return vector


def unitQuad_Edge(lens, N=3):
    """
    lens, list of edge lengths
    N is the number of desired iterations between the left and right degenerate cases, leftDegenerate and rightDegenerate
    
    """
    template = [ np.array([0,0]), np.array([lens[0], 0]), None, None ] #Template from which to generate other Quad Vertex Lists
    leftDegenerate = template.copy()    #Left Limit of quad if you were to rotate edge 3 CCW about the origin until you no longer can
    rightDegenerate = template.copy()   #Right Limit of quad if you were to rotate edge 2 CW about point 1 until you no longer can,
                                        #   or alternatively, how far edge 3 can rotate CW until the quad is degenerate
    try:
        leftDegenerate[3] = np.array( circleIntersection(leftDegenerate[0], lens[3], leftDegenerate[1], lens[1]+lens[2]) )
        leftDegenerate[2] = ( lens[1] / (lens[2]+lens[1]) ) * (leftDegenerate[3]-leftDegenerate[1]) + leftDegenerate[1]
    except: 
        leftDegenerate[3] = np.array([-lens[3],0])
        leftDegenerate[2] = np.array( circleIntersection(leftDegenerate[3], lens[2], leftDegenerate[1], lens[1]) )

    try:
        rightDegenerate[2] = np.array( circleIntersection(rightDegenerate[0], lens[2]+lens[3], rightDegenerate[1], lens[1]) )
        rightDegenerate[3] = ( lens[3] / (lens[3]+lens[2]) ) * rightDegenerate[2]
    except:
        rightDegenerate[2] = np.array([lens[0]+lens[1], 0])
        rightDegenerate[3] = np.array( circleIntersection(rightDegenerate[0], lens[3], rightDegenerate[2], lens[2]))
        
    rightOfOrigin = np.array([1,0]) #Theta = 0 on the Unit Circle
    thetaMin = angle_between(leftDegenerate[3], rightOfOrigin) #Angle of 
    thetaMax = angle_between(rightDegenerate[3], rightOfOrigin)
    pitch = (thetaMax - thetaMin) / (N-1)

    result = []
    result.append(leftDegenerate) 
    for i in range(1, N-1):
        result.append(template.copy())
        result[i][3] = lens[3]*unitCircPt(i*pitch+thetaMin)
        result[i][2] = np.array(circleIntersection( result[i][3], lens[2], result[i][1], lens[1]))
    result.append(rightDegenerate) 

    return listify(result)



def listify(result):
    for quad in result:
        for i, point in enumerate(quad):
            quad[i] = list(point)
    return(result)



########################################
#Main
########################################



if __name__ == "__main__" and len(sys.argv) == 1:
    print("usage")
   
#if __name__ == "__main__" and sys.argv[2] == "E":
#    templates = []
#
#    if len(sys.argv) == 3:
#        with open(sys.argv[2], 'r') as f:
#            lines = f.readlines()
#
#        for i, list_i in lines:
#            lines[i] = list_i.split()
#            for 

    count = 10
    f = randomEdgeLengths
#    count = 1
#    f = [1, 0.9, 0.258, 0.242]

    for i in range(count):
        #template = [1, 1, 1, 1]
        N = 18
        template = f()
        quads = unitQuad_Edge(template, N)
        print("Success with {}".format(str(template)))

        tCopy = template.copy()
        for i, item in enumerate(tCopy):
            tCopy[i] = "{0:0>4}".format(int(item*1000))
        offname = listFilenameFormat(tCopy)
        print("Created {}.off\n\n".format(offname))
        
        RARs = []
        xlist = [j/N for j in range(1, N+1)]

        minRAR = float("inf") 
        minQuad = None
        for quad in quads:
            RARs.append(robinsonAspectRatio(quad))
            #RARs.append(robinsonAspectRatio(quad))
            if RARs[-1] < minRAR:
                minRAR = RARs[-1]
                minQuad = quad
        
        print(xlist)
        print(RARs)
        plt.plot(xlist, RARs)
        plt.axis([0,1,0, 10])
        plt.title(str(template))
        plt.savefig(offname)
        plt.clf() 
        output = [quads[0], minQuad, quads[-1]]

        exportQuadsToOFF( output, offname )

        print("\n"*5)

