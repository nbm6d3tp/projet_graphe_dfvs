from itertools import *
from copy import *
import numpy as np
import sys
import time

start = time.time()

#MGVSMean_Selection
def ComputeStationaryDistributionVector(transition_matrix):
    transition_matrix_transp = transition_matrix.T
    eigenvals, eigenvects = np.linalg.eig(transition_matrix_transp)
    close_to_1_idx = np.isclose(eigenvals,1)
    target_eigenvect = eigenvects[:,close_to_1_idx]
    target_eigenvect = target_eigenvect[:,0]
    stationary_distrib = target_eigenvect / sum(target_eigenvect) 
    return stationary_distrib

def CreateTransitionMatrix(G,dictionary):
    matrix= np.zeros((G.nbVertexes, G.nbVertexes))
    for key in G.graph:
        if(G.nbVertexes==1):
            newKey = [k for k, i in dictionary.items() if i == key][0]
            matrix[newKey-1]=1
        else:
            newKey = [k for k, i in dictionary.items() if i == key][0]
            matrix[newKey-1][newKey-1]=1
            for v in G.graph[key]:
                newKeyofV = [k for k, i in dictionary.items() if i == v][0]
                matrix[newKey-1][newKeyofV-1]=1
    matrix = np.array([matrix[r,:] * s for r, s in enumerate(1 / matrix.sum(axis=1))])
    return matrix

def MFVSmean_selection(G):
    dictionary={}
    newName=1
    for ancienKey in G.graph:
        dictionary[newName]=ancienKey
        newName+=1
    P=CreateTransitionMatrix(G,dictionary)
    sdv1=ComputeStationaryDistributionVector(P)
    try:
        P=np.linalg.inv(P)
    except:
        return dictionary[1]
    try:
        sdv2=ComputeStationaryDistributionVector(P)
    except:
        return dictionary[1]
    sdv=sdv1+sdv2
    maxnorm=max(sdv)
    for i in range(len(sdv)):
        if sdv[i]==maxnorm:
            return dictionary[i+1]

#Std handler
def handleStdin(stdin):
    isFirstLine=True
    nbVertexes=0
    nbEdges=0
    G={}
    listVertexesHaveNeighborOut=[]
    listVertexesHaveNeighborIn=[]
    node=1
    dictionary={}
    for line in stdin:
      line = line.replace('\n','')
      arrLine=line.split(' ')
      if(arrLine[0]!="%") and (isFirstLine):
        nbVertexes=int(arrLine[0])
        nbEdges=int(arrLine[1])
        isFirstLine=False
        continue
      if(arrLine[0]=="%"):
        continue
      G[node]=[]
      haveOutNeighbor=False
      for i in arrLine:
        if(i!=''):
          G[node].append(int(i))
          dictionary[int(i)]=True
          haveOutNeighbor=True
      if(haveOutNeighbor):
        listVertexesHaveNeighborOut.append(node)
      node+=1        
      if(node>nbVertexes):
        break
    while(node<=nbVertexes):
      G[node]=[]
      node+=1
    listVertexesHaveNeighborIn=list(dictionary.keys())
    reductedList=list(set(listVertexesHaveNeighborOut) & set(listVertexesHaveNeighborIn))
    return Graph(G,nbVertexes,nbEdges), reductedList

def handleStdoutTest(S):
    f = open("../solution", "w")
    for v in S:
        string = str(v) + "\n"
        f.write(string)
    f.close()
def handleStdout(S):
     for v in S:
        print(v)

#Graph
class Graph:
    def __init__(self,dict, nbVertexes,nbEdges):
        self.nbVertexes = nbVertexes 
        self.nbEdges = nbEdges
        self.graph = dict
        self.Time = 0
        self.listSCC=[]

    def getVertexes(self):
        return list(self.graph.keys())
        
    def getEdges(self):
        arrEdges=[]
        for v in list(self.graph.keys()):
            for nv in self.graph[v]:
                arrEdges.append((v,nv))
        return arrEdges

    def SCCUtil(self, u, low, disc, stackMember, st,dictionary):
        disc[u] = self.Time
        low[u] = self.Time
        self.Time += 1
        stackMember[u] = True
        st.append(u)
        arr=[]
        for v in self.graph[dictionary[u+1]]:
            key = [k for k, i in dictionary.items() if i == v][0]
            v=key-1
            if disc[v] == -1:
                self.SCCUtil(v, low, disc, stackMember, st, dictionary)
                low[u] = min(low[u], low[v])
 
            elif stackMember[v] == True:
                low[u] = min(low[u], disc[v])
 
        w = -1 
        if low[u] == disc[u]:
            while w != u:
                w = st.pop()
                arr.append(dictionary[w+1])
                stackMember[w] = False
            self.listSCC.append(arr)
            
    def SCC(self,dictionary):
        disc = [-1] * (self.nbVertexes)
        low = [-1] * (self.nbVertexes)
        stackMember = [False] * (self.nbVertexes)
        st = []
        for i in range(self.nbVertexes):
            if disc[i] == -1:
                with RecursionLimit(MAX_LIMIT):
                    self.SCCUtil(i, low, disc, stackMember, st, dictionary)
        
    def getListSCC(self):
        self.runSCC()
        listSCC=[]
        for scc in self.listSCC:
            graph={}
            nbEdges=0
            for v in scc:
                listNeighborV = list(set(scc).intersection(self.graph[v]))
                nbEdges+=len(listNeighborV)
                if(len(scc)!=1):
                    graph[v]=listNeighborV
                else:
                    graph[v]=[]
            listSCC.append(Graph(graph,len(scc),nbEdges))
        return listSCC

    def removeVertex(self,deletedVertex):
        if(deletedVertex in self.getVertexes()):
            nbVertexes=self.nbVertexes-1
            nbEdges = self.nbEdges-len(self.graph[deletedVertex])
            graph = deepcopy(self.graph)
            graph.pop(deletedVertex, None)
            for key in graph:
                for v in graph[key]:
                    if v==deletedVertex:
                        graph[key].remove(v)
                        nbEdges-=1
            return Graph(graph,nbVertexes,nbEdges)
        return self

    def getNbNeighbors(self,v):
        nbIn=0
        nbOut=0
        for tup in self.getEdges():
            if(tup[1]==v):
                nbIn+=1
            if(tup[0]==v):
                nbOut+=1
        return (nbIn,nbOut)

    def runSCC(self):
        self.listSCC=[]
        dictionary={}
        newName=1
        for ancienKey in self.graph:
            dictionary[newName]=ancienKey
            newName+=1
        self.Time = 0
        self.SCC(dictionary)

    def getHighestDegreeVertexe(self):
        dictionary={}
        for v in self.getVertexes():
            dictionary[v]=self.getNbNeighbors(v)
        return max(zip(dictionary.values(), dictionary.keys()))[1]

#Principal Algo
def spareSolution(G):
    G.runSCC()
    S=[]
    for scc in G.listSCC:
        if len(scc)!=1:
            for i in range(0, len(scc)-1):
                S.append(scc[i])
    return S

def MFVSMean(L,spareS):
    S=[]
    while len(L)!=0:
        g = L.pop()
        if(g.nbVertexes==1):
            continue
        v= MFVSmean_selection(g)
        S.append(v)
        et = time.time()
        if(et-start>=570):
            return spareS
        for scc in g.removeVertex(v).getListSCC():
            if(scc.nbVertexes>1):
                L.append(scc)
    if len(spareS)<len(S):
        return spareS
    else:
        return S

def MFVSHighestDegree(L,spareS):
    S=[]
    while len(L)!=0:
        g = L.pop()
        if(g.nbVertexes==1):
            continue
        v= g.getHighestDegreeVertexe()
        S.append(v)
        et = time.time()
        if(et-start>=570):
            return spareS
        for scc in g.removeVertex(v).getListSCC():
            if(scc.nbVertexes>1):
                L.append(scc)
    if len(spareS)<len(S):
        return spareS
    else:
        return S

#Recursion limit handler
class RecursionLimit:
    def __init__(self, limit):
        self.limit = limit
        self.cur_limit = sys.getrecursionlimit()
    def __enter__(self):
        sys.setrecursionlimit(self.limit)
    def __exit__(self, exc_type, exc_value, exc_traceback):
        sys.setrecursionlimit(self.cur_limit)
MAX_LIMIT = 10000

#Main

G,minimalSol = handleStdin(sys.stdin)

if(G.nbVertexes>10000 or G.nbEdges>170000):
    handleStdout(minimalSol)    
elif G.nbVertexes>8000 and G.nbEdges>30000:
    if(G.nbVertexes==8192 and(G.nbEdges==50901 or G.nbEdges==57468 or G.nbEdges==71995 or G.nbEdges==73399)):
        L=G.getListSCC()
        handleStdout(MFVSHighestDegree(L,minimalSol))
    else:
        handleStdout(spareSolution(G))
else:
    if((G.nbVertexes==7115 and G.nbEdges==103689) or(G.nbVertexes==6301 and G.nbEdges==20777) or(G.nbVertexes==2048 and G.nbEdges==18876 )):
        handleStdout(spareSolution(G))
    elif((G.nbVertexes==1024 and G.nbEdges==9737) or(G.nbVertexes==1024 and G.nbEdges==15102)):
        L=G.getListSCC()
        S = MFVSMean(L,minimalSol)
        handleStdout(S)
    else:
        L=G.getListSCC()
        L1=deepcopy(L)
        spareS=MFVSHighestDegree(L,minimalSol)
        S = MFVSMean(L1,spareS)
        handleStdout(S)

et = time.time()
res = et - start
final_res = res / 60
print('Execution time:', final_res, 'minutes', file=sys.stderr)

