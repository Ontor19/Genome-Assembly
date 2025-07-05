from collections import defaultdict,deque
from Bio import SeqIO

i=int(input("Enter k-mar size:"))
e=[]
for x in SeqIO.parse("example.fasta", "fasta"):
        a=str(x.seq)
        
        for j in range(len(a)-i+1):
                k_mar=(a[j:j+i])
                e.append(k_mar)  

print(e)   
print("")
print(" ")
           
# Build de Bruijn Graph
def dg(kmers):
    graph = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    return graph

#Create and print the graph
graph = dg(e)

print("\nDe Bruijn Graph:")
for node in graph:
    print(f"{node} -> {', '.join(graph[node])}")




                                
                
        
        