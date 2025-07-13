import tkinter as tk
from tkinter import scrolledtext
from collections import defaultdict, deque
from Bio import SeqIO

k = 17
kmers = []
for record in SeqIO.parse("file.fasta", "fasta"):
    sequence = str(record.seq)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append(kmer)

for i in SeqIO.parse("file-2.fasta", "fasta"):
    arg = str(i.seq)

# Build Graph
def build_graph(kmers):
    graph = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    return graph

graph = build_graph(kmers)

# In-degree and Out-degree
inn = defaultdict(int)
out = defaultdict(int)
for node in graph:
    out[node] = len(graph[node])
    for neighbor in graph[node]:
        inn[neighbor] += 1

# Start Node
start_node = None
all_nodes = set(inn.keys()).union(set(out.keys()))
for node in all_nodes:
    if out[node] > inn[node]:
        start_node = node
        break
if not start_node:
    start_node = next(iter(graph))

# Eulerian Path
def find_path(graph, start):
    path = []
    stack = [start]
    #local_graph = {node: deque(neighbors) for node, neighbors in graph.items()}
    local_graph = {}
    for node in graph:
         local_graph[node] = deque(graph[node])

    while stack:
        current = stack[-1]
        if current in local_graph and local_graph[current]:
            next_node = local_graph[current].popleft()
            stack.append(next_node)
        else:
            path.append(stack.pop())
    return path[::-1]

path = find_path(graph, start_node)

# Contig Build
def build_contig(path):
    if not path:
        return ""
    contig = path[0]
    a=path[0]
    for node in path[1:]:
        if node != a:
            contig += node[-1]
        a=node
    return contig

contig_seq = build_contig(path)

# Tkinter GUI
root = tk.Tk()
root.title("Assembly Overlap & ARG Visualization")

arg_label = tk.Label(root, text=f"Pieces Of Contig:", font=("Arial", 15))
arg_label.pack(pady=5)

# Path Overlap Visualization (Staircase Style)
path_text = scrolledtext.ScrolledText(root, width=80, height=20, wrap=tk.NONE)
path_text.pack(pady=20)

# Staircase-style overlapping visualization
offset = 0
for i in range(len(path)):
    line = " " * offset + path[i] + "\n"
    path_text.insert(tk.END, line, "normal")
    if i < len(path)-1:  
        if path[i+1] != path[i]:
            offset += 1   
path_text.tag_config("normal", font=("Courier", 20))

arg_label = tk.Label(root, text=f"Scaffold:", font=("Arial", 15))
arg_label.pack(pady=5)

# Contig Visualization with ARG match
contig_text = scrolledtext.ScrolledText(root, width=80, height=5, wrap=tk.WORD)
contig_text.pack(pady=10)
contig_text.tag_config("contig", font=("Courier", 20))
start_idx = contig_seq.find(arg)
if start_idx != -1:
    contig_text.insert(tk.END, contig_seq[:start_idx],"contig")
    contig_text.insert(tk.END, contig_seq[start_idx:start_idx+len(arg)], "match")
    contig_text.insert(tk.END, contig_seq[start_idx+len(arg):],"contig")
else:
    contig_text.insert(tk.END, contig_seq,"contig")


contig_text.tag_config("match", background="yellow", foreground="blue", font=("Courier", 20, "bold"))
contig_text.tag_config("contig",  font=("Courier", 20, "bold"))


# ARG Show
arg_label = tk.Label(root, text=f"ARG Sequence:\n{arg}", font=("Arial", 20))
arg_label.pack(pady=20)

# Result
result = tk.Label(root, text="ARG Found in Scaffold" if start_idx != -1 else "ARG Not Found", fg="green" if start_idx != -1 else "red", font=("Arial", 20, "bold"))
result.pack(pady=20)

root.mainloop()
