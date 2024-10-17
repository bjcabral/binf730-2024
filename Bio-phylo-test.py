# a script to test BioPhylo in BioPython
from Bio import Phylo
from io import StringIO

# Newick tree string
treedata = "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460)"
#treedata = "(((A,B),(C,D)),(E,f,g));"

# Create a StringIO object
handle = StringIO(treedata)

# Read the tree from the StringIO object
tree = Phylo.read(handle, "newick")

# print the tree
print(tree)

# Draw an ASCII representation of the tree
Phylo.draw_ascii(tree)

