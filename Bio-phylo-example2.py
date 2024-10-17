from Bio import Phylo
import urllib.request

try:
    with open('phylo.xml') as handle:
        tree = Phylo.read(handle, 'phyloxml')
    Phylo.draw(tree)
except urllib.error.URLError as e:
    print(f"Failed to fetch the file: {e}")
except Exception as e:
    print(f"An error occurred: {e}")
