import webutils
import os, json
from ete3 import Tree, NCBITaxa, TreeStyle, NodeStyle, faces, AttrFace

# ------------------------------
# TREE LAYOUT FUNCTIONS AND STYLES

# Adds attribute faces for scientific names and ranks; base function to be used in other layouts
def taxaNameplate(node):
    faces.add_face_to_node(AttrFace("sci_name"), node, column=0)
    faces.add_face_to_node(AttrFace("rank"), node, column=0)

# Colors nodes based on rank, using predetermined palette; (Includes taxa nameplate layout)
def rankedLayout(node):
    palette = data['rankedPalette']
    
    taxaNameplate(node)
    for rank in palette.keys():
        if node.rank == rank:
            nstyle = NodeStyle()
            nstyle['bgcolor'] = palette[rank]
            node.set_style(nstyle)

# ------------------------------
# NCBI TREE FUNCTIONS

# Takes an NCBI tree and a given rank, removes all nodes of lower ranks
def pruneToRank(t, rank, unclassified=False, clean=True):
    ranks = data['ranks']
    if clean:
        cleanTree(t)
    
    # converts str rank input to int
    if (type(rank)) == str and rank in ranks:
        rank = ranks.index(rank)
    
    if type(rank) == int and rank <= len(ranks) - 1: 
        # Iterates through all ranks, starting from the lowest
        for r in reversed(ranks):
            # Detaches all nodes with ranks lower than given rank
            if ranks.index(r) > rank:
                for node in t.search_nodes(rank=r):
                    node.detach()
            # Detaches all descendants of nodes of given rank
            if ranks.index(r) == rank:
                for node in t.search_nodes(rank=r):
                    for subnode in node.iter_descendants():
                        subnode.detach()
        
        # If unclassified is not true, remove unclassified clades with no valid descendants
        if not unclassified:
            for node in t.iter_descendants():
                if 'unclassified' in node.sci_name or 'incertae' in node.sci_name:
                    keep = False
                    for subnode in node.iter_descendants():
                        if subnode.rank in ranks:
                            keep = True
                    if not keep:        
                        node.detach()
    else:
        print("Invalid Rank Type - Requires Valid Rank String or Int Between 0 and 26")

# removes 'environmental samples' nodes and taxa that have been replaced
def cleanTree(t):
    for node in t.iter_descendants():
        if 'environmental' in node.sci_name:
            node.detach()
        if node.sci_name in data['oldTaxa']:
            node.detach()
        
# returns parent of given taxa of given rank
def getParent(taxid, rank, mode='taxid'):
    lineage = ncbi.get_lineage(taxid)
    for parent in lineage:
        if ncbi.get_rank([parent])[parent] == rank:
            if mode == 'taxid':
                return parent
            if mode == 'name':
                return ncbi.get_taxid_translator([parent])[parent]

# ------------------------------           
# main function currently for testing only
def main():
    global root
    global data
    
    root = os.getcwd()
    data = json.load(open(root + '/taxa.json', 'r'))
    
    ncbi = NCBITaxa()

    t = ncbi.get_descendant_taxa(58019, return_tree=True)
    pruneToRank(t, 'family')

    wiki = webutils.wikiScraper()
    wiki.getThumbnails(t)

    t.show(layout=rankedLayout)
    
if __name__ == '__main__': main()