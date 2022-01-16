import os, json
from ete3 import Tree, NCBITaxa, TreeStyle, NodeStyle, faces, AttrFace, ImgFace

root = os.getcwd()
data = json.load(open(os.path.join(root, 'taxa.json')))
ncbi = NCBITaxa()

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
def pruneToRank(tree, rank, unclassified=False, clean=True):
    prunedTree = tree
    ranks = data['ranks']
    
    if clean:
        cleanTree(prunedTree)
    
    # converts str rank input to int
    if (type(rank)) == str and rank in ranks:
        rank = ranks.index(rank)
    
    if type(rank) == int and rank <= len(ranks) - 1: 
        # Iterates through all ranks, starting from the lowest
        for r in reversed(ranks):
            # Detaches all nodes with ranks lower than given rank
            if ranks.index(r) > rank:
                for node in prunedTree.search_nodes(rank=r):
                    node.detach()
            # Detaches all descendants of nodes of given rank
            if ranks.index(r) == rank:
                for node in prunedTree.search_nodes(rank=r):
                    for subnode in node.iter_descendants():
                        subnode.detach()
        
        # If unclassified is not true, remove unclassified clades with no valid descendants
        if not unclassified:
            for node in prunedTree.iter_descendants():
                if 'unclassified' in node.sci_name or 'incertae' in node.sci_name:
                    keep = False
                    for subnode in node.iter_descendants():
                        if subnode.rank in ranks:
                            keep = True
                    if not keep:        
                        node.detach()
        return prunedTree
    else:
        print("Invalid Rank Type - Requires Valid Rank String or Int Between 0 and 26")
        return None

# removes 'environmental samples' nodes and taxa that have been replaced
def cleanTree(t):
    for node in t.iter_descendants():
        if 'environmental' in node.sci_name:
            node.detach()
        if node.sci_name in data['oldTaxa']:
            node.detach()

# renders given tree to given format, output in trees directory
def saveTree(tree, _format='svg', layout=rankedLayout):
    directory = os.path.join(root, 'trees', f'{tree.sci_name}.{_format}')
    
    tree.render(directory, layout=layout)

# ------------------------------
# TAXA UTILITIES

# returns parent at specified rank of given taxa
def getParent(taxa, rank, mode='taxid'):
    taxa = getTaxid(taxa)
    
    lineage = ncbi.get_lineage(taxa)
    for parent in lineage:
        if ncbi.get_rank([parent])[parent] == rank:
            if mode == 'taxid':
                return parent
            if mode == 'name':
                return ncbi.get_taxid_translator([parent])[parent]

# takes single scientific name and returns single corresponding taxid
def getTaxid(taxa):
    if type(taxa) == int:
        return taxa
    else:
        for spec_char, char in data['specialChars'].items():
            taxa.replace(spec_char, char)

        taxa_to_taxid = ncbi.get_name_translator([taxa]).get(taxa)
        if taxa_to_taxid is not None:
            taxid = taxa_to_taxid[0]
            return taxid
    
# takes single taxid and returns single corresponding scientific name         
def getName(taxa):
    if type(taxa) == str:
        for spec_char, char in data['specialChars'].items():
            taxa.replace(spec_char, char)
        return taxa
    else:
        taxid_to_taxa = ncbi.get_taxid_translator([taxa]).get(taxa)
        if taxid_to_taxa is not None:
            taxa = taxid_to_taxa
            return taxa

# takes single taxa and returns single corresponding rank
def getRank(taxa):
    if type(taxa) == str:
        taxa = getTaxid(taxa)
    
    rank = ncbi.get_rank([taxa])[taxa]
    return rank