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