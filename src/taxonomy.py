from webutilities import *
from taxutilities import *
from ete3 import Tree, NCBITaxa, TreeStyle, NodeStyle, faces, AttrFace, ImgFace

class Taxonomy:
    def __init__(self):
        self.root = os.getcwd()
        self.data = json.load(open(os.path.join(root, 'taxa.json')))
        self.ranks = self.data['ranks']
        self.ncbi = NCBITaxa()
        self.wiki = wikiScraper()

    ''' One-step thumbnail function, takes ncbi tree,
        uses urls to download all thumbnails that haven't been downloaded yet,
        then attaches them as image face to corresponding nodes on the tree'''
    def getThumbnails(self, t, size=200, column=0, position='aligned'):
        urls = self.thumbnailUrls(t, size)
        thumbnails = {}
        for taxa, url in urls.items():
            filedir = os.path.join(self.root, 'images', url.split('/')[-1])
            
            if os.path.exists(filedir):
                print(f"Using image from memory for {taxa}")
            else:
                print(f"Scraping thumbnail for {taxa}")
                saveImage(url)
                
            thumbnails[taxa] = filedir         
        
        for node in t.iter_leaves():
            if node.sci_name in thumbnails.keys():
                node.add_face(ImgFace(thumbnails[node.sci_name]), column=column, position=position)
    
    ''' Major sub function, takes ncbi tree and returns urls
        from WikiSpecies for thumbnails that correspond
        to each leaf in the tree '''
    def thumbnailUrls(self, t, size):
        # Compile names and ranks of the trees leaves
        taxa = {}
        for node in t.iter_leaves():
            taxa[node.sci_name] = node.rank
        
        # Get list of queries using names and ranks
        queries = self.thumbnailQueries(taxa, size)
        
        # Compiles pages into a single list
        pages = []
        for query in queries:
            for page in query['pages'].values():
                pages.append(page)
                
        # Matches names with urls
        urls = {}
        for name in taxa:
            for page in pages:
                if name in page['title']:
                    try:
                        if any(i in page['thumbnail']['source'] for i in self.data['genericImages']):
                            continue
                        else:
                            urls[name] = page['thumbnail']['source']
                    except:
                        continue
        
        # Attempts to find thumbnails for missed taxa
        for node in t.iter_leaves():
            if node.sci_name not in urls.keys() and node.rank in self.ranks and self.ranks.index(node.rank) < self.ranks.index('species'):
                node_t = ncbi.get_descendant_taxa(node.name, return_tree=True)
                node_taxa = {}
                if type(node_t) is not list:
                    for subnode in node_t.iter_leaves():
                        node_taxa[subnode.sci_name] = subnode.rank

                    node_queries = self.thumbnailQueries(node_taxa, size)
                    node_pages = []
                    for query in node_queries:
                        for page in query['pages'].values():
                            node_pages.append(page)

                    for name in node_taxa:
                        for page in node_pages:
                            if name in page['title']:
                                try:
                                    if any(i in page['thumbnail']['source'] for i in self.data['genericImages']):
                                        continue
                                    else:
                                        urls[node.sci_name] = page['thumbnail']['source']
                                except:
                                    continue
                                            
        return urls
    
    ''' Batch queries taxa names 48/49 at once and returns list of queries '''
    def thumbnailQueries(self, taxa, size):
        params = {'prop': 'pageimages', 'piprop': 'thumbnail', 'pithumbsize': size}
        
        queries = []
        
        titles = []
        for name, rank in taxa.items():
            titles.append(name)
            
            if rank in self.ranks and self.ranks.index(rank) < self.ranks.index('family'):
                fam = getParent(name, rank='family', mode='name')
                titles.append(f'{name} ({fam})')
                
            if len(titles) >= 48:
                params['titles'] = '|'.join(titles)
                query = self.wiki.query(params=params)
                if query is not None:
                    queries.append(query)
                
                titles = []
                params['titles'] = ''
                
        if len(titles) > 0:
            params['titles'] = '|'.join(titles)
            query = self.wiki.query(params=params)
            if query is not None:
                queries.append(query)
                
        return queries
        
    ''' All-in-one function for getting a formatted tree with optional thumbnails for a given taxid '''
    def getTree(self, taxa, rank='family', unclassified=False, clean=True, thumbnails=True):
        taxa = getTaxid(taxa)
        
        if taxa is not None:
            tree = self.ncbi.get_descendant_taxa(taxa, return_tree=True)
            prunedTree = self.pruneToRank(tree, rank)
            if prunedTree is not None:
                if thumbnails:
                    self.getThumbnails(prunedTree)
                return prunedTree
        else:
            print('Taxa invalid, try again')
            return
        
    ''' Using multiple sub functions, 
        takes an NCBI tree and a given rank, 
        removes all nodes of lower ranks '''
    def pruneToRank(self, t, rank, unclassified=False, clean=True):
        if clean:
            self.cleanTree(t)

        prunedTree = self.detachLowerRanks(t, rank)
        if prunedTree is not None:
            if not unclassified:
                self.removeUnclassified(prunedTree)
            return prunedTree
        else:
            return
        
    ''' Takes an ncbi tree, a rank, and a list of lower taxa
        Prunes to given rank while keeping given lower taxa intact
    '''
    def pruneTaxa(self, t, rank, taxa, unclassified=False, clean=True):
        for key, tax in enumerate(taxa):
            if type(tax) == str:
                taxa[key] = getTaxid(tax)
        
        taxaParents = {}
        for tax in taxa:
            lineage = self.ncbi.get_lineage(tax)
            for parent in lineage:
                if getRank(parent) == rank:
                    taxNode = self.pruneToRank(self.ncbi.get_descendant_taxa(tax, return_tree=True), rank=getRank(tax))
                    taxaParents[taxNode] = parent
        
        if clean:
            self.cleanTree(t)
        
        prunedTree = self.detachLowerRanks(t, rank)
        if prunedTree is not None:
            if not unclassified:
                self.removeUnclassified(prunedTree)
            if taxaParents:
                for leaf in prunedTree.iter_leaves():
                    if int(leaf.name) in taxaParents.values():
                        for taxNode, parent in taxaParents.items():
                            if int(leaf.name) == parent:
                                leaf.add_child(taxNode)
                            
            return prunedTree
        else:
            return
    
    ''' Removes all nodes below given rank
        Starts by removing lowest ranks, then moves upward
        Then removes subnodes of taxa at given rank '''
    def detachLowerRanks(self, t, rank):
        prunedTree = t
        
        if (type(rank)) == str and rank in self.ranks:
            rank = self.ranks.index(rank)

        if type(rank) == int and rank <= len(self.ranks) - 1:
            for r in reversed(self.ranks):
                if self.ranks.index(r) > rank:
                    for node in prunedTree.search_nodes(rank=r):
                        node.detach()
                if self.ranks.index(r) == rank:
                    for node in prunedTree.search_nodes(rank=r):
                        for subnode in node.iter_descendants():
                            subnode.detach()
            return prunedTree
        else:
            print("Invalid Rank Type - Requires Valid Rank String or Int Between 0 and 26")
            return None
    
    ''' Removes unclassified and incertae sedis taxa, if they have no valid subtaxa '''
    def removeUnclassified(self, t):
        for node in t.iter_descendants():
            if 'unclassified' in node.sci_name or 'incertae' in node.sci_name:
                keep = False
                for subnode in node.iter_descendants():
                    if subnode.rank in self.ranks:
                        keep = True
                if not keep:
                    node.detach()

    ''' Removes 'environmental samples' nodes and outdated taxa '''
    def cleanTree(self, t):
        for node in t.iter_descendants():
            if 'environmental' in node.sci_name:
                node.detach()
            if node.sci_name in data['oldTaxa']:
                node.detach()

    ''' Renders given tree to given format, outputs to trees subdirectory '''
    def saveTree(self, t, _format='svg', layout=rankedLayout):
        directory = os.path.join(root, 'trees', f'{t.sci_name}.{_format}')
        t.render(directory, layout=layout)

def main():
    tax = Taxonomy()
    
    '''
    taxa = input("Enter a Taxa (Scientific Name or NCBI ID): ")
    rank = input("Rank: ").lower()
    
    savedTaxa = []
    savedTaxaResponse = input("Would you like to keep any lower taxa (Y|N)").lower()
    if savedTaxaResponse == 'y' or savedTaxaResponse == 'yes':
        savedTaxaNames = input("Enter taxa names (i.e. 'taxa, taxa, taxa')").split(", ")
        for name in savedTaxaNames:
            nameToID = getTaxid(name)
            savedTaxa.append(nameToID)
    
    thumbnailsResponse = input("Thumbnails? (Y|N)").lower()
    if thumbnailsResponse == 'y' or thumbnailsResponse == 'yes':
        thumbnails = True
    else:
        thumbnails = False
    
    tree = tax.getTree(taxa, rank=rank, thumbnails=thumbnails)
    if tree is not None:
        tax.saveTree(tree)
    '''
    
if __name__ == '__main__': main()