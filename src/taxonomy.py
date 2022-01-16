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
    
    def thumbnailQueries(self, taxa, size):
        params = {'prop': 'pageimages', 'piprop': 'thumbnail', 'pithumbsize': size}
        
        queries = []
        
        # Creates list of queries for each 48/49 titles
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
        
    # All-in-one function for getting a formatted tree with optional thumbnails for a given taxid
    def getTree(self, taxa, rank='family', unclassified=False, clean=True, thumbnails=True):
        taxa = getTaxid(taxa)
        
        tree = self.ncbi.get_descendant_taxa(taxa, return_tree=True)
        
        prunedTree = pruneToRank(tree, rank)\
        
        if prunedTree != None:
            if thumbnails:
                self.getThumbnails(prunedTree)
            return prunedTree

def main():
    tax = Taxonomy()
    
    taxa = input("Enter a Taxa (Scientific Name or NCBI ID): ")
    rank = input("Rank: ").lower()
    
    response = input("Thumbnails? (Y|N)").lower()
    if response == 'y' or response == 'yes':
        thumbnails = True
    else:
        thumbnails = False
    
    tree = tax.getTree(taxa, rank=rank)
    saveTree(tree)
    
if __name__ == '__main__': main()