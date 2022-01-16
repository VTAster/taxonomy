from taxutilities import *
import os, io, time, json, shutil, csv
import requests, lxml
import folium, imgkit
import pandas as pd
from PIL import Image
from bs4 import BeautifulSoup
from zipfile import *
from ete3 import Tree, NCBITaxa, TreeStyle, NodeStyle, faces, AttrFace, ImgFace

root = os.getcwd()
data = json.load(open(os.path.join(root, 'taxa.json'), 'r'))
headers = data['headers']

# ------------------------------
# GENERAL UTILITIES
def saveImage(url):
    url_content = requests.get(url, headers=headers).content

    filename = url.split('/')[-1]
    directory = os.path.join(root, 'images', filename)

    _bytes = io.BytesIO(url_content)
    image = Image.open(_bytes).convert('RGB')
    with open(directory, 'wb') as f:
        image.save(f, 'JPEG', quality=85)
        
    return directory

# ------------------------------
# WEB SCRAPERS

# WIKIMEDIA API/SITE SCRAPER
class wikiScraper:
    def __init__(self):
        self.root = os.getcwd()
        self.images = os.path.join(self.root, 'images')
        self.data = json.load(open(self.root + '/taxa.json', 'r'))

        self.api = self.data['wikiAPI']
        self.ranks = self.data['ranks']
  
    # function that uses MediaWiki Query action, and returns the json retrieved
    def query(self, _format='json', params={}):
        params['action'] = 'query'
        params['format'] = _format
        
        r = requests.get(self.api, params=params, headers=headers)
        if r.status_code == 200:
            try:
                return r.json()['query']
            except:
                print('No query found for this request')
        else:
            print('Request failed')
            
# GBIF DATASET SCRAPER
class GBIFScraper:
    def __init__(self):
        self.root = os.getcwd()
        self.data = json.load(open(self.root + '/taxa.json', 'r'))
        
        self.datasets = os.path.join(self.root, 'datasets')

        self.api = self.data['gbifAPI']
    
    ''' One-step function to get a GBIF dataset from a given uuid
    Gets url using datasetRequest method, then checks if dataset already exists
    Returns datafile directories if so. If not, calls datasetDownload first, then returns the datafiles' directories
    '''
    def getDataset(self, uuid):
        url = self.datasetRequest(uuid)
        
        if url is not None:
            setname = url.split('/')[-1].replace('.zip', '')
            setdir = os.path.join(self.datasets, setname)

            if os.path.exists(setdir):
                return self.datasetFiles(setdir)
            else:
                setdir = self.datasetDownload(url)
                return self.datasetFiles(setdir)
    
    # Queries GBIF api for the endpoint of a given dataset, returns url
    def datasetRequest(self, uuid):
        try:
            request = requests.get("{}/dataset/{}/endpoint".format(self.api, uuid))
            return request.json()[0]['url']
        except:
            print('Dataset Request Failed')
    
    # Finds and returns the directory of any csv files in the given directory
    def datasetFiles(self, setdir):
        datafiles = []
        for root, dirs, files in os.walk(setdir):
            for file in files:
                if '.csv' in file:
                    datafiles.append(os.path.join(root, file))
        
        if len(datafiles) > 0:
                return datafiles
        else:
            print('No datafiles found')
    
    # Takes dataset url, downloads and extracts zip file into datasets directory, returns directory of resulting folder
    # Doesnt work with big datasets, idk why, fml
    def datasetDownload(self, url):
        filename = url.split('/')[-1]
        setzip = os.path.join(self.datadir, filename)
        setdir = os.path.join(self.datadir, filename.replace('.zip', ''))
        
        try:
            with requests.get(url, stream=True) as r:
                with open(setzip, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)
            with ZipFile(setzip) as archive:
                archive.extractall(setdir)
            return setdir
        except:
            return

    # function for instantiating NameUsages from taxids; filters invalid taxa;
    def nameUsage(self, taxa):
        if type(taxa) == int:
            taxa = getName(taxa)

        r = requests.get(f'http://api.gbif.org/v1/species/match?scientificName={taxa}')
        if r.status_code == 200:
            return r.json()
        
    # function that takes a nameUsage and optionally a dataset key, and returns a json page for occurrences of that taxa
    def occurrences(self, nameUsage, datasetKey=None):
        rank = nameUsage['rank'].lower()

        params = {}
        params[f'{rank}Key'] = nameUsage[f'{rank}Key']
        if datasetKey is not None:
            params['datasetKey'] = datasetKey

        r = requests.get(f'http://api.gbif.org/v1/occurrence/search?', params=params)
        return r.json()

# Handles Plants of the World Online data scraping
class POTWOScraper:
    def __init__(self):
        gbif = GBIFScraper()
        ipniCSV = gbif.getDataset('046bbc50-cae2-47ff-aa43-729fbf53f7c5')
        for csv in ipniCSV:
            self.ipni = pd.read_csv(csv, header=None)
            self.ipni.rename(columns={0:'lsid', 1:'name', 2:'author', 3:'rank', 4:'family', 5:'lsid2', 6:'lsid3', 7:'lsid4', 8:'year', 9:'citation', 10:'?', 11:'link'}, inplace=True)
            
        self.ranks = ['family', 'genus', 'species']
    
    ''' Searches IPNI for LSIDs matching given taxa, then returns as list'''
    def getLSIDS(self, taxa):
        if type(taxa) == int:
            taxa = getName(taxa)
            
        indecies = []
        for index, name in self.ipni.name.items():
            if name == taxa:
                indecies.append(index)
        lsids = []
        for index, lsid in self.ipni.lsid.items():
            if index in indecies:
                lsids.append(lsid)
                
        return lsids
    
    ''' Takes taxa as id or name, returns distribution from POTWO as 
        list of regions or as dictionary of names and distributions
    '''
    def getDistribution(self, taxa):
        if type(taxa) == str:
            taxa = getTaxid(taxa)
        
        if taxa is not None:
            dist = None

            rank = getRank(taxa)
            if rank not in self.ranks:
                print(f'Distributions for {rank} not currently supported')
                return

            dist = self.distribution(taxa)
            if dist is not None:
                return dist
            else:
                dist = self.distributionFromDescendants(taxa, rank)
                return dist
    
            
    # Internal function that finds valid lsid for given taxa, returns native distribution listed on POTWO page
    def distribution(self, taxa):
        lsids = self.getLSIDS(taxa)
    
        for lsid in lsids:
            if lsid is not None:
                request = requests.get(f'http://www.plantsoftheworldonline.org/taxon/{lsid}')
                if request.status_code == 200:
                    soup = BeautifulSoup(request.content, 'lxml')
                    try:
                        distraw = soup.find('div', id='distribution-listing').find('p')
                        distclean = distraw.get_text().strip()

                        dist = []
                        for loc in distclean.split(', \r\n              \r\n              '):
                            dist.append(loc)
                        return dist
                    except:
                        print(f'Distribution not found for {taxa}')
    
    # Takes a taxid with no distribution on POTWO, returns a distribution using its descendants' pages
    def distributionFromDescendants(self, taxa, rank):
        dist = []
        dists = {}
        
        descendants = ncbi.get_descendant_taxa(taxa, return_tree=True)
        if type(descendants) == list:
            return
        
        for node in descendants.iter_descendants():
            if node.rank in self.ranks and self.ranks.index(node.rank) == self.ranks.index(rank) + 1:
                node_dist = self.getDistribution(int(node.name))
                if node_dist is not None:
                    dists[int(node.name)] = node_dist
                elif node.rank in self.ranks and node.rank != 'species':
                    dists[int(node.name)] = self.distributionFromDescendants(int(node.name), node.rank)
        
        for d, regions in dists.items():
            if regions is not None:
                for region in regions:
                    if region not in dist:
                        dist.append(region)
                    
        return dist
        
''' eFloras.org Data Scraper
    Takes flora name or id on initialization
    New instance must be created for each flora
    Provides functions to scrape taxa names and relationships
    Works with all major floras
'''
class eFlora:
    def __init__(self, flora):
        self.ncbi = NCBITaxa()
        self.home = 'http://www.efloras.org/'
        self.floraId = self.getFloraID(flora)
    
    ''' Uses related functions to compile and return a dictionary
        containing taxon ids and scientific names for all taxa in
        the flora, properly organized taxonomically
    '''
    def fullTree(self):
        tree = {}
        
        f = self.families()
        for f_id, f_name in f.items():
            genera = {}
            g = self.genera(f_id)
            for g_id, g_name in g.items():
                genera[g_id] = {g_name: self.species(g_id)}
                
            tree[f_id] = {f_name: genera}
            
        return tree
    
    ''' Calls and returns browseTaxa with no arguments, returning 
        dictionary of taxon ids and corresponding families
    '''
    def families(self):
        return self.browseTaxa()
    
    ''' Returns dictionary of taxon ids and corresponding genus names
        for given taxa, or for the whole flora if no argument given
    '''
    def genera(self, taxonID=None):
        genera = {}
        
        for g_ID, g_name in self.browseTaxa(taxonID).items():
            genera[g_ID] = g_name
        
        return genera
    
    ''' Returns a dictionary of taxon ids and corresponding species names
        for given taxa, or for the whole flora if no argument given
    '''
    def species(self, taxonID=None):
        species = {}
        
        if len(self.browseTaxa(taxonID)) == 0:
            species = self.browseTaxa(taxonID, species=True)  
        else:
            for g_ID, g_name in self.browseTaxa(taxonID).items():
                species[g_ID] = {g_name: self.browseTaxa(g_ID, species=True)}
            
        return species
    
    ''' Returns dictionary of taxon ids and corresponding names for
        every taxon below given taxa ID. When no arguments are passed,
        defaults to the highest ranking, which is family.
    '''
    def browseTaxa(self, taxonID=None, species=False):
        pages = []
        pageNumber = 1
        while pageNumber != 0:
            page = requests.get(self.home + 'browse.aspx?', 
                                 headers=headers, 
                                 params={'flora_id': self.floraId, 'page': pageNumber, 'start_taxon_id': taxonID})
            
            if "No taxa found" in page.text:
                pageNumber = 0
            else:
                pages.append(page)
                pageNumber += 1
        
        taxa = {}
        for page in pages:
            soup = BeautifulSoup(page.text, 'lxml')
            
            taxa_list = soup.find_all('tr', class_='underline')
            for item in taxa_list:
                id_listing = item.find('td', class_='small')
                
                if id_listing is not None:
                    taxa_id = id_listing.get_text()
                    taxa_name = item.find('a').get_text()
                    lower_taxa = item.find('a', title='lower taxa')
                    if lower_taxa is None:
                        lower_taxa = item.find('a', title='lower taxon')
                    
                    if taxa_id != taxonID and 'x\\' not in taxa_name and '×' not in taxa_name:
                        if species and lower_taxa is None and 'subg.' not in taxa_name:
                            taxa[taxa_id] = taxa_name
                        elif not species and lower_taxa is not None and len(taxa_name.split(' ')) == 1:
                            taxa[taxa_id] = taxa_name

        return taxa
    
    # converts flora names to flora ids on initialization
    def getFloraID(self, flora):
        if type(flora) == int:
            return flora
        elif type(flora) == str:    
            for floraId, eFlora in data['eFloras'].items():
                if eFlora.lower() == flora.lower():
                    return floraId
    
# ------------------------------
# UNDER CONSTRUCTION
# ------------------------------

''' JEPSON EFLORA DATA SCRAPER - UNDER CONSTRUCTION
'''
class Jepson:
    def __init__(self):
        self.home = 'http://ucjeps.berkeley.edu/eflora/'
    
    ''' Takes taxon id for taxon in Jepson eFlora
        Returns a dictionary containing all retrievable information
        Including higher taxa, name, status, distribution, etc.
    '''
    def taxon(self, tid):
        r = requests.get(self.home + 'eflora_display.php?', headers=headers, params={'tid': tid})
        soup = BeautifulSoup(r.text, 'lxml')
        
        # Identify and store higher taxa and body tables, for later scraping
        content = soup.find('div', {'id': 'content'})
        tables = content.find_all('table')
        
        higher_taxa = []
        for table in tables:
            if table.find('table', class_='taxonomy_table') is not None and table.find('td', class_='column1') is not None:
                higher_taxa.append(table.find('td', class_='column1').get_text())
                print('Upper taxa found')
            if table.find('div', {'id': 'familydesc'}) is None and table.find('div', {'id': 'genusdesc'}) is None:
                if table.find('div', class_='bodyText') is not None:
                    body = table
        
        # Declares and returns dictionary of relevant info. Stops when it reaches taxon's author.
        taxon = {}
        
        for t in higher_taxa:
            t = t.split(': ')
            taxon[t[0]] = t[1]
        
        if body is not None:
            taxon['Name'] = body.find('div', class_='pageMajorHeading').find('b').get_text()

            for item in body.find('div', class_='bodyText').find_all('b'):
                attr_name = item.get_text().replace(':', '')
                attr_value = item.next_sibling.get_text().strip(' ').capitalize()

                if attr_name.isupper():
                    taxon['Status'] = attr_name
                else:
                    taxon[attr_name] = attr_value

                if 'eflora author' in attr_name.lower():
                    break
        
        return taxon

# FOLIUM CLASS FOR CREATING MAPS FROM POTWO DISTRIBUTIONS - UNDER CONSTRUCTION
class MapMaker:
    def __init__(self):
        self.root = os.getcwd()
        self.level3 = os.path.join(self.root, 'maps', 'wgsrpd', 'level3.geojson')
        
        '''
        self.options = {'format': 'png', 
                        'height': '1000',
                        'width': '3000',
                        'crop-h': '100',
                        'crop-w': '100',
                        'crop-x': '50',
                        'crop-y': '50'}
        '''
        
    ''' Currently the only function of this class. Generates folium map of the 
        world, applies level 3 WGSRPD region overlay, then changes opacity of
        regions not in given distribution to 0. Returns processed map.
    '''
    def distributionMap(self, dist):
        m = folium.Map(location=[0, 0], zoom_start=2.4)
        
        self.dist = dist
        folium.GeoJson(self.level3, 
                       name='WGSRPD Level 3', 
                       style_function=self.dist_function
                      ).add_to(m)
        
        return m
        
        ''' Work-In-Progress method of rendering image from map
        path = os.path.join(self.root, 'temp', 'map.html')
        m.save(path)
        
        img = imgkit.from_file(path, os.path.join(self.root, 'maps', 'map.png'), options=self.options)
        return img
        '''
    
    # style function that sets opacity of regions not within
    def dist_function(self, feature):
        if feature['properties']['LEVEL3_NAM'] in self.dist:
            return { 'fillOpacity': 1,
                     'opacity': 1 }
        else:
            return { 'fillOpacity': 0,
                     'opacity': 0 }