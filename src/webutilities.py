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

# ------------------------------
# GENERAL UTILITIES
def saveImage(url):
    headers = data['headers']
    
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

        self.headers = self.data['headers']
        self.api = self.data['wikiAPI']
        self.ranks = self.data['ranks']
  
    # function that uses MediaWiki Query action, and returns the json retrieved
    def query(self, _format='json', params={}):
        params['action'] = 'query'
        params['format'] = _format
        
        r = requests.get(self.api, params=params, headers=self.headers)
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

        self.headers = self.data['headers']
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
        
        for node in descendants.iter_descendants():
            if node.rank in self.ranks and self.ranks.index(node.rank) == self.ranks.index(rank) + 1:
                node_dist = self.getDistribution(int(node.name))
                if node_dist is not None:
                    dists[int(node.name)] = node_dist
                elif node.rank in self.ranks and node.rank != 'species':
                    dists[int(node.name)] = self.descendantDistribution(int(node.name), node.rank)
        
        for d, regions in dists.items():
            if regions is not None:
                for region in regions:
                    if region not in dist:
                        dist.append(region)
                    
        return dist

    
# ------------------------------
# UNDER CONSTRUCTION
# ------------------------------

# JEPSON EFLORA DATA SCRAPER - UNDER CONSTRUCTION
class Jepson:
    def __init__(self):
        

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