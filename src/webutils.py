import io, time
import requests, lxml
from PIL import Image
from bs4 import BeautifulSoup as bs

class wikiScraper:
    def __init__(self, data):
        self.headers = data['headers']
        self.api = data['wikiAPI']
        self.params = data['thumbnailFetch']
        self.ranks = data['ranks']
    
    def getThumbnails(self, t, size=200, column=0, position='aligned'):
        queries = []
        
        names = {}
        for node in t.iter_leaves():
            if node.rank in self.ranks and self.ranks.index(node.rank) > 13:
                fam = getParent(node.name, 'family', mode='name')
                names[node.sci_name] = fam
            else:
                names[node.sci_name] = None
            if len(names.keys()) == 24:
                queries.append(self.requestThumbnails(names, size))
                names = {}
        if len(names.keys()) < 24:
            queries.append(self.requestThumbnails(names, size))
        
        urls = {}
        for query in queries:
            for node in t.iter_leaves():
                for page in query['pages'].values():
                    if node.sci_name in page['title']:
                        try:
                            urls[node.sci_name] = page['thumbnail']['source']
                        except:
                            continue
        
        dontScrape = []
        for taxa, url in urls.items():
            filename = url.split('/')[-1]
            filedir = '{}/images/{}'.format(root, filename)
            try:
                image = Image.open(filedir)
                for node in t.iter_leaves():
                    if node.sci_name == taxa:
                        print("Using image from memory for {}".format(taxa))
                        node.add_face(ImgFace(filedir), column=column, position=position)
                        dontScrape.append(taxa)
            except:
                continue
        for taxa in dontScrape:
            urls.pop(taxa)
        
        for node in t.iter_leaves():
            for taxa in urls.keys():
                if node.sci_name == taxa:
                    print("Scraping thumbnail for {}".format(taxa))
                    imgDir = self.saveImage(urls[taxa])
                    node.add_face(ImgFace(imgDir), column=column, position=position)
                    
    def requestThumbnails(self, names, size):
        titles = ''
        for name, family in names.items():
            titles += name + '|'
            if names[name] is not None:
                titles += "{} ({})|".format(name, names[name])
        
        self.params['titles'] = titles
        self.params['pithumbsize'] = size
        try:
            request = requests.get(self.api, params=self.params, headers=self.headers)
            query = request.json()['query']
            time.sleep(1)
            return query
        except:
            print("Request failed...")
            return

    def saveImage(self, url):
        try:
            url_content = requests.get(url, headers=self.headers).content
            time.sleep(1)
        except:
            print("Url could not be reached")
        try:
            bytes = io.BytesIO(url_content)

            filename = url.split('/')[-1]
            directory = root + '/images/' + filename

            image = Image.open(bytes).convert('RGB')
            with open(directory, 'wb') as f:
                image.save(f, 'JPEG', quality=85)
            return directory
        except:
            print("Image could not be saved")