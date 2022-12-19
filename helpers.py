import numpy as np
import pandas as pd
from Bio import Entrez
Entrez.email = 'your_email@gmail.com'
def pmids_for_query(query):
    """
    Return PMIDs resulting frmo a query
    
    """
    
    # Search pubmed for the query, returning the PMIDs fo all results (up to 1e4 results)
    handle = Entrez.esearch(db='pubmed', retmax=10000, retmode='xml', term=query)
    searchResults = Entrez.read(handle)
    
    pmids = searchResults['IdList']
    
    return pmids

def pubmed_articles_for_query(query):
    """
    Return a dataframe of articles resulting from a pubmed query.
    """        
    # Search pubmed for the query, returning the PMIDs fo all results (up to 1e4 results)
    handle = Entrez.esearch(db='pubmed', retmax=10000, retmode='xml', term=query)
    searchResults = Entrez.read(handle)
    
    pmids = pmids_for_query(query)
    
    # Get the articles from the PMIDs
    handle = Entrez.efetch(db='pubmed', retmode='xml', id=pmids)
    results = Entrez.read(handle)
    articles = results['PubmedArticle']
    
    df = pd.DataFrame({'PMID': pmids})

    journal_names = []
    titles = []
    dois = []
    abstracts = []
    years = []
    
    results = []

    for i in range(len(articles)):
        # Store the desired information in a dictionary
        result = {}
        
        result['PMID'] = pmids[i]
        
        # MedlineCitation contains all the data of interest.
        article = articles[i]
        citation = article['MedlineCitation']

        # Retrieve the fields of interest. Some have multiple fallback locations.
        result['Journal'] = citation['Article']['Journal']['Title']
        result['Title'] = citation['Article']['ArticleTitle']
        try:
            result['Abstract'] = citation['Article']['Abstract']['AbstractText'][0]
        except:
            result['Abstract'] = None
        try:
            result['Year'] = citation["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]
        except:
            try:
                result['Year'] = citation["Article"]["ArticleDate"][0]["Year"]
            except:
                try:
                    result['Year'] = citation["DateRevised"]["Year"]
                except:
                    result['Year'] = None

        # Get the DOI
        refArray = article['PubmedData']['ArticleIdList']
        doi = None
        for entry in refArray:
            if entry.attributes["IdType"] == "doi":
                doi = entry
        result['doi'] = doi
        
        # Get the article type
        try:
            types = []
            typelist = citation['Article']['PublicationTypeList']
            for t in typelist:
                tt = str(t)
                if "Research Support" not in tt:
                    types.append(tt)
                
            result['Types'] = types
        except:
            result['Types'] = None
        
        results.append(result)
        
    df = pd.DataFrame(results)
    df['Year'] = df['Year'].astype(int)
    return df

def avg_articles_per_year_last5(query):
    """
    Return the average number of articles a year for from 2017-2021 for `query`
    
    """
    
    # Only get last 6 years of data
    query += (" AND \"last 6 years\"[dp]")
    
    try:
        data = pubmed_articles_for_query(query)
        avgs = [len(data[data['Year'] == val]) for val in [2017, 2018, 2019, 2020, 2021]]
        print(np.mean(avgs))
        return np.mean(avgs)
    except Exception as e:
        print('failed on: %s'%query)
        print(e)
        return np.nan
    
