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
    
def truncated_name(name):
    """
    Convert a string name of the form "Firstname Minitial Lastname" or "Firstname Lastname" to Lastname FinitialMinitial
    
    e.g. Hersh K Bhargava -> Bhargava HK; Hersh Bhargava -> Bhargava H
    
    """
    
    split = name.split(' ')
    if len(split) == 3:
        result = '%s %s%s' % (split[2], split[0][0], split[1][0])
    elif len(split) == 2:
        result = '%s %s' % (split[1], split[0][0])
    else:
        print('invalid result from split, suspect invalid name. INptu was : %s'%name)
        result = None
    return result

def author_query(fullname, altnames=None, affiliations=None, author_position='any', orcid=None):
    """
    Generate a pubmed query string to match:
    
        ((ANY name+author_position) AND ((ANY affiliation) OR orcid)
    
    Very annoyingly, Pubmed doesn't support [1au] or [lastau] tags with full author names
    So would need to do:
    
        (Hersh K Bhargava[FAU]) AND (Bhargava HK[1au])
        
    Also can't search by ORCID position
        
    Parameters
    ----------
    name : str
        Primary search term for the name of the author. Must be of the form 'Hersh K Bhargava' or 'Hersh Bhargava'
    altnames : list[str]
        List of alternative names (each matched with OR). Used as-is.
    affiliations : list[str]
        List of affiliations, will be matched for at least one
    position : str
        Identifier for the author position. Can be 'any', 'first', or 'last'
    orcid : str
        ORCID identifier for the author
    
    """
    
    # Assemble author names
    names = [fullname]
    if altnames != None:
        names += altnames
    
    # Figure out which author position
    if author_position == 'any':
        poshandle = '[au]'
    elif author_position == 'first':
        poshandle = '[1au]'
    elif author_position == 'last':
        poshandle = '[lastau]'
    else:
        poshandle = '[au]'
        
    query = ''
        
    # First build the name query (each added with OR)
    for _name in names:
        # Figure out whether it's truncated or full format
        if len(query) > 0:
            query += ' OR '
        else:
            query += '('
            
        query += '(%s[FAU]' % _name
        
        # If searching by author position, must use truncated query format (Bhargava HK)
        if poshandle == '[1au]' or poshandle == '[lastau]':
            trunc_name = truncated_name(_name)
            query += ' AND %s%s)' % (trunc_name, poshandle)
        else:
            query += ')'
        
    # Add affiliation tags and ORCID if present
    if affiliations != None:
        query += ' AND ('
        for i, affil in enumerate(affiliations):
            if i > 0:
                query += ' OR '
            query += '%s[affil]' % affil 
            
        if orcid is not None:
            query += ' OR %s[auid]' % orcid
        
        query += ')'

    else:
        # No affiliations, just ORCID
        if orcid is not None:
            query += ' AND %s[auid]' % orcid
        
    query += ')'
    return query

def author_query_from_row(row):
    author_position = 'last'
    orcid = None
    affiliations = None
    altnames = None
    
    if row['Position'] != '':
        author_position = row['Position']
    
    if row['ORCID'] != '':
        orcid = row['ORCID']
        
    if row['Affiliations'] != '':
        affiliations = row['Affiliations'].split(',')
        affiliations = [a.strip() for a in affiliations]
        
    if row['Alt Names'] != '':
        altnames = row['Alt Names'].split(',')
        altnames = [a.strip() for a in altnames]
    
    return author_query(fullname=row['Name'], altnames=altnames, affiliations=affiliations, author_position=author_position, orcid=orcid)