import numpy as np
import pandas as pd

def load_input_from_spreadsheet(sheet_id: str):
    """
    Load author, keyword, and journal criteria from a pubmed-sieve compatible Google Sheet based on the sheet identifier.
    
    Parameters
    ----------
    sheet_id: str
        Sheet identifier, i.e. the part of the URL after `/d/` and before any trailing arguments.
        For the sheet at `https://docs.google.com/spreadsheets/d/1HD52dXGvVEHbBejDzUAzVhsvi_T3RT99uNhGzHZlQjo/edit#gid=1061598050` it would be `1HD52dXGvVEHbBejDzUAzVhsvi_T3RT99uNhGzHZlQjo.
        See here for formatting example: https://docs.google.com/spreadsheets/d/1HD52dXGvVEHbBejDzUAzVhsvi_T3RT99uNhGzHZlQjo/edit#gid=1061598050.
        
    Returns
    -------
    pd.DataFrame
        Authors dataframe containing the authors criteria
    pd.Series
        Series containing the keywords
    pd.Series
        Series containing the journals
    
    """
    authors_sheet_name = 'authors'
    keywords_and_journals_sheet_name = 'keywords_and_journals'
    
    # Build URLs to get the sheets via the CSV endpoint.
    authors_url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={authors_sheet_name}"
    keywords_and_journals_url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={keywords_and_journals_sheet_name}"
    
    # Load the sheets into Pandas dataframes.
    try:
        authors_df = pd.read_csv(authors_url, keep_default_na=False, on_bad_lines='skip')
        keywords_journals_df = pd.read_csv(keywords_and_journals_url, keep_default_na=False, on_bad_lines='skip')
        
    except Exception as e:
        print('pubmed-sieve: Failed to get the input Google Sheet. Please make sure it has the same sheet names as the template and is set to "anyone with the link can view".')
        print('pubmed-sieve: The error was: %s' % e)
        
    # Postprocess the authors dataframe
    
    # Remove Unnamed columns
    authors_df = authors_df.loc[:, ~authors_df.columns.str.contains('^Unnamed')]

    # Make sure the name column is present (the only one that is actually required.)
    if 'Name' not in authors_df.columns:
        raise Exception('pubmed-sieve: Didnt find required column `Name` in the authors sheet. Please make sure the spreadsheet is set to "Anyone with the link can view."')
    
    # Remove the example entry, if present
    if authors_df.iloc[0]['Name'] == "Name of the author, formatted like 'Hersh K Bhargava'":
        authors_df = authors_df.drop(index=0)
    
    # Postprocess the keywords and journals dataframes
    
    # Remove Unnamed columns
    keywords_journals_df = keywords_journals_df.loc[:, ~keywords_journals_df.columns.str.contains('^Unnamed')]
    
    # Make sure the required columns are present
    if 'List of Keywords' not in keywords_journals_df.columns:
        raise Exception('pubmed-sieve: Didnt find required column `List of Allowed Keywords` in the keywords/journals sheet.')

    if 'List of Allowed Journals' not in keywords_journals_df.columns:
        raise Exception('pubmed-sieve: Didnt find required column `List of Allowed Journals` in the keywords/journals sheet.')

    # Remove the example entry, if present
    if keywords_journals_df.iloc[0]['List of Keywords'] == "Keywords to match within the title and/or abstract.":
        keywords_journals_df = keywords_journals_df.drop(index=0)
            
    # Get the keywords only
    keywords = keywords_journals_df[keywords_journals_df['List of Keywords'] != '']['List of Keywords']
    
    # Get the journals df
    
    journals = keywords_journals_df[keywords_journals_df['List of Allowed Journals'] != '']['List of Allowed Journals']

    return authors_df, keywords, journals

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

def author_query(fullname, altnames=None, affiliations=None, author_position=None, orcid=None):
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
    
    if author_position is None:
        author_position = 'any'

    # Figure out which author position
    if author_position == 'any':
        poshandle = '[au]'
    elif author_position == 'first':
        poshandle = '[1au]'
    elif author_position == 'last':
        poshandle = '[lastau]'
    elif author_position == 'first_or_last':
        poshandle = ['[1au]', '[lastau]']
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
        elif isinstance(poshandle, list):
            # query += ' AND '

            pos_str = ''

            for i, _poshandle in enumerate(poshandle):
                trunc_name = truncated_name(_name)

                if i == 0:
                    pos_str += '%s%s' % (trunc_name, _poshandle)
                else:
                    pos_str += ' OR %s%s' % (trunc_name, _poshandle)
                
            query += ' AND (%s)' % (pos_str)

            query += ')'
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

def author_query_from_row(row: pd.core.series.Series):
    """
    Helper function to call `author_query` on a dataframe row
    
    Parameters
    ----------
    row: pandas.core.series.Series
        Row containing KVPS consistent with `authors_df` for `build_authors_query`.
        
    Returns
    -------
    str
        Query string for the author represented by the `row`.
        
    """
    
    # Initialize arguments for the author with default values
    author_position = None
    orcid = None
    affiliations = None
    altnames = None

    # Parse all args out of the row (beware of empty strings).
    author_position = row.get('Position')
    if author_position == '':
        author_position = None

    orcid = row.get('ORCID')
    if orcid == '':
        orcid = None

    # Affilitaions and altnames may be lists (comma separated)
    affiliations = row.get('Affiliations')
    if isinstance(affiliations, str) and affiliations != '':
        affiliations = affiliations.split(',')
        affiliations = [a.strip() for a in affiliations]
    else:
        affiliations = None

    altnames = row.get('Alt Names')
    if isinstance(altnames, str) and altnames != '':
        altnames = altnames.split(',')
        altnames = [a.strip() for a in altnames]
    else:
        altnames = None

    query = author_query(fullname=row['Name'], altnames=altnames, affiliations=affiliations, author_position=author_position, orcid=orcid)
    return query

def build_authors_query(authors_df: pd.DataFrame, require_hasabstract=True):
    """
    Build a Pubmed query based on information about authors.
    
    The query for each author will have, at most, the following structure:
    
        ((ANY name+author_position) AND ((ANY affiliation) OR orcid)
    
    Parameters
    ----------
    authors_df: pd.DataFrame
        DataFrame containing the following columns:
            * "Name" (required): Name of the author in `Firstname M Lastname` format
            * "Alt Names" (optional): Comma separated alterative names in same format as `Name`
            * "Position" (optional): Author position, can be `first`, `last`, or `any`. Defaults to `any`.
            * "ORCID" (optional): ORCID of the authors
            * "Affiliations": Comma separated list of affiliations.
            
    require_hasabstract: bool
        Whether to include the `hasabstract` keyword, which is helpful in filtering out non-peer-reviewed errata, news, etc.
        
    Returns
    -------
    str
        PubMed query based on the constraints in the dataframe.
        
    """
    # If no authors, return empty string
    if len(authors_df) == 0:
        return ''

    # Build queries for each author row in place.
    authors_df['query'] = authors_df.apply(lambda row: author_query_from_row(row), axis=1)

    # Concatenate all the queries into a final authors query
    authors_query = ''
    i = 0
    for _, row in authors_df.iterrows():
        authors_query += row['query'] 
        if i < len(authors_df)-1:
            authors_query += ' OR '
        i += 1

    # Add the hasabstract flag if desired
    if require_hasabstract and len(authors_df) > 0:
        authors_query += ' AND hasabstract'
    
    return authors_query

def build_keyword_and_journal_query(keywords, journals, require_hasabstract=True):
    """
    Build a Pubmed query for keywords and journals. 
    
    The resulting query will have the logical structure:
    
        `(ANY keyword) AND (ANY journal) AND hasabstract`
    
    Parameters
    ----------
    keywords: list-like
        List of keywords
    journals: list-like
        List of journal names
    require_hasbastract: bool
        Whether ot add the `hasabstract` requirement
    
    Returns
    -------
    str
        The query
        
    """
    # First, add the keywords with [tiab] and OR flags.
    kwquery = ""
    if len(keywords) > 0:
        kwquery += "("
        for i, kw in enumerate(keywords):
            kwquery += '(\"%s\"[tiab])' % kw

            if i < len(keywords) - 1:
                kwquery += ' OR '

        kwquery += ')'
    
    # Now form the journals query
    jquery = ""
    if len(journals) > 0:
        jquery += '('
        for i, jn in enumerate(journals):
            jquery += '(\"%s\"[journal])' % jn

            if i < len(journals) - 1:
                jquery += ' OR '

        jquery += ')'

    # Join the queries and add hasabstract if desired
    if len(kwquery) > 0 and len(jquery) > 0:
        query = kwquery + ' AND ' + jquery
    elif len(kwquery) > 0 and len(jquery) == 0:
        query = kwquery
    elif len(kwquery) == 0 and len(jquery) > 0:
        print('Journal list provided, but no keywords. Journal filters only apply to the keyword query, so will be ignored in this case.')
        print('If you want to filter an authors-only search by journals, you can append the following to the authors query:')
        print(' AND %s' % jquery)
        query = kwquery
    else:
        return ''
    
    if require_hasabstract and len(query) > 0:
        query += ' AND hasabstract'

    return query

def gen_pubmed_rss_link_for_query(feed_name: str, query_string: str):
    """
    Use a headless browser to generate a PubMed RSS feed link given a search string.
    
    Note: This code uses a headless browser and could break unexpectedly if NCBI / PubMed changes their site structure.
    
    Parameters
    ----------
    feed_name: str
        Name of the feed
    query_string: str
        Pubmed query string from which to build the feed
    
    Returns
    -------
    str
        String URL of the RSS feed
        
    """
    illegal_chars = ['&', '=', '<', '>', '/']
    if any([c in feed_name for c in illegal_chars]):
        raise ValueError('Feed name contains illegal characters: %s. Pubmed doesnt allow any of "&=<>/" in RSS feed names.' % illegal_chars)

    try:
    
        # Use selenium to interface with the headless browser
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support.select import Select

        # Set up a headless Chromium driver.
        chrome_options = webdriver.ChromeOptions()
        chrome_options.add_argument('--headless')
        chrome_options.add_argument('--no-sandbox')
        chrome_options.add_argument('--disable-dev-shm-usage')
        chrome_options.add_argument("--window-size=1920,1080")
        user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.50 Safari/537.36'
        chrome_options.add_argument(f'user-agent={user_agent}')
        chrome_options.add_argument('--allow-running-insecure-content')
        
        driver = webdriver.Chrome(options=chrome_options)

        # Navigate to the PubMed website
        driver.get("https://www.ncbi.nlm.nih.gov/pubmed/")

        # Enter the search query in the search box
        search_box = driver.find_element("id", "id_term")
        search_box.send_keys(query_string)

        # Submit the search
        search_box.submit()

        driver.implicitly_wait(3) # implicit wait for 3 seconds to account for loading time

        # Find and click the 'Create RSS button'
        rss_button = driver.find_element("link text", "Create RSS")
        driver.implicitly_wait(1) # gives an implicit wait for 20 seconds
        rss_button.click()

        # Find and fill the RSS feed name field
        name_box = driver.find_element("id", "rss-name")
        name_box.clear()
        name_box.send_keys(feed_name)

        # Set the number of results selector to 100
        limit_input = Select(driver.find_element("id", "rss-limit"))
        limit_input.select_by_visible_text('100')

        # Find and click the Create RSS Button
        driver.find_element(By.XPATH,'//button[normalize-space()="Create RSS"]').click()

        # Wait for the link to appear
        WebDriverWait(driver, 30).until(lambda driver: driver.find_element(By.XPATH, "//input[@id='rss-link']").get_attribute('value').strip() != '')

        # Extract and return the link
        link_output = driver.find_element("id", "rss-link").get_attribute('value')

        return link_output

    except Exception as e:
        print('pubmed-sieve: Exception occurred while generating the Pubmed RSS link. You can generate on manually by navigating to Pubmed and clicking the "Create RSS" button.')
        print(e)

def build_query_from_spreadsheet_url(url: str, sheet_id=None) -> str:
    """
    Build an overall query string from a Google Sheet at a URL.

    Parameters
    ----------
    url: str
        URL of the Google Sheet in question, which should contain an `authors` and `keywords_and_journals` subsheet.

    Returns
    -------
    str
        Query

    """
    if sheet_id is None:
        try:
            sheet_id = url.split('/d/')[1].split('/')[0]
        except Exception as e:
            print('pubmed-sieve: Unable to parse the sheet_id out of the URL. Perhaps try supplying manually.')

    authors_df, keywords, journals = load_input_from_spreadsheet(sheet_id)
    print('pubmed-sieve parsed spreadsheet with id: %s' % sheet_id)
    print('pubmed-sieve building query with %i authors, %i keywords, and %i journals.' % (len(authors_df), len(keywords), len(journals)))

     # Build the authors query
    authors_query = build_authors_query(authors_df=authors_df, require_hasabstract=True)

    # Build the keywords+journals query
    kw_query = build_keyword_and_journal_query(keywords=keywords, journals=journals, require_hasabstract=True)

    # Stitch the queries together, accounting for cases where one query is empty.
    if len(authors_query) > 0 and len(kw_query) > 0:
        final_query = '(%s) OR (%s)' % (authors_query, kw_query)
    elif len(authors_query) > 0 and len(kw_query) == 0:
        final_query = authors_query
    elif len(authors_query) == 0 and len(kw_query) > 0:
        final_query = kw_query
    else:
        raise Exception('No query was generated. THe spreadsheet may have been blank.')

    return final_query