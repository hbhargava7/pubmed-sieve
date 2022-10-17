# follow-the-literature
# TLDR
Build pubmed queries (that can be turned into RSS feeds) from Google sheet describing parameters.

### Type 1: People search – return articles from people based on name/descriptors
### Type 2: Keyword search – return articles from keywords + journals

Example Google Sheet: https://docs.google.com/spreadsheets/d/108U4eq7zkwbqaygDMlCKe1cp0q0BlJ6Z0beI7Q9Egzs/edit#gid=0

# Motivation

Typical approaches for following the literature aren't that good at differentiating the signal from the noise.

Examples of signal:

- Papers in and around primary field (from any level of journal if it's sufficiently close, but needs to be more relevant as source IF decreases)
- Exciting papers of broad relevance to the field (e.g. cell engineering / cell therapies writ large)
- Papers by key people, such as personal connections and famous scientists

Examples of noise:

- Intractable numbers of papers to screen
- Papers in predatory journals or from obscure people/places
- Non research/review articles (e.g. news, notes, etc.)

## Pubmed search helpful tips

* use `hasabstract` to filter out news, notes, errata, etc.
* use `NOT review[pt]` to filter out reviews; see [HERE](https://pubmed.ncbi.nlm.nih.gov/help/#publication-types)
* use `last X years[dp]` to filter by recency; see [HERE](https://pubmed.ncbi.nlm.nih.gov/help/#filter-strategy-pubdate)

Also –

* Wildcards can be used, like `Nat Rev*[Journal]` for all Nat Rev family journals
* `[1au]` and `[lastau]` flags only work with truncated names like Bhargava HK, not with full name


