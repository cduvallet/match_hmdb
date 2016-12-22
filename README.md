# blast_hmdb
blast mz's against HMDB database

much of this code is based off of Isaac Rockafellow's code, in public repo https://github.com/irockafe/search_hmdb.git

# data

`hmdb_metabolites.xml` was download from the HMDB website and contains all of
the HMDB metabolites (at least as of summer 2016).
When you download multiple HMDB metabolites, it just concatenates the individual
xml files into one file. However, this breaks downstream parsers. Isaac wrote
`remove_excess_xml_declarations.py` to clean up the downloaded file into one that
can be parsed by a parser. `hmdb_metabolites_clean.xml` is the cleaned up version,
and should be the default input HMDB file to `blast_hmdb.py`.

Check out the `toy_feature_table.csv` for an example of what your feature table
can look like.

You should gunzip the `hmdb_metabolites_clean.xml.gz` file before trying to use it.

# contributors
- Claire Duvallet (duvallet at mit dot edu)
- Isaac Rockafellow