# blast_hmdb

Script to blast mz's against HMDB database.

In other words, given a feature table with a column labeled `mz`, find all HMDB features within a certain
ppm tolerance for each of the mz's in that column. If a feature doesn't have a mass match in
HMDB, the output file will have one row corresponding to that mz, with the feature name, adduct type,
and corresponding neutral mass.

Most of the parsing I did in this code is based off of Isaac Rockafellow's code, which you
can find on his [github](https://github.com/irockafe/search_hmdb.git).

# Usage

1. Clone this repo: ```git clone https://github.mit.edu/duvallet/blast_hmdb.git```

2. Gunzip the `hmdb_metabolites_clean.xml.gz` file: ```gunzip data/hmdb_metabolites_clean.xml.gz```

3. Run the code:

   ```python blast_hmdb.py -x data/hmdb_metabolites_clean.xml -f data/toy_feature_table.csv -s , -p 2 -o tmp.out```

   Where the `-f` flag corresponds to the feature table with the mz's you want to search against HMDB.

4. See more options for running the code:

   ```python blast_hmdb.py --help```

# Input data

`hmdb_metabolites.xml` was download from the HMDB website and contains all of
the HMDB metabolites (at least as of summer 2016).
When you download multiple HMDB metabolites, it just concatenates the individual
xml files into one file. However, this breaks downstream parsers. Isaac wrote
`remove_excess_xml_declarations.py` to clean up the downloaded file into one that
can be parsed by a parser. You only need to use this if you want to use a newly
downloaded xml file from HMDB. Otherwise, you can just use `hmdb_metabolites_clean.xml`,
which Isaac already downloaded and cleaned up for us.

`hmdb_metabolites_clean.xml` is the cleaned up version containing all HMDB metabolites,
and should be the default input HMDB file to `blast_hmdb.py`.

Check out the `toy_feature_table.csv` for an example of what your feature table
can look like.

You should gunzip the `hmdb_metabolites_clean.xml.gz` file before trying to use it.

# To do

* store the parsed HMDB dict as a pickle file, so users don't have to re-calculate it every time the script is run

# Contributors
- Claire Duvallet (duvallet at mit dot edu)
- Isaac Rockafellow