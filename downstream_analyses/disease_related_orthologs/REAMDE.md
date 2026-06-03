# Analysis of disease related orthologs

## Data source:

###### LncRNADisease.v3.0:

We focuse on v.3.0 as it is the most recent and complete version of the database (check: Nucleic Acids Research, **2024**).
In general, website sux. Well formatted and most complete data can only be downloaded through the Search section (http://www.rnanut.net/lncrnadisease/index.php/home/search).
To download data I used advanced filter with specified settings:

* **ncRNA Category:** lncRNA
* **Species:** Homo sapiens
* **Causality:** Yes + Unknown (downloaded separately)
* **Validated Method:** Nothing selected
* **Score:** From: 0 To: 1

I renamed file to indicate used filters and saved using tsv extension in Raw_data directory. Something was wrong with 5th column. While reading with read_table() R function it was skipping ~ 20% of rows without any warning. Problem was caused by the 5th column (Disease_Name), thus it was removed and tables were subsequently joined together by simply:

`cat LncRNADisease.v3.0_search_lncRNA_hs_causlity-yes_score0-1.tsv | cut -f5 --complement > LncRNADisease.v3.0.tsv`
`tail +2 LncRNADisease.v3.0_search_lncRNA_hs_causlity-unknown_score0-1.tsv | cut -f5 --complement >> LncRNADisease.v3.0.tsv`

There are inconsistent hyphens present in the dataset. I had a good day so I fixed them as well (**in place**):

`sed -i $'s/\u2011/-/g' LncRNADisease.v3.0.tsv`    # 48 occurences
`sed -i $'s/\u2010/-/g' LncRNADisease.v3.0.tsv`    # 108 occurences
`sed -i $'s/\u2013/-/g' LncRNADisease.v3.0.tsv`    # 2 occurences
`sed -i $'s/\u00a0//g' LncRNADisease.v3.0.tsv`      # 2 occurences

Above mumbo-jumbo led to detect 2 more genes.

**Note:**

Forget about numbers on the website. They are not reproducible using any avaiable resources. Numbers of experimentally supported lncRNA/circRNA-disease Associations listed on the website sums up to 25,440 associacions in the `website_alldata.tsv`. But ratio of lncRNA/circRNA do not match. `website_alldata.tsv` dataset is chaotic and after matching with Gencode geneIDs leads fewer genes compared to described extraction method using Search section (2,194 vs 2,443 geneIDs). Forget about 6,066 lncRNAs listed on the website. There are only 5,540 in the `website_alldata.tsv` but more (6,042) if you join all datasets, for all species using Search section. Notice, that while downloading data using Search section you get tabular data with duplicated rows. If you exclude the header there are exacly twice as many. For LncRNADisease.v3.0_search_lncRNA_hs_causlity-yes_score0-1.tsv it's 10,649 vs 5,325.

###### Cancer lncRNA Census 3 (CLC3):

Described in the *Non-Coding RNA* **2022** article (https://doi.org/10.3390/ncrna8060082). Table in xlsx was downloaded from https://zenodo.org/records/7075104#.YyCB0C1Q3T8.
From xlsx table GENCODE_ID column was manually extracted and saved as `CLC3.tsv`.

###### LncBook 2.0:

Described in the Nucleic Acids Res. **2022** article (https://doi.org/10.1093/nar/gkac999).
To download only disease/trait-association genes the Genes section was used (https://ngdc.cncb.ac.cn/lncbook/genes).
Applied settings:

* **Variation (disease/trait-association):** yes

Downloaded `LncBook_GenesFilter.csv` file has problematic formatting. Column content is placed within "" and each column is comma-separated. Unfortunately, comma is also used within columns. Table was reformatted, where double quotation marks were removed and columns are now tab separated:

`cat LncBook_GenesFilter.csv | head -n 1 | sed 's/,/\t/g' > LncBook.v2.0.tsv`
`cat LncBook_GenesFilter.csv | tail +2 | sed 's/\",\"/\"|\"/g' | sed 's/\"//g' | awk -F'|' 'BEGIN {OFS="\t"} {$1=$1; print}' >> LncBook.v2.0.tsv`

**Note:**
Genes will be matched with Gencode geneIDs using Symbol column (when available).
