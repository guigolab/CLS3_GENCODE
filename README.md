# GENCODE CLS Project
## GENCODE: massively expanding the lncRNA catalog through capture long-read RNA sequencing
Tamara Perteghella<sup>1,2,*\*</sup>, Gazaldeep Kaur<sup>1,*\*</sup>,  Sílvia Carbonell-Sala<sup>1,*\*</sup>, Jose Gonzalez-Martinez<sup>3,*\*</sup>, Toby Hunt<sup>3,*\*</sup>, Tomasz Mądry<sup>4\*</sup>, Irwin Jungreis<sup>5,6\*</sup>, Fabien Degalez<sup>1\*</sup>, Carme Arnan<sup>1\*</sup>, Ramil Nurtdinov<sup>1\*</sup>, Julien Lagarde<sup>1,7\*</sup>, Beatrice Borsari<sup>8,9\*</sup>, Cristina Sisu<sup>10\*</sup>, Yunzhe Jiang<sup>8,9\*</sup>, Ruth Bennett<sup>3\*</sup>, Andrew Berry<sup>3\*</sup>, Marta Blangiewicz<sup>4\*</sup>, Daniel Cerdán-Vélez<sup>11\*</sup>, Kelly Cochran<sup>12\*</sup>, Covadonga Vara<sup>13\*</sup>, Claire Davidson<sup>3\*</sup>, Sarah Donaldson<sup>3\*</sup>, Cagatay Dursun<sup>8,9\*</sup>, Silvia González-López<sup>1,2\*</sup>, Sasti Gopal Das<sup>4\*</sup>, Kathryn Lawrence<sup>14\*</sup>, Daniel Nachun<sup>14\*</sup>, Matthew Hardy<sup>3\*</sup>, Zoe Hollis<sup>3\*</sup>, Mike Kay<sup>3\*</sup>, José Carlos Montañés<sup>13\*</sup>, Pengyu Ni<sup>8,9\*</sup>, Emilio Palumbo<sup>1\*</sup>, Carlos Pulido-Quetglas<sup>15,16\*</sup>, Marie-Marthe Suner<sup>3\*</sup>, Xuezhu Yu<sup>8,9\*</sup>, Dingyao Zhang<sup>8,9\*</sup>, Francois Aguet<sup>6\*</sup>, Kristin Ardlie<sup>6\*</sup>, Stephen B. Montgomery<sup>14,17,18\*</sup>, Jane E. Loveland<sup>3\*</sup>, M. Mar Albà<sup>13,19\*</sup>, Mark Diekhans<sup>20\*</sup>, Andrea Tanzer<sup>21\*</sup>, Jonathan M. Mudge<sup>3\*</sup>, Paul Flicek<sup>3\*</sup>, Fergal J Martin<sup>3\*</sup>, Mark Gerstein<sup>8,9\*</sup>, Manolis Kellis<sup>5,6\*</sup>, Anshul Kundaje<sup>12,14\*</sup>, Benedict Paten<sup>20\*</sup>, Michael L. Tress<sup>11\*</sup>, Rory Johnson<sup>15,16\*</sup>, Barbara Uszczynska-Ratajczak<sup>4\*</sup>, Adam Frankish<sup>3\*</sup>, Roderic Guigó<sup>1,2\*</sup>
 
    1. Centre for Genomic Regulation (CRG), The Barcelona Institute of Science and Technology, Dr. Aiguader 88, Barcelona 08003, Catalonia, Spain.
    2. Departament de Ciències Experimentals i de la Salut, Universitat Pompeu Fabra (UPF).
    3. European Molecular Biology Laboratory, European Bioinformatics Institute, Wellcome Genome Campus, Hinxton, Cambridge CB10 1SD, UK.
    4. Department of Computational Biology of Noncoding RNA, Institute of Bioorganic Chemistry, Polish Academy of Sciences, Noskowskiego 12/14, 61-704 Poznan, Poland.
    5. Computer Science and Artificial Intelligence Lab, Massachusetts Institute of Technology, 32 Vassar St, Cambridge, MA 02139, USA.
    6. The Broad Institute of MIT and Harvard, 415 Main Street, Cambridge, MA 02142, USA.
    7. Flomics Biotech, SL, Carrer de Roc Boronat 31, 08005 Barcelona, Catalonia, Spain.
    8. Program in Computational Biology and Bioinformatics, Yale University, New Haven, Connecticut 06520, USA.
    9. Department of Molecular Biophysics and Biochemistry, Yale University, New Haven, Connecticut 06520, USA.
    10. Department of Life Sciences, Brunel University London, Uxbridge, London, UB8 3PH, UK.
    11. Bioinformatics Unit, Spanish National Cancer Research Centre (CNIO), Calle Melchor Fernandez Almagro, 3, 28029 Madrid, Spain.
    12. Department of Computer Science, Stanford University, Stanford, CA, USA.
    13. Hospital del Mar Research Institute, Dr. Aiguader 88, Barcelona 08003, Spain.
    14. Department of Genetics, Stanford University School of Medicine, Stanford, CA, USA.
    15. Department of Medical Oncology, Bern University Hospital, Murtenstrasse 35, 3008 Bern, Switzerland.
    16. School of Biology and Environmental Science, University College Dublin, University College Dublin, Belfield, Dublin 4, D04 V1W8, Ireland.
    17. Department of Pathology, Stanford University School of Medicine, Stanford, CA, USA
    18. Department of Biomedical Data Science, Stanford University School of Medicine, Stanford, CA, USA
    19. Catalan Institute for Research and Advanced Studies (ICREA), Barcelona, Spain.
    20. UC Santa Cruz Genomics Institute, 2300 Delaware Avenue, University of California, Santa Cruz, CA 95060, USA.
    21. University of Vienna, Department of Biochemistry and Cell Biology, Vienna, Austria
    
    * Equal contribution
    Correspondence should be addressed to R.G. (roderic.guigo@crg.cat)

GENCODE is a 20-year international project focused on producing high-quality annotations for human and mouse genomes, crucial for understanding gene function. While the human gene catalog for protein-coding genes is nearly complete, long non-coding RNA (lncRNA) annotations have remained inconsistent across different catalogs. To address this, GENCODE used targeted RNA sequencing to unify and expand lncRNA annotations in human and mouse, employing full-length sequencing across diverse tissues. This effort resulted in 16,817 new human genes and 22,210 new mouse genes, significantly increasing the lncRNA catalog and improving orthology mapping between species. These new annotations enhance the functional interpretation of genome data, linking previously unannotated regions to biological functions.

In this repository:

### [Data Preprocessing and Quality Assessment](https://github.com/guigolab/CLS3_GENCODE/tree/main/data_preprocessing)
Summary of the steps taken to process long-read data, upon sequencing but prior to LyRic. Measures undertaken to assess the quality of the data prior to downstream processing are also detailed here.
### [Data Release](https://github.com/guigolab/CLS3_GENCODE/tree/main/data_release)
List of the files used in this work, complemented with descriptions of the steps taken to generate them, links to direct download, and detailed information about formats and tags.
### [Complementary Data](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data)
Datasets used in this work, complemented with useful information regarding the files and their processing prior to analysis.
### [Downstream Analyses](https://github.com/guigolab/CLS3_GENCODE/tree/main/downstream_analyses)
Codes used in various downstream analyses.
