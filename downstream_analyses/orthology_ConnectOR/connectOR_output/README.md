# ConnectOR_output

#### Description

This directory contains output files from two separate runs of ConnectOR
([https://github.com/cobRNA/ConnectOR-optimized](https://github.com/cobRNA/ConnectOR-optimized)).

1. Analysis performed using human GENCODE v27 and mouse GENCODE vM16 annotations.
   Only main chromosomes were considered, excluding `chrM`.
2. Analysis performed using human GENCODE v47 and mouse GENCODE vM36 annotations.
   Only main chromosomes were considered, excluding `chrM`.

Additionally, three analysis-specific configuration files are included:

1. `dictionaries.json` — contains links to the required chain files;
2. `config_content_for_v27-vM16_analysis` — snapshot of the ConnectOR `config` file used for the analysis based on GENCODE v27 and vM16 annotations. Before use, the annotation paths in the `annotation` column must be updated
3. `config_content_for_v47-vM36_analysis` — snapshot of the ConnectOR `config` file used for the analysis based on GENCODE v47 and vM36 annotations. Before use, the annotation paths in the `annotation` column must be updated
