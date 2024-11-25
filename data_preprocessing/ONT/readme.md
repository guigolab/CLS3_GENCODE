# Data Preprocessing and Quality Assessment

This preprocessing step uses the [split_on_adapter]([url](https://github.com/nanoporetech/duplex-tools/blob/master/fillet.md)) utility from [duplex_tools](https://github.com/nanoporetech/duplex-tools), to evaluate the presence of and further process the concatamers in the ONT sequencing data.

This utility was adjusted in accordance with the CapTrap-Seq adapter and primer sequences and is forked [here](https://github.com/Gazal90/duplex-tools). 

In addition to these changes, two rounds of “read splitting” were performed, followed by a final quality control step.

#### Round 1: 
In the first step, the complete set including the ONT adapter linked to the CapTrap-Seq primer was used to detect concatemers and split them into sub-reads. Due to presence of incomplete ONT adapter within some reads, these concatemers were left undetected. 
#### Round 2: 
Therefore, a second round of splitting was performed wherein the splitting was based on presence of only the CapTrap-Seq adapter within the reads. 
#### Round 3: 
Later, as a final step of this quality control, the multi-split reads either reported by the utility or split in both the consecutive rounds of splitting were discarded. 
Overall, the sequencing data was preprocessed to meet the necessary quality standards for reliable downstream analysis.

## 
<img width="1165" alt="image" src="https://github.com/user-attachments/assets/f6842bff-d4d1-43a2-a90c-52b9ac4485f6">

