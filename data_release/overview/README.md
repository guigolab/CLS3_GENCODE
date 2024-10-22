# CLS Transcripts Overview

The separate projection of spliced and monoexonic transcript models resulted in a set of 468,598 unique intron chains and 57,709 monoexonic regions in human, and a total of 407,194 unique intron chains and 76,231 unspliced regions in mouse. with the term “CLS transcript models” we refer to the union of unique intron chains and monoexonic transcripts (526,307 
in human and 483,425 in mouse). The set of models yield through PacBio is almost a subset of those gained through ONT. More precisely, in human, 25.5% of the transcript models are shared between the two technologies (28% and 76% of the entire ONT and PacBio sets, respectively), while 28.3% are shared in mouse (31% and 74% of the entire ONT and PacBio sets, respectively). 
The overlap also increases when comparing the transcript models detected pre- and post-capture; in human, 76.4% of the models are detected post-capture, 19.3% are shared across design (45% and 25% of the total number of models obtained pre- and post-capture, respectively). Similarly, 67.5% of the transcript models are detected post-capture in mouse, 
and overall 24% are shared across pre- and post-capture samples (42% and 35% of the total number obtained pre- and post-capture, respectively). About 20% of the transcripts are detected in three or more tissues.

See a graphical representation for human:
![image](https://github.com/user-attachments/assets/932e9f24-f94d-4c98-b488-f12031deb1b5)

See a graphical representation for mouse:
![image](https://github.com/user-attachments/assets/7fe6fe57-7991-4d1f-8dc5-0b4bbcbda78f)

The plots have been produced with focus on intersection between ONT and PacBio, see script [here](https://github.com/guigolab/CLS3_GENCODE/blob/main/data_release/overview/get_upset.by_technology.R).
The same plot with focus on the intersection between developmental stages can be found [here](https://github.com/guigolab/CLS3_GENCODE/blob/main/data_release/overview/get_upset.by_stage.R).
