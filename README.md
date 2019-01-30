## Code for Balancing sensitivity and specificity in distinguishing TCR groups by CDR sequence similarity 

Link to preprint: https://www.biorxiv.org/content/10.1101/526467v1

### General notes and dependencies
- All of the code is run using python 2.7
- You will need the [swalign package](https://pypi.python.org/pypi/swalign), which can be acquired with `pip install swalign` 

### Data pre-processing
- All data pre-processing scripts are in the directory data_preprocessing, which has further instructions on how to acquire and process the data
- Once you run the data pre-processing scripts, you will need to produce the outputs in a data directory and know the full path to the directory. Our example of this directory is called `datasets/`
- Some of the data (ex. McPAS) does not need to be pre-processed, and the script to generate results is called directly on the data

### Generating figure 1
To create figure 1, we picked some CDR sequences and called the function `getDistanceSW(cdr1, cdr2)` in `sw_scoring.py` to get the distances between them, which we then approximately visualized.

### Generating figure 2
To generate the density estimate plots, as well as cumulative density estimates, call:

`python structural_cluster_histograms.py <output_dir> <file_title> <list of repertoire files with full path>`

We used these commands to generate the alpha and beta plots:

```
python structural_cluster_histograms.py results/stcr_results/ STCR_global_beta "datasets/stcr_data/B3-10_11_12_13-A,datasets/stcr_data/B3-10_11-A,datasets/stcr_data/B3-12_A,datasets/stcr_data/B3-12-B,datasets/stcr_data/B3-13_14-A,datasets/stcr_data/B3-14-A"
```

```
python structural_cluster_histograms.py results/stcr_results/ STCR_global_alpha "datasets/stcr_data/A3-10-A,datasets/stcr_data/A3-13-B,datasets/stcr_data/A3-13-A,datasets/stcr_data/A3-10_11_12-A_only_10"
```


To generate threshold plots:
`python run_threshold_analysis.py  <output_dir> <file_title> <list of repertoire files with full path>`

We used these commands to generate the alpha and beta plots:

```
python run_threshold_analysis.py results/stcr_results/ STCR_global_alpha "datasets/stcr_data/A3-10-A,datasets/stcr_data/A3-13-B,datasets/stcr_data/A3-13-A,datasets/stcr_data/A3-10_11_12-A_only_10"
```

```
python run_threshold_analysis.py results/stcr_results/ STCR_global_beta "datasets/stcr_data/B3-10_11_12_13-A,datasets/stcr_data/B3-10_11-A,datasets/stcr_data/B3-12_A,datasets/stcr_data/B3-12-B,datasets/stcr_data/B3-13_14-A,datasets/stcr_data/B3-14-A"
```


### Generating figure 3

To generate the within and between cluster closest pairs, output as a csv file: 
`python stcr_pairs.py <output_dir> <a3 or b3> <list of repertoire files with full path>`

We used these commands to generate the csv files:

```
python stcr_pairs.py results/stcr_results/ a3 "datasets/stcr_data/A3-10-A,datasets/stcr_data/A3-13-B,datasets/stcr_data/A3-13-A,datasets/stcr_data/A3-10_11_12-A_only_10"
```

```
python stcr_pairs.py results/stcr_results/ b3 "datasets/stcr_data/B3-10_11_12_13-A,datasets/stcr_data/B3-10_11-A,datasets/stcr_data/B3-12_A,datasets/stcr_data/B3-12-B,datasets/stcr_data/B3-13_14-A,datasets/stcr_data/B3-14-A"
```

These results were plotted with `seq_vs_struct_similarity.py`.


### Generating figure 4
To generate threshold plots for all of the pairs of twins, call:
`python twins_threshold_analysis.py <output_dir> <list of all of the twin repertoire files with full path>`

We used these commands to generate the twin threshold plots, for the most frequent 1000 sequences in the alpha and beta repertoires:

```
python twins_threshold_analysis.py results/twin_results/ "datasets/twins_data/TwA1_B_top_1000.txt,datasets/twins_data/TwA2_B_top_1000.txt,datasets/twins_data/TwC1_B_top_1000.txt,datasets/twins_data/TwC2_B_top_1000.txt,datasets/twins_data/TwD1_B_top_1000.txt,datasets/twins_data/TwD2_B_top_1000.txt"
```

```
python twins_threshold_analysis.py results/twin_results/ "datasets/twins_data/TwA1_A_top_1000.txt,datasets/twins_data/TwA2_A_top_1000.txt,datasets/twins_data/TwC1_A_top_1000.txt,datasets/twins_data/TwC2_A_top_1000.txt,datasets/twins_data/TwD1_A_top_1000.txt,datasets/twins_data/TwD2_A_top_1000.txt"
```

### Generating figure 5

To extract all of the clusters in one individual, generate motifs from the clusters, and score those motifs against the clusters in the other individuals using our distance method, call:
`python twins_distance_scores.py <output_dir> <clustering_distance> <individual: TwA1, TwA2, TwC1, TwC2, TwD1, TwD2> <individuals's twin> <list of all of the twin repertoire files with full path>`

To extract all of the clusters from individual A1's CDR3b sequences, clustered at a distance of 0.3, we used:

```
python twins_distance_scores.py results/twin_results/motifs_distance_scoring/ 0.3 TwA1 TwA2 "datasets/twins_data/TwA1_A_top_1000.txt,datasets/twins_data/TwA2_A_top_1000.txt,datasets/twins_data/TwC1_A_top_1000.txt,datasets/twins_data/TwC2_A_top_1000.txt,datasets/twins_data/TwD1_A_top_1000.txt,datasets/twins_data/TwD2_A_top_1000.txt"
```

For this script to work, the files must be named similarly to ours, with the identifier "TwA1, TwA2, TwC1, TwC2, TwD1, or TwD2" in the title of the file.

To layout the motifs as we did in figure 5, use `layout-twin-motifs`. 


### Generating figure 6

To generate threshold plots:
`python run_threshold_analysis.py  <output_dir> <file_title> <list of repertoire files with full path>`

We used the commands:

```
python run_threshold_analysis.py results/dash_results/ Dash_beta_mice "datasets/dash_data/NP_beta.txt,datasets/dash_data/PA_beta.txt,datasets/dash_data/PB1_beta.txt,datasets/dash_data/PB1-F2_beta.txt"
```

```
python run_threshold_analysis.py results/dash_results/ Dash_beta_human "datasets/dash_data/BMLF1_beta.txt,datasets/dash_data/M1_beta.txt,datasets/dash_data/M38_beta.txt,datasets/dash_data/M45_beta.txt,datasets/dash_data/m139_beta.txt,datasets/dash_data/p65_beta.txt"
```

```
python run_threshold_analysis.py results/dash_results/ Dash_alpha_mice "datasets/dash_data/NP_alpha.txt,datasets/dash_data/PA_alpha.txt,datasets/dash_data/PB1_alpha.txt,datasets/dash_data/PB1-F2_alpha.txt"
```

```
python run_threshold_analysis.py results/dash_results/ Dash_alpha_human "datasets/dash_data/BMLF1_alpha.txt,datasets/dash_data/M1_alpha.txt,datasets/dash_data/M38_alpha.txt,datasets/dash_data/M45_alpha.txt,datasets/dash_data/m139_alpha.txt,datasets/dash_data/p65_alpha.txt"
```

### Generating figures 7 and 8

To generate the motifs, with distance scores against all other human/mice specific repertoires and cluster progressions as the threshold is varied from 0.2-0.4:
`python dash_distance_scores.py <output_dir> <distance to cluster to> <True if human, False if mice> <epitope name>`

We used:

```
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 True BMLF1    
```  
                                                                                           
```                                                                                             
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 True p65      
```          
                                                                                   
```                                                                                             
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 True m139     
```   
                                                                                          
```                                                                                             
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 True M1       
``` 
                                                                                                                                                                                                                                                                                 
```                                                                                             
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False M38     
```  
                                                                                           
```                                                                                             
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False M45     
```                                                                                             
```                                                                                             
python2 dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False NP      
```                                                                                             
```                                                                                             
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False PA      
```  
                                                                                           
```                                                                                             
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False PB1     
```   
                                                                                          
```                                                                                             
python dash_distance_scores.py results/dash_results/motifs_distance_scoring/ 0.2 False PB1-F2  
```                                                                                                                                                                                          
                                                                                                
### Generating figure 9

To generate the results spreadsheets, use `glanville_threshold_analysis.py`. Note that the spreadsheet needs to be converted to CSV before running our script.


We used:

```
python glanville_threshold_analysis.py datasets/glanville_data/glanville_s2_curated_data.csv results/glanville_results/ glanville_antigen
```

To plot the data, run the script `plot-glanville-antigen.R`.


### Generating figure 10

You can acquire the raw data file at http://friedmanlab.weizmann.ac.il/McPAS-TCR/ -> download the complete database.

To generate the results spreadsheets, use: 

`python mcpas_threshold_analysis.py <McPAS raw data file path> <output_dir> <small results filename> <large results filename>`

We used:

```
python mcpas_threshold_analysis.py datasets/mcpas_data/McPAS-TCR.csv results/mcpas_results/ small large
```

To plot the data, run the script `plot-mcpas.R`.



