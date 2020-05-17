The "expected_out" folder show the expected output folder from running the following command:

```
perl concur.pl --input demo/dom34_ncs2_elp6.YPD.rep1_100k.bam --genome sc3 --out demo/out
```

The demo bam files has been shortened to contain 100,000 reads. Analyzing it using CONCUR 
should take about 1.5 minutes to run using a single core (@ 3.60GHz) and would print the 
following message to the terminal:

```
Running CONCUR v0.9

### Parameters ###
Input file: demo/dom34_ncs2_elp6.YPD.rep1_100k.bam
Genome: sc3
Fragment size range tested: 20-50
Output folder: demo/out
Output file name: dom34_ncs2_elp6.YPD.rep1_100k
Run without R: FALSE
##################
[Step 1/10] Mapping genomic reads to transcripts...
[Step 2/10] Calculating periodicity...
[Step 3/10] Predicting frame per read length...
[Step 4/10] Calculating codon frequency per read length...
[Step 5/10] Calculating codon frequency (step 2) per read length...
[Step 6/10] Plotting codon per read length correlation...
null device 
          1 
null device 
          1 
[Step 7/10] Calculating correlation between read lengths...
[Step 8/10] Plotting codon per read length correlation...
null device 
          1 
null device 
          1 
[Step 9/10] Calculate final codon frequency...
[Step 10/10] Make final figures for validation...
null device 
          1 
```
