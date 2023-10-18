## 1. Raw Fastq treatment

for data pre-treatment we follow [VTAM](https://vtam.readthedocs.io/en/latest/) pipeline

## Merge fastq files
~~~
conda activate vtam

vtam merge --fastqinfo ses1_16S1/fastqinfo_SES1_16S1.txt --fastqdir ses1_16S1/fastq --fastainfo ses1_16S1/run_16S1/fastainfo_16S1_Output.tsv --fastadir ses1_16S1/run_16S1/merged -v --log ses1_16S1/vtam.log
~~~

## random_seq: Create a smaller randomized dataset from the main dataset (using Pearl script)
~~~
perl Pearl/random_seq.pl -fastainfo run_16S1/fastainfo_Cesar.tsv -fastadir run_16S1/
merged -random_seqdir run_16S1/randomizedPearl -random_seqinfo  run_16S1/randomSeqInfoPearl.tsv -samplesize 10000000

$VAR1 = {
          'fastadir' => 'run_16S1/merged',
          'fastainfo' => 'run_16S1/fastainfo_Cesar.tsv',
          'random_seqdir' => 'run_16S1/randomizedPearl',
          'samplesize' => '10000000',
          'random_seqinfo' => 'run_16S1/randomSeqInfoPearl.tsv'
        };
~~~

## Sort reads: Demultiplex and trim the reads
~~~
vtam sortreads --fastainfo run_16S1/randomSeqInfoPearl.tsv --fastadir run_16S1/randomizedPearl --sorteddir ses1-16S1/run_16S1/sortedPearl -v --log ses1-16S1/vtam.log
~~~

## Filter 1 - Filter variants and create the ASV table (using lfn_variant_replicate)

**params_filter.yml**
~~~
global_read_count_cutoff: 10
skip_filter_indel: 1
skip_filter_codon_stop: 1
~~~
~~~
vvtam filter --db db_ses_var_repl.sqlite --sortedinfo sortedPearl/sortedinfo.tsv --sorteddir sortedPearl --asvtable filter3/asvtable_default_variant_repl.tsv --params user_input/params_filter.yml -v --lfn_variant_replicate
~~~

## Make known occurrences - Create file containing the known_occurences.tsv to be used as an inut for *optimize*
~~~
vtam make_known_occurrences --asvtable filter3/asvtable_default_variant_repl.tsv --sample_types user_input/sample_types.tsv --mock_composition user_input/mock_ses.tsv --known_occurrences filter1/known_occurrences.tsv --missing_occurrences filter3/missing_occurrences3.tsv -v 
~~~

## Optimize - Compute optimal filter parameters based on mock and negative samples

**params_optimize.yml**
~~~
global_read_count_cutoff: 10
lfn_sample_replicate_cutoff: 0.001
pcr_error_var_prop: 0.1
lfn_variant_replicate_cutoff: 0.001
lfn_read_count_cutoff: 10
~~~

~~~
vtam optimize --db db_ses_var_repl.sqlite --sortedinfo sortedPearl/sortedinfo.tsv --sorteddir sortedPearl --known_occurrences filter3/known_occurrences.tsv --outdir optimize3 -v --params user_input/params_optimize.yml --lfn_variant_replicate
~~~

**params_optimize_b.yml**

~~~
global_read_count_cutoff: 10
lfn_sample_replicate_cutoff: 0.002
pcr_error_var_prop: 0.1
lfn_variant_replicate_cutoff: 0.001
lfn_read_count_cutoff: 10
~~~

Delete Variant 130248 from known_occurrences.tsv => known_occurrences_WO_conta.tsv

~~~
nohup vtam optimize --db db_ses_var_repl.sqlite --sortedinfo sortedPearl/sortedinfo.tsv --sorteddir sortedPearl --known_occurrences filter3/known_occurrences_WO_conta.tsv --outdir optimize4 -v --params user_input/params_optimize_b.yml --lfn_variant_replicate --until OptimizeLFNreadCountAndLFNvariant &
~~~
## Filter 2 - Create an ASV table with optimal parameters and assign variants to taxa

**params_filter2.yml**
~~~
global_read_count_cutoff: 10
skip_filter_indel: 1
skip_filter_codon_stop: 1
lfn_sample_replicate_cutoff: 0.002
pcr_error_var_prop: 0.1
lfn_variant_replicate_cutoff: 0.016
lfn_read_count_cutoff: 70
~~~

~~~
nohup vtam filter --db db_ses_var_repl.sqlite --sortedinfo sortedPearl/sortedinfo.tsv --sorteddir sortedPearl --asvtable filter4/asvtable_optimized_variant_repl.tsv --params user_input/params_filter4.yml -v --lfn_variant_replicate 
~~~
make_known_occurrences to count the nb of False Positives
~~~
vtam make_known_occurrences --asvtable filter4/asvtable_optimized_variant_repl.tsv --sample_types user_input/sample_types.tsv --mock_composition user_input/mock_ses.tsv --known_occurrences filter4/known_occurrences4.tsv --missing_occurrences filter4/missing_occurrences4.tsv -v
~~~
Contamination in negative controls and fals positive where eliminated from de dataset


## Taxonomic assignation
Starting with (asvtable_optimized_variant_repl.tsv)
Convert to fastafile -> asvtable_optimizedFASTA
Taxa assign with online RPD classifier (http://rdp.cme.msu.edu/classifier/class_help.jsp#conf) accesed in November 2022
