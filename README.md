# The holobiont *Agelas oroides* 


[![Sesam_img](https://www.imbe.fr/local/cache-vignettes/L400xH186/d46a112bebd61c35-0c5b6.png?1668533164)](https://sesam-anr.imbe.fr/)


**Sponges** are ancient sesile animals inhabiting almost all aquatic ecosystems in our planet. Sponges are well known to produce bioactive compounds that can be produced by sponge itself, by their associated microorngaims or by both. The concept of sponge holobiont gathers the host animal, its associated microorganisms and the biological processed occurring between this entities. Our objectif is tho describe the *Agelas oroides* holobiont by describing their cytological components and the microorganism assocaited with this specie. In the present study we  present the bioinformatic pipeline used to describe the microbial community associated to the spong *Agelas oroides*.

[![Agelas_img](https://inpn.mnhn.fr/photos/uploads/webtofs/inpn/3/139323.jpg){width=50% height=50%}](https://inpn.mnhn.fr/espece/cd_nom/71479)

**Figure 1.** *A. oroides* is a Mediterranean sponge inhabiting coraligenous formations between 12-60 m depth.



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

## Add your files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.osupytheas.fr/sesam/agelas-oroides-microbiome.git
git branch -M main
git push -uf origin main
```

## Integrate with your tools

- [ ] [Set up project integrations](https://gitlab.osupytheas.fr/sesam/agelas-oroides-microbiome/-/settings/integrations)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Set auto-merge](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

# Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!). Thank you to [makeareadme.com](https://www.makeareadme.com/) for this template.

## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
