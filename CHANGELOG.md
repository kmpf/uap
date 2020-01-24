## 1.2 (24.01.2020)

**Fixes**
 * completing slurm jobs are not considered running (#100)
 * invoce uap from anywhere (#114)
 * avoid depricated pip wrapper (#122)
 * all step input connections are used (#125)

**Features**
 * steps now run within their temp directory (#29, #31)
 * use relative paths; uap output can now be moved (#115)
 * slurm jobs as array per step (#105)
 * error handling and --debugging option (#118)
 * platform info in annotation file (#62)
 * output hash now sensitive to tool versions (#116)
 * shorter lmod config (#104)
 * uap tools musst not be referenced in config (#119)
 * uap path in PATH environmental variable not required (#107)
 * default job quota is 0 -> no quota (#105)
 * `_depends` amended implicitly through `_connect` (#127)
 * reference assembly is optional for `stringtieMerge` and `stringtie` (#124)
 * object `ConnectionsCollector` for handling input from multiple steps (#47)
 * pass reference assembly through `in/reference` connection to `stringtieMerge` or `stringtie`
 * `stringtieMerge` option `run_id` changed to `output_prefix`
 * introduce --profiling option to analyse uap runtime (#132)
 * introduce --legacy option for submit-to-cluster to use none array jobs

## 1.1 (20.01.2020)

**Fixed**
 * tests in uap_test repo (#113)
 * CI pipeline (#110)
 * automatic volatilization (#98)
 * fastqscreen: move html output files (#95)
 * removed option --optional in patched fastq_screen version (#94)
 * display correct uap version (#93)
 * deprecated warning from python package PyYAML (#91)
 * _cluster_job_quota is not read on slurm (#40)
 * fastq_screen: forgot to modify nohits option (#30)
 * fixed fastqscreen and rseqc file path issues (#120)

**Features**
 * tools sections defaults (#103)

**additional stuff**
 * fastq_screen is not running on ribnode018 (#97)
 * slurm cluster gives finished for failed runs (#60)
 * released documentation with gitlab pages (#117)
 * added uap_test as git submodule and modify gitlab ci process (#108)
