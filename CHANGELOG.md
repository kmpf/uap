## 1.2 (24.01.2020)

**Fixes**
 * completing slurm jobs are not considered running (#100)
 * invoke uap from anywhere (#114)
 * avoid deprecated pip wrapper (#122)
 * one parent step can now be used for multiple inputs (#125)
 * revised connection checks (#37, #138, #142)
 * update from broken sha1 to sha256 (#141)
 * required options stay required even if a default is set (#144)
 * fail cluster job if UAP fails (#60)

**Features**
 * steps now run within their temp directory (#29, #31)
 * use relative paths; uap output can now be moved (#115)
 * slurm jobs as array job per step (#105)
 * introduce --legacy option for submit-to-cluster to use none array jobs
 * error handling and --debugging option (#118)
 * platform info in annotation file (#62)
 * output hash now sensitive to tool versions (#116)
 * output hash robust against tool location and order
 * tool versions can optionally be ignored with `ignore_version`
 * shorter lmod config (#104)
 * uap tools must not be referenced in config (#119)
 * uap path in PATH environmental variable not required (#107)
 * default job quota is now 0, which is no quota (#105)
 * dependencies are now completed through selected connections (#127)
 * reference assembly `-G` is optional for `stringtieMerge` and `stringtie` (#124)
 * input connection for reference assembly in `stringtieMerge` or `stringtie`
 * object `ConnectionsCollector` for handling input from multiple steps (#47)
 * introduce --profiling option to analyze uap runtime (#132)
 * introduce --path option to retrieve the UAP installation path
 * optional step connections (#35)
 * improved single end support and sensitivity to respective connections (#38, #139)
 * forward None values for options to step declaration (#140)
 * raise exeption if an unknown configuration key is used
 * configure cluster default options (#76)
 * introduce --job-ids as option for status to report config specifc jobs (#58)
 * introduce --first-error option for submit-to-cluster for faster debugging (#61)
 * parallel tool check for major performance gain (#82)
 * revised run-locally message output and regulation with verbosity level
 * print stderr of failed processes in UAP log
 * watcher report with proc names printed on verbosity level 3 (-vv)
 * propagate verbosity level to cluster jobs

**additional stuff**
 * updated documentation and resolved sphinx warnings
 * `stringtieMerge` option `run_id` changed to `output_prefix`
 * step connection documentation (#137)
 * use common names for stringtie gtf connections `in/features`
 * hisat2 now takes `library_type` option
 * complete merge with https://github.com/yigbt/uap
 * chronological naming for log files
 * move ping files on job error or interruption (#147)

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
