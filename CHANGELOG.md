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
 * no costly tool check if not needed (#82)
 * more stable `status` while jobs are running
 * make signal traps work on cluster
 * do not mix up lines in stderr tail of failed tasks
 * run-info uses list2cmdline for accurate command representation
 * run-locally checks parent task states befor running
 * include output file list to detect output changes

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
 * new connection information in command line `steps` doc
 * improved single end support and sensitivity to respective connections (#38, #139)
 * forward None values for options to step declaration (#140)
 * raise exeption if an unknown configuration key is used
 * configure cluster default options (#76)
 * introduce --job-ids as option for status to report config specifc jobs (#58)
 * introduce --first-error option for fix-problems for faster debugging (#61)
 * parallel tool check for major performance gain (#82)
 * revised message output and regulation with verbosity level
 * print stderr of failed processes in UAP log
 * watcher report with proc names printed on verbosity level 3 (-vv)
 * allow empty step options to set None values/unset defaults
 * io and net stats in annotation file (#34)
 * if step config has wrong options pass -v to display all available options
 * log comprehensible command structure in annotation file
 * new CHANGED state for tasks based on commands and tool versions
 * sha256 blockchain to determine CHANGED state
 * new BAD state for tasks when UAP caught an error
 * new VILATILIZED state if a step has been volatilized
 * status --details to view changes of tasks and errors of bad tasks
 * --force to overwrite changed tasks
 * --ignore to ignore changes in tasks and consider them finished
 * show status of some uap processes with tqdm
 * config can be changed without effecting submitted tasks
 * task error logged in annotation and reported in status --details
 * `base_working_directory` can be configured and paths set relatively to it
 * status bar during tool check and task status
 * improved resource monitoring report in annotation file
 * introduce status --hash to validate sha256sum of output files
 * introduce fix-problems --file-modification-date to fix output after copy

**additional stuff**
 * updated documentation and resolved sphinx warnings
 * `stringtieMerge` option `run_id` changed to `output_prefix`
 * step connection documentation (#137)
 * use common names for stringtie gtf connections `in/features`
 * hisat2 now takes `library_type` option
 * complete merge with https://github.com/yigbt/uap
 * chronological naming for log files
 * move ping files on job error or interruption (#147)
 * no hash in output directories
 * nothing is written outside of `destination_path` (no `config.yaml-out` link)
 * only one annotation and no symlinks in output directory
 * multiple executions do not accumulate files in output directory
 * log uap version in annotation
 * always use highest verbosity level on cluster jobs
 * run states overhaul
 * caching is now run specific

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
