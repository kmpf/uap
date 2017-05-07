..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: Tested Platforms

.. _platforms:

****************
Tested Platforms
****************

So far **uap** has been tested on several operating systems (OS) and cluster
engines.
The tables below list the combinations of OS and cluster engine we successfully
tested. 

+---------+------------+------------+--------------+---------------+
| ``download_human_gencode_release.yaml``                          |
+=========+============+============+==============+===============+
|         | |uge_link| | |oge_link| | |slurm_link| | local machine |
+---------+------------+------------+--------------+---------------+
| Cent OS | |check|    | untested   | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Fedora  | untested   | |check|    | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Ubuntu  | untested   | untested   | |check|      | |check|       |
+---------+------------+------------+--------------+---------------+

+---------+------------+------------+--------------+---------------+
| ``index_homo_sapiens_hg19_genome.yaml``                          |
+=========+============+============+==============+===============+
|         | |uge_link| | |oge_link| | |slurm_link| | local machine |
+---------+------------+------------+--------------+---------------+
| Cent OS | |check|    | untested   | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Fedora  | untested   | |check|    | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Ubuntu  | untested   | untested   | |check|      | |check|       |
+---------+------------+------------+--------------+---------------+

+---------+------------+------------+--------------+---------------+
| ``index_mycoplasma_genitalium_ASM2732v1_genome.yaml``            |
+=========+============+============+==============+===============+
|         | |uge_link| | |oge_link| | |slurm_link| | local machine |
+---------+------------+------------+--------------+---------------+
| Cent OS | |check|    | untested   | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Fedora  | untested   | |check|    | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Ubuntu  | untested   | untested   | |check|      | |check|       |
+---------+------------+------------+--------------+---------------+

+---------+------------+------------+--------------+---------------+
| ``2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml``         |
+=========+============+============+==============+===============+
|         | |uge_link| | |oge_link| | |slurm_link| | local machine |
+---------+------------+------------+--------------+---------------+
| Cent OS | |check|    | untested   | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Fedora  | untested   | |check|    | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Ubuntu  | untested   | untested   | |check|      | |check|       |
+---------+------------+------------+--------------+---------------+

+---------+------------+------------+--------------+---------------+
| ``2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml``                   |
+=========+============+============+==============+===============+
|         | |uge_link| | |oge_link| | |slurm_link| | local machine |
+---------+------------+------------+--------------+---------------+
| Cent OS | |check|    | untested   | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Fedora  | untested   | |check|    | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Ubuntu  | untested   | untested   | |check|      | |check|       |
+---------+------------+------------+--------------+---------------+

+---------+------------+------------+--------------+---------------+
| ``2014-RNA_CaptureSeq-Mercer_et_al_download.yaml``               |
+=========+============+============+==============+===============+
|         | |uge_link| | |oge_link| | |slurm_link| | local machine |
+---------+------------+------------+--------------+---------------+
| Cent OS | |check|    | untested   | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Fedora  | untested   | untested   | untested     | untested      |
+---------+------------+------------+--------------+---------------+
| Ubuntu  | untested   | untested   | |check|      | |check|       |
+---------+------------+------------+--------------+---------------+

+---------+------------+------------+--------------+---------------+
| ``2014-RNA_CaptureSeq-Mercer_et_al.yaml``                        |
+=========+============+============+==============+===============+
|         | |uge_link| | |oge_link| | |slurm_link| | local machine |
+---------+------------+------------+--------------+---------------+
| Cent OS | |check|    | untested   | untested     | |check|       |
+---------+------------+------------+--------------+---------------+
| Fedora  | untested   | untested   | untested     | untested      |
+---------+------------+------------+--------------+---------------+
| Ubuntu  | untested   | untested   | |check|      | |check|       |
+---------+------------+------------+--------------+---------------+

.. |check| unicode:: U+2713

.. |uge_link| raw:: html
 
   <a href="http://www.univa.com/products/" target="_blank">UGE</a>

.. |oge_link| raw:: html

   <a href="http://www.univa.com/oracle" target="_blank">OGE/SGE</a>

.. |slurm_link| raw:: html
      
   <a href="http://slurm.schedmd.com/" target="_blank">SLURM</a>
