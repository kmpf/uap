..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: Tested Platforms

.. _platforms:

****************
Tested Platforms
****************

So far **uap** has been tested on several operating systems and cluster engines.
The table below lists the combinations we successfully tested the example
configurations on.

+-----------------------------+----------------------------------------+
| Operating Systems           | Cluster Engines                        |
+---------+---------+---------+------------+------------+--------------+
| Cent OS | Fedora  | Ubuntu  | |uge_link| | |oge_link| | |slurm_link| |
+=========+=========+=========+============+============+==============+
| download_human_gencode_release.yaml                                  |
| :ref:`example_download_gencode`                                      |
+---------+---------+---------+------------+------------+--------------+
| |check| | |check| | |check| | |check|    | |check|    | |check|      |
+---------+---------+---------+------------+------------+--------------+
| index_homo_sapiens_hg19_genome.yaml                                  |
| :ref:`example_index_hg19`                                            |
+---------+---------+---------+------------+------------+--------------+
| |check| | |check| | |check| | |check|    | |check|    | |check|      |
+---------+---------+---------+------------+------------+--------------+
| index_mycoplasma_genitalium_ASM2732v1_genome.yaml                    |
| :ref:`example_index_mycoplasma`                                      |
+---------+---------+---------+------------+------------+--------------+
| |check| | |check| | |check| | |check|    | |check|    | |check|      |
+---------+---------+---------+------------+------------+--------------+
| 2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml                  |
| :ref:`example_barski_download`                                       |
+---------+---------+---------+------------+------------+--------------+
| |check| | |check| | |check| | |check|    | |check|    | |check|      |
+---------+---------+---------+------------+------------+--------------+
| 2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml                           |
| :ref:`example_barski`                                                |
+---------+---------+---------+------------+------------+--------------+
| |check| | |check| | |check| | |check|    | |check|    | |check|      |
+---------+---------+---------+------------+------------+--------------+
| 2014-RNA_CaptureSeq-Mercer_et_al_download.yaml                       |
| :ref:`example_mercer_download`                                       |
+---------+---------+---------+------------+------------+--------------+
| |check| | |check| | |check| | |check|    | |check|    | |check|      |
+---------+---------+---------+------------+------------+--------------+
| 2014-RNA_CaptureSeq-Mercer_et_al.yaml                                |
| :ref:`example_mercer`                                                |
+---------+---------+---------+------------+------------+--------------+
| |check| | |check| | |check| | |check|    | |check|    | |check|      |
+---------+---------+---------+------------+------------+--------------+


.. |check| unicode:: U+2713

.. |uge_link| raw:: html
 
   <a href="http://www.univa.com/products/" target="_blank">UGE</a>

.. |oge_link| raw:: html

   <a href="http://www.univa.com/oracle" target="_blank">OGE/SGE</a>

.. |slurm_link| raw:: html
      
   <a href="http://slurm.schedmd.com/" target="_blank">SLURM</a>
