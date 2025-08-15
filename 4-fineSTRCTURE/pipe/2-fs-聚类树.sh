#!/bin/bash
#todo 请修改4-fineSTRCTURE/python/Fat_Initial_Fs.R 中的source进行修改


##Visulization of fs
/usr/bin/Rscript \
    /mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/python/Fat_Initial_Fs.R \
    /mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/整体1369_HP_linked

# /usr/bin/Rscript \
#     /mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/python/2-fs_xml→nwk.R \
#     /mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/整体1369_HP_linked_tree.xml \
#     /mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/整体1369_HP_linked_tree.nwk


# ###样本量过大
# ulimit -s unlimited 
# /usr/bin/Rscript /home/biosoftware/ppgv1/R_packages/Fat_Initial_Fs.R $2