#for f in TCR_Project/Results/*/TCRdist/cluster_metrics_R30.0/sample_level_cluster_metrics.tsv
# do
#    batch=$(basename $(dirname $(dirname $(dirname "$f"))))
#    scp "$f" xuliu@192.168.1.36:/hdd/xuliu/file36_home/TCRdata/TCR_item_full/new_characteristic/${batch}_cluster_metrics_R30.tsv
# done

# CI结果
for f in TCR_Project/Results/*/TCRdist/TCR_CI_radius_30.0.csv
 do
    batch=$(basename $(dirname  $(dirname "$f")))
    scp "$f" xuliu@192.168.1.36:/hdd/xuliu/file36_home/TCRdata/TCR_item_full/TCR_CI/${batch}_R30.csv
 done