# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


# curl -s http://togows.dbcls.jp/entry/pathway/hsa00010/genes.json \
#     | awk 'BEGIN{OFS="\t"}$0~/KO:/{
#         match($0,/"([^"]+)": "([^;[]+);?(.+?) \[(KO:[^\]]+)/,a);
#         a[3]=a[3]==""?"-":a[3];
#         print a[1],a[2],a[3],a[4]
#         }' \
#     | head -10

