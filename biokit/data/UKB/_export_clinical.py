# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


def build_fieldnames(field_df, fields):
    fieldnames = []
    for field in fields:
        field_info = field_df.loc[int(field)]
        instance_min, instance_max = field_info[['instance_min', 'instance_max']]
        array_min, array_max = field_info[['array_min', 'array_max']]

        if instance_max > 0:
            for instance in range(instance_min, instance_max + 1):
                if array_max > 0:
                    for array in range(array_min, array_max + 1):
                        fieldnames.append(f'p{field}_i{instance}_a{array}')
                else:
                    fieldnames.append(f'p{field}_i{instance}')
        else:
            fieldnames.append(f'p{field}')
    return fieldnames


def export_clinical_data(field_df, clinical_df, fields, output_path):
    """

    :param field_df:
    :param clinical_df:
    :param fields:
    :param output_path:
    :return:
    """

    cmd_templete = """
    #!/bin/bash

mkdir -p olink_single_protein

fail_count=0
total_count=0

> failed_batches.txt
for instance in 0 1 2; do
  echo "Processing instance $instance"
  for f in batch_*; do
    outfile="olink_single_protein/${f}_instance${instance}.tsv"

    # ======================
    # 如果文件存在则跳过
    # ======================
    if [ -f "$outfile" ]; then
        echo "Skip existing $outfile"
        continue
    fi

    tmp_fields=$(mktemp)

    echo "olink_instance_${instance}.eid" > "$tmp_fields"

    while read protein; do
        echo "olink_instance_${instance}.${protein}" >> "$tmp_fields"
    done < "$f"

    echo "Exporting $f (instance $instance)"

    ((total_count++))

    dx extract_dataset \
      app1005423_20251217031229.dataset \
      --entities olink_instance_${instance} \
      --fields-file "$tmp_fields" \
      -o "$outfile"

    status=$?

    if [ $status -ne 0 ]; then
        echo "FAILED $f instance $instance"
        echo "$f instance $instance" >> failed_batches.txt
        ((fail_count++))
    fi

    rm "$tmp_fields"

  done

done

echo "=============================="
echo "Total jobs: $total_count"
echo "Failed jobs: $fail_count"
echo "Failed list saved to failed_batches.txt"
echo "=============================="
    
    """
