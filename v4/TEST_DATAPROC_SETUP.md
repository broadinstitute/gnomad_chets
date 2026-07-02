# Test Intervals and Dataproc Commands

## Test Intervals

```python
TEST_INTERVALS = {
    "PCNT": "chr21:46324141-46445769",
    "COL18A1": "chr21:45405123-45513721",
    "AHNAK2": "chr14:104937244-104978374",
    "TTN": "chr2:178525989-178830802",
    "FLG": "chr1:152302165-152325239",
    "OBSCN": "chr1:228208044-228378876",
    "HRNR": "chr1:152212076-152224193",
    "NBPF10": "chr1:146064711-146229000",
    "PLEC": "chr8:143915153-143976734",
    "PDE4DIP": "chr1:148808181-149048286",
    "FCGBP": "chr19:39863323-39934626",
    "NEB": "chr2:151485336-151734487",
    "LAMA5": "chr20:62307955-62367312",
    "SYNE1": "chr6:152121687-152637801",
}
```

## Execution Time Summary

| Step | Description | Variants/Variant Pairs | Partitions | Time (seconds) | Time (minutes) | Time (hours) |
|------|-------------|------------------------|------------|----------------|----------------|--------------|
| 1 | Create Variant Filter HT and Filtered VMT | 195,937 variants | 202 | 602.26 | ~10.0 | ~0.17 |
| 2 | Create Variant Pair List HT | 21,210,652 variant pairs | 3,559 | 5412.37 | ~90.2 | ~1.50 |
| 3 | Create Dense Filtered MT | 195,628 variants | 368 | 4812.08 | ~80.2 | ~1.34 |
| 4 | Create Variant Pair Genotype HT and Counts HT | - | - | - | - | - |
| **Total** | Steps 1-3 | - | - | **10826.71** | **~180.4** | **~3.01** |

## Cluster Setup

### Hail 0.2.120 Cluster (test-chets-hail120)

```bash
conda activate hail-2.120

hailctl dataproc start test-chets-hail120 \
    --requester-pays-allow-all \
    --packages="git+https://github.com/broadinstitute/gnomad_methods.git@main","git+https://github.com/broadinstitute/gnomad_qc.git@main" \
    --no-off-heap-memory \
    --autoscaling-policy=max-20 \
    --max-idle 540m \
    --properties='dataproc:dataproc.logging.stackdriver.enable=false,dataproc:diagnostic.capture.enabled=false,dataproc:dataproc.logging.syslog.enabled=false,dataproc:dataproc.logging.extended.enabled=false' \
    --project broad-mpg-gnomad \
    --labels gnomad_release=gnomad_v4_1,chets_run=chets_test_large_vp_genes_01_27_2026_3_cluster_split
```

### Hail 0.2.134 High-Memory Cluster (test-chets-highmem)

```bash
conda activate hail-2.134

hailctl dataproc start test-chets-highmem \
    --requester-pays-allow-all \
    --packages="git+https://github.com/broadinstitute/gnomad_methods.git@main","git+https://github.com/broadinstitute/gnomad_qc.git@main" \
    --no-off-heap-memory \
    --autoscaling-policy=max-20 \
    --max-idle 540m \
    --master-machine-type=n1-highmem-16 \
    --worker-machine-type=n1-highmem-8 \
    --properties='dataproc:dataproc.logging.stackdriver.enable=false,dataproc:diagnostic.capture.enabled=false,dataproc:dataproc.logging.syslog.enabled=false,dataproc:dataproc.logging.extended.enabled=false' \
    --project broad-mpg-gnomad \
    --labels gnomad_release=gnomad_v4_1,chets_run=chets_test_large_vp_genes_01_27_2026_3_cluster_splitm
```

## Pipeline Commands (Execution Order)

### Step 1: Create Variant Filter HT and Filtered VMT

**Cluster:** `test-chets-hail120` (Hail 0.2.120)

```bash
hailctl dataproc submit test-chets-hail120 gnomad_chets/v4/create_vp_matrix.py \
    --create-variant-filter-ht \
    --filter-vmt \
    --test \
    --output-postfix julia_lg_vp_test \
    --overwrite \
    --pyfiles gnomad_chets
```

**Output:**

```
Submitting to cluster 'test-chets-hail120'...
gcloud command:
gcloud dataproc jobs submit pyspark gnomad_chets/v4/create_vp_matrix.py \
    --files= \
    --py-files=/var/folders/r8/f581hggx4_n6rmtz01r7nljm0000gq/T/pyscripts_k7yntrue.zip \
    --properties= \
    -- \
    --create-variant-filter-ht \
    --filter-vmt \
    --test \
    --output-postfix \
    julia_lg_vp_test \
    --overwrite
Job [a5c94a41012c48dda388667cda4da3c9] submitted.
Waiting for job output...
SLF4J: No SLF4J providers were found.
SLF4J: Defaulting to no-operation (NOP) logger implementation
SLF4J: See https://www.slf4j.org/codes.html#noProviders for further details.
SLF4J: Class path contains SLF4J bindings targeting slf4j-api versions 1.7.x or earlier.
SLF4J: Ignoring binding found at [jar:file:/usr/lib/spark/jars/log4j-slf4j-impl-2.17.2.jar!/org/slf4j/impl/StaticLoggerBinder.class]
SLF4J: See https://www.slf4j.org/codes.html#ignoredBindings for an explanation.
0.2.120
/opt/conda/default/lib/python3.10/site-packages/hailtop/aiocloud/aiogoogle/user_config.py:43: UserWarning: Reading spark-defaults.conf to determine GCS requester pays configuration. This is deprecated. Please use `hailctl config set gcs_requester_pays/project` and `hailctl config set gcs_requester_pays/buckets`.
  warnings.warn(
Running on Apache Spark version 3.3.0
SparkUI available at http://test-chets-hail120-m.c.broad-mpg-gnomad.internal:32913
Welcome to
     __  __     <>__
    / /_/ /__  __/ /
   / __  / _ `/ / /
  /_/ /_/\_,_/_/_/   version 0.2.120-f00f916faf78
LOGGING: writing to /create_vp_matrix.log
01/28/2026 04:31:44 AM (create_vp_matrix 634): 
        Running script with the following parameters:

            Data type: exomes
            Test: True
            Output postfix: julia_lg_vp_test
            Overwrite: True
            Tmp dir: gs://gnomad-tmp-4day
            Least consequence: 3_prime_UTR_variant
            Max freq: 0.05
        
01/28/2026 04:31:44 AM (create_vp_matrix 661): Creating variant filter Table...
01/28/2026 04:31:48 AM (create_vp_matrix 671): Filtering filter_ht, freq_ht, and vep_ht to test interval...
01/28/2026 04:31:52 AM (gnomad.utils.vep 959): Filtering to Ensembl transcripts...
01/28/2026 04:31:52 AM (gnomad.utils.vep 962): Filtering to protein coding transcripts...
01/28/2026 04:31:52 AM (gnomad.utils.vep 981): Filtering to variants with additional criteria...
2026-01-28 04:32:37.823 Hail: INFO: wrote table with 195937 rows in 202 partitions to gs://gnomad-tmp-4day/exomes.variant_filter.julia_lg_vp_test.ht
01/28/2026 04:32:38 AM (create_vp_matrix 683): Number of variants in the VEP Table that pass QC, have a consequence at least as severe as 3_prime_UTR_variant, and have a gnomAD AF <= 0.05: 195937
01/28/2026 04:32:38 AM (create_vp_matrix 690): Filtering gnomAD v4 exomes variant data MatrixTable...
01/28/2026 04:32:43 AM (basic_resources 177): Filtering to 13 intervals...
01/28/2026 04:32:44 AM (basic_resources 189): Dropping excessively multi-allelic site at chr19:5787204...
01/28/2026 04:32:44 AM (basic_resources 203): Removing 27 duplicate UKB samples by column index...
01/28/2026 04:32:48 AM (basic_resources 272): Total number of UKB samples removed from the VDS: 27
01/28/2026 04:32:48 AM (basic_resources 352): Filtering VDS to release samples only...
01/28/2026 04:32:55 AM (basic_resources 530): Splitting multiallelics...5) / 16]
2026-01-28 04:41:35.112 Hail: INFO: wrote matrix table with 195937 rows and 730947 columns in 202 partitions to gs://gnomad-tmp-4day/exomes.filtered_vmt.julia_lg_vp_test.mt
01/28/2026 04:41:36 AM (create_vp_matrix 708): The filtered VDS has been written...
01/28/2026 04:41:36 AM (create_vp_matrix 775): Time taken to run the script is 602.2579773039997 seconds.
Job [a5c94a41012c48dda388667cda4da3c9] finished successfully.
done: true
driverControlFilesUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/f9607727-6d6a-4f75-bcc3-9a489601b2c5/jobs/a5c94a41012c48dda388667cda4da3c9/
driverOutputResourceUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/f9607727-6d6a-4f75-bcc3-9a489601b2c5/jobs/a5c94a41012c48dda388667cda4da3c9/driveroutput
jobUuid: 01ee1ac6-3252-3d87-9a33-e2767ac50ab7
placement:
  clusterName: test-chets-hail120
  clusterUuid: f9607727-6d6a-4f75-bcc3-9a489601b2c5
pysparkJob:
  args:
  - --create-variant-filter-ht
  - --filter-vmt
  - --test
  - --output-postfix
  - julia_lg_vp_test
  - --overwrite
  mainPythonFileUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/f9607727-6d6a-4f75-bcc3-9a489601b2c5/jobs/a5c94a41012c48dda388667cda4da3c9/staging/create_vp_matrix.py
  pythonFileUris:
  - gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/f9607727-6d6a-4f75-bcc3-9a489601b2c5/jobs/a5c94a41012c48dda388667cda4da3c9/staging/pyscripts_k7yntrue.zip
reference:
  jobId: a5c94a41012c48dda388667cda4da3c9
  projectId: broad-mpg-gnomad
status:
  state: DONE
  stateStartTime: '2026-01-28T04:41:39.949248Z'
statusHistory:
- state: PENDING
  stateStartTime: '2026-01-28T04:31:26.548211Z'
- state: SETUP_DONE
  stateStartTime: '2026-01-28T04:31:26.565052Z'
- details: Agent reported job success
  state: RUNNING
  stateStartTime: '2026-01-28T04:31:26.770522Z'
yarnApplications:
- name: Hail
  progress: 1.0
  state: FINISHED
  trackingUrl: http://test-chets-hail120-m.c.broad-mpg-gnomad.internal.:8088/proxy/application_1769572196642_0005/

Time taken to run the script is 602.2579773039997 seconds.
```

### Step 2: Create Variant Pair List HT

**Cluster:** `test-chets-highmem` (Hail 0.2.134)

```bash
hailctl dataproc submit test-chets-highmem gnomad_chets/v4/create_vp_matrix.py \
    --create-variant-pair-list-ht \
    --test \
    --output-postfix julia_lg_vp_test \
    --overwrite \
    --pyfiles gnomad_chets
```

**Output:**

```
Submitting to cluster 'test-chets-highmem'...
gcloud command:
gcloud dataproc jobs submit pyspark gnomad_chets/v4/create_vp_matrix.py \
    --files= \
    --py-files=/var/folders/r8/f581hggx4_n6rmtz01r7nljm0000gq/T/pyscripts_q4bp0e2h.zip \
    --properties= \
    -- \
    --create-variant-pair-list-ht \
    --test \
    --output-postfix \
    julia_lg_vp_test \
    --overwrite
Job [7768ed392ad242d6b2c72314235572cd] submitted.
Waiting for job output...
0.2.134
01/28/2026 05:00:53 AM (create_vp_matrix 624): WARNING: Using Hail version 0.2.134-952ae203dbbe which is greater than 0.2.120. This will cause issues in create_variant_pair_ht, please use Hail version 0.2.120 or lower.
/opt/conda/default/lib/python3.11/site-packages/hailtop/aiocloud/aiogoogle/user_config.py:43: UserWarning: Reading spark-defaults.conf to determine GCS requester pays configuration. This is deprecated. Please use `hailctl config set gcs_requester_pays/project` and `hailctl config set gcs_requester_pays/buckets`.
  warnings.warn(
Running on Apache Spark version 3.5.0
SparkUI available at http://test-chets-highmem-m.c.broad-mpg-gnomad.internal:39259
Welcome to
     __  __     <>__
    / /_/ /__  __/ /
   / __  / _ `/ / /
  /_/ /_/\_,_/_/_/   version 0.2.134-952ae203dbbe
LOGGING: writing to /create_vp_matrix.log
01/28/2026 05:01:06 AM (create_vp_matrix 634): 
        Running script with the following parameters:
            Data type: exomes
            Test: True
            Output postfix: julia_lg_vp_test
            Overwrite: True
            Tmp dir: gs://gnomad-tmp-4day
            Least consequence: 3_prime_UTR_variant
            Max freq: 0.05
        
01/28/2026 05:01:06 AM (create_vp_matrix 711): Creating variant pair list Table...
2026-01-28 05:01:08.715 Hail: WARN: entries(): Resulting entries table is sorted by '(row_key, col_key)'.
    To preserve row-major matrix table order, first unkey columns with 'key_cols_by()'
2026-01-28 05:09:06.262 Hail: INFO: Ordering unsorted dataset with network shuffle
WARNING: An illegal reflective access operation has occurred
WARNING: Illegal reflective access by org.apache.spark.util.SizeEstimator$ (file:/usr/lib/spark/jars/spark-core_2.12-3.5.0.jar) to field java.lang.ref.Reference.referent
WARNING: Please consider reporting this to the maintainers of org.apache.spark.util.SizeEstimator$
WARNING: Use --illegal-access=warn to enable warnings of further illegal reflective access operations
WARNING: All illegal access operations will be denied in a future release
2026-01-28 05:15:06.048 Hail: INFO: wrote table with 10159581 rows in 202 partitions to gs://gnomad-tmp-4day/create_variant_pair_ht.gene_sample_grouped-yZwLTCzZH8BruDcTxb2Iux.ht
2026-01-28 05:22:08.634 Hail: INFO: wrote table with 4863611021 rows in 202 partitions to gs://gnomad-tmp-4day/persist_TablebfW8CLxC1s
2026-01-28 06:25:58.875 Hail: INFO: wrote table with 4863611021 rows in 3562 partitions to gs://gnomad-tmp-4day/persist_Table5mOXNVBJmN
2026-01-28 06:28:48.358 Hail: INFO: Coerced sorted dataset 3:(202 + 14) / 202]
2026-01-28 06:31:05.741 Hail: INFO: wrote table with 21210652 rows in 3559 partitions to gs://gnomad-tmp-4day/exomes.variant_pairs.julia_lg_vp_test.ht
01/28/2026 06:31:06 AM (create_vp_matrix 717): The variant pair list Table has been written...
The number of unique variant pairs is 21210652
01/28/2026 06:31:06 AM (create_vp_matrix 775): Time taken to run the script is 5412.365701876 seconds.
Job [7768ed392ad242d6b2c72314235572cd] finished successfully.
done: true
driverControlFilesUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/585b459a-8149-4d24-94c2-0f21cf29485e/jobs/7768ed392ad242d6b2c72314235572cd/
driverOutputResourceUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/585b459a-8149-4d24-94c2-0f21cf29485e/jobs/7768ed392ad242d6b2c72314235572cd/driveroutput
jobUuid: 171f7c99-8df6-3f8d-981c-c3931cf90046
placement:
  clusterName: test-chets-highmem
  clusterUuid: 585b459a-8149-4d24-94c2-0f21cf29485e
pysparkJob:
  args:
  - --create-variant-pair-list-ht
  - --test
  - --output-postfix
  - julia_lg_vp_test
  - --overwrite
  mainPythonFileUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/585b459a-8149-4d24-94c2-0f21cf29485e/jobs/7768ed392ad242d6b2c72314235572cd/staging/create_vp_matrix.py
  pythonFileUris:
  - gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/585b459a-8149-4d24-94c2-0f21cf29485e/jobs/7768ed392ad242d6b2c72314235572cd/staging/pyscripts_q4bp0e2h.zip
reference:
  jobId: 7768ed392ad242d6b2c72314235572cd
  projectId: broad-mpg-gnomad
status:
  state: DONE
  stateStartTime: '2026-01-28T06:31:10.017248Z'
statusHistory:
- state: PENDING
  stateStartTime: '2026-01-28T05:00:42.819946Z'
- state: SETUP_DONE
  stateStartTime: '2026-01-28T05:00:42.842685Z'
- details: Agent reported job success
  state: RUNNING
  stateStartTime: '2026-01-28T05:00:43.178882Z'
yarnApplications:
- name: Hail
  progress: 1.0
  state: FINISHED
  trackingUrl: http://test-chets-highmem-m.c.broad-mpg-gnomad.internal.:8088/proxy/application_1769575879460_0001/

Time taken to run the script is 5412.365701876 seconds.
```

### Step 3: Create Dense Filtered MT

**Cluster:** `test-chets-hail120` (Hail 0.2.120)

```bash
hailctl dataproc submit test-chets-hail120 gnomad_chets/v4/create_vp_matrix.py \
    --create-dense-filtered-mt \
    --test \
    --output-postfix julia_lg_vp_test \
    --overwrite \
    --pyfiles gnomad_chets
```

**Output:**

```
Submitting to cluster 'test-chets-hail120'...
gcloud command:
gcloud dataproc jobs submit pyspark gnomad_chets/v4/create_vp_matrix.py \
    --files= \
    --py-files=/var/folders/r8/f581hggx4_n6rmtz01r7nljm0000gq/T/pyscripts_2rb7uio4.zip \
    --properties= \
    -- \
    --create-dense-filtered-mt \
    --test \
    --output-postfix \
    julia_lg_vp_test \
    --overwrite
Job [90066708d0184a28aa9d5695f2c7ab74] submitted.
Waiting for job output...
SLF4J: No SLF4J providers were found.
SLF4J: Defaulting to no-operation (NOP) logger implementation
SLF4J: See https://www.slf4j.org/codes.html#noProviders for further details.
SLF4J: Class path contains SLF4J bindings targeting slf4j-api versions 1.7.x or earlier.
SLF4J: Ignoring binding found at [jar:file:/usr/lib/spark/jars/log4j-slf4j-impl-2.17.2.jar!/org/slf4j/impl/StaticLoggerBinder.class]
SLF4J: See https://www.slf4j.org/codes.html#ignoredBindings for an explanation.
/opt/conda/default/lib/python3.10/site-packages/hailtop/aiocloud/aiogoogle/user_config.py:43: UserWarning: Reading spark-defaults.conf to determine GCS requester pays configuration. This is deprecated. Please use `hailctl config set gcs_requester_pays/project` and `hailctl config set gcs_requester_pays/buckets`.
  warnings.warn(
Running on Apache Spark version 3.3.0
SparkUI available at http://test-chets-hail120-m.c.broad-mpg-gnomad.internal:36931
Welcome to
     __  __     <>__
    / /_/ /__  __/ /
   / __  / _ `/ / /
  /_/ /_/\_,_/_/_/   version 0.2.120-f00f916faf78
LOGGING: writing to /create_vp_matrix.log
01/28/2026 02:28:47 PM (create_vp_matrix 633): 
        Running script with the following parameters:

            Data type: exomes
            Test: True
            Output postfix: julia_lg_vp_test
            Overwrite: True
            Tmp dir: gs://gnomad-tmp-4day
            Least consequence: 3_prime_UTR_variant
            Max freq: 0.05
        
01/28/2026 02:28:47 PM (create_vp_matrix 722): Creating dense filtered MatrixTable...
2026-01-28 14:29:29.567 Hail: INFO: Coerced sorted dataset===>(3555 + 4) / 3559]
2026-01-28 14:29:54.092 Hail: INFO: Ordering unsorted dataset with network shuffle
2026-01-28 14:32:11.204 Hail: INFO: wrote table with 195628 rows in 4144 partitions to gs://gnomad-tmp-4day/create_dense_filtered_mt.variants-nWIsV0Ex2DWrspoGUpLTCo.ht
01/28/2026 02:32:18 PM (basic_resources 177): Filtering to 13 intervals...
01/28/2026 02:32:19 PM (basic_resources 189): Dropping excessively multi-allelic site at chr19:5787204...
01/28/2026 02:32:20 PM (basic_resources 203): Removing 27 duplicate UKB samples by column index...
01/28/2026 02:32:24 PM (basic_resources 272): Total number of UKB samples removed from the VDS: 27
01/28/2026 02:32:24 PM (basic_resources 352): Filtering VDS to release samples only...
01/28/2026 02:32:30 PM (basic_resources 530): Splitting multiallelics...2) / 16]
2026-01-28 15:48:43.301 Hail: INFO: wrote matrix table with 195628 rows and 730947 columns in 368 partitions to gs://gnomad-tmp-4day/exomes.filtered.dense.julia_lg_vp_test.mt
01/28/2026 03:48:46 PM (create_vp_matrix 742): The dense filtered MatrixTable has been written...
The number of rows in the dense filtered MatrixTable is 195628
01/28/2026 03:48:46 PM (create_vp_matrix 774): Time taken to run the script is 4812.0827281540005 seconds.
Job [90066708d0184a28aa9d5695f2c7ab74] finished successfully.
done: true
driverControlFilesUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/b1b0f564-e1cb-4ec0-b3ca-333bb29e7bdc/jobs/90066708d0184a28aa9d5695f2c7ab74/
driverOutputResourceUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/b1b0f564-e1cb-4ec0-b3ca-333bb29e7bdc/jobs/90066708d0184a28aa9d5695f2c7ab74/driveroutput
jobUuid: afc01b88-afc7-3bcd-87a5-3dd0ed340a0f
placement:
  clusterName: test-chets-hail120
  clusterUuid: b1b0f564-e1cb-4ec0-b3ca-333bb29e7bdc
pysparkJob:
  args:
  - --create-dense-filtered-mt
  - --test
  - --output-postfix
  - julia_lg_vp_test
  - --overwrite
  mainPythonFileUri: gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/b1b0f564-e1cb-4ec0-b3ca-333bb29e7bdc/jobs/90066708d0184a28aa9d5695f2c7ab74/staging/create_vp_matrix.py
  pythonFileUris:
  - gs://dataproc-faa46220-ec08-4f5b-92bd-9722e1963047-us-central1/google-cloud-dataproc-metainfo/b1b0f564-e1cb-4ec0-b3ca-333bb29e7bdc/jobs/90066708d0184a28aa9d5695f2c7ab74/staging/pyscripts_2rb7uio4.zip
reference:
  jobId: 90066708d0184a28aa9d5695f2c7ab74
  projectId: broad-mpg-gnomad
status:
  state: DONE
  stateStartTime: '2026-01-28T15:48:49.510928Z'
statusHistory:
- state: PENDING
  stateStartTime: '2026-01-28T14:28:24.108714Z'
- state: SETUP_DONE
  stateStartTime: '2026-01-28T14:28:24.138260Z'
- details: Agent reported job success
  state: RUNNING
  stateStartTime: '2026-01-28T14:28:24.436297Z'
yarnApplications:
- name: Hail
  progress: 1.0
  state: FINISHED
  trackingUrl: http://test-chets-hail120-m.c.broad-mpg-gnomad.internal.:8088/proxy/application_1769610172905_0001/
```

### Step 4: Create Variant Pair Genotype HT and Counts HT

**Cluster:** `test-chets-highmem` (Hail 0.2.134)

```bash
hailctl dataproc submit test-chets-highmem gnomad_chets/v4/create_vp_matrix.py \
    --create-variant-pair-genotype-ht \
    --create-variant-pair-genotype-counts-ht \
    --test \
    --output-postfix julia_lg_vp_test \
    --overwrite \
    --pyfiles gnomad_chets
```

**Output:** (Pending)
