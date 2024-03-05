
# Setting up snpEff database
By following these steps, you'll ensure that your snpEff database is correctly configured for subsequent usage with snpkit. 


## 1. Create a directory

> First, create a directory named after your reference genome within the designated *data* directory.
For instance, if your reference genome is identified as KPNIH1, create a directory named **KPNIH1** within *data/KPNIH1*.
```

mkdir data/KPNIH1

```

## 2. Upload Genbank and GFF files

  Upload the corresponding genbank and GFF files into the newly created directory.
  

## 3. Rename Genbank and GFF files
After uploading the files, it's crucial to rename them for compatibility, otherwise snpEff will not run.

> Rename the genbank file to _**genes.gbk**_. For example, if your genbank file is named *KPNIH1.gbk*, rename it to ***genes.gbk***.
```

mv KPNIH1.gbk genes.gbk

```
    
> Similarly, rename the GFF file to ***genes.gff***. For instance, if your GFF file is named *KPNIH1.gff*, rename it to ***genes.gff***.
```

mv KPNIH1.gff genes.gff

```