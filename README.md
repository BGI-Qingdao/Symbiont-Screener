# ContamFilter
Filter contamination from filial tgs reads by parental kmer sets.


## USAGE

```
Usage    :
    ./ContamFilter.sh [OPTION]

Filter contamination from filial tgs reads by parental kmer sets.

Options  :
        --filial            filial TGS reads file in FASTA/Q format.
                            file in gzip format can be accepted, but filename must end by .gz.

        --format            fasta/fastq . set the format of --filial.
                            [ optional, default fasta. ]

        --thread            thread num.
                            [ optional, default 8 threads. ]

        --memory            x (GB) of memory to used by meryl.
                            [ optional, default 50GB. ]

        --mer               mer-size
                            [ optional, default 21. ]

        --use_existing_lib  on/off
                            [ optional, default off]

        --paternal          paternal NGS reads file in FASTQ format.
                            ( note : gzip format is NOT supported. )
                            [ optional, needed when --use_existing_libs off ]

        --maternal          maternal NGS reads file in FASTQ format.
                            ( note : gzip format is NOT supported. )
                            [ optional, needed when --use_existing_libs off ]

        --paternal_mer      existing paternal specific kmer lib file.
                            [ optional, needed when --use_existing_libs on ]

        --maternal_mer      existing maternal specific kmer lib file.
                            [ optional, needed when --use_existing_libs on ]

        --common_mer        existing parental common kmer lib file.
                            [ optional, needed when --use_existing_libs on ]

        --density_common    the minimum density of common_kmer
                            [ optional, default 0.01]

        --density_specific  the minimum density of parental specific kmer
                            [ optional, default 0.00005]

        --count_specific    the minimum count of parental specific kmer
                            [ optional, default 2]

        --help              print this usage message.

Examples :

    ./ContamFilter.sh --filial son.fasta  --paternal father.fastq --maternal mater.fastq

    # if the filial read file follow fastq format :
    ./ContamFilter.sh --paternal father.fastq --maternal mater.fastq --filial son.fastq --format fastq

    # if there are more than one filial read files :
    ./ContamFilter.sh --paternal father.fastq --maternal mater.fastq --filial son.L01.fasta --filial son.L02.fasta

    # use existing libs
    ./ContamFilter.sh --filial son.fasta  --use_existing_lib on  --paternal_mer p.kmer\
                                                                 --maternal_mer m.kmer\
                                                                 --common_mer   c.kmer


```


