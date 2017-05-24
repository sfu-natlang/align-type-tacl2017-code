# align-type-tacl2017-code
Joint prediction of word alignment with alignment types

## Usage
    USAGE: java -jar type-aligner.jar [OPTIONS]

    Option                      description
    
    --train <trainFile>         [required] parallel training corpus prefix

    --test <testFile>           [required] parallel test corpus prefix
    
    --source <sourceLang>       [required] source language suffix (suffix of the source train and test data)
    
    --target <targetLang>       [required] target language suffix (suffix of the target train and test data)
    
    --model <model>             [required] alignment model name. Choose among the implemented models (ibm1, hmm, ibm1+type+gen, hmm+type+gen)
    --trainSize <trainSize>     [required] size of the training data (number of the sentences from the beginning of the parallel train data to consider for training).
    --testSize <testSize>       [required] size of the test data (number of sentences from the beginning of the parallel test data).
    
    --reference <referenceFile> reference file for the test data which is annotated with gold alignment
    
## Citation

If you use this code in your research, please cite the following paper:

Anahita Mansouri Bigvand, Te Bu, and Anoop Sarkar, "Joint prediction of word alignment with alignment types", *Transactions of the Association for Computational Linguistics* (2017)

## Contacts
* Anahita Mansouri Bigvand <amansour@cs.sfu.ca>
* Anoop Sarkar <anoop@cs.sfu.ca>

## References


