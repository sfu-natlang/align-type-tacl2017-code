# align-type-tacl2017-code
Joint prediction of word alignment with alignment types

## Usage
    USAGE: java -jar type-aligner.jar [OPTIONS]

    Option                      description
    
    --train <trainFilePrefix>   [required] parallel training corpus prefix

    --test <testFilePrefix>     [required] parallel test corpus prefix
    
    --source <sourceLang>       [required] source language suffix (suffix of the source train and test data)
    
    --target <targetLang>       [required] target language suffix (suffix of the target train and test data)
    
    --model <model>             [required] alignment model name. Choose among the implemented models (ibm1, hmm, ibm1+type+gen, hmm+type+gen)
    --trainSize <trainSize>     [required] size of the training data (number of the sentences from the beginning of the parallel
                                train data to consider for training).
    --testSize <testSize>       [required] size of the test data (number of sentences from the beginning of the parallel test data).
    
    --reference <referenceFile> reference file for the test data which is annotated with gold alignment. If provided, precision, recall and F1-score will be computed.
    --output <outputFile>       the prefix of the predicted alignment files for the test data (default: <sourceLang>-<targetLang>.wa)
    
    --augmentedTrain <augmentTrainFilePrefix>   augmented training file prefix. This data is used to augment the training data.
    --augmentedTrainSize <augmentedTrainSize>   augmented training file size. This is a required option if --augmentedTrain is set.
    
    

## Examples:
We can train a baseline HMM-based word alignment model (Och and Ney, 2003) on a parallel Chinese-English corpus `{<trainFilePrefix>.cn, <trainFilePrefix>.en}`, and then using this trained model compute the Viterbi alignment for a test corpus `{<testFilePrefix>.cn, <testFilePrefix>.en}`. The output alignment file will be saved in `<outputFile>.wa`. To do this, run the following command:
```
java -jar type-aligner.jar --train <trainFilePrefix> --test testFilePrefix --source cn --target en --trainSize <trainSize> --testSize <testSize> --model ibm1 --reference <referenceFile> --output <outputFile>
```

## Citation

If you use this code in your research, please cite the following paper:

Anahita Mansouri Bigvand, Te Bu, and Anoop Sarkar, "Joint prediction of word alignment with alignment types", *Transactions of the Association for Computational Linguistics* (2017)

## Contacts
* Anahita Mansouri Bigvand <amansour@cs.sfu.ca>
* Anoop Sarkar <anoop@cs.sfu.ca>

## References


