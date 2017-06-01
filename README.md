# align-type-tacl2017-code

This is the implementation of our TACL paper "Joint Prediction of Word Alignment with Alignment Types".

## Features
- Baseline IBM1
- Baseline HMM-based word alignment (Och and Ney, 2003)
- Generative models (IBM1+Type+Gen and HMM+Type+Gen) for the joint prediction of word alignment with alignment types

## Citation
If you use this code in your research, please cite the following paper:

Anahita Mansouri Bigvand, Te Bu, and Anoop Sarkar, "Joint prediction of word alignment with alignment types", *Transactions of the Association for Computational Linguistics* (2017)

## Building Instructions
To build the code, make sure you have Apache Ant installed on your machine. Set your current working directory to the root directory of the distribution, where the file "build.xml" is located, and then build with the following command:
```
ant
```

## Usage
    USAGE: java -jar dist/type-aligner.jar [OPTIONS]
| Option            | description |
|---                   |---|
| --train `<trainFilePrefix>` | [required] parallel training corpus prefix |
| --test `<testFilePrefix>`    | [required] parallel test corpus prefix |
| --source `<sourceLang>`      | [required] source language suffix (suffix of the source train and test data)|
| --target `<targetLang>`      | [required] target language suffix (suffix of the target train and test data) |
| --model `<model>`            |  [required] alignment model name. Choose among the implemented models (ibm1, hmm, ibm1+type+gen, hmm+type+gen)|
| --reference `<referenceFile>`  | reference file for the test data which is annotated with gold alignment. If provided, precision, recall and F1-score will be computed. |
| --output `<outputFile>`       | the prefix of the predicted alignment files for the test data (default: `<sourceLang>-<targetLang>.wa`) |
| --augmentedTrain `<augmentTrainFilePrefix>` | augmented training file prefix. This data is used to augment the training data. |
| --augmentedTrainSize `<augmentedTrainSize>`   | augmented training file size. This is a required option if --augmentedTrain is set. |
| --goldTrain `<goldTrainFile>`  | gold training file. The file with gold alignment and alignment types. |
| --POSTaggedTrain `<POSTaggedTrainFilePrefix>` | POS-tagged training file prefix |
| --POSTaggedTest `<POSTaggedTestFilePrefix>`  | POS-tagged test file prefix |
| --POSTaggedAugmentedTrain `<POSTaggedAugmentedTrainFilePrefix>` | POS-tagged augmented training file prefix |
| --referenceAlignmentAndType `<referenceAlignmentAndTypeFile>` | reference file for the test data which is annotated with gold alignment and alignment type |
   

#### Baselines (IBM1 and HMM)
We can train a baseline HMM-based word alignment model (

and Ney, 2003) on a parallel Chinese-English corpus `{<trainFilePrefix>.cn, <trainFilePrefix>.en}`, and then using this trained model compute the Viterbi alignment for a test corpus `{<testFilePrefix>.cn, <testFilePrefix>.en}`. The output alignment file will be saved in `<outputFile>.wa`. To do this, run the following command:
```
java -jar dist/type-aligner.jar --train <trainFilePrefix> --test <testFilePrefix> --source cn --target en --trainSize <trainSize> --testSize <testSize> --model hmm --reference <referenceFile> --output <outputFile>
```
For baseline IBM1, use `--model ibm1` option instead. 

#### Generative Models (IBM1+Type+Gen and HMM+Type+Gen)

To jointly predict alignment and alignment type, you need a parallel training data which is annotated with gold 
word alignment and alignment types. The GALE Chinese-English Word Alignment and Tagging Corpus is an example of such a corpus. 

To train a generative model, you need to provide `--goldTrain <goldTrainFile>` which points to the file with gold alignment and alignment types for the corresponding training corpus `{<trainFilePrefix>.cn, <trainFilePrefix>.en}`. 
For each line (sentence pair) in the training corpus, there is a corresponding line including the gold alignment and alignment types in the `<goldTrainFile>` with a special format. For example, a line from the "GALE Chinese-English Word Alignment and Tagging Corpus" has the following format:

1-2(SEM) 2-3(SEM) 3-5(SEM) 4-7(SEM) 5-6[TEN],8(GIS) 6-9(SEM) 7-10(FUN) 8-11(SEM) 9-12(SEM) 10-13(SEM) 11-14(SEM) 12-15(FUN) 13-16(SEM) 14-17(SEM)

This program works for the above format. Apart from the `--goldTrain <goldTrainFile>` option, `--POSTaggedTrain` and `--POSTaggedTest` are required options for joint prediction of word alignment with alignment types.

where `--POSTaggedTrain <POSTaggedTrainFilePrefix>` and `--POSTaggedTest <POSTaggedTestFilePrefix>` are the prefixes for the POS-tagged train and test corpus files as explained in the paper.

Also, to obtain the precision, recall and F1-score for the joint prediction of word alignment and alignment types, provide the program with `--referenceAlignmentAndType <referenceAlignmentAndTypeFile>` option which indicates the reference file with gold alignment and alignment types. Note that you can still set the `--reference <referenceFile>` option to get the alignment task results. For example to train a HMM+Type+Gen model on the train corpus, and test on the test corpus, and get the results for word alignment task and the joint prediction of word alignment and alignment types task, run the following command:
```
java -jar dist/type-aligner.jar --train <trainFilePrefix> --test <testFilePrefix> --source cn --target en --trainSize <trainSize> --testSize <testSize> --model hmm+type+gen --reference <referenceFile> --output <outputFile> --goldTrain <goldTrainFile> --POSTaggedTrain <POSTaggedTrainFilePrefix> --POSTaggedTest <POSTaggedTestFilePrefix> --referenceAlignmentAndType <referenceAlignmentAndTypeFile>
```
#### Augmented Models:
Since we only have a small amount of training data which is annotated with gold alignment and alignment types, we can augment our training data with a parallel data with no annotation.  
The only extra parameter we have to set apart from `--augmentedTrain <augmentTrainFilePrefix>` and `--augmentedTrainSize <augmentedTrainSize>` is `--POSTaggedAugmentedTrain <POSTaggedAugmentedTrainFilePrefix>` which
specifies the prefix of the POS-tagged augmented data. Note that the augmented training corpus `{<augmentTrainFilePrefix>.cn, <augmentTrainFilePrefix>.en}` include the `{<trainFilePrefix>.cn, <trainFilePrefix>.en}` followed by the rest of data (for which we have no annotation).

For example, to replicate our best results for word alignment task as well as the joint prediction of word alignment with alignment types task, run the following command:
```
java -jar dist/type-aligner.jar --train <trainFilePrefix> --test <testFilePrefix> --source cn --target en --trainSize <trainSize> --testSize <testSize> --model hmm+type+gen --reference <referenceFile> --output <outputFile> --goldTrain <goldTrainFile> --POSTaggedTrain <POSTaggedTrainFilePrefix> --POSTaggedTest <POSTaggedTestFilePrefix> --referenceAlignmentAndType <referenceAlignmentAndTypeFile> --augmentedTrain <augmentTrainFilePrefix> --augmentedTrainSize <augmentedTrainSize> --POSTaggedAugmentedTrain <POSTaggedAugmentedTrainFilePrefix>
```
In order to use the generative model for a machine translation experiment, simply test on the augmented data ( `--test <augmentTrainFilePrefix>`) to find the alignment for the `{<augmentTrainFilePrefix>.cn, <augmentTrainFilePrefix>.en}` corpus. The `<testSize>` and `<augmentedTrainSize>` are also the same accordingly.

## Contacts
* Anahita Mansouri Bigvand <amansour@cs.sfu.ca>
* Anoop Sarkar <anoop@cs.sfu.ca>

## References
- Peter F Brown, Vincent J Della Pietra, Stephen A Della Pietra, and Robert L Mercer. 1993. The mathematics of statistical machine translation: Parameter estimation. Computational linguistics, 19(2):263–311.
- Franz Josef Och and Hermann Ney. 2003. A systematic comparison of various statistical alignment models. Computational linguistics, 29(1):19–51
