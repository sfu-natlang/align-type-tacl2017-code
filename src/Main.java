import java.io.IOException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.*;

public class Main {
	public static void main(String[] args) throws IOException, ClassNotFoundException {
		for (String s : args)
			System.out.println(s);
		Options options = new Options();
		CommandLineParser parser = new DefaultParser();
		
		CommandLine commandLine;
	/*	Option option_A = Option.builder("A").required(true).desc("The A option").longOpt("opt3").hasArgs().build();
		Option option_r = Option.builder("r").required(true).desc("The r option").longOpt("opt1").build();
		Option option_test = Option.builder().required(true).desc("The test option").longOpt("test").build();
	*/	
		Option optionTrain = new Option("train", "train", true, "training file prefix");
		optionTrain.setRequired(true);
        options.addOption(optionTrain);
        
        Option optionSource = new Option("f", "source", true, "source language suffix");
        optionSource.setRequired(true);
        options.addOption(optionSource);
        
        Option optionTarget = new Option("e", "target", true, "target language suffix");
        optionTarget.setRequired(true);
        options.addOption(optionTarget);
        
        Option optionTest = new Option("", "test", true, "test file prefix");
        optionTest.setRequired(true);
        options.addOption(optionTest);
        
        Option modelOption = new Option("m", "model", true, "model name: could be any of the following models \n ibm1, hmm, ibm1+type+gen, ibm1+type+disc, hmm+type+gen, hmm+type+disc");
        modelOption.setRequired(true);
        options.addOption(modelOption);
        
        Option trainingSizeOption = new Option("", "trainSize", true, "size of the training data or number of lines from the training files to do the training");
        trainingSizeOption.setRequired(true);
        options.addOption(trainingSizeOption);
        
        Option testSizeOption = new Option("", "testSize", true, "size of the test data or number of lines from the test files to do the test");
        testSizeOption.setRequired(true);
        options.addOption(testSizeOption);
        
        Option referenceOption = new Option("r", "reference", true, "reference file for the test data which is annotated with gold alignment for the test data");
        options.addOption(referenceOption);
       
        
        Option augmentOption = new Option("at", "augmentedTrain", true, "augmented training file prefix. This data is used to augment the training data.");
        options.addOption(augmentOption);
        
        Option augmentSizeOption = new Option("as", "augmentedTrainSize", true, "augmented training file size. This is a required option if augmentedTrain is set.");
        options.addOption(augmentSizeOption);
       
        Option goldTrainOption = new Option("gt", "goldTrain", true, "gold training file. The file with gold alignment and alignment types.");
        options.addOption(goldTrainOption);
        
        Option POSTaggedTrainOption = new Option("", "POSTaggedTrain", true, "POS-tagged training file prefix");
        options.addOption(POSTaggedTrainOption);

        Option POSTaggedTestOption = new Option("", "POSTaggedTest", true, "POS-tagged test file prefix");
        options.addOption(POSTaggedTestOption);
        
        Option POSTaggedAugmentedTrainOption = new Option("", "POSTaggedAugmentedTrain", true, "POS-tagged augmented training file prefix");
        options.addOption(POSTaggedAugmentedTrainOption);
        
        Option referenceAlignmentAndTypeOption = new Option("", "referenceAlignmentAndType", true, "reference file for the test data which is annotated with gold alignment and alignment type for the test data");
        options.addOption(referenceAlignmentAndTypeOption);
        
        Option outputOption = new Option("o", "output", true, "the prefix of the alignment files for the test data");
        options.addOption(outputOption);
        //java -cp .:commons-cli-1.4.jar Main --train data-te/train.20k.seg.cln --test data-te/test.seg.cln --source cn --target en --trainSize 203 --testSize 19 --model ibm1 --reference data-te/test.seg.cln.gold.wa
        //java -jar type-aligner.jar --train data-te/train.20k.seg.cln --test data-te/test.seg.cln --source cn --target en --trainSize 203 --testSize 19 --model ibm1 --reference data-te/test.seg.cln.gold.wa
		//String[] testArgs = {"--train","asd","-r","opt1","--test","-A","opt3",};
        // For ant run: java -cp bin/.:lib/commons-cli-1.4.jar Main --train data-te/train.20k ...
		String[] testArgs = {"--train", "data-te/train.20k.seg.cln", 
				"--source", "cn",
				"--target", "en",
				"--test","data-te/test.seg.cln",
				"--trainSize", "203",
				"--testSize", "19",
				"--model","ibm1",
				"--reference","data-te/test.seg.cln.gold.wa",
				"--augmentedTrain", "data-te/train.20k+hk",
				"--augmentedTrainSize","303",
			//	"--POSTaggedAugmentedTrain", "data-te/train.20k+hk.tags",
			//	"--goldTrain", "data-te/train.20k.seg.cln.wa",
			//	"--POSTaggedTrain","data-te/train.20k.tags",
			//	"--POSTaggedTest","data-te/test.tags",
			//	"--referenceAlignmentAndType","data-te/test.seg.cln.wa",
				};
				
	/*	options.addOption(option_A);
		options.addOption(option_r);
		options.addOption(option_test);
	*/	int trainingSize = 0;
		int testSize = 0;
		int augmentedTrainSize = 0;
		String augmentedTrain = "";
		String goldTrain = "";
		String POSTaggedTrain = "";
		String POSTaggedTest = "";
		//String refere
		String referenceAlignmentAndType = "";
		String POSTaggedAugmentedTrain = "";
		
		try{
			commandLine = parser.parse(options, args);
			if (commandLine.hasOption("A"))
				System.out.println("option A is present and the value is " + commandLine.getOptionValue("A"));
			if (commandLine.hasOption("test"))
				System.out.println("option test is present and this is a flag option");
			if (commandLine.hasOption("i"))
				System.out.println("with i " + commandLine.getOptionValue("i"));
			if (commandLine.hasOption("trainSize"))
				trainingSize = Integer.parseInt(commandLine.getOptionValue("trainSize"));
			if (commandLine.hasOption("testSize"))
				testSize = Integer.parseInt(commandLine.getOptionValue("testSize"));
			if (commandLine.hasOption("model") && commandLine.getOptionValue("model").indexOf("+type+gen") != -1){
				goldTrain = commandLine.getOptionValue("goldTrain");
				POSTaggedTrain = commandLine.getOptionValue("POSTaggedTrain");
				POSTaggedTest = commandLine.getOptionValue("POSTaggedTest");
				POSTaggedAugmentedTrain = commandLine.getOptionValue("POSTaggedAugmentedTrain");
				System.out.println(goldTrain + " " + POSTaggedTrain + " " + POSTaggedTest+" " + POSTaggedAugmentedTrain);
				if (goldTrain.equals("") || POSTaggedTrain.equals("") || POSTaggedTest.equals("")){
					System.err.println("To train an alignment model with alignment type, the following arguments are required: \n" +
							"goldTrain, POSTaggedTrain and POSTaggedTest.");
					throw new ParseException("AGAIN");
				}
				referenceAlignmentAndType = commandLine.getOptionValue("referenceAlignmentAndType");
			}
			
			String train = commandLine.getOptionValue("train");
			String sourceLang = commandLine.getOptionValue("source");
			String targetLang = commandLine.getOptionValue("target");
			String test = commandLine.getOptionValue("test");
			String output = sourceLang + "-" + targetLang + ".wa";
			if (commandLine.hasOption("output"))
				output = commandLine.getOptionValue("output");
			//write a code to fix for null ref
			String referenceFile = commandLine.getOptionValue("reference");
			if (commandLine.hasOption("augmentedTrainSize"))
				augmentedTrainSize = Integer.parseInt(commandLine.getOptionValue("augmentedTrainSize"));
			if (commandLine.hasOption("augmentedTrain"))
				augmentedTrain = commandLine.getOptionValue("augmentedTrain");
			
			if (commandLine.hasOption("model") && commandLine.getOptionValue("model").equals("ibm1")){
				if (!commandLine.hasOption("augmentedTrain"))
					IBM_Model1.mainIBM1(trainingSize, testSize, train, sourceLang, targetLang, test, referenceFile, output);
				else
					IBM1WithTestTrainHK.mainIBM1(trainingSize, testSize, train, sourceLang, targetLang, test, referenceFile, augmentedTrain, augmentedTrainSize, output);
			}
			if (commandLine.hasOption("model") && commandLine.getOptionValue("model").equals("hmm")){
				if (!commandLine.hasOption("augmentedTrain"))
					HMMWithNullForTuningTe.mainHMM(trainingSize, testSize, train, sourceLang, targetLang, test, referenceFile, output);
				else
					HMMWithNullTestTrainHKTe.mainHMM(trainingSize, testSize, train, sourceLang, targetLang, test, referenceFile, output, augmentedTrain, augmentedTrainSize);	
			}
			if (commandLine.hasOption("model") && commandLine.getOptionValue("model").equals("ibm1+type+gen")){
				if (!commandLine.hasOption("augmentedTrain"))
					IBM1WithAlignmentTypeCorrect.fScoreExperiment(trainingSize, testSize, train, sourceLang, targetLang, test, referenceFile, output, 
						 goldTrain, POSTaggedTrain, POSTaggedTest, referenceAlignmentAndType);
				else
					IBM1WithAlignmentTypeCorrect.augmentedModelExperimentMain(trainingSize, testSize, train, sourceLang, targetLang, test, referenceFile, output, goldTrain, POSTaggedTrain, POSTaggedTest, referenceAlignmentAndType, augmentedTrain, POSTaggedAugmentedTrain, augmentedTrainSize);
				
			}
			if (commandLine.hasOption("model") && commandLine.getOptionValue("model").equals("hmm+type+gen")){
				if (!commandLine.hasOption("augmentedTrain"))
					HMMWithAlignmentTypeNullCorrect.FScoreExperiment1AlignmentTypeMain(trainingSize, testSize, train, sourceLang, targetLang, test, referenceFile, output, goldTrain, POSTaggedTrain, POSTaggedTest, referenceAlignmentAndType);
				//TODO
				else{
					HMMWithAlignmentTypeNullCorrect.augmentedModelExperiment(trainingSize, testSize, train, sourceLang, targetLang, test, referenceFile, output, goldTrain, POSTaggedTrain, POSTaggedTest, referenceAlignmentAndType, augmentedTrain, POSTaggedAugmentedTrain, augmentedTrainSize);
				}
			}
			
			
			for(String argument : commandLine.getArgs()){
				System.out.println(argument);
			}
		}  catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("Main", options);
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
		
	}
	
	
}
