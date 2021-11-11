import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.zip.GZIPInputStream;

import org.apache.commons.cli.*;

public class ApertureMain {
	public final static String APERTURE_VERSION="1.3";
	public static int lCode,lPos,rCode,rPos;
	public static boolean debug=false;
	//public static Integer sync=new Integer(0);
	
	private void readClassify(String r1,String r2,String index,String workDir,String projectName,
			                 int r1BarStart,int r2BarStart,int r1BarLen,int r2BarLen,int r1InsStart,int r2InsStart,int nthreads,int mergeK,
			                 int minIndelLen,double similarity,double maxBpScore,String[] commands)
			                throws FileNotFoundException, IllegalBioFileException, IllegalThreadNumberException, 
			                IOException, InterruptedException, ExecutionException, ClassNotFoundException, IllegalPositionException 
	{
		
		//TranslationTable table=readTranslationTable(Paths.get(index+".tt"));
		//System.out.print(table);
		
		Path r1Path=Paths.get(r1);
		ensureFastqFileExist(r1Path);
		Path r2Path=Paths.get(r2);
		ensureFastqFileExist(r1Path);
		
		Path ttPath=Paths.get(index+".tt");
		Path kcPath=Paths.get(index+".km");
		Path longkcPath=Paths.get(index+".long.km");
		Path spacedkcPath=Paths.get(index+".spaced.km");
		Path ciPath=Paths.get(index+".ci");
		ensureIndexFileExist(ttPath);
		ensureIndexFileExist(kcPath);
		ensureIndexFileExist(longkcPath);
		ensureIndexFileExist(spacedkcPath);
		ensureIndexFileExist(ciPath);
		
		Path outDir=Paths.get(workDir);
		tryCreateWorkDir(outDir);
		
		Path outVcfPath=Paths.get(workDir,projectName+".sv.vcf.gz");
		Path mergedDataPath=Paths.get(workDir,projectName+".bp.dat");
		
		Path[] tempDataList=createTempDataPaths(nthreads,workDir,projectName);
		Path[] tempIndexList=createTempIndexPaths(nthreads,workDir,projectName);

		long time0=System.currentTimeMillis();  
		
		System.out.print("Loading Reference...");
		ReadClassifier classifier=new ReadClassifier(r1Path,r2Path,kcPath,longkcPath,spacedkcPath,nthreads,tempDataList,tempIndexList,
				                                     r1BarStart,r2BarStart,r1BarLen,r2BarLen,r1InsStart,r2InsStart);
		classifier.loadRef();
		System.out.println("Done!");
		
		long time1=System.currentTimeMillis();  
		long interval1=(time1-time0)/1000;  
		System.out.println("Elapsed time: "+interval1+"s");
    	
		System.out.print("K-mer Based Searching...");
    	int[] blockCntList=classifier.classify();
    	classifier.clean();
    	classifier=null;
    	System.out.println("Done!");
    	
    	long time2=System.currentTimeMillis();  
		long interval2=(time2-time1)/1000;  
		System.out.println("Elapsed time: "+interval2+"s");
 
    	System.out.print("Sorting Candidates...");
    	ResultsSorter sorter=new ResultsSorter(tempDataList,tempIndexList,nthreads,mergedDataPath,blockCntList);
    	BpIndexList finalIndexlist=null;
    	try {
        	finalIndexlist=sorter.sortResults();
    	}finally {
    		tryDeleteFiles(tempDataList);
        	tryDeleteFiles(tempIndexList);
    	}
    	
    	tempDataList=null;
    	tempIndexList=null;
    	sorter=null;	
    	System.out.println("Done!");
    	
    	long time3=System.currentTimeMillis();  
		long interval3=(time3-time2)/1000;  
		System.out.println("Elapsed time: "+interval3+"s");
    	
    	System.out.print("Clustering Breakpoint Candidates...");
    	BpClusterManager mergeManager=new BpClusterManager(ttPath,workDir,projectName,nthreads,mergeK,similarity,mergedDataPath);
    	SVCollection svCollection=null;
    	
    	try {
        	svCollection=mergeManager.mergeBreakpoints(finalIndexlist);
    	}finally {
    		tryDeleteAFile(mergedDataPath);
    	}
    	
    	mergedDataPath=null;
    	mergeManager=null;
    	
    	svCollection.merge();
    	Object[] vcfRecordMatrix=svCollection.transAndFilter(r1BarLen>0&&r2BarLen>0,minIndelLen,maxBpScore);
    	
    	svCollection.clean();
    	svCollection=null;
    	System.out.println("Done!");
    	
    	long time4=System.currentTimeMillis();  
		long interval4=(time4-time3)/1000;  
		System.out.println("Elapsed time: "+interval4+"s");
    	
    	System.out.print("Filtering And Saving...");
    	
    //	fw.write("Loading Reference:    "+Long.toString(interval1)+"s\n");
    //	fw.write("Classifying Reads:    "+Long.toString(interval2)+"s\n");
    //	fw.write("Merging Reads:    "+Long.toString(interval3)+"s\n");
    //	fw.write("Merging Breakpoints:    "+Long.toString(interval4)+"s\n");
    //	fw.write("threads:    "+Integer.toString(nthreads)+"s\n");

    	
    	VCFWriter vcfWriter=new VCFWriter(ciPath,outVcfPath,projectName,commands,vcfRecordMatrix);
    	try {
    		vcfWriter.sortAndWriteVcf();
    		vcfWriter.flush();
    	}finally {
    		vcfWriter.close();
    	}
    	vcfWriter=null;
    	for(int i=0,len=vcfRecordMatrix.length;i<len;++i) {
			ArrayList<VCFRecord> vcfList=(ArrayList<VCFRecord>)vcfRecordMatrix[i];
			vcfList.clear();
			vcfRecordMatrix[i]=null;
		}
    	vcfRecordMatrix=null;
    	System.out.println("Done!");
    	
	}
	
	private void databaseBuild(String fasta,String vcf,String save,int nWorkers,int k,int longk,int spacedk,int jump,int minSegLen,int maxSegLen,int geneCodeLen) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException {
		Path fastaPath=Paths.get(fasta);
		Path faiPath=Paths.get(fasta+".fai");
		Path vcfPath=Paths.get(vcf);
		ensureIndexFileExist(fastaPath);
		ensureIndexFileExist(faiPath);
		ensureIndexFileExist(vcfPath);
		
		Path ttPath=Paths.get(save+".tt");
		Path kcPath=Paths.get(save+".km");
		Path longkcPath=Paths.get(save+".long.km");
		Path spacedkcPath=Paths.get(save+".spaced.km");
		Path ciPath=Paths.get(save+".ci");
		
		DatabaseBuilder builder=new DatabaseBuilder(fastaPath,faiPath,vcfPath,ttPath,kcPath,longkcPath,spacedkcPath,ciPath,nWorkers,k,longk,spacedk,jump,minSegLen,maxSegLen,geneCodeLen);
		builder.buildDatabase();
	}
	
	public static void main(String[] args) {
		System.out.println("Aperture (Version: "+APERTURE_VERSION+")");
		ApertureMain aperture=new ApertureMain();
		if(args.length==0) {
			aperture.showApertureHelp();
		}else if(args[0].equals("index")) {
			try {
				aperture.parseArgForIndex(args,aperture);
			}catch(Exception e){
				e.printStackTrace();
			}
			aperture.showArgs(args);
		}else if(args[0].equals("call")) {
			try {
				aperture.parseArgForCall(args,aperture);
			}catch(Exception e) {
				e.printStackTrace();
			}
			aperture.showArgs(args);
		}else {
			System.out.println("Invalid command: " + args[0]+" !");
			aperture.showApertureHelp();
		}
		
	}
	
	private void showApertureHelp() {
		System.out.println("Aperture Help");
		System.out.println("Description: Alignment-free detection of structural variations and viral integrations in circulating tumor DNA");
		System.out.println("Contact: Hongchao Liu <i8q8r9@live.com>");
		System.out.println("");
		System.out.println("Usage: java -jar aperture.jar <command> <arguments>");
		System.out.println("");
		System.out.println("Commands:");
		System.out.println("");
		System.out.println("    index        Build index for Aperture");
		System.out.println("    call         Discover SVs and viral integrations");
		System.out.println("");
	}
	
	private void parseArgForIndex(String[] args,ApertureMain aperture) throws ParseException {
		Options options=new Options();
		
		Option opt = new Option("h", "help", false, "Show help message");
        opt.setRequired(false);
        options.addOption(opt);
        
        opt = new Option("R", "reference", true, "Genome FASTA file with fai index");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("V", "vcf", true, "Common SNPs from dbSNP database in VCF format");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("O", "out", true, "Output path");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("T", "threads", true, "Number of threads");
        opt.setType(Integer.TYPE);
        opt.setRequired(true);
        options.addOption(opt);
        
        CommandLine commandLine = null;
        CommandLineParser parser = new DefaultParser();
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(150);
        
        try {
        	commandLine = parser.parse(options, args);    	
        	if (commandLine.hasOption('h') || commandLine.hasOption("help")) {
                hf.printHelp("java -jar aperture.jar index -R <genome.fa> -V <snp.vcf> -O <out> -T <threads>", options, false);
            }else {
            	try {
            		long time1=System.currentTimeMillis();
    				aperture.databaseBuild(commandLine.getOptionValue("R"),commandLine.getOptionValue("V"),commandLine.getOptionValue("O"),Integer.parseInt(commandLine.getOptionValue("T")),23,41,83,3,30000,65000,5);
    				long time2=System.currentTimeMillis();  
    		    	long interval=(time2-time1)/1000;  
    		    	System.out.println("Index building workflow successfully completed.");
    		    	System.out.println("Elapsed time: "+interval+"s");			
    			}catch(Exception e){
    				e.printStackTrace();
    			}	
            }
        	
        }catch(ParseException e) {
        	System.out.println( "Unexpected exception:" + e.getMessage() );
        	hf.printHelp("java -jar aperture.jar index -R <genome.fa> -V <snp.vcf> -O <out> -T <threads>", options, false);
        }
	}
	
	
	private void parseArgForCall(String[] args,ApertureMain aperture) throws ParseException {
		Options options=new Options();
		
		Option opt = new Option("H", "help", false, "Show help message");
        opt.setRequired(false);
        options.addOption(opt);
        
        if(debug) {
        	opt = new Option("LC",  true, "debug:left code");
            opt.setRequired(false);
            options.addOption(opt);
            opt = new Option("LP",  true, "debug:left pos");
            opt.setRequired(false);
            options.addOption(opt);
            opt = new Option("RC", true, "debug:right code");
            opt.setRequired(false);
            options.addOption(opt);
            opt = new Option("RP", true, "debug:right pos");
            opt.setRequired(false);
            options.addOption(opt);
        }
        
        opt = new Option("1", "r1", true, "Path of R1.fq.gz");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("2", "r2", true, "Path of R2.fq.gz");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("I", "index", true, "Path of Aperture index files");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("D", "dir", true, "Output path");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("P", "project", true, "Project name");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("1BS", "r1BarStart", true, "Barcode start index in R1 (0-based)");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("2BS", "r2BarStart", true, "Barcode start index in R2 (0-based)");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("1BL", "r1BarLen", true, "Length of barcode in R1");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("2BL", "r2BarLen", true, "Length of barcode in R2");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("1S", "r1InsStart", true, "ctDNA fragment start index in R1 (0-based)");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("2S", "r2InsStart", true, "ctDNA fragment start index in R2 (0-based)");
        opt.setRequired(true);
        options.addOption(opt);
        
        opt = new Option("T", "threads", true, "Number of threads");
        opt.setType(Integer.TYPE);
        opt.setRequired(true);
        options.addOption(opt);
        
        CommandLine commandLine = null;
        CommandLineParser parser = new DefaultParser();
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(150);
        
        try {
        	commandLine = parser.parse(options, args);    	
        	if (commandLine.hasOption('h') || commandLine.hasOption("help")) {
                hf.printHelp("java -jar aperture.jar call", options, true);
            }else {
            	if(debug) {
                	if (commandLine.hasOption("LC")  && commandLine.hasOption("RC")) {
                		lCode=Integer.parseInt(commandLine.getOptionValue("LC"));
                        rCode=Integer.parseInt(commandLine.getOptionValue("RC"));
                        lPos=Integer.parseInt(commandLine.getOptionValue("LP"));
                        rPos=Integer.parseInt(commandLine.getOptionValue("RP"));
                    }
            	}
            	try {
    			    String indexPath=commandLine.getOptionValue("I");
    				int r1BarStart=Integer.parseInt(commandLine.getOptionValue("1BS"));
    				int r2BarStart=Integer.parseInt(commandLine.getOptionValue("2BS"));
    				int r1BarLen=Integer.parseInt(commandLine.getOptionValue("1BL"));
    				int r2BarLen=Integer.parseInt(commandLine.getOptionValue("2BL"));
    				int r1InsStart=Integer.parseInt(commandLine.getOptionValue("1S"));
    				int r2InsStart=Integer.parseInt(commandLine.getOptionValue("2S"));
    				int threads=Integer.parseInt(commandLine.getOptionValue("T"));
    				
    				long time1=System.currentTimeMillis();		
    				aperture.readClassify(commandLine.getOptionValue("1"),commandLine.getOptionValue("2"),indexPath,
    						commandLine.getOptionValue("D"),commandLine.getOptionValue("P"),
    						r1BarStart,r2BarStart,r1BarLen,r2BarLen,r1InsStart,r2InsStart,threads,11,50,0.25,3.0,args);
    				long time2=System.currentTimeMillis();  
    		    	long interval=(time2-time1)/1000;   	
    				System.out.println("SV calling workflow successfully completed.");
    				System.out.println("Elapsed time: "+interval+"s");
    			}catch(Exception e){
    				e.printStackTrace();
    			}	
            }
        	
        }catch(ParseException e) {
        	System.out.println( "Unexpected exception:" + e.getMessage() );
        	hf.printHelp("java -jar aperture.jar call", options, true);
        }
	}
	
	private void showArgs(String[] args) {
		System.out.print("Args: java -jar aperture.jar ");
		for(String arg:args) {
    		System.out.print(arg+" ");
    	}
    	System.out.println();
	}
	
	private void ensureIndexFileExist(Path path) throws FileNotFoundException {
		if(!Files.exists(path)) {
			throw new FileNotFoundException("Index file: "+path.toString()+" cannot be found!");
		}
	}
	
	private void ensureFastqFileExist(Path path) throws FileNotFoundException {
		if(!Files.exists(path)) {
			throw new FileNotFoundException("FastQ file: "+path.toString()+" cannot be found!");
		}
	}
	
	private Path[] createTempDataPaths(int threads,String workDir,String projectName) {
		Path[] pathList=new Path[threads];
		for(int i=0;i<threads;++i) {
			pathList[i]=Paths.get(workDir, projectName + "." + i + ".tmp");
		}
		return pathList;
	}
	
	private Path[] createTempIndexPaths(int threads,String workDir,String projectName) {
		Path[] pathList=new Path[threads];
		for(int i=0;i<threads;++i) {
			pathList[i]=Paths.get(workDir, projectName + "." + i + ".tmp.idx");
		}
		return pathList;
	}
	
	private void tryCreateWorkDir(Path dir) throws IOException {

		try{
			Files.createDirectory(dir);
			System.out.println("Working directory: "+dir.toString()+" is created");
		}catch(FileAlreadyExistsException e) {
			if(Files.isDirectory(dir)) {
				System.out.println("Working directory: "+dir.toString()+" already exists");
			}else {
				throw e;
			}
		}
	}
	
	private void tryDeleteFiles(Path[] pathList) {
		for(int i=0,len=pathList.length;i<len;++i) {
			if(pathList[i]!=null) {
				try {
					Files.delete(pathList[i]);
				}catch(IOException e) {
					e.printStackTrace();
				}finally {
					pathList[i]=null;
				}
			}
		}
	}
	
	private void tryDeleteAFile(Path path) {
		try {
			Files.delete(path);
		}catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	private TranslationTable readTranslationTable(Path ttPath) throws IOException, ClassNotFoundException {
		ObjectInputStream objIn=null;
		try {
			objIn=new ObjectInputStream(new GZIPInputStream(new FileInputStream(ttPath.toFile())));
			return (TranslationTable)objIn.readObject();
		}finally {
			if(objIn!=null) {
				try{
					objIn.close();
				}catch(IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
}
