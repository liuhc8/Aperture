import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPOutputStream;



public class DatabaseBuilder {
	public final static long HUMAN_GENOME_LEN=3000000000L;
	
	private final int threads,jump,k,longk,spacedk,minSegLen,maxSegLen,geneCodeLen;
    private final Path fastaPath,faiPath,vcfPath,ttPath,kcPath,longkcPath,spacedkcPath,ciPath;
    
    DatabaseBuilder(Path fastaPath,Path faiPath,Path vcfPath,Path ttPath,Path kcPath,Path longkcPath,Path spacedkcPath,Path ciPath,int threads,int k,int longk,int spacedk,int jump,int minSegLen,int maxSegLen,int geneCodeLen){
    	this.fastaPath=fastaPath;
    	this.faiPath=faiPath;
    	this.vcfPath=vcfPath;
    	this.ttPath=ttPath;
    	this.kcPath=kcPath;
    	this.longkcPath=longkcPath;
    	this.spacedkcPath=spacedkcPath;
    	this.ciPath=ciPath;
    	this.threads=threads;
    	this.k=k;
    	this.longk=longk;
    	this.spacedk=spacedk;
    	this.jump=jump;
    	this.minSegLen=minSegLen;
    	this.maxSegLen=maxSegLen;
    	this.geneCodeLen=geneCodeLen;
    }
	
	private FastaReader loadFasta() throws IllegalPositionException, IllegalBioFileException, IOException  {
		return new FastaReader(this.fastaPath,this.faiPath);
	}
	
	private VCFReader loadVCF() throws IOException {
		return new VCFReader(this.vcfPath);
	}
	
	public void buildDatabase() throws IllegalPositionException, IllegalBioFileException,IOException, InterruptedException {
		FastaReader faReader=loadFasta();
		TranslationTable trTable=null;
		FileOutputStream kmerOut=null;
		
		try {
    		System.out.print("Detecting repetitive k-mers...");
    		SimpleKmerCollection repeatKmer=getRepeatKmer(faReader);
    		System.out.println("Done!");
 		
	    	System.out.print("Cutting genome into segements...");
	    	trTable=cutGenome(faReader,repeatKmer);
	    	repeatKmer.clean();
	    	repeatKmer=null;
    		System.out.println("Done!");
    		
    		ObjectOutputStream trTableOut=null;
    		try {
    			trTableOut=new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(this.ttPath.toFile())));		
    			System.out.print("Saving translation table...");
    			trTableOut.writeObject(trTable);	
    			trTable.cleanRef();
    		
    		}finally {
    			if(trTableOut!=null) {
    				try {
    					trTableOut.close();
    				}catch(Exception e) {
    					e.printStackTrace();
    				}
    			}
    		}
		
    		ObjectOutputStream chromInfoOut=null;
    		try {
    			chromInfoOut=new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(this.ciPath.toFile())));
    			ChromInfo ci=faReader.getChromInfo();
    			chromInfoOut.writeObject(ci);	
    			System.out.println("Done!");
    		
    		}finally {
    			if(chromInfoOut!=null) {
    				try {
    					chromInfoOut.close();
    				}catch(Exception e) {
    					e.printStackTrace();
    				}
    			}
    		}
    		
    		System.out.print("Building k-mer library...");
	    	KmerCollection kc=insertKmer(trTable,faReader);
	    	
	    	VCFReader vcfReader=null;
	    	try {
	    		vcfReader=loadVCF();
	    		addSNPs(trTable,kc,faReader,vcfReader);
	    	}finally {
	    		if(vcfReader!=null) {
	    			vcfReader.close();
	    		}
	    	}
	    	
	    	compact(kc);
    		updateKmerPos(kc,trTable,faReader);
    		
    		try {
    			vcfReader=loadVCF();
    			updateSNPsPos(trTable,kc,faReader,vcfReader);
    		}finally {
    			if(vcfReader!=null) {
	    			vcfReader.close();
	    		}
    		}

    		setKmerScore(kc,trTable,faReader);
	    	System.out.println("Done!");
    		
	    	System.out.print("Building long k-mer library...");
	    	LongKmerCollection longkc=addLongKmer(trTable,kc,faReader);
	    	updateLongKmerPos(longkc,trTable,faReader);
	    	
	    	try {
	    		vcfReader=loadVCF();
	    		updateLongKmerSNPsPos(longkc,trTable,faReader, vcfReader);
	    	}finally {
    			if(vcfReader!=null) {
	    			vcfReader.close();
	    		}
    		}
	    	
	    	System.out.println("Done!");
	    	
	    	try {
	    		kmerOut=new FileOutputStream(this.longkcPath.toFile());
	    		FileChannel kmerFc=kmerOut.getChannel();
			
	    		System.out.print("Saving long k-mer library...");
	    		longkc.saveToDisk(kmerFc);
	    		System.out.println("Done!");
			
    		}finally {
	    		if(kmerOut!=null) {
		    		try {
		    			kmerOut.close();
		    		}catch(Exception e) {
		    			e.printStackTrace();
		    		}
	    		}
    		}
		
     		longkc.clean();
     		longkc=null;
	    	   	
	    	System.out.print("Building spaced seed library...");
	    	LongKmerCollection spacedkc=addSpacedSeed(trTable,kc,faReader);
	    	updateSpacedSeedPos(spacedkc,trTable,faReader);
	    	
	    	try {
	    		vcfReader=loadVCF();
	    		updateSpacedSeedSNPsPos(spacedkc,trTable,faReader, vcfReader);
	    	}finally {
    			if(vcfReader!=null) {
	    			vcfReader.close();
	    		}
    		}
	    	
	    	System.out.println("Done!");
	    	
	    	try {
	    		kmerOut=new FileOutputStream(this.spacedkcPath.toFile());
	    		FileChannel kmerFc=kmerOut.getChannel();
			
	    		System.out.print("Saving spaced seed library...");
	    		spacedkc.saveToDisk(kmerFc);
	    		System.out.println("Done!");
			
    		}finally {
	    		if(kmerOut!=null) {
		    		try {
		    			kmerOut.close();
		    		}catch(Exception e) {
		    			e.printStackTrace();
		    		}
	    		}
    		}
		
	    	spacedkc.clean();
     		longkc=null;
		
    		try {
	    		kmerOut=new FileOutputStream(this.kcPath.toFile());
	    		FileChannel kmerFc=kmerOut.getChannel();
			
	    		System.out.print("Saving k-mer library...");
	    		kc.saveToDisk(kmerFc);
	    		System.out.println("Done!");
			
    		}finally {
	    		if(kmerOut!=null) {
		    		try {
		    			kmerOut.close();
		    		}catch(Exception e) {
		    			e.printStackTrace();
		    		}
	    		}
    		}
		
     		kc.clean();
     		kc=null;
		
    	}finally {
    		if(faReader!=null) {
	        	faReader.close();
	    	}
    	}
     	
	}
	
	private List<Integer> getNPos(byte[] seq){
		List<Integer> array=new ArrayList<Integer>();
		int len=seq.length;
		int start=-1,end=-1;
		for(int i=0;i<len;++i) {
			if(seq[i]==78) {
				if(start!=-1) {
					end=i;
					array.add(end);
					start=-1;
					end=-1;
				}
			}else {
				if(start==-1) {
					start=i;
					array.add(start);
				}
			}
		}
		if((array.size() & 1)==1) {
			array.add(len);
		}
		return array;
	}
	
	private SimpleKmerCollection getRepeatKmer(FastaReader faReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException {
		final int k=this.k;
		final int jump=this.jump;
		final ArrayList<String> chromList=faReader.getChromList();
		final DNASeqWindow refWindow=new DNASeqWindow(k);
		
		SimpleKmerCollection kmerCollect=new SimpleKmerCollection(k);
		
		for(String chrom:chromList) {
			byte[] seq=faReader.getSequence(chrom).getSequence();
			List<Integer> array=getNPos(seq);
			for(int i=0;i<array.size();i=i+2) {
				refWindow.set(seq, array.get(i), jump,array.get(i+1));
				while(refWindow.hasNextKmer()) {
					kmerCollect.insert(refWindow.nextKmer());
				}
			}
		}
		kmerCollect.compact(this.threads);
		
		
		return kmerCollect;
	}
	
	private TranslationTable cutGenome(FastaReader faReader,SimpleKmerCollection repeatKmer) throws IllegalPositionException, IllegalBioFileException, IOException{
		final int k=this.k;
		final int longk=this.longk;
		final int minSegLen=this.minSegLen;
		final int maxSegLen=this.maxSegLen;
		final int geneCodeLen=this.geneCodeLen;
		final ArrayList<String> chromList=faReader.getChromList();
		final DNASeqWindow refWindow=new DNASeqWindow(k);
		final GeneCodeGenerator codeGenerator=new GeneCodeGenerator(geneCodeLen);
		
		TranslationTable trTable=new TranslationTable(300000,k,longk,geneCodeLen);

		for(String chrom:chromList) {
			byte[] seq=faReader.getSequence(chrom).getSequence();
			List<Integer> array=getNPos(seq);
			for(int i=0;i<array.size();i=i+2) {
				int start=array.get(i);
				int end=array.get(i+1);
		inner:	while(end-start>maxSegLen) {
					refWindow.set(seq,start+minSegLen,1,end);
					
					int cnt=0;
					for(int move=0;move<maxSegLen-minSegLen && refWindow.hasNextKmer();++move) {
						if(!repeatKmer.find(refWindow.nextKmer())) {
							++cnt;
							if(cnt>=2) {
						    	addSeperator(codeGenerator,trTable,faReader,chrom,start,start+minSegLen+move+k-1);
						    	start+=minSegLen+move;
						    	continue inner;
							}
						}else {
							cnt=0;
						}
					}
					addSeperator(codeGenerator,trTable,faReader,chrom,start,start+maxSegLen+k-1);
					start+=maxSegLen;
				}
				addSeperator(codeGenerator,trTable,faReader,chrom,start,end);
			}
		}
		trTable.compact();
		
		return trTable;
	}
	
	private void addSeperator(GeneCodeGenerator codeGenerator,TranslationTable trTable,FastaReader faReader,String chrom,int start,int end) throws IllegalPositionException, IllegalBioFileException, IOException {		
		RefSequence ref=faReader.getSequence(chrom,start+1,end+300);
		long[] binSeq=DNASequence.compressSeq(ref);
		int code=codeGenerator.getValidCode();
		trTable.addSeperator(chrom,start,end,code,binSeq);
	}
	
	private KmerCollection insertKmer(TranslationTable trTable,FastaReader faReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException  {
		final int k=this.k;
		final int jump=this.jump;
		final int geneCodeLen=this.geneCodeLen;
		final KmerCollection kc=new KmerCollection(HUMAN_GENOME_LEN/jump,k,geneCodeLen);
		final DNASeqWindow refWindow=new DNASeqWindow(k);
		
		int len=trTable.size();
		for(int i=0;i<len;++i) {
			Segment seg=trTable.getSegment(i);
			int geneCode=seg.code;
			RefSequence ref=faReader.getSequence(seg.chrom,seg.start+1,seg.end);

			refWindow.set(ref,0,jump);
			
			for(int pos=0;refWindow.hasNextKmer();pos+=jump) {
				kc.insert(refWindow.nextKmer(),geneCode,pos);
			}
		}

		//kc.compact(this.nWorkers);     
		//System.out.println("size::"+kc.size());
		return kc;
	}
	
	private void addSNPs(TranslationTable trTable,KmerCollection kc,FastaReader faReader,VCFReader vcfReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException  {
		final int k=this.k;
		final int jump=this.jump;
		final int threads=this.threads;
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new AddSNPsHelper(k,jump,faReader,vcfReader,trTable,kc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
		//kc.compact(this.nWorkers); 	
		//System.out.println("size::"+kc.size());
	}
	
	private void compact(KmerCollection kc) throws InterruptedException {
		kc.compact(threads);
	}
	
	private class AddSNPsHelper implements Runnable{
		private int k;
		private int jump;
		private DNASeqWindow refWindow;
		private FastaReader faReader;
		private VCFReader vcfReader;
		private TranslationTable trTable;
		private KmerCollection kc;
		
		public AddSNPsHelper(int k,int jump,FastaReader faReader,VCFReader vcfReader,TranslationTable trTable,KmerCollection kc) throws IllegalPositionException, IllegalBioFileException, IOException{
			this.k=k;
			this.jump=jump;
			this.faReader=faReader;
			this.vcfReader=vcfReader;
			this.trTable=trTable;
			this.refWindow=new DNASeqWindow(k);
			this.kc=kc;
		}
		
		@Override
		public void run() {
			try {
	     		final TranslationTable trTable=this.trTable;
	     		final KmerCollection kc=this.kc;
	    		final DNASeqWindow refWindow=this.refWindow;
	    		final FastaReader faReader=this.faReader;
	    		final VCFReader vcfReader=this.vcfReader;
	    		final int k=this.k;
			    final int jump=this.jump;
	    		
			    byte[] seq=new byte[220];
			    long[] kmerList=new long[200];
			    RefSequence ref=null;
			    SNP[] snpArray=null;
			    Segment lastSeg=null;
			    while((snpArray=vcfReader.readVCF())!=null) {
			    	for(SNP snp:snpArray) {
			    		if(snp==null) {
			    			break;
			    		}
			    		Segment seg=trTable.findSegment(snp.chrom, snp.pos);
			    		if(seg==null) {
			    			continue;
			    		}
			    		
			    		if(seg!=lastSeg) {
			    	    	ref=faReader.getSequence(seg.chrom,seg.start+1,seg.end);
			    		}
			    		int pos=snp.pos-seg.start-1;
			    		
			    		int refLen=snp.ref.length;
			    		int altLen=snp.alt.length;
			    		if(pos<k-1 || snp.pos+refLen+k-1>seg.end || refLen>150 || altLen>150) {
			    			continue;
			    		}
			    			
			    		ref.copyBases(seq, 0, pos-k+1, k-1);
			    		System.arraycopy(snp.alt, 0, seq, k-1, altLen);
			    		ref.copyBases(seq, altLen+k-1, pos+refLen, k-1);
			    		
			    		refWindow.set(seq,0,1,altLen+2*k-2);
			    		for(int i=0;refWindow.hasNextKmer();++i) {
			    			kmerList[i]=refWindow.nextKmer();
			    		}
			    		
			    		int recentCode=seg.code;
			    		int kmerListLen=altLen+k-1;
			    		int leftStartPos=pos-k+1;
			    		int leftStartRemain=leftStartPos%jump;
			    		int leftStart=leftStartRemain!=0?jump-leftStartRemain:0;
			    		
			    		if(refLen==altLen) {
				    		for(int i=leftStart;i<kmerListLen;i+=jump) {
				    			kc.insert(kmerList[i],recentCode,leftStartPos+i);
				    		}
			    		}else {
			    			for(int i=leftStart;i<kmerListLen;i+=jump) {
				    			kc.insert(kmerList[i],recentCode,-1);
				    		}
			    			kc.insert(kmerList[kmerListLen-1],recentCode,-1);
			    		}
			    		
			    		lastSeg=seg;
			    	}
			    }
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	private LongKmerCollection addLongKmer(TranslationTable trTable,KmerCollection kc,FastaReader faReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException  {
		final int k=this.k,longk=this.longk;
		final int geneCodeLen=this.geneCodeLen;
		final int nWorkers=this.threads;
		final int jump=this.jump;
		
		LongKmerCollection longkc=new LongKmerCollection(100000000,9,longk,geneCodeLen);
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(nWorkers);
		for(int i=0;i<nWorkers;++i) {
			fixedThreadPool.execute(new AddLongKmerHelper(i,nWorkers,k,longk,jump,faReader,trTable,kc,longkc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
		longkc.compact(nWorkers);
		
		//System.out.println("longkc size::"+longkc.size());
		return longkc;
	}
	
	private LongKmerCollection addSpacedSeed(TranslationTable trTable,KmerCollection kc,FastaReader faReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException  {
		final int k=this.k,spacedk=this.spacedk;
		final int geneCodeLen=this.geneCodeLen;
		final int threads=this.threads;
		final int jump=this.jump;
		
		LongKmerCollection spacedkc=new LongKmerCollection(100000000,9,spacedk,geneCodeLen);
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new AddLongKmerHelper(i,threads,k,spacedk,jump,faReader,trTable,kc,spacedkc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
		spacedkc.compact(threads);
		
		//System.out.println("spacedkc size::"+spacedkc.size());
		return spacedkc;
	}
	
	private class AddLongKmerHelper implements Runnable{
		private int no,threads,k,longk,jump,kmerLMove1,kmerLMove2;
		private DNASeqWindow refWindowL,refWindowR;
		private FastaReader faReader;
		private TranslationTable trTable;
		private KmerCollection kc;
		private LongKmerCollection longkc;
		
		public AddLongKmerHelper(int no,int threads,int k,int longk,int jump,FastaReader faReader,TranslationTable trTable,KmerCollection kc,LongKmerCollection longkc) throws IllegalPositionException, IllegalBioFileException, IOException{
			this.no=no;
			this.threads=threads;
			this.k=k;
			this.longk=longk;
			this.jump=jump;
			this.faReader=faReader;
			this.trTable=trTable;
			this.refWindowL=new DNASeqWindow(k);
			this.refWindowR=new DNASeqWindow(k);
			this.kc=kc;
			this.longkc=longkc;
			if(longk<k*2) {
				this.kmerLMove1=(k-(longk-32))*2;
	    		this.kmerLMove2=(longk-k)*2;
			}else {
				this.kmerLMove1=(32-k)*2;
	    		this.kmerLMove2=k*2;
			}
		}
		
		@Override
		public void run() {
			try {
	     		final TranslationTable trTable=this.trTable;
	     		final KmerCollection kc=this.kc;
	     		final LongKmerCollection longkc=this.longkc;
	    		final DNASeqWindow refWindowL=this.refWindowL,refWindowR=this.refWindowR;
	    		final FastaReader faReader=this.faReader;
	    		final int threads=this.threads;
	    		final int k=this.k,longk=this.longk,jump=this.jump,kmerLMove1=this.kmerLMove1,kmerLMove2=this.kmerLMove2;
	    		
	    		int len=trTable.size();
	    		
	    		for(int i=this.no;i<len;i+=threads) {
	    			Segment seg=trTable.getSegment(i);
	    			int geneCode=seg.code;
	    			RefSequence ref=faReader.getSequence(seg.chrom,seg.start+1,seg.end);

	    			refWindowL.set(ref,0,jump);
	    			refWindowR.set(ref,longk-k,jump);
	    			for(int pos=0;refWindowR.hasNextKmer();pos+=jump) {
	    				long kmerL=refWindowL.nextKmer();
	    				long kmerR=refWindowR.nextKmer();
	    				
	    				long findresL=kc.find(kmerL);
	    				int kmerCodeL=(int)(findresL>>>32);
						int kmerPosL=(int)findresL;
						
						long findresR=kc.find(kmerR);
	    				int kmerCodeR=(int)(findresR>>>32);
						int kmerPosR=(int)findresR;
						
						if(kmerPosL==0xFFFF && kmerCodeL!=0 && kmerPosR==0xFFFF && kmerCodeR!=0) {
							longkc.insert((int)(kmerL>>>kmerLMove1), (kmerL<<kmerLMove2)|kmerR, geneCode, pos);
	    				}
	    			}
	    		}
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	
	private void updateKmerPos(KmerCollection kc,TranslationTable trTable,FastaReader faReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException {
		final int k=this.k;
		final int threads=this.threads;
		
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new UpdateHelper(i,threads,k,faReader,trTable,kc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
	}
	
	private class UpdateHelper implements Runnable{
		private int no,threads;
		private DNASeqWindow refWindow;
		private FastaReader faReader;
		private TranslationTable trTable;
		private KmerCollection kc;
		
		public UpdateHelper(int no,int threads,int k,FastaReader faReader,TranslationTable trTable,KmerCollection kc) throws IllegalPositionException, IllegalBioFileException, IOException{
			this.no=no;
			this.threads=threads;
			this.faReader=faReader;
			this.trTable=trTable;
			this.refWindow=new DNASeqWindow(k);
			this.kc=kc;
		}
		
		@Override
		public void run() {
			try {
	     		final TranslationTable trTable=this.trTable;
	     		final KmerCollection kc=this.kc;
	    		final DNASeqWindow refWindow=this.refWindow;
	    		final FastaReader faReader=this.faReader;
	    		final int threads=this.threads;
			
	    		int len=trTable.size();
	    		
		    	for(int i=this.no;i<len;i+=threads) {
	    			Segment seg=trTable.getSegment(i);
	    			int geneCode=seg.code;
	    			RefSequence ref=faReader.getSequence(seg.chrom,seg.start+1,seg.end);;
	    			
	    			refWindow.set(ref,0);
	    			
		    		for(int j=0;refWindow.hasNextKmer();++j) {
		    			long kmer=refWindow.nextKmer();
		    			long rev=refWindow.showThisAsRevComp();
		    			if(j%3==0) {
			    			kc.update(rev,geneCode);
		    			}else {
		    				kc.update(kmer,geneCode);
		    				kc.update(rev,geneCode);
		    			}
	    			}
	    		}
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	private void updateSNPsPos(TranslationTable trTable,KmerCollection kc,FastaReader faReader,VCFReader vcfReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException  {
		final int k=this.k;
		final int jump=this.jump;
		final int threads=this.threads;
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new UpdateSNPsHelper(k,jump,faReader,vcfReader,trTable,kc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }

	}
	
	private class UpdateSNPsHelper implements Runnable{
		private int k,jump;
		private DNASeqWindow refWindow;
		private FastaReader faReader;
		private VCFReader vcfReader;
		private TranslationTable trTable;
		private KmerCollection kc;
		
		public UpdateSNPsHelper(int k,int jump,FastaReader faReader,VCFReader vcfReader,TranslationTable trTable,KmerCollection kc) throws IllegalPositionException, IllegalBioFileException, IOException{
			this.k=k;
			this.jump=jump;
			this.faReader=faReader;
			this.vcfReader=vcfReader;
			this.trTable=trTable;
			this.refWindow=new DNASeqWindow(k);
			this.kc=kc;
		}
		
		@Override
		public void run() {
			try {
	     		final TranslationTable trTable=this.trTable;
	     		final KmerCollection kc=this.kc;
	    		final DNASeqWindow refWindow=this.refWindow;
	    		final FastaReader faReader=this.faReader;
	    		final VCFReader vcfReader=this.vcfReader;
	    		final int k=this.k,jump=this.jump;
	    		
			    byte[] seq=new byte[220];
			    long[] kmerList=new long[200];
			    long[] revKmerList=new long[200];
			    RefSequence ref=null;
			    SNP[] snpArray=null;
			    Segment lastSeg=null;
			    while((snpArray=vcfReader.readVCF())!=null) {
			    	for(SNP snp:snpArray) {
			    		if(snp==null) {
			    			break;
			    		}
			    		Segment seg=trTable.findSegment(snp.chrom, snp.pos);
			    		if(seg==null) {
			    			continue;
			    		}
			    		
			    		if(seg!=lastSeg) {
			    	    	ref=faReader.getSequence(seg.chrom,seg.start+1,seg.end);
			    		}
			    		int pos=snp.pos-seg.start-1;
			    		
			    		int refLen=snp.ref.length;
			    		int altLen=snp.alt.length;
			    		if(pos<k-1 || snp.pos+refLen+k-1>seg.end || refLen>150 || altLen>150) {
			    			continue;
			    		}
			    			
			    		ref.copyBases(seq, 0, pos-k+1, k-1);
			    		System.arraycopy(snp.alt, 0, seq, k-1, altLen);
			    		ref.copyBases(seq, altLen+k-1, pos+refLen, k-1);
			    		
			    		refWindow.set(seq,0,1,altLen+2*k-2);
			    		for(int i=0;refWindow.hasNextKmer();++i) {
			    			kmerList[i]=refWindow.nextKmer();
			    			revKmerList[i]=refWindow.showThisAsRevComp();
			    		}
			    		
			    		int recentCode=seg.code;
			    		int kmerListLen=altLen+k-1;
			    		int leftStartPos=pos-k+1;
			    		int leftStartRemain=leftStartPos%jump;
			    		int leftStart=leftStartRemain!=0?jump-leftStartRemain:0;
			    		

				    	for(int i=0;i+leftStart<kmerListLen;++i) {
				    		if(i%3==0) {
				    			/*if(snp.chrom.equals("chr2") && snp.pos==42526854) {
				    				System.out.println("CXXXXXX");
				    				System.out.println(DNASequence.decompressKmer(revKmerList[i+leftStart],23));
				    				System.out.println("AAA"+"    "+Long.toHexString(kc.find(revKmerList[i+leftStart])));
				    			}*/
					    		kc.update(revKmerList[i+leftStart],recentCode);
					    		/*if(snp.chrom.equals("chr2") && snp.pos==42526854) {
				    				System.out.println(DNASequence.decompressKmer(revKmerList[i+leftStart],23));
				    				System.out.println("BBB"+"    "+Long.toHexString(kc.find(revKmerList[i+leftStart])));
				    			}*/
				    		}else {
				    			kc.update(kmerList[i+leftStart],recentCode);
				    			kc.update(revKmerList[i+leftStart],recentCode);
				    		}
				    	}
			    		
			    		lastSeg=seg;
			    	}
			    }
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	private void setKmerScore(KmerCollection kc,TranslationTable trTable,FastaReader faReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException {
		final int k=this.k;
		final int threads=this.threads;
		final int jump=this.jump;
		final int geneCodeLen=this.geneCodeLen;
		
		kc.prepareScoreMatrix();
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new SetScoreHelper(i,threads,k,jump,geneCodeLen,faReader,trTable,kc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
	}
	
	private class SetScoreHelper implements Runnable{
		private int no,nWorkers,k,jump,geneCodeLen;
		private DNASeqWindow refWindow;
		private FastaReader faReader;
		private TranslationTable trTable;
		private KmerCollection kc;
		
		public SetScoreHelper(int no,int nWorkers,int k,int jump,int geneCodeLen,FastaReader faReader,TranslationTable trTable,KmerCollection kc) throws IllegalPositionException, IllegalBioFileException, IOException{
			this.no=no;
			this.nWorkers=nWorkers;
			this.k=k;
			this.jump=jump;
			this.geneCodeLen=geneCodeLen;
			this.faReader=faReader;
			this.trTable=trTable;
			this.refWindow=new DNASeqWindow(k);
			this.kc=kc;
		}
		
		@Override
		public void run() {
			try {
	     		final TranslationTable trTable=this.trTable;
	     		final KmerCollection kc=this.kc;
	    		final DNASeqWindow refWindow=this.refWindow;
	    		final FastaReader faReader=this.faReader;
	    		final int nWorkers=this.nWorkers,k=this.k,jump=this.jump,geneCodeLen=this.geneCodeLen;
			
	    		int len=trTable.size();
	    		
		    	for(int i=this.no;i<len;i+=nWorkers) {
	    			Segment seg=trTable.getSegment(i);
	    			RefSequence ref=faReader.getSequence(seg.chrom,seg.start+1,seg.end);
	    			
	    			refWindow.set(ref,0,jump);
	    			
	    			long lastKmer=0;
	    			int lastRepeatCategory=0;
	    			while(refWindow.hasNextKmer()) {
	    				long kmer=refWindow.nextKmer();
	    				int repeatCategory=0;
	    				int code=kc.findGeneCode(kmer);
	    				if(code==0xFFFFFFFF) {
	    					repeatCategory=KmerCollection.FULL_REPETITIVE;
	    				}else if(GeneCodeGenerator.testGeneCode(code,geneCodeLen)>0) {
	    					repeatCategory=KmerCollection.LESS_REPETITIVE;
	    				}
	    				if(repeatCategory!=0) {
	    					kc.setScore(kmer, repeatCategory);
	    					long rev=refWindow.showThisAsRevComp();
	    					kc.setScore(rev, repeatCategory);
	    					
	    					if(lastRepeatCategory!=0) {
	    						long combinedKmer=(lastKmer<<(jump*2))|kmer;
	    						int combinedRepeatCategory=(repeatCategory==KmerCollection.FULL_REPETITIVE && lastRepeatCategory==KmerCollection.FULL_REPETITIVE)?
	    								KmerCollection.FULL_REPETITIVE:KmerCollection.LESS_REPETITIVE;
	    						for(int j=0;j<jump-1;++j) {
	    							combinedKmer>>=2;
	    							kc.setScore(combinedKmer,combinedRepeatCategory);
	    							rev=DNASequence.getRevComp(combinedKmer,k);
	    							kc.setScore(rev,combinedRepeatCategory);
	    						}
	    					}
	    				}
	    				lastKmer=kmer;
	    				lastRepeatCategory=repeatCategory;
	    			}
	    		}
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	private void updateLongKmerPos(LongKmerCollection longkc,TranslationTable trTable,FastaReader faReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException {
		final int k=this.k,longk=this.longk,threads=this.threads;
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new UpdateLongHelper(i,threads,k,longk,faReader,trTable,longkc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
	}
	
	private void updateSpacedSeedPos(LongKmerCollection spacedkc,TranslationTable trTable,FastaReader faReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException {
		final int k=this.k,spacedk=this.spacedk,threads=this.threads;
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new UpdateLongHelper(i,threads,k,spacedk,faReader,trTable,spacedkc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
	}
	
	private class UpdateLongHelper implements Runnable{
		private int no,nWorkers,k,longk,kmerLMove1,kmerLMove2;
		private DNASeqWindow refWindowL,refWindowR;
		private FastaReader faReader;
		private TranslationTable trTable;
		private LongKmerCollection longkc;
		
		public UpdateLongHelper(int no,int nWorkers,int k,int longk,FastaReader faReader,TranslationTable trTable,LongKmerCollection longkc) throws IllegalPositionException, IllegalBioFileException, IOException{
			this.no=no;
			this.nWorkers=nWorkers;
			this.k=k;
			this.longk=longk;
			this.faReader=faReader;
			this.trTable=trTable;
			this.refWindowL=new DNASeqWindow(k);
			this.refWindowR=new DNASeqWindow(k);
			this.longkc=longkc;
			if(longk<k*2) {
				this.kmerLMove1=(k-(longk-32))*2;
	    		this.kmerLMove2=(longk-k)*2;
			}else {
				this.kmerLMove1=(32-k)*2;
	    		this.kmerLMove2=k*2;
			}
		}
		
		@Override
		public void run() {
			try {
	     		final TranslationTable trTable=this.trTable;
	     		final LongKmerCollection longkc=this.longkc;
	    		final DNASeqWindow refWindowL=this.refWindowL,refWindowR=this.refWindowR;
	    		final FastaReader faReader=this.faReader;
	    		final int nWorkers=this.nWorkers,k=this.k,longk=this.longk,kmerLMove1=this.kmerLMove1,kmerLMove2=this.kmerLMove2;
			
	    		int len=trTable.size();
	    		
		    	for(int i=this.no;i<len;i+=nWorkers) {
	    			Segment seg=trTable.getSegment(i);
	    			int geneCode=seg.code;
	    			RefSequence ref=faReader.getSequence(seg.chrom,seg.start+1,seg.end);;
	    			
	    			refWindowL.set(ref,0);
	    			refWindowR.set(ref,longk-k);
		    		
	    			for(int j=0;refWindowR.hasNextKmer();++j) {
	    				long kmerL=refWindowL.nextKmer();
	    				long kmerR=refWindowR.nextKmer();
	    				long kmerLRev=refWindowL.showThisAsRevComp();
	    				long kmerRRev=refWindowR.showThisAsRevComp();
	    				
	    				
	    				if(j%3==0) {
	    					longkc.update((int)(kmerRRev>>>kmerLMove1), (kmerRRev<<kmerLMove2)|kmerLRev, geneCode);
		    			}else {
		    				longkc.update((int)(kmerL>>>kmerLMove1), (kmerL<<kmerLMove2)|kmerR, geneCode);
		    				longkc.update((int)(kmerRRev>>>kmerLMove1), (kmerRRev<<kmerLMove2)|kmerLRev, geneCode);
		    			}
	    			}
	    		}
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	
	private void updateLongKmerSNPsPos(LongKmerCollection longkc,TranslationTable trTable,FastaReader faReader,VCFReader vcfReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException {
		final int k=this.k,longk=this.longk,threads=this.threads,jump=this.jump;
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new UpdateSNPsLongHelper(k,longk,jump,faReader,vcfReader,trTable,longkc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
	}
	
	private void updateSpacedSeedSNPsPos(LongKmerCollection spacedkc,TranslationTable trTable,FastaReader faReader,VCFReader vcfReader) throws IllegalPositionException, IllegalBioFileException, IOException, InterruptedException {
		final int k=this.k,spacedk=this.spacedk,threads=this.threads,jump=this.jump;
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new UpdateSNPsLongHelper(k,spacedk,jump,faReader,vcfReader,trTable,spacedkc));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
	}
	
	private class UpdateSNPsLongHelper implements Runnable{
		private int k,longk,jump,kmerLMove1,kmerLMove2;;
		private DNASeqWindow refWindowL,refWindowR;
		private FastaReader faReader;
		private VCFReader vcfReader;
		private TranslationTable trTable;
		private LongKmerCollection longkc;
		
		public UpdateSNPsLongHelper(int k,int longk,int jump,FastaReader faReader,VCFReader vcfReader,TranslationTable trTable,LongKmerCollection longkc) throws IllegalPositionException, IllegalBioFileException, IOException{
			this.k=k;
			this.longk=longk;
			this.jump=jump;
			this.faReader=faReader;
			this.vcfReader=vcfReader;
			this.trTable=trTable;
			this.refWindowL=new DNASeqWindow(k);
			this.refWindowR=new DNASeqWindow(k);
			this.longkc=longkc;
			if(longk<k*2) {
				this.kmerLMove1=(k-(longk-32))*2;
	    		this.kmerLMove2=(longk-k)*2;
			}else {
				this.kmerLMove1=(32-k)*2;
	    		this.kmerLMove2=k*2;
			}
		}
		
		@Override
		public void run() {
			try {
	     		final TranslationTable trTable=this.trTable;
	     		final LongKmerCollection longkc=this.longkc;
	    		final DNASeqWindow refWindowL=this.refWindowL;
	    		final DNASeqWindow refWindowR=this.refWindowR;
	    		final FastaReader faReader=this.faReader;
	    		final VCFReader vcfReader=this.vcfReader;
	    		final int k=this.k,longk=this.longk,jump=this.jump,kmerLMove1=this.kmerLMove1,kmerLMove2=this.kmerLMove2;
	    		
			    byte[] seq=new byte[350];
			    long[] kmerListL=new long[330];
			    long[] revKmerListL=new long[330];
			    long[] kmerListR=new long[330];
			    long[] revKmerListR=new long[330];
			    RefSequence ref=null;
			    SNP[] snpArray=null;
			    Segment lastSeg=null;
			    while((snpArray=vcfReader.readVCF())!=null) {
			    	for(SNP snp:snpArray) {
			    		if(snp==null) {
			    			break;
			    		}
			    		Segment seg=trTable.findSegment(snp.chrom, snp.pos);
			    		if(seg==null) {
			    			continue;
			    		}
			    		
			    		if(seg!=lastSeg) {
			    	    	ref=faReader.getSequence(seg.chrom,seg.start+1,seg.end);
			    		}
			    		int pos=snp.pos-seg.start-1;
			    		
			    		int refLen=snp.ref.length;
			    		int altLen=snp.alt.length;
			    		if(pos<longk-1 || snp.pos+refLen+longk-1>seg.end || refLen>150 || altLen>150) {
			    			continue;
			    		}
			    			
			    		ref.copyBases(seq, 0, pos-longk+1, longk-1);
			    		System.arraycopy(snp.alt, 0, seq, longk-1, altLen);
			    		ref.copyBases(seq, altLen+longk-1, pos+refLen, longk-1);
			    		
			    		refWindowL.set(seq,0,1,altLen+2*longk-2);
			    		refWindowR.set(seq,longk-k,1,altLen+2*longk-2);
			    		for(int i=0;refWindowR.hasNextKmer();++i) {
			    			kmerListL[i]=refWindowL.nextKmer();
			    			revKmerListL[i]=refWindowL.showThisAsRevComp();
			    			kmerListR[i]=refWindowR.nextKmer();
			    			revKmerListR[i]=refWindowR.showThisAsRevComp();
			    		}
			    		
			    		int recentCode=seg.code;
			    		int kmerListLen=altLen+longk-1;
			    		int leftStartPos=pos-longk+1;
			    		int leftStartRemain=leftStartPos%jump;
			    		int leftStart=leftStartRemain!=0?jump-leftStartRemain:0;
			    		

				    	for(int i=0;i+leftStart<kmerListLen;++i) {
				    		if(i%3==0) {
					    		longkc.update((int)(revKmerListR[i+leftStart]>>>kmerLMove1), (revKmerListR[i+leftStart]<<kmerLMove2)|revKmerListL[i+leftStart], recentCode);
				    		}else {
				    			longkc.update((int)(kmerListL[i+leftStart]>>>kmerLMove1), (kmerListL[i+leftStart]<<kmerLMove2)|kmerListR[i+leftStart], recentCode);
				    			longkc.update((int)(revKmerListR[i+leftStart]>>>kmerLMove1), (revKmerListR[i+leftStart]<<kmerLMove2)|revKmerListL[i+leftStart], recentCode);
				    		}
				    	}
			    		
				    	lastSeg=seg;
			    	}
			    }
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
	}

}


