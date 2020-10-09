import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.nio.file.Paths;


public class ReadClassifier {
	private Path[] tempDataList,tempIndexList;
	private Path r1Path,r2Path,kcPath,longkcPath,spacedkcPath;
	private int nthreads;
	private int r1BarStart,r2BarStart,r1BarLen,r2BarLen,r1InsStart,r2InsStart;
	
	private KmerCollection kc;
	private LongKmerCollection longkc,spacedkc;
	
	ReadClassifier(Path r1Path,Path r2Path,Path kcPath,Path longkcPath,Path spacedkcPath,int nthreads,
			       Path[] tempDataList,Path[] tempIndexList,int r1BarStart,int r2BarStart,
	               int r1BarLen,int r2BarLen,int r1InsStart,int r2InsStart){
		this.r1Path=r1Path;
		this.r2Path=r2Path;
		this.kcPath=kcPath;
		this.longkcPath=longkcPath;
		this.spacedkcPath=spacedkcPath;
		this.nthreads=nthreads;
		this.tempDataList=tempDataList;
		this.tempIndexList=tempIndexList;
		this.r1BarStart=r1BarStart;
		this.r2BarStart=r2BarStart;
		this.r1BarLen=r1BarLen;
		this.r2BarLen=r2BarLen;
		this.r1InsStart=r1InsStart;
		this.r2InsStart=r2InsStart;
	}
	
	private FastqReader loadR1() throws FileNotFoundException, IOException {
		return new FastqReader(this.r1Path.toFile());
	}
	
	private FastqReader loadR2() throws FileNotFoundException, IOException {
		return new FastqReader(this.r2Path.toFile());
	}
	
	private KmerCollection loadDatabase(Path kcPath) throws IllegalBioFileException, IOException {
		FileInputStream kcIn=null;
		try {
			kcIn=new FileInputStream(kcPath.toFile());
			FileChannel kcFc=kcIn.getChannel();
			return KmerCollection.loadFromDisk(kcFc);
		}catch(FileNotFoundException e){
			throw new IllegalBioFileException("Cannot Find Database File "+kcPath);
		}catch(IllegalBioFileException e){
			throw new IllegalBioFileException("Database File "+kcPath+" is damaged !");
		}finally {
			if(kcIn!=null) {
				try {
					kcIn.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	private LongKmerCollection loadLongDatabase(Path longkcPath) throws IllegalBioFileException, IOException {
		FileInputStream kcIn=null;
		try {
			kcIn=new FileInputStream(longkcPath.toFile());
			FileChannel kcFc=kcIn.getChannel();
			return LongKmerCollection.loadFromDisk(kcFc);
		}catch(FileNotFoundException e){
			throw new IllegalBioFileException("Cannot Find Database File "+longkcPath);
		}catch(IllegalBioFileException e){
			throw new IllegalBioFileException("Database File "+longkcPath+" is damaged !");
		}finally {
			if(kcIn!=null) {
				try {
					kcIn.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	
	public void temp(KmerCollection kc) throws IllegalPositionException, IllegalBioFileException, IOException {
		FastaReader fa=new FastaReader(Paths.get("/home/lhc/fusion_test/ref/hg38HBV.fa"),Paths.get("/home/lhc/fusion_test/ref/hg38HBV.fasta.fai"));
		String[] chArray= {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"};
		DNASeqWindow refWindowL=new DNASeqWindow(23);
		//DNASeqWindow refWindowR=new DNASeqWindow(23);
		for(String chr:chArray) {
			RefSequence ref=fa.getSequence(chr);

			refWindowL.set(ref,0);
			//refWindowR.set(ref,22);
			
			for(int pos=0;refWindowL.hasNextKmer();++pos) {
				long kmerL=refWindowL.nextKmer();
				//long kmerR=refWindowR.nextKmer();
				long res=kc.find(kmerL);
				//if(res!=0xFFFFL) {
					int code=(int)(res>>>32);
			    	int po=(int)(res&0xFFFFL);
					System.out.println(kmerL+":"+chr+":"+pos+":"+code+":"+po);
				//}
			}
		}
		System.out.println("XXXXXXX");
		System.exit(0);
	}
	
	public void temp2(KmerCollection kc) throws IllegalPositionException, IllegalBioFileException, IOException {
		
		long kmer=18776233498198L;
				long res=kc.find(kmer);
				
					int code=(int)(res>>>32);
			    	int po=(int)(res&0xFFFFL);
					System.out.println(DNASequence.decompressKmer(kmer,23)+":"+code+":"+po);

		System.out.println("XXXXXXX");
		System.exit(0);
	}
	
	public void loadRef() throws IllegalBioFileException, IOException, IllegalPositionException {
		this.kc=loadDatabase(this.kcPath);
		this.longkc=loadLongDatabase(this.longkcPath);
		this.spacedkc=loadLongDatabase(this.spacedkcPath);
		//kc.info();
		//temp2(kc);
	}
	
	public int[] classify() throws FileNotFoundException, IOException, IllegalBioFileException, InterruptedException, IllegalPositionException {
		final KmerCollection kc=this.kc;
		final LongKmerCollection longkc=this.longkc,spacedkc=this.spacedkc;
		
		FastqReader r1=loadR1();
		FastqReader r2=loadR2();
		
		ClassifyManager classifier=new ClassifyManager(kc,longkc,spacedkc,r1,r2,this.nthreads,this.tempDataList,this.tempIndexList,
				                                      this.r1BarStart,this.r2BarStart,this.r1BarLen,this.r2BarLen,this.r1InsStart,this.r2InsStart);
	
		int[] blockcntList=classifier.classifyReads();

		return blockcntList;
	}
	
	public void clean() {
		this.kc.clean();
		this.longkc.clean();
		this.spacedkc.clean();
		this.kc=null;
		this.longkc=null;
		this.spacedkc=null;
	}
	
}
