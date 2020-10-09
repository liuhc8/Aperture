import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

public final class ClassifyWorker implements Runnable{
	//private static int workercnt=0;
	
	/*  Public Resources  */
	private final int no,k,longk,spacedk,geneCodeLen,r1BarStart,r2BarStart,r1BarLen,r2BarLen,r1InsStart,r2InsStart;
	//private final DoubleBufferList dbfList;
	private final DoubleBlockingQueue dbQueue;
	private final KmerCollection kc;
	private final LongKmerCollection longkc,spacedkc;
	//private final String workDir,projectName;
	private final int[] blockcntList;
	//private final Lock lock;
	//private final Condition lastWorkerFinish,writeFinish,readFinish;
	
	/*  Private Resources */
	private byte[] r1Seq,r2Seq,r1Qua,r2Qua;
	private int r1InsLen,r2InsLen,r1r2KmerBorderIdx;
	//private long barcode;
	private int r1Bar,r2Bar;
	private long[] r1Kmers,r2Kmers,r1Rev,r2Rev,bpSeq;
	private int[] forCodeList,revCodeList,forPosList,revPosList,forScoreList,revScoreList;
	private int[] forCodeList2,revCodeList2,forPosList2,revPosList2,forScoreList2,revScoreList2;
	private BpSaver bpSaver;
	private DNASeqWindow seqWindow;
	private ReadKmerCollection rkcTrimer;
	private HashMap<Integer,Integer> posCountMap;
	
	//########### TEMP #########
	Read r1,r2;
	
	public ClassifyWorker(int no,int k,int longk,int spacedk,int geneCodeLen,Path tempDataPath,Path tempIndexPath,
			              int r1BarStart,int r2BarStart,int r1BarLen,int r2BarLen,int r1InsStart,int r2InsStart,DoubleBlockingQueue dbQueue,
			              KmerCollection kc,LongKmerCollection longkc,LongKmerCollection spacedkc,int[] blockcntList) throws IOException{
		this.geneCodeLen=geneCodeLen;
		this.k=k;
		this.longk=longk;
		this.spacedk=spacedk;
		this.no=no;
		this.r1BarStart=r1BarStart;
		this.r2BarStart=r2BarStart;
		this.r1BarLen=r1BarLen;
		this.r2BarLen=r2BarLen;
		this.r1InsStart=r1InsStart;
		this.r2InsStart=r2InsStart;
		this.dbQueue=dbQueue;
		this.kc=kc;
		this.longkc=longkc;
		this.spacedkc=spacedkc;
		this.blockcntList=blockcntList;

		this.bpSeq=new long[15];
		this.r1Kmers=new long[300];
		this.r2Kmers=new long[150];
		this.r1Rev=new long[300];
		this.r2Rev=new long[150];
		
		this.forCodeList=new int[300];
		this.revCodeList=new int[300];
		this.forPosList=new int[300];
		this.revPosList=new int[300];
		this.forScoreList=new int[300];
		this.revScoreList=new int[300];
		
		this.forCodeList2=new int[150];
		this.revCodeList2=new int[150];
		this.forPosList2=new int[150];
		this.revPosList2=new int[150];
		this.forScoreList2=new int[150];
		this.revScoreList2=new int[150];
		
		this.bpSaver=new BpSaver(tempDataPath,tempIndexPath);
		this.seqWindow=new DNASeqWindow(k);
		this.rkcTrimer=new ReadKmerCollection(150);
		
		this.posCountMap=new HashMap<Integer,Integer>();
	}

	@Override
	public void run() {
		bpDetect();
	}

	
	private void bpDetect() {
		final int no=this.no;
		final int[] blockcntList=this.blockcntList;
		final BpSaver bpSaver=this.bpSaver;
		final DoubleBlockingQueue dbQueue=this.dbQueue;
		
		try {
			
			ReadBox lastBox=new ReadBox(ClassifyManager.Read_BOX_CAP);
			
			while(!Thread.interrupted()) {
				
				ReadBox currentBox=dbQueue.tradeIn(lastBox);	
				if(currentBox.isNull()) {
					blockcntList[no]=bpSaver.finalizeWrite();
					return;
				}else {
					analyseReads(currentBox);
				}
				lastBox=currentBox;
			}
		}catch(InterruptedException | IllegalBioFileException | IOException e) {
			System.out.println("Breakpoint Detection Worker Error!");
			e.printStackTrace();
		}finally {
			bpSaver.closeAndClean();
			clean();
		}
	}
	

	private void analyseReads(ReadBox readBox) throws IllegalBioFileException, IOException{
		final int r1InsStart=this.r1InsStart;
		final int r2InsStart=this.r2InsStart;
		//final int remainLimit=this.reaminLimit;
		final int k=this.k;
		//final BreakpointSaver bpSaver=this.bpSaver;
		
		//final Read[] readList=dbfList.getReadList();
		//final int capacity=dbfList.getCapacity();
		
		//boolean requestSave=false;

	//	checkFastq(readList[no],readList[no+capacity]);       //check once in each cycle
		boolean fqChecked=false;
		Iterator<Read> r1Itr=readBox.getR1Itr();
		Iterator<Read> r2Itr=readBox.getR2Itr();
    	
		while(r1Itr.hasNext()) {
	    	Read r1=r1Itr.next();
	    	Read r2=r2Itr.next();
	   
	    	if(r1.isNull() || r2.isNull()) {
	    		return;
	    	}
	    	
	    	if(!fqChecked) {
	    		checkFastq(r1,r2);
	    		fqChecked=true;
	    	}
	    	
	    	this.r1=r1;
	    	this.r2=r2;
	    	
	    	this.r1Seq=r1.getSequence();
	    	this.r2Seq=r2.getSequence();
	    	this.r1Qua=r1.getQua();
	    	this.r2Qua=r2.getQua();
	    	
	    	
			this.r1InsLen=r1.getLength()-r1InsStart;
			this.r2InsLen=r2.getLength()-r2InsStart;
			
			this.r1r2KmerBorderIdx=this.r1InsLen-k+1;
			
			this.r1Bar=0;
			this.r2Bar=0;
			
			
			if(filterLowQua()) {
				continue;
			}
			
			
			cutIntoKmers();
			
		
			//System.out.println(r1.toString(this.r1len+this.r1start));
			//System.out.println(r2.toString(this.r2len+this.r2start));
			
			
			overlapAnalysis();
			
			
			//System.out.println(r1.toString(this.r1len+this.r1start));
			//System.out.println(r2.toString(this.r2len+this.r2start));
			
			
			if(this.r1InsLen-k+1<=50) {
				continue;
			}
			
			searchKmers();
			
			//String a="....";
		//	for(int m=0;m<r1KmersLen+r2KmersLen;++m) {
		//		System.out.print(r1Kmers[i]+"  ::  ");
				
				//System.out.println(a+DNASequence.decompressKmer(r1Kmers[i],k));
				//a=a+".";
		//		System.out.println(r1Rev[i]);
		//		System.out.print(kmersCode[m]+"  ::  ");
		//		System.out.println(revCode[m]);
				
		//		if(m==r1KmersLen-1) {
		//			System.out.println("####################################");
		//		}
		//	}
			
			breakpointAnalysis();
			
	    }

	}
	
	
	private void checkFastq(Read r1,Read r2) throws IllegalBioFileException {
    	final byte[] r1Name=r1.getName();
    	final byte[] r2Name=r2.getName();
    	final int r1Len=r1Name.length;
    	for(int i=0;i<r1Len;++i) {
			if(r1Name[i]==32 || r1Name[i]==47) {
				break;
			}
			if(r1Name[i]!=r2Name[i]) {
				throw new IllegalBioFileException("R1 and R2 not matched!");
			}
		}
	}
	
	/*private int writeRes(FileChannel seqfc,FileChannel idxfc) throws IOException {
		return this.bpSaver.writeRes(seqfc, idxfc);
	}*/
	
	private boolean filterLowQua() {
		final int r1start=this.r1InsStart,r1len=this.r1InsLen,r1BarStart=this.r1BarStart,r1InsStart=this.r1InsStart,r1BarLen=this.r1BarLen;
		final int r2start=this.r2InsStart,r2len=this.r2InsLen,r2BarStart=this.r2BarStart,r2InsStart=this.r2InsStart,r2BarLen=this.r2BarLen;
		final byte[] r1qua=this.r1Qua,r1seq=this.r1Seq;
		final byte[] r2qua=this.r2Qua,r2seq=this.r2Seq;
		
		try {
			if(r1BarLen<=0) {
				this.r1Bar=DNASequence.compressBarcode(r1seq, r1InsStart+10, 8);
			}else {
				this.r1Bar=DNASequence.compressBarcode(r1seq, r1BarStart, r1BarLen)^DNASequence.compressBarcode(r1seq, r1InsStart+20, r1BarLen);
			}
			if(r2BarLen<=0) {
				this.r2Bar=DNASequence.compressBarcode(r2seq, r2InsStart+10, 8);
			}else {
				this.r2Bar=DNASequence.compressBarcode(r2seq, r2BarStart, r2BarLen)^DNASequence.compressBarcode(r2seq, r2InsStart+20, r2BarLen);
			}
		}catch(IllegalBarcodeException e) {
		    return true;
		}
		
		int cntN1=0;
		int cntL1=0;
		for(int i=0;i<r1len;++i) {
			if(r1seq[i+r1start]=='N') {
				//r1seq[i+r1start]=DNASequence.getRandomBase();
				++cntN1;
			}
			if(r1qua[i+r1start]<=25 + 33) {
				++cntL1;
			}
		}
		
		if(cntL1/r1len>=0.4 || cntN1>0) {
			return true;
		}
		
		int cntN2=0;
		int cntL2=0;
		for(int i=0;i<r2len;++i) {
			if(r2seq[i+r2start]=='N') {
				//r2seq[i+r2start]=DNASequence.getRandomBase();
				++cntN2;
			}
			if(r2qua[i+r2start]<=25 + 33) {
				++cntL2;
			}
		}
		
		if(cntL2/r2len>=0.4 || cntN2>0) {
			return true;
		}
		
		return false;
	}

	private void cutIntoKmers() {
		final DNASeqWindow seqWindow=this.seqWindow;
		final int k=this.k;
		final int r1KmersLen=this.r1InsLen>=k?this.r1InsLen-k+1:0;
		final int r2KmersLen=this.r2InsLen>=k?this.r2InsLen-k+1:0;
		
		final int r1InsStart=this.r1InsStart;
		final int r1InsLen=this.r1InsLen;
		final byte[] r1Seq=this.r1Seq;
		
		final int r2InsStart=this.r2InsStart;
		final int r2InsLen=this.r2InsLen;
		final byte[] r2Seq=this.r2Seq;
		
		final long[] r1Kmers=this.r1Kmers;
		final long[] r2Kmers=this.r2Kmers;
		final long[] r1Rev=this.r1Rev;
		final long[] r2Rev=this.r2Rev;

		
		if(r1KmersLen>0) {
			seqWindow.set(r1Seq,r1InsStart,1,r1InsLen+r1InsStart);
			for(int i=0;seqWindow.hasNextKmer();++i) {
				r1Kmers[i]=seqWindow.nextKmer();
				r1Rev[i]=seqWindow.showThisAsRevComp();
			}
		}
		
		if(r2KmersLen>0) {
			seqWindow.set(r2Seq,r2InsStart,1,r2InsLen+r2InsStart);
			for(int i=0;seqWindow.hasNextKmer();++i) {
				r2Kmers[i]=seqWindow.nextKmer();
				r2Rev[i]=seqWindow.showThisAsRevComp();
			}
		}
		
		
	}
	
	private void overlapAnalysis() {
		final int k=this.k;            
		
		final ReadKmerCollection rkcTrimer=this.rkcTrimer;


		final int r1InsLen=this.r1InsLen;
		final int r2InsStart=this.r2InsStart;
		final int r2InsLen=this.r2InsLen;
		
		final byte[] r2Seq=this.r2Seq;
		
		final long[] r1Kmers=this.r1Kmers;
		final long[] r1Rev=this.r1Rev;
		final long[] r2Rev=this.r2Rev;
		final int r1KmersLen=r1InsLen>=k?r1InsLen-k+1:0;
		final int r2KmersLen=r2InsLen>=k?r2InsLen-k+1:0;
		final byte[] nt4BinRev=DNASequence.NT4_UPPER_REV_BIN_TABLE;

		
		for(short i=0;i<r1KmersLen;++i) {
			rkcTrimer.insert(r1Kmers[i],i);
		}

		
		try {
			
	     	int pos=rkcTrimer.find(r2Rev[0]);     //first k-mer in reads
	    	if(pos>=0) {                             //Cut Adapter
		    	this.r1InsLen=pos+k;
		    	this.r2InsLen=0;
		    	this.r1r2KmerBorderIdx=this.r1InsLen-k+1;
		    	//System.out.println("AAAAA");
		    	//System.out.println(pos);
		    	return;
	    	}

		
	    	for(int i=r2KmersLen-1;i>=0;i-=5) {
	    		pos=rkcTrimer.find(r2Rev[i]);
		    	if(pos>=0) {
		    		int insertSize=pos+i+k;
			    	if(insertSize>r1InsLen) {
			    		//System.out.println(DNASequence.decompressKmer(r1Kmers[pos],23));
			    		//System.out.println(DNASequence.decompressKmer(r2Rev[i],23));
			    		//System.out.println(DNASequence.decompressKmer(r1Kmers[r1KmersLen-1],23));
				    	//System.out.println(DNASequence.decompressKmer(r2Rev[pos+i-r1KmersLen+1],23));
			    		if(!DNASequence.isSimilar(r1Kmers[r1KmersLen-1], r2Rev[pos+i-r1KmersLen+1])) {
			    		//if(r1Kmers[r1KmersLen-1] != r2Rev[pos+i-r1KmersLen+1]) {
			    			//System.out.println("DDDDD");
					    	//System.out.println(pos);
			    			return;
			    		}
			    		this.r2InsLen=insertSize-r1InsLen+k-1;
				    	copyKmersR2ToR1();
				    	//System.out.println("BBBB");
				    	//System.out.println(i+":::"+pos);
				    	
				    	return;
			    	}else {
			     		this.r1InsLen=insertSize;
			    		this.r2InsLen=0;
			    		this.r1r2KmerBorderIdx=this.r1InsLen-k+1;
			    		//System.out.println("CCCC");
				    	//System.out.println(pos);
			    		return;
			    	}
		    	}
	     	}

		
	    	/*long r1LastKmer=r1Kmers[r1KmersLen-1];
	     	long r2LastKmer=r2Rev[r2KmersLen-1];
		
	    	long flag=(1L<<(k*2))-1;
	    	long flag2=(1L<<(k*2))-1;
	    	for(int i=0;i<k-15;++i) {
		    	r2LastKmer>>>=2;
		        flag>>>=2;
	     		if(((r1LastKmer^r2LastKmer)&flag)==0) {
				
		    		for(int j=0;j<i;++j) {
		    			r1LastKmer<<=2;
		    			r1LastKmer|=nt4BinRev[r2Seq[r2InsStart+r2InsLen-k+i-j]];
		    	    	r1LastKmer&=flag2;
		    			r1Kmers[r1KmersLen+j]=r1LastKmer;
			    		r1Rev[r1KmersLen+j]=DNASequence.getRevComp(r1LastKmer,k);
		    		}
		    		this.r1InsLen+=(i);
		    		copyKmersR2ToR1();
		    		//System.out.println("DDDDD");
			    	//System.out.println(i);
	       			return;
	     		}
     		}*/

		}finally {
			rkcTrimer.clean();
		}
	}
	
	private void copyKmersR2ToR1() {
		final int k=this.k;
		final int r1InsLen=this.r1InsLen;
		final int r2InsLen=this.r2InsLen;
		final long[] r1Kmers=this.r1Kmers;
		final long[] r2Kmers=this.r2Kmers;
		final long[] r1Rev=this.r1Rev;
		final long[] r2Rev=this.r2Rev;
		final int r1KmersLen=r1InsLen>=k?r1InsLen-k+1:0;
		final int r2KmersLen=r2InsLen>=k?r2InsLen-k+1:0;
		
	    //System.out.println(this.r1len+" :: "+this.r2len);
		
		int totalLen=r1KmersLen+r2KmersLen-1;
     	for(int i=0;i<r2KmersLen;++i) {
     		int idx=totalLen-i;
     		r1Rev[idx]=r2Kmers[i];
     		r1Kmers[idx]=r2Rev[i];
     		
     		
    	}
     	this.r1InsLen=r1InsLen+r2InsLen-k+1;
     	this.r2InsLen=0;
	}
	
	private void searchKmers() {
		final int k=this.k,longk=this.longk,spacedk=this.spacedk;
		final KmerCollection kc=this.kc;
		final LongKmerCollection longkc=this.longkc,spacedkc=this.spacedkc;
		final int r1KmersLen=this.r1InsLen>=k?this.r1InsLen-k+1:0;
		final int r2KmersLen=this.r2InsLen>=k?this.r2InsLen-k+1:0;

		//System.out.println("r1KmersLen::"+r1KmersLen);
		// System.out.println("r2KmersLen::"+r2KmersLen);
		
		searchKmers0(kc,r1KmersLen,this.r1Kmers,this.forCodeList,this.forPosList,this.forScoreList);
		searchKmers0(kc,r1KmersLen,this.r1Rev,this.revCodeList,this.revPosList,this.revScoreList);
		
		searchKmers0(kc,r2KmersLen,this.r2Kmers,this.forCodeList2,this.forPosList2,this.forScoreList2);
		searchKmers0(kc,r2KmersLen,this.r2Rev,this.revCodeList2,this.revPosList2,this.revScoreList2);
		
		searchLongKmers0(k,longk,longkc,r1KmersLen,this.r1Kmers,this.r1Rev,this.forCodeList,this.revCodeList,this.forPosList,this.revPosList,this.forScoreList,this.revScoreList);
		searchLongKmers0(k,longk,longkc,r2KmersLen,this.r2Kmers,this.r2Rev,this.forCodeList2,this.revCodeList2,this.forPosList2,this.revPosList2,this.forScoreList2,this.revScoreList2);
		
		searchLongKmers0(k,spacedk,spacedkc,r1KmersLen,this.r1Kmers,this.r1Rev,this.forCodeList,this.revCodeList,this.forPosList,this.revPosList,this.forScoreList,this.revScoreList);
		searchLongKmers0(k,spacedk,spacedkc,r2KmersLen,this.r2Kmers,this.r2Rev,this.forCodeList2,this.revCodeList2,this.forPosList2,this.revPosList2,this.forScoreList2,this.revScoreList2);
		
		copyFindResR2ToR1();
	}
	
	private void searchKmers0(KmerCollection kc,int kmersLen,long[] kmerList,int[] codeList,int[] posList,int[] scoreList) {
		
	//	System.out.println("kmerlen::"+kmersLen);
		for(int i=0;i<kmersLen;i+=2) {	
			long findres=kc.find(kmerList[i]);
			decodeFindRes(i,findres,codeList,posList,scoreList);
		}

		for(int i=1;i<kmersLen;i+=2) {			
			if(i<2 || i>kmersLen-3) {
				long findres=kc.find(kmerList[i]);
				decodeFindRes(i,findres,codeList,posList,scoreList);
			}else {
				if(codeList[i-3]==codeList[i+3]) {
					if(Math.abs(posList[i-3]-posList[i+3])==6) {
						codeList[i]=codeList[i-3];
						posList[i]=(posList[i-3]+posList[i+3])/2;
						scoreList[i]=Math.max(scoreList[i-3],scoreList[i+3]);
					}else if(posList[i-3]==0xFFFF && posList[i+3]==0xFFFF) {
						codeList[i]=codeList[i-3];
						posList[i]=0xFFFF;
						if(codeList[i]==0) {
							scoreList[i]=0;
						}else {
							scoreList[i]=1;
						}
					}else {
						long findres=kc.find(kmerList[i]);
						decodeFindRes(i,findres,codeList,posList,scoreList);	
					}
				}else {
					long findres=kc.find(kmerList[i]);
					decodeFindRes(i,findres,codeList,posList,scoreList);	
				}
			}	
		}
	}
	
	private void decodeFindRes(int idx,long findres,int[] codeList,int[] posList,int[] scoreList) {
		if(findres==0xFFFFL) {
			codeList[idx]=0;
			posList[idx]=0xFFFF;
			scoreList[idx]=0;
		}else {
    		codeList[idx]=(int)(findres>>>32);
	    	int pos=(int)(findres&0xFFFFL);
    		posList[idx]=pos;
	    	if(pos==0xFFFF) {
	    		scoreList[idx]=1;
    		}else {
    			int repeatCategory=(int)((findres&0x30000L)>>>16);
    			if(repeatCategory==KmerCollection.UNIQUE) {
    				scoreList[idx]=5;
    			}else {
    				scoreList[idx]=3;
    			}
    		}
		}
	}
	
	private void searchLongKmers0(int k,int longk,LongKmerCollection longkc,int kmersLen,long[] kmerList,long[] revList,int[] kmersCodeList,int[] revCodeList,int[] kmersPosList,int[] revPosList,int[] kmersScoreList,int[] revScoreList) {
		
		final int gap=longk-k;
		int kmerLMove1=0,kmerLMove2=0;
		
		if(longk<k*2) {
			kmerLMove1=(k-(longk-32))*2;
    		kmerLMove2=(longk-k)*2;
		}else {
			kmerLMove1=(32-k)*2;
    		kmerLMove2=k*2;
		}
		
		
		for(int i=0;i<kmersLen-gap;++i) {
			if(kmersCodeList[i]!=0 && kmersPosList[i]==0xFFFF && kmersCodeList[i+gap]!=0 && kmersPosList[i+gap]==0xFFFF) {					
				long findres=longkc.find((int)(kmerList[i]>>>kmerLMove1),(kmerList[i]<<kmerLMove2)|kmerList[i+gap]);		
					//long findres=longkc.find((kmers[i+gap]^(kmers[i]*3)));
				int code=(int)(findres>>>32);
				int pos=(int)findres;
				if(code!=0) {
				    kmersCodeList[i]=code;
				    if(pos!=0xFFFF) {
				    	kmersPosList[i]=pos;
		    		    kmersScoreList[i]=4;
				    }
				}
			}
		}
		for(int i=gap;i<kmersLen;++i) {
			if(revCodeList[i]!=0 && revPosList[i]==0xFFFF && revCodeList[i-gap]!=0 && revPosList[i-gap]==0xFFFF) {
				long findres=longkc.find((int)(revList[i]>>>kmerLMove1),(revList[i]<<kmerLMove2)|revList[i-gap]);
				//long findres=longkc.find((rev[i-gap]^(rev[i]*3)));
				int code=(int)(findres>>>32);
				int pos=(int)findres;
				if(code!=0) {
					revCodeList[i]=code;
					if(pos!=0xFFFF) {
						revPosList[i]=pos;
		    			revScoreList[i]=4;
				    }
				}
			}
		}
		for(int i=0;i<kmersLen-gap;++i) {
			if(kmersScoreList[i]==4) {
				if(kmersCodeList[i]==kmersCodeList[i+3] && kmersCodeList[i]==kmersCodeList[i+6] && 
						kmersPosList[i+3]*2==kmersPosList[i]+kmersPosList[i+6]) {
					continue;
				}
				if(i>=3) {
			    	if(kmersCodeList[i]==kmersCodeList[i-3] && kmersCodeList[i]==kmersCodeList[i+3] &&
			    			kmersPosList[i]*2==kmersPosList[i-3]+kmersPosList[i+3]) {
			    		continue;
			    	}
				}
				if(i>=6) {
					if(kmersCodeList[i]==kmersCodeList[i-3] && kmersCodeList[i]==kmersCodeList[i-6] &&
							kmersPosList[i-3]*2==kmersPosList[i]+kmersPosList[i-6]) {
						//kmersScoreList[i]+=gap;
						for(int j=3;j<gap;j+=3) {
							if(kmersScoreList[i+j]==1) {
					     		kmersScoreList[i+j]=3;
							}
						}
						continue;
					}
				}
				kmersCodeList[i]=0xFFFFFFFF;
		    	kmersPosList[i]=0xFFFF;
		    	kmersScoreList[i]=1;
			}
		}
		for(int i=gap;i<kmersLen;++i) {
			if(revScoreList[i]==4) {
				if(revCodeList[i]==revCodeList[i-3] && revCodeList[i]==revCodeList[i-6] &&
						revPosList[i-3]*2==revPosList[i]+revPosList[i-6]) {
					continue;
				}
				if(i<kmersLen-3) {
			    	if(revCodeList[i]==revCodeList[i-3] && revCodeList[i]==revCodeList[i+3] &&
			    			revPosList[i]*2==revPosList[i-3]+revPosList[i+3]) {
			    		continue;
			    	}
				}
				if(i<kmersLen-6) {
			    	if(revCodeList[i]==revCodeList[i+3] && revCodeList[i]==revCodeList[i+6] &&
			    			revPosList[i+3]*2==revPosList[i]+revPosList[i+6]) {
			    		//revScoreList[i]+=gap;
			    		for(int j=3;j<gap;j+=3) {
			    			if(revScoreList[i-j]==1) {
			    		    	revScoreList[i-j]=3;
			    			}
						}
			    		continue;
			    	}
				}
				revCodeList[i]=0xFFFFFFFF;
		    	revPosList[i]=0xFFFF;
		    	revScoreList[i]=1;
			}
		}
	}
	
	private void copyFindResR2ToR1() {
		final int k=this.k;
		final int r1KmersLen=this.r1InsLen>=k?this.r1InsLen-k+1:0;
		final int r2KmersLen=this.r2InsLen>=k?this.r2InsLen-k+1:0;
		
		final int[] kmersCode=this.forCodeList;
		final int[] revCode=this.revCodeList;
		final int[] kmersPos=this.forPosList;
		final int[] revPos=this.revPosList;
		final int[] kmersScore=this.forScoreList;
		final int[] revScore=this.revScoreList;
		
		final int[] kmersCode2=this.forCodeList2;
		final int[] revCode2=this.revCodeList2;
		final int[] kmersPos2=this.forPosList2;
		final int[] revPos2=this.revPosList2;
		final int[] kmersScore2=this.forScoreList2;
		final int[] revScore2=this.revScoreList2;
		
		int wholelen=r1KmersLen+r2KmersLen-1;
		for(int i=0;i<r2KmersLen;++i) {
			int idx=wholelen-i;
			
			revCode[idx]=kmersCode2[i];
			revPos[idx]=kmersPos2[i];
			revScore[idx]=kmersScore2[i];
			
			kmersCode[idx]=revCode2[i];
			kmersPos[idx]=revPos2[i];
			kmersScore[idx]=revScore2[i];
		}
	}
  
    
	private void breakpointAnalysis() throws IOException {
		final int k=this.k;
		final int geneCodeLen=this.geneCodeLen;
		final int r1KmersLen=this.r1InsLen>=k?this.r1InsLen-k+1:0;
		final int r2KmersLen=this.r2InsLen>=k?this.r2InsLen-k+1:0;
		final int totalLen=r1KmersLen+r2KmersLen;
		final int r1r2KmerBorderIdx=this.r1r2KmerBorderIdx;
		final int[] forCodeList=this.forCodeList;
		final int[] revCodeList=this.revCodeList;
		final int[] forPosList=this.forPosList;
		final int[] revPosList=this.revPosList;
		final int[] forScoreList=this.forScoreList;
		final int[] revScoreList=this.revScoreList;
		
		int currentForIdx=0,currentRevIdx=0;
		int leftForCode=0xFFFFFFFF,leftRevCode=0xFFFFFFFF;
		int leftForPos=0xFFFF,leftRevPos=0xFFFF;
		int leftForPosIdx=-1,leftRevPosIdx=-1;	
		boolean leftForPosNeedCorrect=false,leftRevPosNeedCorrect=false;
		

		for(int forGapNum=0;currentForIdx<totalLen;++currentForIdx) {
			int currentForCode=forCodeList[currentForIdx];
			int currentForPos=forPosList[currentForIdx];
			
			if(currentForIdx==r1KmersLen) {
				forGapNum=0;
			}
			
	    	if(currentForCode!=0) {
	    		forGapNum=0;
	    		if(GeneCodeGenerator.testGeneCode(leftForCode&currentForCode,geneCodeLen)<0) {
	    			break;
	    		}else {
	    			//++leftShunKmers;
	    			leftForCode&=currentForCode;
	        		if(currentForPos!=0xFFFF) {
	        			if(Math.abs(currentForPos-leftForPos-(currentForIdx-leftForPosIdx))>50 && leftForPos!=0xFFFF && !(leftForPosIdx<r1KmersLen && currentForIdx>=r1KmersLen)) {
	        				/*synchronized(kc) {
	        					System.out.println(forPos+":"+leftForPos+":"+forIdx+":"+leftForPosIdx);
	        				}*/
	        				if(leftForPosIdx>=r1r2KmerBorderIdx || currentForIdx<r1r2KmerBorderIdx) {
	        					leftForPosNeedCorrect=true;
	        				}
	        			}
		        		leftForPos=currentForPos;
		         		leftForPosIdx=currentForIdx;
	        		}
	    		}
	    	}else {
		    	++forGapNum;
		     	if(forGapNum>=3) {
			    	break;
		    	}
	    	}
	    	
		}
		
		for(int revGapNum=0;currentRevIdx<totalLen;++currentRevIdx) {
			int currentRevCode=revCodeList[currentRevIdx];
			int currentRevPos=revPosList[currentRevIdx];
			
			if(currentRevIdx==r1KmersLen) {
				revGapNum=0;
			}
			
	      	if(currentRevCode!=0) {
     			revGapNum=0;
     			if(GeneCodeGenerator.testGeneCode(leftRevCode&currentRevCode,geneCodeLen)<0) {
     				break;
     			}else {
     				//++leftFanKmers;
     				leftRevCode&=currentRevCode;
    				if(currentRevPos!=0xFFFF) {
    					if(Math.abs(leftRevPos-currentRevPos-(currentRevIdx-leftRevPosIdx))>50 && leftRevPos!=0xFFFF && !(leftRevPosIdx<r1KmersLen && currentRevIdx>=r1KmersLen)) {
    						/*synchronized(kc) {
	        					System.out.println(revPos+":"+leftRevPos+":"+revIdx+":"+leftRevPosIdx);
	        				}*/
    						if(leftRevPosIdx>=r1r2KmerBorderIdx || currentRevIdx<r1r2KmerBorderIdx) {
    							leftRevPosNeedCorrect=true;
	        				}
	        			}
	    				leftRevPos=currentRevPos;
	    				leftRevPosIdx=currentRevIdx;
    				}
     			}
     		}else {
	    		++revGapNum;
     			if(revGapNum>=3) {
	    			break;
	    		}
	    	}
		}
		
		if(currentForIdx>=totalLen-3||currentRevIdx>=totalLen-3) {
			return;
		}
		
		int leftForKmers=0,leftRevKmers=0,leftConfi=0;
		int destIdx=currentForIdx>currentRevIdx?currentForIdx:currentRevIdx;
		int score=0;
		for(int i=0;i<destIdx;++i) {
			score=forScoreList[i]-revScoreList[i];
			if(score>0) {
				++leftForKmers;
			}else if(score<0) {
				++leftRevKmers;
			}
			leftConfi+=score;
		}
		
		if(leftConfi>=0) {
			if(leftForPos!=0xFFFF) {
				int sameCodeKmer=0;
				for(int i=0;i<currentForIdx;++i) {
					if(forCodeList[i]==leftForCode && forPosList[i]!=0xFFFF) {
						++sameCodeKmer;
					}
				}
				leftForPos=sameCodeKmer>=2?leftForPos:0xFFFF;
			}
		}else {
			if(leftRevPos!=0xFFFF) {
				int sameCodeKmer=0;
				for(int i=0;i<currentRevIdx;++i) {
					if(revCodeList[i]==leftRevCode && revPosList[i]!=0xFFFF) {
						++sameCodeKmer;
					}
				}
				leftRevPos=sameCodeKmer>=2?leftRevPos:0xFFFF;
			}
		}
		
		if(leftConfi>=0) {
			if(leftForPos==0xFFFF) {
				findRightEnd(currentForIdx,0,0xFFFF,-1,1,0,0);
			}else {
				leftForPos=leftForPosNeedCorrect?posCorrect(true,0,leftForPosIdx):leftForPos;
				if(leftForPos<0) {
					return;
				}
				findRightEnd(currentForIdx,leftForCode,leftForPos,leftForPosIdx,leftConfi>0?leftConfi:1,leftForKmers,0);
			}
		}else {
			if(leftRevPos==0xFFFF) {
				findRightEnd(currentRevIdx,0,0xFFFF,-1,-1,0,0);
			}else {
				leftRevPos=leftRevPosNeedCorrect?posCorrect(false,0,leftRevPosIdx):leftRevPos;
				if(leftRevPos<0) {
					return;
				}
				findRightEnd(currentRevIdx,leftRevCode,leftRevPos,leftRevPosIdx,leftConfi,leftRevKmers,0);
			}	
		}
		
	}
	
	private int posCorrect(boolean isFor,int startIdx,int currentIdx) {
		HashMap<Integer,Integer> posCountMap=this.posCountMap;
		int[] codeList=isFor?this.forCodeList:this.revCodeList;
		int[] posList=isFor?this.forPosList:this.revPosList;
		int sign=isFor?1:-1;
		int code=codeList[currentIdx];
		for(int i=currentIdx;i>=startIdx;i-=3) {
			if(codeList[i]==code && posList[i]!=0xFFFF) {
		    	int pos=posList[i]-(i*sign);
		    	if(posCountMap.containsKey(pos)) {
			    	posCountMap.put(pos, posCountMap.get(pos)+1);
		    	}else {
			     	posCountMap.put(pos, 1);
		    	}
			}
		}
		int maxCount=0,correctPos=0;
		for(Entry<Integer,Integer> entry:posCountMap.entrySet()) {
			int cnt=entry.getValue();
			if(cnt>maxCount) {
				maxCount=cnt;
				correctPos=entry.getKey();
			}
		}
		posCountMap.clear();
		
		/*synchronized(kc) {
			if(correctPos+(currentIdx*sign)>65500 || correctPos+(currentIdx*sign)<0) {
			System.out.println(r1.toString());
	     	System.out.println(r2.toString());

	     	final int r1KmersLen=this.r1InsLen>=k?this.r1InsLen-k+1:0;
			final int r2KmersLen=this.r2InsLen>=k?this.r2InsLen-k+1:0;
		    for(int m=0;m<r1KmersLen+r2KmersLen;++m) {
		    		System.out.print(m+"  ::  ");
		    		System.out.print(forCodeList[m]+"  ::  ");
		    		System.out.print(forPosList[m]+"  ::  ");
		    		System.out.print(forScoreList[m]+"  ::  ");
		    		System.out.print(revCodeList[m]+"  ::  ");
		    		System.out.print(revPosList[m]+"  ::  ");
		    		System.out.print(revScoreList[m]+"  ::  ");
			    	System.out.println(DNASequence.decompressKmer(r1Kmers[m],23)+" :: "+DNASequence.decompressKmer(r1Rev[m],23)+" :: ");
			    	//System.out.println(r1Kmers[m]+" :: "+r1Rev[m]);
			 	
		         if(m==r1KmersLen-1) {
		     			System.out.println("####################################");
		    	}
	    	}
			
     		System.out.println(r1KmersLen+"::"+r2KmersLen+"::"+startIdx+"::"+currentIdx+"::"+isFor+"::"+(correctPos+(currentIdx*sign)));
			}
    	}*/
		
		return correctPos+(currentIdx*sign);
	}
	
	private boolean findRightEnd(int idx,int leftCode,int leftPos,int lefti,int leftConfi,int leftKmers,int bpCnt) throws IOException {
		
		if(++bpCnt>=4) {
			return false;
		}
		
		final int k=this.k;
		final int geneCodeLen=this.geneCodeLen;
		final int r1KmersLen=this.r1InsLen>=k?this.r1InsLen-k+1:0;
		final int r2KmersLen=this.r2InsLen>=k?this.r2InsLen-k+1:0;
		final int totalLen=r1KmersLen+r2KmersLen;
		final int[] forCodeList=this.forCodeList;
		final int[] revCodeList=this.revCodeList;
		final int[] forPosList=this.forPosList;
		final int[] revPosList=this.revPosList;	
		final int[] forScore=this.forScoreList;
		final int[] revScore=this.revScoreList;
		final int r1r2KmerBorderIdx=this.r1r2KmerBorderIdx;
		
		int leftBpIdx=idx;
		int rightCode=0xFFFFFFFF,rightPos=0xFFFF;
		int bpStart=lefti+1,bpEnd=totalLen;
		
		int[] selectedCode=null;
		int[] selectedPos=null;

		boolean getSameRightCode=false;
		if(leftPos!=0xFFFF) {
			selectedCode=leftConfi>0?forCodeList:revCodeList;
			selectedPos=leftConfi>0?forPosList:revPosList;
	    	for(int i=idx+1;i<totalLen;++i) {
		    	if(selectedCode[i]==leftCode && selectedPos[i]!=0xFFFF && Math.abs(Math.abs(leftPos-selectedPos[i])-(i-lefti))<15) {
		    		bpEnd=i;
		    		rightCode=selectedCode[i];
			    	rightPos=selectedPos[i];
			    	getSameRightCode=true;
		    		break;
	     		}
	    	}
		}
		
		if(!getSameRightCode) {		
		    int confi=0;
		    for(int i=idx+1;i<totalLen;++i) {
		    	confi+=(forScore[i]-revScore[i]);
	    		if(Math.abs(confi)>20) {
		    		break;
		    	}
	    	}
	    	if(confi==0) {
	    		confi=leftConfi>=0?1:-1;
	    	}
	    	selectedCode=confi>=0?forCodeList:revCodeList;
	    	selectedPos=confi>=0?forPosList:revPosList;  
		
    		for(int i=idx+1;i<totalLen;++i) {
		     	if(selectedPos[i]!=0xFFFF) {
		    		if(i<r1KmersLen-3 || (i>=r1KmersLen && i<totalLen-3)) {
			    		if(selectedCode[i+3]==selectedCode[i] && selectedPos[i+3]!=0xFFFF && Math.abs(selectedPos[i+3]-selectedPos[i])==3) {
			    			bpEnd=i;
				     		rightCode=selectedCode[i];
				    		rightPos=selectedPos[i];
			    			break;
		    			}
	    			}
	    		}
	    	}
		
		}
		
		
		int rightKmers=0,rightConfi=0,rightBpIdx=bpEnd;
		int nextBpIdx=0,nextLeftCode=rightCode,nextLeftPos=0xFFFF,nextLeftIdx=-1;
		boolean hasNextBp=false,rightPosNeedCorrect=false;
		
		if(rightPos!=0xFFFF) {	
	     	for(int i=bpEnd-1,gapNum=0,tempRightCode=rightCode;i>0;--i) {
		    	if(i==r1KmersLen-1) {
		    		gapNum=0;
	    		}
		    	if(selectedCode[i]!=0) {
		    		gapNum=0;
		     		tempRightCode&=selectedCode[i];
		    	    if(GeneCodeGenerator.testGeneCode(tempRightCode,geneCodeLen)<0) {
			        	rightBpIdx=i;
			        	break;
		    	    }
	    		}else {
		    		++gapNum;
		    		if(gapNum>=3) {
		    			rightBpIdx=i;
		    			break;
		    		}
	    		}
	    	}
		

	    	for(int i=bpEnd,gapNum=0;i<totalLen;++i) {
		    	if(i==r1KmersLen) {
			     	gapNum=0;
		    	}
		    	if(selectedCode[i]!=0) {
			    	gapNum=0;
			    	//++rightKmers;
			    	if(GeneCodeGenerator.testGeneCode(selectedCode[i]&rightCode,geneCodeLen)<0) {
			    		nextBpIdx=i;
			    		hasNextBp=true;
			    		break;
		    		}
			    	if(selectedPos[i]!=0xFFFF) {
			    		if(Math.abs(Math.abs(selectedPos[i]-nextLeftPos)-(i-nextLeftIdx))>50 && nextLeftPos!=0xFFFF && !(nextLeftIdx<r1KmersLen && i>=r1KmersLen)) {
			    			if(rightBpIdx>=r1r2KmerBorderIdx || i<r1r2KmerBorderIdx) {
			    		    	rightPosNeedCorrect=true;
			    			}
	        			}
			    		nextLeftPos=selectedPos[i];
			    		nextLeftIdx=i;
		     		}
		    	}else {
			    	++gapNum;
			    	if(gapNum>=3) {
			    		nextBpIdx=i;
			    		hasNextBp=true;
			    		break;
	     			}
	    		}
			
	    		int score=forScore[i]-revScore[i];
		    	rightConfi+=score;
	    		if(selectedCode==forCodeList) {
		    		if(score>0) {
		    			++rightKmers;
	    			}
	    		}else {
		    		if(score<0) {
		    			++rightKmers;
	    			}
	    		}
	    	}
		
	    	if(selectedCode==forCodeList) {
	    		rightConfi=rightConfi>0?rightConfi:1;
	    	}else {
	    		rightConfi=rightConfi<0?rightConfi:-1;
	    	}
		//rightConfi=(selectedCode==kmersCode)?Math.abs(rightConfi):-Math.abs(rightConfi);
		}else{
			rightCode=0;
			rightConfi=leftConfi>0?1:-1;
			rightKmers=0;
		}
		
		int leftBpLen=leftBpIdx-bpStart;
		int rightBpLen=bpEnd-rightBpIdx;
		
		if(rightPosNeedCorrect) {
			nextLeftPos=posCorrect(rightConfi>=0,bpEnd,nextLeftIdx);
			rightPos=rightConfi>=0?nextLeftPos-(nextLeftIdx-bpEnd):nextLeftPos+(nextLeftIdx-bpEnd);
			if(nextLeftPos<0 || rightPos<0) {
				return false;
			}
			/*synchronized(kc) {
				if(rightPos>65500 || rightPos<0) {
				System.out.println(r1.toString());
		     	System.out.println(r2.toString());

			    for(int m=0;m<r1KmersLen+r2KmersLen;++m) {
			    		System.out.print(m+"  ::  ");
			    		System.out.print(forCodeList[m]+"  ::  ");
			    		System.out.print(forPosList[m]+"  ::  ");
			    		System.out.print(forScoreList[m]+"  ::  ");
			    		System.out.print(revCodeList[m]+"  ::  ");
			    		System.out.print(revPosList[m]+"  ::  ");
			    		System.out.print(revScoreList[m]+"  ::  ");
				    	System.out.println(DNASequence.decompressKmer(r1Kmers[m],23)+" :: "+DNASequence.decompressKmer(r1Rev[m],23)+" :: ");
				    	//System.out.println(r1Kmers[m]+" :: "+r1Rev[m]);
				 	
			         if(m==r1KmersLen-1) {
			     			System.out.println("####################################");
			    	}
		    	}
				
	     		System.out.println(r1KmersLen+"::"+r2KmersLen+"::"+nextLeftIdx+"::"+bpEnd);
				}
	    	}*/
		}
		
		if(hasNextBp) {		
			boolean savApproval=findRightEnd(nextBpIdx,nextLeftCode,nextLeftPos,nextLeftIdx,rightConfi,rightKmers,bpCnt);
			
	       if(!savApproval) {
	    		return false;
	    	}
	    	
		}
		
		if(ApertureMain.debug) {
		synchronized(kc) {
			if((leftCode==ApertureMain.lCode && rightCode==ApertureMain.rCode)||(leftCode==ApertureMain.rCode && rightCode==ApertureMain.lCode)) {
		
		    //if((leftCode==ApertureMain.lCode && rightCode==ApertureMain.rCode && leftPos<ApertureMain.lPos+200 && leftPos>ApertureMain.lPos-200 && rightPos<ApertureMain.rPos+200 && rightPos>ApertureMain.rPos-200)
    		//		||(leftCode==ApertureMain.rCode && rightCode==ApertureMain.lCode && rightPos<ApertureMain.lPos+200 && rightPos>ApertureMain.lPos-200 && leftPos<ApertureMain.rPos+200 && leftPos>ApertureMain.rPos-200)) {
			//  if((leftCode==ApertureMain.lCode && rightCode==ApertureMain.rCode && leftPos==ApertureMain.lPos && rightPos==ApertureMain.rPos)
			//		||(leftCode==ApertureMain.rCode && rightCode==ApertureMain.lCode && rightPos==ApertureMain.lPos && leftPos==ApertureMain.rPos)) {
			//if(true) {
	       		System.out.println(r1.toString());
	     	    System.out.println(r2.toString());

		    	for(int m=0;m<r1KmersLen+r2KmersLen;++m) {
		    		System.out.print(m+"  ::  ");
		    		System.out.print(forCodeList[m]+"  ::  ");
		    		System.out.print(forPosList[m]+"  ::  ");
		    		System.out.print(forScore[m]+"  ::  ");
		    		System.out.print(revCodeList[m]+"  ::  ");
		    		System.out.print(revPosList[m]+"  ::  ");
		    		System.out.print(revScore[m]+"  ::  ");
			    	System.out.println(DNASequence.decompressKmer(r1Kmers[m],23)+" :: "+DNASequence.decompressKmer(r1Rev[m],23)+" :: ");
			    	//System.out.println(r1Kmers[m]+" :: "+r1Rev[m]);
			 	
		         	if(m==r1KmersLen-1) {
		     			System.out.println("####################################");
		    		}
	    		}
			
     			System.out.println(r1KmersLen+"::"+r2KmersLen+"::"+bpStart+"::"+bpEnd+"::"+leftCode+"::"+leftPos+"::"+rightCode+"::"+rightPos+"::"+leftConfi+"::"+leftKmers+"::"+rightConfi+"::"+rightKmers);
    		}
     	} 
		}
		//}
		//System.out.println(leftCode+":"+leftPos+":"+leftConfi+":"+rightCode+":"+rightPos+":"+rightConfi);

		if(leftPos!=0xFFFF && rightPos!=0xFFFF) {
	    	if(bpStart<=r1KmersLen && bpEnd>=r1KmersLen) {	
	    		if(leftKmers>=10 && rightKmers>=10 && !(leftCode==rightCode && Math.abs(rightPos-leftPos)<50)) {
	    	    	if(bpStart>=r1KmersLen-10 && bpEnd<r1KmersLen+10) {
		             	saveBreakpoint(true,false,false,leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,0,0,0,0);
	    	        }else {
	    	        	if(leftBpIdx<r1KmersLen && rightBpIdx<r1KmersLen) {
	    	        		if(bpStart<r1KmersLen-10) {
			             		saveBreakpoint(false,false,true,leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,leftBpLen,0,bpStart,r1KmersLen);
			         		}
	    	        	}else if(leftBpIdx>=r1KmersLen && rightBpIdx>=r1KmersLen) {
	    	        		if(bpEnd>=r1KmersLen+10) {
			         			saveBreakpoint(false,true,false,leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,0,rightBpLen,r1KmersLen,bpEnd);
			        		}
	    	        	}else {
	    	        		if(r1KmersLen-bpStart>bpEnd-r1KmersLen) {
	    	        			if(bpStart<r1KmersLen-10) {
	    		             		saveBreakpoint(false,false,true,leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,leftBpLen,0,bpStart,r1KmersLen);
	    		         		}
	    	        		}else {
	    	        			if(bpEnd>=r1KmersLen+10) {
	    		         			saveBreakpoint(false,true,false,leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,0,rightBpLen,r1KmersLen,bpEnd);
	    		        		}
	    	        		}
	    	        	}
	        	    }
	    		}
	     	}else {
	     		if(bpEnd-bpStart>=15 && !(leftCode==rightCode && bpEnd-bpStart==Math.abs(rightPos-leftPos)-1)) {    
		         	saveBreakpoint(false,false,false,leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,leftBpLen,rightBpLen,bpStart,bpEnd);
	     		}
	    	}
		}else if(leftPos!=0xFFFF && rightPos==0xFFFF && leftKmers>=15) {			
			if(!(bpStart<r1KmersLen && bpEnd==totalLen)) {
				if(bpEnd-bpStart>10) {
			        saveBreakpoint(false,false,false,leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,0,0,leftBpLen,rightBpLen,bpStart,bpEnd);
				}
			}
		}else if(rightPos!=0xFFFF && leftPos==0xFFFF && rightKmers>=15) {
			if(!(bpStart==0 && bpEnd>=r1KmersLen)) {
				if(bpEnd-bpStart>10) {
			        saveBreakpoint(false,false,false,leftCode,leftPos,0,0,rightCode,rightPos,rightConfi,rightKmers,leftBpLen,rightBpLen,bpStart,bpEnd);
				}
			}
		}
		return true;
		
	}
	
	
	private void saveBreakpoint(boolean isPE,boolean leftVague,boolean rightVague,int leftCode,int leftPos,
			                     int leftConfi,int leftKmers,int rightCode,int rightPos,int rightConfi,int rightKmers,
			                     int leftBpLen,int rightBpLen,int bpStart,int bpEnd) throws IOException {
		
		final int k=this.k;
		final int r1KmersLen=this.r1InsLen>=k?this.r1InsLen-k+1:0;
		final int r2KmersLen=this.r2InsLen>=k?this.r2InsLen-k+1:0;
		final int len=r1KmersLen+r2KmersLen;
		final long[] bpSeq=this.bpSeq;
		final BpSaver bpSaver=this.bpSaver;

		boolean flip,needSecondary;
		if(isPE) {
			needSecondary=false;
			flip=Math.abs(rightConfi) > Math.abs(leftConfi);
		}else if(leftVague || leftPos==0xFFFF) {
			flip=true;
			needSecondary=false;
		}else if(rightVague || rightPos==0xFFFF) {
			flip=false;
			needSecondary=false;
		}else {
			needSecondary=true;
			flip=Math.abs(rightConfi) > Math.abs(leftConfi);
		} 

		int seqLen=0,bpScore=0;
		if(!isPE) {
			seqLen=getBreakpointSeq(flip,bpStart,bpEnd,bpStart<r1KmersLen);
			bpScore=getBreakpointScore(bpStart,bpEnd);
		}
			
    	if(flip) {	
	        bpSaver.insert(isPE,leftVague,leftPos!=0xFFFF,
	        			   rightCode,rightPos,-rightConfi,rightKmers,leftCode,leftPos,-leftConfi,leftKmers,
			               bpSeq,bpEnd-bpStart+k-1,seqLen,rightBpLen,leftBpLen,bpScore,r2Bar,r1Bar,flip,false);
	    }else {
	    	bpSaver.insert(isPE,rightVague,rightPos!=0xFFFF,
	    				   leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,
			               bpSeq,bpEnd-bpStart+k-1,seqLen,leftBpLen,rightBpLen,bpScore,r1Bar,r2Bar,flip,false);
	    }
		
    	if(needSecondary) {
    		boolean secondaryFlip=!flip;
    		seqLen=getBreakpointSeq(secondaryFlip,bpStart,bpEnd,bpStart<r1KmersLen);
    		if(secondaryFlip) {
    			 bpSaver.insert(isPE,leftVague,leftPos!=0xFFFF,
	        			   rightCode,rightPos,-rightConfi,rightKmers,leftCode,leftPos,-leftConfi,leftKmers,
			               bpSeq,bpEnd-bpStart+k-1,seqLen,rightBpLen,leftBpLen,bpScore,r2Bar,r1Bar,secondaryFlip,true);
    		}else {
    			bpSaver.insert(isPE,rightVague,rightPos!=0xFFFF,
	    				   leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,
			               bpSeq,bpEnd-bpStart+k-1,seqLen,leftBpLen,rightBpLen,bpScore,r1Bar,r2Bar,secondaryFlip,true);
    		}
    	}
		
	}
	
	private int getBreakpointScore(int bpStart,int bpEnd) {
		final int[] kmersScore=this.forScoreList;
		final int[] revScore=this.revScoreList;
		
		int bpScore=0;
		int bpScoreRev=0;
		for(int i=bpStart;i<bpEnd;++i) {
			if(kmersScore[i]!=0) {
				bpScore+=(6-kmersScore[i]);
			}
			if(revScore[i]!=0) {
				bpScoreRev+=(6-revScore[i]);
			}
		}
		return Math.max(bpScore, bpScoreRev);
	}
	
	private int getBreakpointSeq(boolean flip,int breakStart,int breakEnd,boolean inR1) {
		final int k=this.k;
		final int r1KmersLen=this.r1InsLen>=k?this.r1InsLen-k+1:0;
		final int r2KmersLen=this.r2InsLen>=k?this.r2InsLen-k+1:0;
		final int totalLen=r1KmersLen+r2KmersLen;
		final long[] r1Kmers=this.r1Kmers;
		final long[] r2Kmers=this.r2Kmers;
		final long[] r1Rev=this.r1Rev;
		final long[] r2Rev=this.r2Rev;
		final long[] bpSeq=this.bpSeq;
		
		//final int len=breakEnd-breakStart;
		//int start=turnAround?breakEnd-1:breakStart;
		//System.out.println(r1KmersLen+"::"+r2KmersLen+"::"+start+"::"+len);

		int len=0;
		int leftShift=64-2*k,nextMove=32-k;
		if(!flip) {
			if(inR1) {
				int start=breakStart-k,end=r1KmersLen,idx=0;
				len=end-start-1;
				//len=r1KmersLen-breakStart+k-1;
				while(end-start>32) {
					long seq=r1Kmers[start+k]<<leftShift;
					seq |= r1Kmers[start+32];
					bpSeq[idx++]=seq;
					start+=32;
				}
				if(end-start>k) {
					long seq=r1Kmers[start+k]<<leftShift;
					int remain=end-start-k-1;
					seq |= (r1Kmers[end-1]<<((nextMove-remain)*2));
					bpSeq[idx++]=seq;
				}else{
					int remain=end-start-1;
					long seq=r1Kmers[end-1]<<((32-remain)*2);
					bpSeq[idx++]=seq;
				}
			}else {
				int start=totalLen-1-breakStart+k,end=-1,idx=0;
				len=start-end-1;
				//len=totalLen-breakStart+k-1;
				while(start-end>32) {
					long seq=r2Rev[start-k]<<leftShift;
					seq |= r2Rev[start-32];
					bpSeq[idx++]=seq;
					start-=32;
				}
				if(start-end>k) {
					long seq=r2Rev[start-k]<<leftShift;
					int remain=start-end-k-1;
					seq |= (r2Rev[end+1]<<((nextMove-remain)*2));
					bpSeq[idx++]=seq;
				}else{
					int remain=start-end-1;
					long seq=r2Rev[end+1]<<((32-remain)*2);
					bpSeq[idx++]=seq;
				}
			}
		}else {
			if(inR1) {
				int start=breakEnd-1+k,end=-1,idx=0;
				len=start-end-1;
				//len=breakEnd+k-1;
				while(start-end>32) {
					long seq=r1Rev[start-k]<<leftShift;
					seq |= r1Rev[start-32];
					bpSeq[idx++]=seq;
					start-=32;
				}
				if(start-end>k) {
					long seq=r1Rev[start-k]<<leftShift;
					int remain=start-end-k-1;
					seq |= (r1Rev[end+1]<<((nextMove-remain)*2));
					bpSeq[idx++]=seq;
				}else{
					int remain=start-end-1;
					long seq=r1Rev[end+1]<<((32-remain)*2);
					bpSeq[idx++]=seq;
				}
			}else {
				int start=totalLen-breakEnd-k,end=r2KmersLen,idx=0;
				len=end-start-1;
				//len=breakEnd-r1KmersLen+k-1;
				while(end-start>32) {
					long seq=r2Kmers[start+k]<<leftShift;
					seq |= r2Kmers[start+32];
					bpSeq[idx++]=seq;
					start+=32;
				}
				if(end-start>k) {
					long seq=r2Kmers[start+k]<<leftShift;
					int remain=end-start-k-1;
					seq |= (r2Kmers[end-1]<<((nextMove-remain)*2));
					bpSeq[idx++]=seq;
				}else{
					int remain=end-start-1;
					long seq=r2Kmers[end-1]<<((32-remain)*2);
					bpSeq[idx++]=seq;
				}
			}
		}
		
		return len;
	}

	private void clean() {
		this.bpSeq=null;
		this.r1Kmers=null;
		this.r2Kmers=null;
		this.r1Rev=null;
		this.r2Rev=null;
		
		this.forCodeList=null;
		this.revCodeList=null;
		this.forPosList=null;
		this.revPosList=null;
		this.forScoreList=null;
		this.revScoreList=null;
		
		this.forCodeList2=null;
		this.revCodeList2=null;
		this.forPosList2=null;
		this.revPosList2=null;
		this.forScoreList2=null;
		this.revScoreList2=null;
		
		this.bpSaver=null;
		this.seqWindow=null;
		this.rkcTrimer=null;
		
		this.posCountMap=null;
	}
	
}
