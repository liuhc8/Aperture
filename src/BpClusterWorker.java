import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.BufferUnderflowException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicInteger;

public class BpClusterWorker implements Runnable{

	private final BpIndexList.SynchronizedGetItr synGetItr;
	private final int mergeK;
	private final Path mergedDataPath;
	private final double similarity;
//	private boolean hasBar;
	private final TranslationTable trTable;
	
	private int code;
	private RefCountList forRefCntList,revRefCntList;
	private boolean[] availableList;
	private ArrayList<SplitSupport> bpListA,bpListB,bpListC;
	private ArrayList<PESupport> peList,peListAfterMerge;
	private ArrayList<BpSeq> seqListA,seqListC,seqListE,seqListB,seqListD;
	
	private SVCollection svCollection;
	
	public BpClusterWorker(BpIndexList.SynchronizedGetItr synGetItr,Path mergedDataPath,int k,int mergeK,
			                     double similarity,SVCollection svCollection,TranslationTable trTable) {
		this.synGetItr=synGetItr;
		this.mergeK=mergeK;
		this.mergedDataPath=mergedDataPath;
		this.similarity=similarity;
		this.svCollection=svCollection;
		this.trTable=trTable;
		
		this.forRefCntList=new RefCountList(mergeK,k,true);
		this.revRefCntList=new RefCountList(mergeK,k,false);
		
		this.availableList=new boolean[100];
		this.bpListA=new ArrayList<SplitSupport>(100);
		this.bpListB=new ArrayList<SplitSupport>(100);
		this.bpListC=new ArrayList<SplitSupport>(100);
		this.peList=new ArrayList<PESupport>(100);
		this.peListAfterMerge=new ArrayList<PESupport>(100);
		this.seqListA=new ArrayList<BpSeq>(100);
		this.seqListC=new ArrayList<BpSeq>(100);
		this.seqListE=new ArrayList<BpSeq>(100);
		this.seqListB=new ArrayList<BpSeq>(100);
		this.seqListD=new ArrayList<BpSeq>(100);
	}

	private void fillKmerList(Segment segment) {
		final int mergeK=this.mergeK;
		
		final int[] refKmerList=new int[65536];
		final int[] refRevList=new int[65536];
		final long[] seq=segment.ref;
		
		long[] middle=new long[seq.length];
		for(int i=0;i<seq.length-1;++i) {
			middle[i]= (seq[i]<<32) | (seq[i+1]>>>32);
		}
		middle[seq.length-1]=seq[seq.length-1]<<32;
		
		int kmerListLen=segment.end-segment.start-mergeK;
		
		long mask=(1L<<(2*mergeK))-1;
		int start=64-2*mergeK,end=34-2*mergeK;
		outer:
		for(int i=0,j=0;;++j) {
			for(int p=start;p>=end;p-=2) {
				refKmerList[i++]=(int)((seq[j]>>>p)&mask);
				if(i>=kmerListLen) {
					break outer;
				}
			}
			for(int p=start;p>=end;p-=2) {
				refKmerList[i++]=(int)((middle[j]>>>p)&mask);
				if(i>=kmerListLen) {
					break outer;
				}
			}
		}
		
		for(int i=0;i<kmerListLen;++i) {
			refRevList[i]=DNASequence.getRevComp(refKmerList[i], mergeK);
		}
		
		this.forRefCntList.loadRefKmerList(refKmerList, kmerListLen);
		this.revRefCntList.loadRefKmerList(refRevList, kmerListLen);
	}

	private void readBreakpoints(long fileStart,int fileLen) throws IOException, IllegalBioFileException {
		FileInputStream in=null;
		ByteBuffer bb=ByteBuffer.allocate(fileLen);
		try {
			in=new FileInputStream(this.mergedDataPath.toFile());
			FileChannel fc=in.getChannel();
			fc.read(bb, fileStart);
			bb.flip();
			//mbb=fc.map(FileChannel.MapMode.READ_ONLY, fileStart, fileLen);
			travelBreakpoints(bb,fileStart,fileLen);
			bb=null;
		}finally {
			//BufferCleaner.clean(mbb);
			if(in!=null) {
				try {
					in.close();
				}catch(IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	private void travelBreakpoints(ByteBuffer mbb,long fileStart,int fileLen) throws IllegalBioFileException{
		final ArrayList<PESupport> peList=this.peList;
		final ArrayList<BpSeq> seqListA=this.seqListA;
		final int mergeK=this.mergeK;
	//	final boolean hasBar=this.hasBar;
		
		while(mbb.hasRemaining()) {
			byte flag=mbb.get();
			
			boolean isSecondary=(flag&64)!=0;
			boolean flip=(flag&32)!=0;
			boolean isPE=(flag&8)!=0;
			boolean rightVague=(flag&4)!=0;
			boolean hasRightPos=(flag&2)!=0;
			
			int leftCode,leftPos,rightCode,rightPos,r1Bar,r2Bar;
			short leftConfi,leftKmers,rightConfi,rightKmers;
			int bpLen=-1,seqLen=-1,bpScore=0;
			short leftBpLen=-1,rightBpLen=-1;
			long[] seq=null;
			try {
			
			 leftCode=mbb.getInt();	
			 if(leftCode!=this.code) {
				 throw new IllegalBioFileException("Wrong code!!");
			 }
			 leftPos=mbb.getShort()&0xFFFF;
			 leftConfi=mbb.getShort();
			 leftKmers=isSecondary?0:mbb.getShort();
			
			
			 rightCode=mbb.getInt(); 
			 rightPos=hasRightPos?(mbb.getShort()&0xFFFF):-1;
			 rightConfi=mbb.getShort();
			 rightKmers=isSecondary?0:mbb.getShort();
			
			
			int checkpoint1=mbb.getInt();
			if(checkpoint1!=777) {
				throw new IllegalBioFileException("Breakpoint file is damaged! Start:"+fileStart+"  Len:"+fileLen+"  "+"isSecondary:"+isSecondary+" flip:"+flip+" isPE:"+isPE+" rightVague:"+rightVague+" hasRightPos:"+hasRightPos);
			}
			
			r1Bar=isSecondary?0:mbb.getInt();
			r2Bar=isSecondary?0:mbb.getInt();
			
			bpLen=-1;seqLen=-1;bpScore=0;
			leftBpLen=-1;rightBpLen=-1;
			seq=null;
			if(!isPE) {
				bpLen=mbb.getShort()&0xFFFF;
				seqLen=mbb.getShort()&0xFFFF;

				int longLen=(seqLen-1)/32+1;
				seq=new long[longLen];
				for(int i=0;i<longLen;++i) {
					seq[i]=mbb.getLong();
				}
				
				int checkpoint2=mbb.getInt();
				
				if(checkpoint2!=999) {
					throw new IllegalBioFileException("Breakpoint file is damaged! Start:"+fileStart+"  Len:"+fileLen+"  "+"isSecondary:"+isSecondary+" flip:"+flip+" isPE:"+isPE+" rightVague:"+rightVague+" hasRightPos:"+hasRightPos);
				}
	         	
	     		leftBpLen=(short)mbb.get();
	     		rightBpLen=(short)mbb.get();
	     		
	     		bpScore=mbb.getShort();
			}
			
			}catch(BufferUnderflowException r) {
				throw new IllegalBioFileException("Breakpoint file is damaged! Start:"+fileStart+"  Len:"+fileLen+"  "+"isSecondary:"+isSecondary+" flip:"+flip+" isPE:"+isPE+" rightVague:"+rightVague+" hasRightPos:"+hasRightPos);
			}
			
			if(ApertureMain.debug) {
				synchronized(svCollection) {
					if((leftCode==ApertureMain.lCode && rightCode==ApertureMain.rCode)||(leftCode==ApertureMain.rCode && rightCode==ApertureMain.lCode)){
						System.out.println("Prepare: "+leftPos+" + "+rightPos+" + "+(leftConfi>0)+" + "+(rightConfi>0));
					}
				}
			}
			
			
			if(isPE) {
				if( (!(leftCode==rightCode && Math.abs(leftPos-rightPos)<500)) || leftConfi*rightConfi<0 ) {
			    	peList.add(new PESupport(leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,r1Bar,r2Bar,flip));
				}
			}else{
				//System.out.println(seq.length+":"+bpLen+":"+seqLen);
				SplitSupport bp=new SplitSupport(rightVague,leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,bpScore,bpLen,r1Bar,r2Bar,flip,isSecondary);
				BpSeq bpSeq=new BpSeq(seq,bpLen,seqLen,mergeK,leftBpLen,rightBpLen,bp);
				bp.setBpSeq(bpSeq);
				seqListA.add(bpSeq);
				
				/*if((breakpoint.getLeftCode()==809501704 && breakpoint.getRightCode()==270549520)||(breakpoint.getLeftCode()==270549520 && breakpoint.getRightCode()==809501704)){
					System.out.print("AAAA "+breakpoint.getLeftCode()+" + "+breakpoint.getRightCode());
				}*/
				
			}
			
		}
		
	}	
	
	private void prepareMergeRight() {
		Collections.sort(this.seqListA);     //Sort in descending order
	}
	
	private void mergeEqual() {
		final ArrayList<BpSeq> seqListA=this.seqListA;
		final ArrayList<BpSeq> seqListB=this.seqListB;
		
		int len=seqListA.size();
		
		BpSeq bpSeqJ=null;
		SplitSupport bpJ=null;
		
		bpSeqJ=seqListA.get(0);
		bpJ=bpSeqJ.getSVSupport();
		seqListB.add(bpSeqJ);
		for(int i=1;i<len;++i) {
			BpSeq bpSeqI=seqListA.get(i);
			SplitSupport bpI=bpSeqI.getSVSupport();
			if(bpI.equals(bpJ)){
				bpJ.mergeWith(bpI);
		  	}else {
				seqListB.add(bpSeqI);
				bpSeqJ=bpSeqI;
				bpJ=bpI;
			}
		}
		
		seqListA.clear();
		
		if(ApertureMain.debug) {
			synchronized(svCollection) {
				for(BpSeq seq:seqListB) {
					SplitSupport bp=seq.getSVSupport();
					if((bp.getLeftCode()==ApertureMain.lCode && bp.getRightCode()==ApertureMain.rCode)||(bp.getLeftCode()==ApertureMain.rCode && bp.getRightCode()==ApertureMain.lCode)){
						System.out.println("After merging equal: "+bp.getLeftPos()+" + "+bp.getRightPos()+" + "+bp.getLeftCnt()+" + "+bp.getRightCnt()+" + "+seq.toStringBp()+" + "+bp.isLeftForward()+" + "+bp.isRightForward());
					}
				}
			}
		}
	}
	
	private void prepareMergePEEqual() {
		Collections.sort(this.peList);     //Sort in descending order
	}
	
	private void mergePEEqual() {
		final ArrayList<PESupport> peList=this.peList;
		final ArrayList<PESupport> peListAfterMerge=this.peListAfterMerge;
		
		int len=peList.size();
		
		PESupport peJ=null;
		
		peJ=peList.get(0);
		peListAfterMerge.add(peJ);
		for(int i=1;i<len;++i) {
			PESupport peI=peList.get(i);
			if(peI.equals(peJ)){
				peJ.mergeWith(peI);
		  	}else {
		  		peListAfterMerge.add(peI);
				peJ=peI;
			}
		}
		
		peList.clear();

	}
	
	private void mergeRightEnd() {
		final ArrayList<BpSeq> seqListB=this.seqListB;
		final ArrayList<BpSeq> seqListC=this.seqListC;
		final ArrayList<PESupport> peListAfterMerge=this.peListAfterMerge;

		Iterator<BpSeq> iter=seqListB.iterator();
		while(iter.hasNext()) {
			BpSeq bpSeq0=iter.next();
			SplitSupport bp0=bpSeq0.getSVSupport();
	    	if(bp0.hasRightPos() && !bp0.isRightVague()) {
	         	seqListC.add(bpSeq0);
	         	break;
	    	}
		}
		
		while(iter.hasNext()) {
			BpSeq bpSeqI=iter.next();
			SplitSupport bpI=bpSeqI.getSVSupport();
			
			if(!bpI.hasRightPos()) {
				SplitSupport bestMergeBp=null;
				int maxSplit=-1;
				for(int j=seqListC.size()-1;j>=0;--j) {
					BpSeq bpSeqJ=seqListC.get(j);
					if(bpSeqJ.isNullSeq()) {
						continue;
					}
					SplitSupport bpJ=bpSeqJ.getSVSupport();
					if(bpSeqI.isPrefix(bpSeqJ) && bpI.getLeftPos()==bpJ.getLeftPos()){
						if(bpJ.getSplitSupport()>maxSplit) {
							maxSplit=bpJ.getSplitSupport();
							bestMergeBp=bpJ;
						}
					}else {
						break;
					}
				}
				if(bestMergeBp!=null) {
					bestMergeBp.mergeWith(bpI);					
				}
			}else if(bpI.isRightVague()) {
				SplitSupport bestMergeBp=null;
				int maxSplit=-1;
				for(int j=seqListC.size()-1;j>=0;--j) {
					BpSeq bpSeqJ=seqListC.get(j);
					if(bpSeqJ.isNullSeq()) {
						continue;
					}
					SplitSupport bpJ=bpSeqJ.getSVSupport();
					if(bpSeqI.isPrefix(bpSeqJ) && bpI.getLeftPos()==bpJ.getLeftPos()){
						if(bpI.getRightCode()==bpJ.getRightCode() && Math.abs(bpI.getRightPos()-bpJ.getRightPos())<100) {
				        	if(bpJ.getSplitSupport()>maxSplit) {
					     		maxSplit=bpJ.getSplitSupport();
						    	bestMergeBp=bpJ;
				    		}
						}
					}else {
						break;
					}
				}
				if(bestMergeBp!=null) {
					bestMergeBp.mergeWith(bpI);
				}else {
					if((!(bpI.getLeftCode()==bpI.getRightCode() && Math.abs(bpI.getLeftPos()-bpI.getRightPos())<500))||(bpI.isLeftForward()!=bpI.isRightForward())) {
				    	peListAfterMerge.add(new PESupport(bpI));     //######## convert split to PE
					}
				}
			}else {
				for(int j=seqListC.size()-1;j>=0;--j) {
					BpSeq bpSeqJ=seqListC.get(j);
					if(bpSeqJ.isNullSeq()) {
						continue;
					}
					SplitSupport bpJ=bpSeqJ.getSVSupport();
					if(bpSeqI.isPrefix(bpSeqJ) && bpI.getLeftPos()==bpJ.getLeftPos()){
						if(bpI.getRightCode()==bpJ.getRightCode() && Math.abs(bpI.getRightPos()-bpJ.getRightPos())<50) {
				     		if(bpSeqI.getBpLen()<bpSeqJ.getBpLen()) {
					    		bpI.mergeWith(bpJ);
						    	seqListC.set(j, BpSeq.getNullSeq());
					     	}
						}
					}else {
						break;
					}
				}
				seqListC.add(bpSeqI);
			}
		}
		
		seqListB.clear();
		
		if(ApertureMain.debug) {
			synchronized(svCollection) {
				for(BpSeq seq:seqListC) {
					if(seq.isNullSeq()) {
						continue;
					}
					SplitSupport bp=seq.getSVSupport();
					if((bp.getLeftCode()==ApertureMain.lCode && bp.getRightCode()==ApertureMain.rCode)||(bp.getLeftCode()==ApertureMain.rCode && bp.getRightCode()==ApertureMain.lCode)){
						System.out.println("After merging right: "+bp.getLeftPos()+" + "+bp.getRightPos()+" + "+bp.getLeftCnt()+" + "+bp.getRightCnt()+" + "+seq.toStringBp()+" + "+bp.isLeftForward()+" + "+bp.isRightForward());
				    }
				}
			}
		}
		
	}
	
	
	
	private void prepareMergeLeftEnd() {
		final ArrayList<BpSeq> seqListC=this.seqListC;
		final ArrayList<BpSeq> seqListD=this.seqListD;
		
		for(int i=0,len=seqListC.size();i<len;++i) {
			BpSeq bpSeq=seqListC.get(i);
			if(!bpSeq.isNullSeq()) {
				bpSeq.reverse();
				seqListD.add(bpSeq);
			}
		}
		seqListC.clear();
		Collections.sort(seqListD);     //Sort in descending order
	}
	
	
	private void mergeLeftEnd() {
		final ArrayList<BpSeq> seqListD=this.seqListD;
		final ArrayList<BpSeq> seqListE=this.seqListE;
		final ArrayList<SplitSupport> bpListA=this.bpListA;

		Iterator<BpSeq> iter=seqListD.iterator();
		if(iter.hasNext()) {
			BpSeq bpSeq0=iter.next();
			seqListE.add(bpSeq0);
		}

		while(iter.hasNext()) {
			BpSeq bpSeqI=iter.next();
			SplitSupport bpI=bpSeqI.getSVSupport();
			
			for(int j=seqListE.size()-1;j>=0;--j) {
				BpSeq bpSeqJ=seqListE.get(j);			
				if(bpSeqJ.isNullSeq()) {
					continue;
				}
				SplitSupport bpJ=bpSeqJ.getSVSupport();
				
				if(bpSeqI.isPrefix(bpSeqJ) && Math.abs(bpI.getLeftPos()-bpJ.getLeftPos())<50) {
					if(bpI.getRightCode()==bpJ.getRightCode() && bpI.getLeftCode()==bpJ.getLeftCode()) {
						if(bpSeqI.getBpLen()<bpSeqJ.getBpLen()) {
							bpI.mergeWith(bpJ);
							seqListE.set(j, BpSeq.getNullSeq());
						}
					}
				}else {
					break;
				}
			}
			seqListE.add(bpSeqI);
		}
		
		
		if(ApertureMain.debug) {
			synchronized(svCollection) {
				for(BpSeq seq:seqListE) {
					if(seq.isNullSeq()) {
						continue;
					}
					SplitSupport bp=seq.getSVSupport();
					if((bp.getLeftCode()==ApertureMain.lCode && bp.getRightCode()==ApertureMain.rCode)||(bp.getLeftCode()==ApertureMain.rCode && bp.getRightCode()==ApertureMain.lCode)){
						System.out.println("After merging left: "+bp.getLeftPos()+" + "+bp.getRightPos()+" + "+bp.getLeftCnt()+" + "+bp.getRightCnt()+" + "+seq.toStringBp()+" + "+bp.isLeftForward()+" + "+bp.isRightForward());
				    }
				}
			}
		}
		
		for(int i=0,len=seqListE.size();i<len;++i) {
			BpSeq bpSeq=seqListE.get(i);
			if(!bpSeq.isNullSeq()) {
				bpSeq.reverse();
				bpListA.add(bpSeq.getSVSupport());
			}
		}
		
		seqListD.clear();
		seqListE.clear();

	}
	
	
	private void prepareMergeCluster() {
		final ArrayList<SplitSupport> bpListA=this.bpListA;
		final int len=bpListA.size();
		for(int i=0;i<len;++i) {
			BpSeq seq=bpListA.get(i).getBpSeq();
			seq.prepareKmer();
			//seq.prepareCntList();
		}
		Collections.sort(bpListA);
	}
	
	private void mergeCluster() {
		final ArrayList<SplitSupport> bpListA=this.bpListA;
		final int len=bpListA.size();
		
		int start=0;
		while(start<len) {
			int leftPos=bpListA.get(start).getLeftPos();
			int end=start+1;
			while(end<len) {
				if(leftPos+50>bpListA.get(end).getLeftPos()) {
					++end;
				}else {
					break;
				}
			}
			mergeCluster0(start,end);
			start=end;
		}
		
		bpListA.clear();
		
		if(ApertureMain.debug) {
			synchronized(svCollection) {
				for(SplitSupport bp:bpListB) {
					if((bp.getLeftCode()==ApertureMain.lCode && bp.getRightCode()==ApertureMain.rCode)||(bp.getLeftCode()==ApertureMain.rCode && bp.getRightCode()==ApertureMain.lCode)){
						System.out.println("After merging cluster: "+bp.getLeftPos()+" + "+bp.getRightPos()+" + "+bp.getLeftCnt()+" + "+bp.getRightCnt()+" + "+bp.isLeftForward()+" + "+bp.isRightForward());
				    }
		        }
			}
		}
		
	}
	
	
	private void mergeCluster0(int start,int end) {
		final RefCountList forRefCntList=this.forRefCntList,revRefCntList=this.revRefCntList;
		final ArrayList<SplitSupport> bpListA=this.bpListA,bpListB=this.bpListB;
		final int len=end-start;
		//final int mergeK=this.mergeK;
		final double similarity=this.similarity;
		
		if(len>this.availableList.length) {
			this.availableList=new boolean[len*2];
		}
		final boolean[] availableList=this.availableList;	

			
		for(int i=0;i<len;++i) {
			double dis=0;
			SplitSupport bp=bpListA.get(i+start);
			BpSeq bpSeq=bp.getBpSeq();
			
			//synchronized(this.indexList) {
				
			if(bp.isLeftForward()) {
				forRefCntList.updateRefCntList(bp);
				dis=bpSeq.getDistance(forRefCntList);
			}else {
				revRefCntList.updateRefCntList(bp);
				dis=bpSeq.getDistance(revRefCntList);
			}
			
			if(ApertureMain.debug) {
			synchronized(svCollection) {
				if((bp.getLeftCode()==ApertureMain.lCode && bp.getRightCode()==ApertureMain.rCode )||(bp.getLeftCode()==ApertureMain.rCode && bp.getRightCode()==ApertureMain.lCode)){
					System.out.println(bp.isLeftForward()+":::"+dis);
					if(bp.isLeftForward()) {
						System.out.println(forRefCntList);
					}else {
						System.out.println(revRefCntList);
					}
					System.out.println(bpSeq.toStringAll());
				}	
				
			}
			}

			if(dis<similarity) {
				availableList[i]=false;
			}else {
				availableList[i]=true;
			}
		}
		
		int maxSupIdx=-1;
		int maxSup=-1;
		SplitSupport maxSupBp=null,lastMaxSupBp=null;
		BpSeq lastMaxSupSeq=null;
		for(int i=0;i<len;++i) {
			if(availableList[i]==false) {
				continue;
			}
			SplitSupport bp=bpListA.get(i+start);
			if(bp.getSplitSupport()>maxSup && bp.hasRightPos() && bp.getRightCode()==bp.getLeftCode()) {
				maxSup=bp.getSplitSupport();
				maxSupBp=bp;
				maxSupIdx=i;
			}
		}
		
		if(maxSupBp==null) {
			maxSup=-1;
			for(int i=0;i<len;++i) {
				if(availableList[i]==false) {
					continue;
				}
				SplitSupport bp=bpListA.get(i+start);
				if(bp.getSplitSupport()>maxSup && bp.hasRightPos()) {
					maxSup=bp.getSplitSupport();
					maxSupBp=bp;
					maxSupIdx=i;
				}
			}
		}
		
		while(maxSupBp!=null) {
			availableList[maxSupIdx]=false;
			bpListB.add(maxSupBp); 
			
			lastMaxSupBp=maxSupBp;
			lastMaxSupSeq=lastMaxSupBp.getBpSeq();
			
			maxSup=-1;
			maxSupIdx=-1;
			maxSupBp=null;
			
			for(int i=0;i<len;++i) {
				if(!availableList[i]) {
					continue;
				}
				SplitSupport bp=bpListA.get(i+start);
				BpSeq bpSeq=bp.getBpSeq();
				double dis=lastMaxSupSeq.getDistance(bpSeq); //-Math.abs(seq.getSeqLen()-lastMaxSupSeq.getSeqLen());
				//int minLen=Math.min(seq.getSeqLen(),lastMaxSupSeq.getSeqLen());
				/*synchronized(fbpList) {
					int disssss=dis-Math.abs(seq.getSeqLen()-lastMaxSupSeq.getSeqLen());
		     		System.out.println(disssss+":"+Math.min(seq.getSeqLen(),lastMaxSupSeq.getSeqLen())+":"+(double)disssss/Math.min(seq.getSeqLen(),lastMaxSupSeq.getSeqLen()));
				}*/
				//System.out.println(dis+":"+minLen+":"+dis/minLen);
				if(dis<similarity && (lastMaxSupBp.getRightCode()==bp.getRightCode() && Math.abs(lastMaxSupBp.getRightPos()-bp.getRightPos())<100)) {
					availableList[i]=false;
					//lastMaxSupBp.updateRight(bp);
					lastMaxSupBp.mergeWith(bp);
				}else {
					if(bp.getSplitSupport()>maxSup && bp.hasRightPos()) {
						maxSup=bp.getSplitSupport();
						maxSupBp=bp;
						maxSupIdx=i;
					}
				}
			}
		}
	
	}
	
	
	private void mergeClusterAgain() {
		final ArrayList<SplitSupport> bpListB=this.bpListB;
		final int len=bpListB.size();
		
		int start=0;
		while(start<len) {
			int leftPos=bpListB.get(start).getLeftPos();
			int end=start+1;
			while(end<len) {
				if(leftPos+100>bpListB.get(end).getLeftPos()) {
					++end;
				}else {
					break;
				}
			}
			mergeClusterAgain0(start,end);
			start=end;
		}
		
		bpListB.clear();
		
		if(ApertureMain.debug) {
			synchronized(svCollection) {
				for(SplitSupport bp:bpListC) {
					if((bp.getLeftCode()==ApertureMain.lCode && bp.getRightCode()==ApertureMain.rCode)||(bp.getLeftCode()==ApertureMain.rCode && bp.getRightCode()==ApertureMain.lCode)){
						System.out.println("After merging cluster again: "+bp.getLeftPos()+" + "+bp.getRightPos()+" + "+bp.getLeftCnt()+" + "+bp.getRightCnt()+" + "+bp.isLeftForward()+" + "+bp.isRightForward());
				    }
		        }
			}
		}
		
	}
	
	private void mergeClusterAgain0(int start,int end) {
		final ArrayList<SplitSupport> bpListB=this.bpListB;
		final ArrayList<SplitSupport> bpListC=this.bpListC;
		final int len=end-start;
		final double similarity=this.similarity;
		
		if(len>this.availableList.length) {
			this.availableList=new boolean[len*2];
		}
		final boolean[] availableList=this.availableList;	
		
		for(int i=0;i<len;++i) {
			availableList[i]=true;
		}
		
		int maxSupIdx=-1;
		int maxSup=-1;
		SplitSupport maxSupBp=null,lastMaxSupBp=null;
		BpSeq lastMaxSupSeq=null;
		for(int i=0;i<len;++i) {
			if(availableList[i]==false) {
				continue;
			}
			SplitSupport bp=bpListB.get(i+start);
			if(bp.getSplitSupport()>maxSup && bp.hasRightPos() && bp.getRightCode()==bp.getLeftCode()) {
				maxSup=bp.getSplitSupport();
				maxSupBp=bp;
				maxSupIdx=i;
			}
		}
		
		if(maxSupBp==null) {
			maxSup=-1;
			for(int i=0;i<len;++i) {
				if(availableList[i]==false) {
					continue;
				}
				SplitSupport bp=bpListB.get(i+start);
				if(bp.getSplitSupport()>maxSup && bp.hasRightPos()) {
					maxSup=bp.getSplitSupport();
					maxSupBp=bp;
					maxSupIdx=i;
				}
			}
		}
		
		while(maxSupBp!=null) {
			availableList[maxSupIdx]=false;
			bpListC.add(maxSupBp); 
			
			lastMaxSupBp=maxSupBp;
			lastMaxSupSeq=lastMaxSupBp.getBpSeq();
			
			maxSup=-1;
			maxSupIdx=-1;
			maxSupBp=null;
			
			for(int i=0;i<len;++i) {
				if(!availableList[i]) {
					continue;
				}
				SplitSupport bp=bpListB.get(i+start);
				BpSeq seq=bp.getBpSeq();
				double dis=lastMaxSupSeq.getDistance(seq); //-Math.abs(seq.getSeqLen()-lastMaxSupSeq.getSeqLen());
				//int minLen=Math.min(seq.getSeqLen(),lastMaxSupSeq.getSeqLen());
				/*synchronized(fbpList) {
					int disssss=dis-Math.abs(seq.getSeqLen()-lastMaxSupSeq.getSeqLen());
		     		System.out.println(disssss+":"+Math.min(seq.getSeqLen(),lastMaxSupSeq.getSeqLen())+":"+(double)disssss/Math.min(seq.getSeqLen(),lastMaxSupSeq.getSeqLen()));
				}*/
				//System.out.println(dis+":"+minLen+":"+dis/minLen);
				if(dis<similarity && (lastMaxSupBp.getRightCode()==bp.getRightCode() && Math.abs(lastMaxSupBp.getRightPos()-bp.getRightPos())<100)) {
					availableList[i]=false;
					//lastMaxSupBp.updateRight(bp);
					lastMaxSupBp.mergeWith(bp);
				}else {
					if(bp.getSplitSupport()>maxSup && bp.hasRightPos()) {
						maxSup=bp.getSplitSupport();
						maxSupBp=bp;
						maxSupIdx=i;
					}
				}
			}
		}
	}
	
	private boolean isLargeOrINVPE(PESupport pes,Segment leftSeg) {
		final TranslationTable trTable=this.trTable;
		if(pes.isLeftForward()!=pes.isRightForward()) {
			return true;
		}
		if(pes.getLeftCode()==pes.getRightCode()) {
			if(Math.abs(pes.getLeftPos()-pes.getRightPos())>1000) {
				return true;
			}else {
				return false;
			}
		}
		Segment rightSeg=trTable.findSegment(pes.getRightCode());
		if(leftSeg.chrom.equals(rightSeg.chrom)) {
			long leftPos=leftSeg.start+pes.getLeftPos();
			long rightPos=rightSeg.start+pes.getRightPos();
			if(Math.abs(leftPos-rightPos)>1000) {
				return true;
			}else {
				return false;
			}
		}else{
			return true;
		}
	}
	
	private void updatePEAndAddPE(Segment leftSeg) {
		final ArrayList<PESupport> peListAfterMerge=this.peListAfterMerge;
		final ArrayList<SplitSupport> bpListC=this.bpListC;
		final SVCollection svCollection=this.svCollection;
		
		final int len=bpListC.size();

		
		for(PESupport pes:peListAfterMerge) {
			int idx=Collections.binarySearch(bpListC, pes);
			idx=idx<0?(-idx-1):idx;
			if(idx>=len) {
				continue;
			}
			if(pes.isLeftForward()) {
				for(int i=idx;i<len;++i) {
					SplitSupport bp=bpListC.get(i);
					if(bp.getLeftPos()-pes.getLeftPos()>200) {
						if(isLargeOrINVPE(pes,leftSeg)) {
							svCollection.adjustAndAdd(pes);
							break;
						}else {
							break;
						}
					}
					if(bp.similarWith(pes)) {
						bp.mergeWith(pes);
						break;
					}
				}
			}else {
				for(int i=idx;i>=0;--i) {
					SplitSupport bp=bpListC.get(i);
					if(pes.getLeftPos()-bp.getLeftPos()>200) {
						if(isLargeOrINVPE(pes,leftSeg)) {
							svCollection.adjustAndAdd(pes);
							break;
						}else {
							break;
						}
					}
					if(bp.similarWith(pes)) {
						bp.mergeWith(pes);
						break;
					}
				}
			}
		}
		peListAfterMerge.clear();
	}

	
	private void addSplit() {
		final ArrayList<SplitSupport> bpListC=this.bpListC;
		final SVCollection svCollection=this.svCollection;;
		//final int minIndelLen=this.minIndelLen;
		//final int minSupport=this.minSupport;
		
		for(SplitSupport bp:bpListC) {
			/*if((bp.getLeftCode()==809501704 && bp.getRightCode()==270549520)||(bp.getLeftCode()==270549520 && bp.getRightCode()==809501704)){
				System.out.print("FFFF "+bp.getLeftCode()+" + "+bp.getRightCode());
			}*/
			svCollection.adjustAndAdd(bp);
		}
		
		bpListC.clear();
	}

	
	private void clean() {
		this.forRefCntList.clean();
		this.revRefCntList.clean();
		
		this.availableList=null;
		this.bpListA=null;
		this.bpListB=null;
		this.bpListC=null;
		this.peList=null;
		this.seqListA=null;
		this.seqListC=null;
		this.seqListE=null;
		this.seqListB=null;
		this.forRefCntList=null;
		this.revRefCntList=null;
	}
	
	@Override
	public void run() {
		final TranslationTable trTable=this.trTable;
		final BpIndexList.SynchronizedGetItr synGetItr=this.synGetItr;

		CodeStartLen pack=new CodeStartLen();
		
		try {
			while((pack=synGetItr.getNext(pack))!=null) {
				int code=pack.code;
				long start=pack.start;
				int len=pack.len;
				
				this.code=code;
				Segment segment=trTable.findSegment(code);
				
				if(segment!=null) {
					fillKmerList(segment);
	    		    readBreakpoints(start,len);
	    		    if(!peList.isEmpty()) {
	    		    	prepareMergePEEqual();
		    		    mergePEEqual();
	    		    }	    	
	    		    if(!seqListA.isEmpty()) {
	                    prepareMergeRight();
	    		        mergeEqual();
	         		    mergeRightEnd();
	        		    prepareMergeLeftEnd();
	        		    mergeLeftEnd();
	        		    prepareMergeCluster();
	        		    mergeCluster();
	        		    mergeClusterAgain();
	    		    }
	    		    updatePEAndAddPE(segment);
	    		    addSplit();
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}finally {
			clean();
		}
		
	}
}

class RefCountList{
	private int[] refKmerList;
	private byte[] cntList;
	private int start,end,mergeK,k,refKmerListLen;
	private boolean isForward,init;
	
	public RefCountList(int mergeK,int k,boolean isForward){
		this.cntList=new byte[1<<(mergeK*2)];
		this.isForward=isForward;
		this.init=false;
		this.mergeK=mergeK;
		this.k=k;
	}
	
	public String toString() {
		StringBuilder str=new StringBuilder(1000);
		String tab="\t";
		/*for(int i:cntList) {
	     	str.append(i);
	     	str.append(tab);
		}
		str.append("\n");*/
		str.append(this.isForward);
		if(this.isForward) {
			for(int i=start;i<end;++i) {
		    	str.append(DNASequence.decompressKmer(refKmerList[i],this.mergeK).toString());
		    	//str.append(find(refKmerList[i]));
		    	str.append(tab);
			}
		}else {
			for(int i=start;i>end;--i) {
				str.append(DNASequence.decompressKmer(refKmerList[i],this.mergeK).toString());
				//str.append(find(refKmerList[i]));
		     	str.append(tab);
			}
		}
		str.append("\n");
		return str.toString();
	}
	
	public void loadRefKmerList(int[] refKmerList,int refKmerListLen) {
		this.refKmerList=refKmerList;
		this.refKmerListLen=refKmerListLen;
		this.init=false;
	}
	
	private void initRefCntList(SplitSupport bp0) {
		final int[] refKmerList=this.refKmerList;
		final byte[] cntList=this.cntList;
		final BpSeq bpSeq0=bp0.getBpSeq();
		
		for(int i=0,len=cntList.length;i<len;++i) {
			cntList[i]=0;
		}
		
		int start=0,end=0;
		
		if(this.isForward) {
			start=bp0.getLeftPos()+1;
			end=start+bpSeq0.getSeqLen()-this.mergeK+10;
			end=end<this.refKmerListLen?end:this.refKmerListLen;
			for(int i=start;i<end;++i) {
				//System.out.print(DNASequence.decompressKmer(refKmerList[i],this.mergeK).toString()+" ");
				cntList[refKmerList[i]]++;
			}
		}else {
			start=bp0.getLeftPos()+this.k-this.mergeK-1;
			end=start-bpSeq0.getSeqLen()+this.mergeK-10;
			end=end>=0?end:0;
			for(int i=start;i>end;--i) {
				//System.out.print(DNASequence.decompressKmer(refKmerList[i],this.mergeK).toString()+" ");
				cntList[refKmerList[i]]++;
			}
		}
		
		//System.out.println(start+" "+end);
		
		this.start=start;
		this.end=end;
	}
	
	public void updateRefCntList(SplitSupport bp) {
		if(!init) {
			initRefCntList(bp);
			this.init=true;
			return;
		}
		
		final int[] refKmerList=this.refKmerList;
		final byte[] cntList=this.cntList;
		final BpSeq bpSeq=bp.getBpSeq();
		
		if(this.isForward) {
			int newStart=bp.getLeftPos()+1;
			int newEnd=newStart+bpSeq.getSeqLen()-this.mergeK+10;
			newEnd=newEnd<this.refKmerListLen?newEnd:this.refKmerListLen;
			
			if(newStart<end-this.mergeK) {
				for(int i=start;i<newStart;++i) {
					--cntList[refKmerList[i]];
				}
				if(end<newEnd) {
					for(int i=end;i<newEnd;++i) {
						++cntList[refKmerList[i]];
					}
				}else{
					for(int i=newEnd;i<end;++i) {
						--cntList[refKmerList[i]];
					}
				}	
				this.start=newStart;
				this.end=newEnd;
				
			}else {
				initRefCntList(bp);
			}
			
		}else {	
			int newStart=bp.getLeftPos()+this.k-this.mergeK-1;
			int newEnd=newStart-bpSeq.getSeqLen()+this.mergeK-10;
			newEnd=newEnd>=0?newEnd:0;
			
			if(newEnd<start-this.mergeK) {
				for(int i=newStart;i>start;--i) {
					++cntList[refKmerList[i]];
				}
				if(newEnd>end) {
					for(int i=newEnd;i>end;--i) {
						--cntList[refKmerList[i]];
					}
				}else {
					for(int i=end;i>newEnd;--i) {
						++cntList[refKmerList[i]];
					}
				}
				this.start=newStart;
				this.end=newEnd;
				
			}else {
				initRefCntList(bp);
			}
		}
	}
	
	public int getLen() {
		return this.isForward?this.end-this.start:this.start-this.end;
	}
	
	public boolean find(int kmer) {
		//System.out.print(this.cntList[kmer]);
		return this.cntList[kmer]>0;
	}
	
	public void clean() {
		this.refKmerList=null;
		this.cntList=null;
	}
}

class KmerInfo{
	public int support;
	public int[] source;
	private int idx;
	
	KmerInfo(int support,int source){
		this.support=support;
		this.source=new int[10];
		this.source[0]=source;
		this.idx=1;
	}
	
	void update(int support,int source) {
		this.support+=support;
		if(this.idx>=this.source.length) {
			int[] sourceTmp=new int[this.idx*2];
			System.arraycopy(this.source, 0, sourceTmp, 0, this.idx);
			this.source=sourceTmp;
		}
		this.source[this.idx++]=source;
	}
	
	boolean available() {
		return this.idx>0;
	}
	
	int getLen() {
		return this.idx;
	}
	
	void reset() {
		this.support=0;
		this.source=null;
		this.idx=0;
	}
}

