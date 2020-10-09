import java.util.Arrays;

public abstract class SVSupport implements Comparable<SVSupport>{
	protected byte flag;    // 0x1:Left breakend has smaller genomic position; 0x2:Left breakend has greater genomic position
	protected boolean rightVague,leftForward,rightForward;
	protected int leftCode,leftPos,leftKmers;
	protected int rightCode,rightPos,rightKmers;
	protected double leftScoreSum,rightScoreSum;
	protected int[] r1BarList,r2BarList;
	protected boolean[] flipList;
	protected int barCnt;
	protected Barcode[] longBarList;
	
	public SVSupport(int leftCode,int leftPos,short leftConfi,short leftKmers,
            int rightCode,int rightPos,short rightConfi,short rightKmers,boolean rightVague,
            int r1Bar,int r2Bar,boolean flip,boolean isSecondary) {	
		this.leftCode=leftCode;
		this.leftPos=leftPos;
		this.leftKmers=leftKmers;
		this.leftScoreSum=isSecondary?0:Math.abs((double)leftConfi/leftKmers);
		this.leftForward=leftConfi>0;

		this.rightCode=rightCode;
		this.rightPos=rightPos;
		this.rightKmers=rightKmers;
		this.rightScoreSum=(rightPos==-1||isSecondary)?0:Math.abs((double)rightConfi/rightKmers);
		this.rightForward=rightConfi>0;
		
		this.rightVague=rightVague;
		this.flag=0;	

	    this.r1BarList=new int[5];
	    this.r2BarList=new int[5];
        this.flipList=new boolean[5];
        	
        if(isSecondary) {
        	this.barCnt=0;
        }else {
    	    this.r1BarList[0]=r1Bar;
     	    this.r2BarList[0]=r2Bar;
     	     this.flipList[0]=flip;
     	    this.barCnt=1;
        }
	}

	public void setLeftSmallFlag() {
		this.flag=1;
	}
	
	public void setLeftHighFlag() {
		this.flag=2;
	}
	
	public boolean hasSameDirect(SVSupport sv) {
		return this.leftForward==sv.leftForward && this.rightForward==sv.rightForward;
	}
	
	public boolean isLeftForward() {
		return this.leftForward;
	}
	
	public boolean isRightForward() {
		return this.rightForward;
	}
	
	public int getLeftPos() {
		return this.leftPos;
	}

	public int getLeftCode() {
		return this.leftCode;
	}
	
	public int getRightCode() {
		return this.rightCode;
	}
	
	public int getRightPos() {
		return this.rightPos;
	}
	
	public boolean isRightVague() {
		return this.rightVague;
	}
	
	public boolean hasRightPos() {
		return this.rightPos!=-1;
	}
	
	private void prepareLongBarList() {
		if(this.longBarList!=null) {
			return;
		}
		final int[] r1BarList=this.r1BarList;
		final int[] r2BarList=this.r2BarList;
		final boolean[] flipList=this.flipList;
		final int barCnt=this.barCnt;
		
		Barcode[] longBarList=new Barcode[barCnt];
		for(int i=0;i<barCnt;++i) {
			longBarList[i]=new Barcode(r1BarList[i],r2BarList[i],flipList[i]);
		}
	    Arrays.sort(longBarList,0,barCnt);
		
	    this.longBarList=longBarList;
	}
	
	public double molecularCnt() {
		prepareLongBarList();
		final Barcode[] longBarList=this.longBarList;
		
		double molecularCnt=1;
		for(int i=1,len=longBarList.length;i<len;++i) {
			if(longBarList[i].getBar()!=longBarList[i-1].getBar()) {
				++molecularCnt;
			}else {
				if(longBarList[i].getFlip()!=longBarList[i-1].getFlip()) {
					molecularCnt+=0.5;
				}
			}
		}
		return molecularCnt;
	}
	
	public int uniqueBarCnt() {
		prepareLongBarList();
		final Barcode[] longBarList=this.longBarList;
		
		int uniqueBarCnt=0,len=longBarList.length;
		if(len==1) {
			uniqueBarCnt=1;
		}else {
	    	for(int i=0;i<len;++i) {
	     		if(i==0) {
		    		if(longBarList[i].getBar()!=longBarList[i+1].getBar()) {    
			    		++uniqueBarCnt;
			    	}
		     	}else if(i==len-1) {
			    	if(longBarList[i].getBar()!=longBarList[i-1].getBar()) {
			    		++uniqueBarCnt;
		    		}
		    	}else {
		    		if(longBarList[i].getBar()!=longBarList[i-1].getBar() && longBarList[i].getBar()!=longBarList[i+1].getBar()) {
		    			++uniqueBarCnt;
		    		}
	     		}
	    	}
		}
		return uniqueBarCnt;
	}
	
	public void updateBarList(SVSupport sv) {
		if(this.barCnt+sv.barCnt>this.r1BarList.length) {
			int newCapacity=(this.barCnt+sv.barCnt)*2;
			int[] newR1BarList=new int[newCapacity];
			System.arraycopy(this.r1BarList, 0, newR1BarList, 0, this.barCnt);
			this.r1BarList=newR1BarList;
			
			int[] newR2BarList=new int[newCapacity];
			System.arraycopy(this.r2BarList, 0, newR2BarList, 0, this.barCnt);
			this.r2BarList=newR2BarList;
			
			boolean[] newFlipList=new boolean[newCapacity];
			System.arraycopy(this.flipList, 0, newFlipList, 0, this.barCnt);
			this.flipList=newFlipList;
		}
		
		System.arraycopy(sv.r1BarList, 0, this.r1BarList, this.barCnt, sv.barCnt);
		System.arraycopy(sv.r2BarList, 0, this.r2BarList, this.barCnt, sv.barCnt);
		System.arraycopy(sv.flipList, 0, this.flipList, this.barCnt, sv.barCnt);
		this.barCnt+=sv.barCnt;
	}

	public boolean similarWith(SVSupport sv) {
		if(this.leftCode!=sv.leftCode || this.rightCode!=sv.rightCode || !this.hasSameDirect(sv)) {
			return false;
		}else {
			if(this instanceof SplitSupport && sv instanceof SplitSupport) {
				if(this.leftCode==this.rightCode && Math.abs(this.leftPos-this.rightPos)<200) {
					return Math.abs(this.leftPos-sv.leftPos)<15 && Math.abs(this.rightPos-sv.rightPos)<15;
				}else {
					return Math.abs(this.leftPos-sv.leftPos)<100 && Math.abs(this.rightPos-sv.rightPos)<100;
				}
			}else {
				return Math.abs(this.leftPos-sv.leftPos)<200 && Math.abs(this.rightPos-sv.rightPos)<200;
			}
		}
	}
	
	public boolean mergeWith(SVSupport sv) {
		if(sv instanceof SplitSupport) {
			return mergeWith((SplitSupport)sv);
		}else {
			return mergeWith((PESupport)sv);
		}
	}
	
	 public boolean hasShorterJunction(SVSupport sv) {
		 return this.getBpLen()<=sv.getBpLen();
	 }
	
	@Override
	public int compareTo(SVSupport o) {
		return this.leftPos-o.leftPos;
	}
	
	@Override
	public abstract boolean equals(Object anObject);	
	public abstract boolean mergeWith(SplitSupport bp) ;
	public abstract boolean mergeWith(PESupport pes);
	public abstract void flip();
	public abstract int getBpLen();
	public abstract double getLeftScoreAvg();
	public abstract double getRightScoreAvg();
	public abstract SVCandidate translateAndFilter(TranslationTable trTable,SVCandidate svc);
}

class SplitSupport extends SVSupport{
	private double bpScoreSum;
	private BpSeq seq;
	private short splitCnt,peCnt,leftCnt,rightCnt;

	
	public SplitSupport(boolean rightVague,int leftCode,int leftPos,short leftConfi,short leftKmers,
	                   int rightCode,int rightPos,short rightConfi,short rightKmers,int bpScore,int bpLen,
	                   int r1Bar,int r2Bar,boolean flip,boolean isSecondary) {	
		super(leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,rightVague,r1Bar,r2Bar,flip,isSecondary);
		this.bpScoreSum=(double)bpScore/bpLen;
		this.peCnt=0;
		if(isSecondary) {
			this.splitCnt=0;
			this.leftCnt=0;
			this.rightCnt=0;
		}else {
			this.splitCnt=1;
			this.leftCnt=1;
			this.rightCnt=rightPos==-1?(short)0:(short)1;
		}
	}
	
	@Override
	public boolean mergeWith(SplitSupport bp) {	
	    updateBarList(bp);
		
    	this.leftKmers+=bp.leftKmers;
    	this.rightKmers+=bp.rightKmers;
		
	    this.leftCnt+=bp.leftCnt;
	    this.rightCnt+=bp.rightCnt;
	    this.splitCnt+=bp.splitCnt;
	    this.peCnt+=bp.peCnt;
		
	   	this.leftScoreSum+=bp.leftScoreSum;
	   	this.rightScoreSum+=bp.rightScoreSum;
		
	    this.bpScoreSum+=bp.bpScoreSum;
	    
	    this.flag|=bp.flag;
	    
	    return true;
	}

	
	@Override
	public boolean mergeWith(PESupport pes) {	
		updateBarList(pes);
		
		this.leftKmers+=pes.leftKmers;
		this.rightKmers+=pes.rightKmers;
		
		this.peCnt+=pes.getPECnt();
		
		this.leftScoreSum+=pes.leftScoreSum;
		this.rightScoreSum+=pes.rightScoreSum;

		this.flag|=pes.flag;
		
		return true;
	}
	
	public void setBpSeq(BpSeq bpSeq) {
		this.seq=bpSeq;
	}
	
	public BpSeq getBpSeq() {
		return this.seq;
	}

	@Override
	public int getBpLen() {
		return this.seq.getBpLen();
	}
	
	
	@Override
	public double getLeftScoreAvg() {
		return Math.abs(this.leftScoreSum)/(double)(this.leftCnt+this.peCnt);
	}
	
	@Override
	public double getRightScoreAvg() {
		return Math.abs(this.rightScoreSum)/(double)(this.rightCnt+this.peCnt);
	}
	
	public int getSupportCnt() {
		return this.peCnt+this.splitCnt;
	}
	
	public int getSplitSupport() {
		return this.splitCnt;
	}
	
	public int getLeftCnt() {
		return this.leftCnt;
	}
	
	public int getRightCnt() {
		return this.rightCnt;
	}

	@Override
	public void flip() {
		int temp=this.leftCode;
		this.leftCode=this.rightCode;
		this.rightCode=temp;
		
		temp=this.leftPos;
		this.leftPos=this.rightPos;
		this.rightPos=temp;
		
		temp=this.leftKmers;
		this.leftKmers=this.rightKmers;
		this.rightKmers=temp;
		
		double temp2=this.leftScoreSum;
		this.leftScoreSum=this.rightScoreSum;
		this.rightScoreSum=temp2;
	
		short temp3=this.leftCnt;
		this.leftCnt=this.rightCnt;
		this.rightCnt=temp3;
		
		boolean temp4=this.leftForward;
		this.leftForward=!this.rightForward;
		this.rightForward=!temp4;
		
		
	    int [] tempList=this.r1BarList;
	    this.r1BarList=this.r2BarList;
	    this.r2BarList=tempList;
		
	    final int barCnt=this.barCnt;
	    final boolean[] flipList=this.flipList;
	    for(int i=0;i<barCnt;++i) {
	    	flipList[i]=!flipList[i];
	    }
	    
		this.seq.reverse();
	}
	
	@Override
	public SVCandidate translateAndFilter(TranslationTable trTable,SVCandidate svc) {
		
		if(this.splitCnt+this.peCnt==0) {
			return null;
		}
		
		Segment leftSeg=trTable.findSegment(this.leftCode);
		Segment rightSeg=trTable.findSegment(this.rightCode);
		if(rightSeg==null) {
			return null;
		}
		int k=trTable.getK();
		
		svc.setPassFakeBpFilter(flag);
		svc.setLeft(leftSeg,this.leftPos,this.leftKmers, getLeftScoreAvg(),this.leftForward,k);
		
		svc.setRight(rightSeg,this.rightPos,this.rightKmers, getRightScoreAvg(),this.rightForward,k);

		/*int leftIdentical,rightIdentical;
		synchronized(ApertureMain.syn) {
		System.out.println("AAAAAAA");
		leftIdentical=seqAlignmentLeft(leftSeg,k);
		System.out.println();
		System.out.println(seq.toStringBp());
		System.out.println(leftIdentical);
		System.out.println("BBBBBB");
		rightIdentical=seqAlignmentRight(rightSeg,k);
		System.out.println();
		System.out.println(seq.toStringRevBp());
		System.out.println(rightIdentical);
		}*/
		
		int leftIdentical=seqAlignmentLeft(leftSeg,k);
		int rightIdentical=seqAlignmentRight(rightSeg,k);
		
		int bpLen=this.seq.getBpLen();
		if(leftIdentical==bpLen || leftIdentical<=0 || rightIdentical==bpLen || rightIdentical<=0) {
			return null;
		}
		
		svc.setSplit(seq, bpScoreSum/splitCnt, splitCnt, peCnt, leftCnt, rightCnt,leftIdentical,rightIdentical);
		svc.setUMI(molecularCnt(), uniqueBarCnt());
		
		svc.setCode(this.leftCode, this.leftPos, this.rightCode, this.rightPos);
		return svc;
	}


	private int seqAlignmentLeft(Segment leftSeg,int k) {
		int bpLen=this.seq.getBpLen();
		int leftBpLen=this.seq.getLeftBpLen();
		int maxLen=leftBpLen>=3?leftBpLen+k:3+k;
		maxLen=maxLen>bpLen?bpLen:maxLen;
		if(leftForward) {
			//System.out.println("leftForward "+leftBpLen);
			
			int identicalCnt=countIdenticalFor(leftSeg.ref,this.leftPos+1,maxLen,this.seq.getBpSeq());
			return identicalCnt>bpLen?bpLen:identicalCnt;
		}else {
			//System.out.println("leftRev "+leftBpLen);
			int identicalCnt=countIdenticalRev(leftSeg.ref,this.leftPos+k-1,maxLen,this.seq.getBpSeq());
			return identicalCnt>bpLen?bpLen:identicalCnt;
		}
	}
	
	private int seqAlignmentRight(Segment rightSeg,int k) {
		int bpLen=this.seq.getBpLen();
		int rightBpLen=this.seq.getRightBpLen();
		int maxLen=rightBpLen>=3?rightBpLen+k:3+k;
		maxLen=maxLen>bpLen?bpLen:maxLen;
		if(rightForward) {
			//System.out.println("rightForward "+rightBpLen);
			int identicalCnt=countIdenticalRev(rightSeg.ref,this.rightPos+k-1,maxLen,this.seq.getRevBpSeq());
			return identicalCnt>bpLen?bpLen:identicalCnt;
		}else {
			//System.out.println("rightRev "+rightBpLen);
			int identicalCnt=countIdenticalFor(rightSeg.ref,this.rightPos+1,maxLen,this.seq.getRevBpSeq());
			return identicalCnt>bpLen?bpLen:identicalCnt;
		}
	}
	
	private int countIdenticalFor(long[] ref,int start,int len,long[] bpSeq) {
		int sameBaseCnt=0;
		int blockPos=start>>>5;   // divide by 32
        int offset=start&0x1F;
        long sub;
        boolean end=false;
        for(int i=0;len>0;++i) {
        	if(blockPos+1==ref.length) {
        		sub=(ref[blockPos]<<(offset<<1))|(((1<<(offset<<1))-1)&(~bpSeq[i]));
        		end=true;
        	}else if(blockPos+1<ref.length){
        		sub=(ref[blockPos]<<(offset<<1))|(ref[blockPos+1]>>>((32-offset)<<1));
        	}else {
        		return sameBaseCnt;
        	}
	    	//System.out.print(DNASequence.decompressKmer(sub,32).toString()+" ");
	        int cnt=Long.numberOfLeadingZeros(sub^bpSeq[i])>>>1;
        	if(cnt<=10 && sameBaseCnt==0) {
        		int move=(cnt+1)<<1;
        		cnt=Long.numberOfLeadingZeros(((sub^bpSeq[i])<<move)>>>move)>>>1;
        	}
	        sameBaseCnt+=cnt;
	        if(cnt<32||end) {
	        	return sameBaseCnt;
	        }
	        ++blockPos;
	        len-=32;
        }
		return sameBaseCnt;
	}
	
	private int countIdenticalRev(long[] ref,int start,int len,long[] bpSeq) {
		int sameBaseCnt=0;
		int blockPos=start>>>5;   // divide by 32
        int offset=start&0x1F;
        long sub;
        boolean end=false;
        for(int i=0;len>0;++i) {
        	if(blockPos==0) {
        		sub=ref[blockPos]>>>((32-offset)<<1);
        		sub=DNASequence.getRevComp(sub, 32);
        		sub|=(((1<<((32-offset)<<1))-1)&(~bpSeq[i]));
        		end=true;
        	}else if(blockPos>0){
        		sub=(ref[blockPos-1]<<(offset<<1))|(ref[blockPos]>>>((32-offset)<<1));
        		sub=DNASequence.getRevComp(sub, 32);
        	}else {
        		return sameBaseCnt;
        	}
	        //System.out.print(DNASequence.decompressKmer(sub,32).toString()+" ");
	        int cnt=Long.numberOfLeadingZeros(sub^bpSeq[i])>>>1;
        	if(cnt<=10 && sameBaseCnt==0) {
            	int move=(cnt+1)<<1;
            	cnt=Long.numberOfLeadingZeros(((sub^bpSeq[i])<<move)>>>move)>>>1;
            }
	        sameBaseCnt+=cnt;
	        if(cnt<32 ||end) {
	        	return sameBaseCnt;
	        }
	        --blockPos;
	        len-=32;
        }
		return sameBaseCnt;
	}
	
	@Override
	public boolean equals(Object anObject) {
		if (this == anObject) {
            return true;
        }
        if (anObject instanceof SplitSupport) {
        	SplitSupport anotherSplit = (SplitSupport)anObject;;
            if (seq.equals(anotherSplit.seq) && this.leftCode==anotherSplit.leftCode && this.rightCode==anotherSplit.rightCode 
            		&& this.leftPos==anotherSplit.leftPos && this.rightPos==anotherSplit.rightPos) {
                return true;
            }
        }
        return false;
	}

}


class PESupport extends SVSupport{	

	private short peCnt;
	
	public PESupport(int leftCode,int leftPos,short leftConfi,short leftKmers,
			         int rightCode,int rightPos,short rightConfi,short rightKmers,int r1Bar,int r2Bar,boolean flip,boolean isSecondary) {
		super(leftCode,leftPos,leftConfi,leftKmers,rightCode,rightPos,rightConfi,rightKmers,true,r1Bar,r2Bar,flip,isSecondary);
		
		this.peCnt=isSecondary?(short)0:(short)1;
	}
	
	public int getPECnt() {
		return this.peCnt;
	}

	@Override
	public boolean mergeWith(SplitSupport bp) {
		return false;
	}

	@Override
	public boolean mergeWith(PESupport pes) {
		updateBarList(pes);
		
		this.leftKmers+=pes.leftKmers;
		this.rightKmers+=pes.rightKmers;
		
		this.peCnt+=pes.peCnt;
		
		this.leftScoreSum+=pes.leftScoreSum;
		this.rightScoreSum+=pes.rightScoreSum;

		this.flag|=pes.flag;
		
		return true;
	}
	
	@Override
	public void flip() {
		int temp=this.leftCode;
		this.leftCode=this.rightCode;
		this.rightCode=temp;
		
		temp=this.leftPos;
		this.leftPos=this.rightPos;
		this.rightPos=temp;
		
		temp=this.leftKmers;
		this.leftKmers=this.rightKmers;
		this.rightKmers=temp;
		
		double temp2=this.leftScoreSum;
		this.leftScoreSum=this.rightScoreSum;
		this.rightScoreSum=temp2;
		
		boolean temp4=this.leftForward;
		this.leftForward=!this.rightForward;
		this.rightForward=!temp4;
	
		
	    int [] tempList=this.r1BarList;
	    this.r1BarList=this.r2BarList;
	    this.r2BarList=tempList;
		
	    final int barCnt=this.barCnt;
	    final boolean[] flipList=this.flipList;
	    for(int i=0;i<barCnt;++i) {
	    	flipList[i]=!flipList[i];
	    }
	    
	}
	
	@Override
	public int getBpLen() {
		return Integer.MAX_VALUE;
	}

	@Override
	public double getLeftScoreAvg() {
		return Math.abs(this.leftScoreSum)/this.peCnt;
	}
	
	@Override
	public double getRightScoreAvg() {
		return Math.abs(this.rightScoreSum)/this.peCnt;
	}

	//@Override
	public SVCandidate translateAndFilter(TranslationTable trTable,SVCandidate svc) {
		
		if(this.peCnt==0) {
			return null;
		}
		
		Segment leftSeg=trTable.findSegment(this.leftCode);
		Segment rightseg=trTable.findSegment(this.rightCode);
		if(rightseg==null) {
			return null;
		}
		
		final int k=trTable.getK();
		
		svc.setPassFakeBpFilter(flag);
		svc.setLeft(leftSeg,this.leftPos, this.leftKmers, getLeftScoreAvg(),this.leftForward,k);
		svc.setRight(rightseg,this.rightPos, this.rightKmers, getRightScoreAvg(),this.rightForward,k);
		
		svc.setPE(peCnt);
		svc.setUMI(molecularCnt(), uniqueBarCnt());
		
		svc.setCode(this.leftCode, this.leftPos, this.rightCode, this.rightPos);
		return svc;
	}
	
	@Override
	public boolean equals(Object anObject) {
		if (this == anObject) {
            return true;
        }
        if (anObject instanceof PESupport) {
        	PESupport anotherPE = (PESupport)anObject;;
            if (this.leftCode==anotherPE.leftCode && this.rightCode==anotherPE.rightCode 
            		&& this.leftPos==anotherPE.leftPos && this.rightPos==anotherPE.rightPos) {
                return true;
            }
        }
        return false;
	}
}

class Barcode implements Comparable<Barcode>{
	private long bar;
	private boolean flip;
	
	Barcode(int r1Bar,int r2Bar,boolean flip){
		this.bar=((long)r2Bar & 0xFFFFFFFFL) | (((long)r1Bar << 32) & 0xFFFFFFFF00000000L);
		this.flip=flip;
	}
	@Override
	public int compareTo(Barcode o) {
		long res=this.bar-o.bar;
		if(res<0L) {
			return -1;
		}else if(res>0L) {
			return 1;
		}else {
     		return (this.flip == o.flip) ? 0 : (this.flip ? 1 : -1);
		}
	}
	
	public long getBar(){
		return this.bar;
	}
	
	public boolean getFlip() {
		return this.flip;
	}
	
	public String toString() {
		StringBuilder str=new StringBuilder(100);
		String tab="\t";
		str.append(bar);
		str.append(tab);
		str.append(flip);
		str.append(tab);
		return str.toString();
	}
}

