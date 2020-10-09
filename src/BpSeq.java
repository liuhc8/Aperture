import java.util.HashSet;

public class BpSeq implements Comparable<BpSeq>{
	private final static BpSeq NULLSEQ;
	
	private long[] seq,bpSeq,bpSeqRev;
	private int[] kmerList;
	//private byte[] kmerCntList;
	private HashSet<Integer> kmerSet;
	private int k,bpLen,seqLen;
	private short leftBpLen,rightBpLen;
	private SplitSupport bp;
	
	static {
		NULLSEQ=new BpSeq();
	}
	
	public BpSeq(long[] seq,int bpLen,int seqLen,int k,short leftBpLen,short rightBpLen,SplitSupport bp) {	
		this.seq=seq;
		
		int longBpSeqLen=(bpLen-1)/32+1;
		long[] bpSeq=new long[longBpSeqLen];
		System.arraycopy(seq, 0, bpSeq, 0, longBpSeqLen);
		
		int notZeroCnt=(bpLen&0x1F)*2;
		int zeroCnt=64-notZeroCnt;
		bpSeq[longBpSeqLen-1] >>>= zeroCnt;
		bpSeq[longBpSeqLen-1] <<= zeroCnt;
		
		this.bpSeq=bpSeq;
		this.bpSeqRev=null;
		
		this.k=k;
		
		this.bpLen=bpLen;
		this.seqLen=seqLen;
		
		this.leftBpLen=leftBpLen;
		this.rightBpLen=rightBpLen;
		this.bp=bp;
	}
	
	private BpSeq() {}
	
	public static BpSeq getNullSeq() {
		return NULLSEQ;
	}
	
	public boolean isNullSeq() {
		return NULLSEQ==this;
	}
	
	public void reverse() {
		if(this.bpSeqRev==null) {
			this.bpSeqRev=DNASequence.getRevComp(bpSeq, bpLen);
		}
		
		long[] temp=this.bpSeq;
		this.bpSeq=this.bpSeqRev;
		this.bpSeqRev=temp;
		
		short tmp=this.leftBpLen;
		this.leftBpLen=(short)this.rightBpLen;
		this.rightBpLen=(short)tmp;
	}
	
	public void prepareKmer() {
		final int k=this.k;
		final long[] seq=this.seq;
		
		int middleLen=seq.length-1;
		long[] middle=new long[seq.length];
		for(int i=0;i<middleLen;++i) {
			middle[i]= (seq[i]<<32) | (seq[i+1]>>>32);
		}
		middle[middleLen] = (seq[middleLen]<<32);
		
		int kmerListLen=this.seqLen-k+1;
		long mask=(1L<<(2*k))-1;
		int[] kmerList=new int[kmerListLen];
		
		int start=64-2*k,end=34-2*k;
		outer:
		for(int i=0,j=0;;++j) {
			for(int p=start;p>=end;p-=2) {
				kmerList[i++]=(int)((seq[j]>>>p)&mask);
				if(i>=kmerListLen) {
					break outer;
				}
			}
			for(int p=start;p>=end;p-=2) {
				kmerList[i++]=(int)((middle[j]>>>p)&mask);
				if(i>=kmerListLen) {
					break outer;
				}
			}
		}
		
		this.kmerList=kmerList;
	}
	
	/*private void prepareCntList() {
		this.kmerCntList=new byte[1<<(this.k*2)];
		for(int i=0,len=kmerList.length;i<len;++i) {
			++kmerCntList[kmerList[i]];
		}
	}*/
	
	private void prepareKmerSet() {
		HashSet<Integer> kmerSet=new HashSet<Integer>(kmerList.length);
		for(int i=0,len=kmerList.length;i<len;++i) {
			kmerSet.add(kmerList[i]);
		}
		this.kmerSet=kmerSet;
	}
	
	public double getDistance(BpSeq bpX) {
		if(this.kmerSet==null) {
			prepareKmerSet();
		}
		
		final int[] bpXKmerList=bpX.kmerList;
		final HashSet<Integer> kmerSet=this.kmerSet;
		
		int overlap=0,j=0;
		for(int i=0,len=bpXKmerList.length;i<len;++i) {
			if(kmerSet.contains(bpXKmerList[i])) {
				++j;
			}else {
				if(j>0) {
					overlap+=j+this.k-1;
					j=0;
				}
			}
		}
		if(j>0) {
			overlap+=j+k-1;
		}
		if(bpX.seqLen>this.seqLen) {
			overlap-=(bpX.seqLen-this.seqLen);
		}
		return 1-(((double)overlap)/bpX.seqLen);
	}
	
	/*public double getDistance(BreakpointSeq bpX) {
		if(this.kmerCntList==null) {
			prepareCntList();
		}
		
		final int[] bpXKmerList=bpX.kmerList;
		final byte[] kmerCntList=this.kmerCntList;
		
		int overlap=0,j=0;
		for(int i=0,len=bpXKmerList.length;i<len;++i) {
			if(kmerCntList[bpXKmerList[i]]>0) {
				++j;
			}else {
				if(j>0) {
					overlap+=j+this.k-1;
					j=0;
				}
			}
		}
		if(j>0) {
			overlap+=j+k-1;
		}
		if(bpX.seqLen>this.seqLen) {
			overlap-=(bpX.seqLen-this.seqLen);
		}
		return 1-(((double)overlap)/bpX.seqLen);
	}*/
	
	public double getDistance(RefCountList refCntList) {
		final int[] kmerList=this.kmerList;
		final int k=this.k;
		
		int len=Math.min(kmerList.length,refCntList.getLen()+5);
		int overlap=0,j=0;
		for(int i=0;i<len;++i) {
			boolean res=refCntList.find(kmerList[i]);
			if(res) {
				++j;
			}else {
				if(j>0) {
					overlap+=j+k-1;
					j=0;
					//i+=k-1;
				}
			}
		}
		if(j>0) {
			overlap+=j+k-1;
		}
		//System.out.print(" "+overlap+"   "+this.seqLen+" ");
		//int len=Math.min(this.seqLen, refCntList.getLen());
		return 1-(((double)overlap)/(len+k-1));
	}
	
	public boolean shareBpWith(BpSeq o) {
		final long[] seq1=this.bpSeq;
		final long[] seq2=o.bpSeq;
		final int len1=seq1.length;
		final int len2=seq2.length;
		
		if(len1!=len2) {
			return false;
		}
		
		for(int i=0;i<len1;++i) {
			if(seq1[i]!=seq2[i]) {
				return false;
			}
		}
		
		return true;
	}
	
	public boolean isPrefix(BpSeq o) {
		final long[] seq1=this.bpSeq;
		final long[] seq2=o.bpSeq;
		final int len1=seq1.length;
		final int len2=seq2.length;
		final int minLen=Math.min(len1, len2);
		
		int consensus=0;
		for(int i=0;i<minLen;++i) {
			long xor=seq1[i]^seq2[i];
			if(xor==0L) {
				consensus+=64;
			}else {
				consensus+=Long.numberOfLeadingZeros(xor);
				break;
			}
		}
		
		if(Math.min(this.bpLen,o.bpLen)<=consensus/2) {
			return true;
		}else {
			return false;
		}
	}
	
	
	/*public void update(BreakpointSeq o) {
		this.seq=o.seq;
		this.bpSeq=o.bpSeq;
		this.bpLen=o.bpLen;
		this.seqLen=o.seqLen;
		this.leftBpLen=o.leftBpLen;
		this.rightBpLen=o.rightBpLen;
	}*/
	
	public long[] getBpSeq() {
		return this.bpSeq;
	}
	
	public long[] getRevBpSeq() {
		return this.bpSeqRev;
	}
	
	public short getLeftBpLen() {
		return this.leftBpLen;
	}
	
	public short getRightBpLen() {
		return this.rightBpLen;
	}
	
	public int[] getKmerList() {
		return this.kmerList;
	}
	
	public int getBpLen() { 
		return this.bpLen;
	}
	
	public int getSeqLen() {
		return this.seqLen;
	}
	
	public SplitSupport getSVSupport() {
		return this.bp;
	}
	
	public void setSVSupport(SplitSupport ss) {
		this.bp=ss;
	}

	public String toStringBp() {
		byte[] byteseq=new byte[this.bpLen];
		DNASequence.decompressSeq(byteseq,this.bpSeq,this.bpLen);
		return new String(byteseq);
	}
	
	public String toStringBp(int offset,int length) {
		byte[] byteseq=new byte[this.bpLen];
		DNASequence.decompressSeq(byteseq,this.bpSeq,this.bpLen);
		
		return new String(byteseq,offset,length);
	}
	
	public String toStringRevBp() {
		byte[] byteseq=new byte[this.bpLen];
		DNASequence.decompressSeq(byteseq,this.bpSeqRev,this.bpLen);
		return new String(byteseq);
	}
	
	public String toStringRevBp(int offset,int length) {
		byte[] byteseq=new byte[this.bpLen];
		DNASequence.decompressSeq(byteseq,this.bpSeqRev,this.bpLen);
		return new String(byteseq,offset,length);
	}

	public String toStringAll() {
		byte[] byteseq=new byte[this.seqLen];
		DNASequence.decompressSeq(byteseq,this.seq,this.seqLen);
		return new String(byteseq);
	}
	
	public void clean() {
		this.seq=null;
		this.kmerSet.clear();
		this.kmerSet=null;
		this.kmerList=null;
		this.bp=null;
		this.seq=null;
		this.bpSeq=null;
		this.bpSeqRev=null;
	}
	

	@Override
	public int compareTo(BpSeq o) {    //descending order
		final long[] seq1=this.bpSeq;
		final long[] seq2=o.bpSeq;
		final int len1=seq1.length;
		final int len2=seq2.length;
		final int minLen=Math.min(len1, len2);
		
		for(int i=0;i<minLen;++i) {
			if(seq1[i]<seq2[i]) {
				return 1;
			}else if(seq1[i]>seq2[i]) {
				return -1;
			}
		}
		if(this.bpLen<o.bpLen) {
			return 1;
		}else if(this.bpLen>o.bpLen) {
			return -1;
		}else {
			return 0;
		}
	}

	@Override
	public boolean equals(Object anObject) {
		if (this == anObject) {
            return true;
        }
        if (anObject instanceof BpSeq) {
        	BpSeq anotherBpSeq = (BpSeq)anObject;;
            if (this.compareTo(anotherBpSeq)==0) {
                return true;
            }
        }
        return false;
	}
}


