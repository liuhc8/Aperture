import java.nio.ByteBuffer;
import java.util.Random;

class DNASequence {
	protected final static byte[] NT4_UPPER_REV_TABLE= {
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0,84, 0,71,  0, 0, 0,67,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 0, 0, 0, 65, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	};
	
	protected final static byte[] NT4_UPPER_REV_BIN_TABLE= {
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 3, 0, 2,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, -1, 0,
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	};
	
	protected final static byte[] NT4_BIN_TABLE= {
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
			0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0,-1, 0,
			0, 0, 0, 0,  3, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	};
	
	protected final static byte[] NT4_ASCII_TABLE= {
			65, 67, 71, 84,
	};
	
	protected final static byte[] NT4_BIN_REV16_TABLE= {
		    15,11, 7, 3,  14,10, 6, 2,  13, 9, 5, 1,  12, 8, 4, 0,
	};

	protected final static byte[] NT4_BIN_REV4_TABLE= {
		    3, 2, 1, 0,
	};
	
	protected final static Random rand=new Random();
	
	protected byte[] sequence;
	protected int sequenceLen;

	protected DNASequence(byte[] seq){
		this.sequence=seq;
		this.sequenceLen=seq.length;
	}
	
	protected DNASequence(byte[] seq,int start,int length){
		byte[] newseq=new byte[length];
		System.arraycopy(seq, start, newseq, 0, length);
		this.sequence=newseq;
		this.sequenceLen=length;
	}
	
	protected DNASequence(int len){           //create empty object
		this.sequence=new byte[len];
		this.sequenceLen=0;
	}
	
	static DNASequence createDNASequence(byte[] seq) {
		return new DNASequence(seq);
	}
	
	byte getBaseAt(int pos0) {
		return this.sequence[pos0];
	}
	
	byte[] getSequence() {
		return this.sequence;
	}
	
	int getLength() {
		return this.sequenceLen;
	}
	
	DNASequence getReverseComplementarity() {
		final int sequenceLen=this.sequenceLen;
		final int cycle=(sequenceLen+1)/2;
		final byte[] sequence=this.sequence;
		final byte[] nt4Table=DNASequence.NT4_UPPER_REV_TABLE;
		
		byte tmp=0;
		for(int i=0;i<cycle;++i) {
			tmp=nt4Table[sequence[sequenceLen-1-i]];
			sequence[sequenceLen-1-i]=nt4Table[sequence[i]];
			sequence[i]=tmp;
		}
		return this;
	}
	
	@Override
	public int hashCode() {
		int h=0;
		byte[] val = this.sequence;
		int len=val.length;
		
        for (int i = 0; i < len; i++) {
            h = 31 * h + val[i];
        }
        return h;
	}
	
	@Override
	public boolean equals(Object anObject) {
		if (this == anObject) {
            return true;
        }
        if (anObject instanceof DNASequence) {
        	DNASequence anotherSeq = (DNASequence)anObject;
            int n = this.sequenceLen;
            if (n == anotherSeq.sequenceLen) {
                byte[] v1 = this.sequence;
                byte[] v2 = anotherSeq.sequence;
                int i = 0;
                while (n-- != 0) {
                    if (v1[i] != v2[i])
                        return false;
                    i++;
                }
                return true;
            }
        }
        return false;
	}
	
	@Override
	 public String toString(){  
	        return new String(this.sequence); 
	 } 
	 
	 public void copyBases(byte[] dest,int destStart,int srcStart,int len) {
		 System.arraycopy(this.sequence, srcStart, dest, destStart, len);
	 }
	
	public static byte getRandomBase() {
		return (byte)(rand.nextInt()&3);
	}
	
	private static long reverseBase(long i) {
        i = (i & 0x3333333333333333L) << 2 | (i >>> 2) & 0x3333333333333333L;
        i = (i & 0x0f0f0f0f0f0f0f0fL) << 4 | (i >>> 4) & 0x0f0f0f0f0f0f0f0fL;
        i = (i & 0x00ff00ff00ff00ffL) << 8 | (i >>> 8) & 0x00ff00ff00ff00ffL;
        i = (i << 48) | ((i & 0xffff0000L) << 16) |
            ((i >>> 16) & 0xffff0000L) | (i >>> 48);
        return i;
    }
	
	private static int reverseBase(int i) {
	    i = (i & 0x33333333) << 2 | (i >>> 2) & 0x33333333;
	    i = (i & 0x0f0f0f0f) << 4 | (i >>> 4) & 0x0f0f0f0f;
	    i = (i << 24) | ((i & 0xff00) << 8) |
	        ((i >>> 8) & 0xff00) | (i >>> 24);
	    return i;
	}
	
	public static long getRevComp(long seq,int k) {
		seq=reverseBase(seq);
		seq=~seq;
		seq>>>=((32-k)*2);
		return seq;
	}
	
    public static int getRevComp(int seq,int k) {	
		seq=reverseBase(seq);
		seq=~seq;
		seq>>>=((16-k)*2);
		return seq;
	}
    
    public static long[] getRevComp(long[] seq,int len) {
		int longLen=(len-1)/32+1;
		int nonZeroCnt=(len&0x1F)*2;
		int zeroCnt=64-nonZeroCnt;
		long[] seqRevComp=new long[longLen];
		seqRevComp[0]=(~reverseBase(seq[longLen-1]))<<zeroCnt;
		for(int i=1;i<longLen;++i) {
			long s=~reverseBase(seq[longLen-1-i]);
			seqRevComp[i-1] |= (s>>>nonZeroCnt);
			seqRevComp[i] = (s<<zeroCnt);
		}
		return seqRevComp;
	}

	public static DNASequence decompressKmer(long seq,int k) {
		final byte[] nt4_ascii=NT4_ASCII_TABLE;
		byte[] byteseq=new byte[k];
		long flag=3;
		for(int i=k-1;i>=0;--i) {
			byteseq[i]=nt4_ascii[(int)(seq&flag)];
			seq>>>=2;
		}
		return new DNASequence(byteseq);
	}
	
	public static void decompressSeq(byte[] byteSeq,byte[] binSeq,int len) {
		final byte[] nt4=NT4_ASCII_TABLE;
		for(int i=0;true;++i) {
			byte base=binSeq[i];
			for(int j=0;j<4;++j) {
				int idx=i*4+j;
				byteSeq[idx]=nt4[(base>>>(6-j*2))&3];
				if(idx==len-1) {
					return;
				}
			}
		}
	}
	
	public static void decompressSeq(byte[] byteSeq,long[] binSeq,int len) {
		final byte[] nt4=NT4_ASCII_TABLE;
		for(int i=0;true;++i) {
			long bases=binSeq[i];
			for(int j=0;j<32;++j) {
				int idx=i*32+j;
				byteSeq[idx]=nt4[(int)((bases>>>(62-j*2))&3)];
				if(idx>=len-1) {
					return;
				}
			}
		}
	}
	
	public static void decompressSeq(byte[] byteSeq,ByteBuffer bb,int len) {
		final byte[] nt4=NT4_ASCII_TABLE;
		for(int i=0;true;++i) {
			byte base=bb.get();
			for(int j=0;j<4;++j) {
				int idx=i*4+j;
				byteSeq[idx]=nt4[(base>>>(6-j*2))&3];
				if(idx>=len-1) {
					return;
				}
			}
		}
	}
	
	public static long[] compressSeq(DNASequence dnaSeq) {
		return compressSeq(dnaSeq.getSequence(),dnaSeq.getLength());
	}
	
	public static long[] compressSeq(byte[] byteSeq,int len) {
		final byte[] nt4_bin=DNASequence.NT4_BIN_TABLE;
		
		int binSeqLen=((len-1)/32)+1;
		long[] binSeq=new long[binSeqLen];
		outer:
		for(int i=0;i<binSeqLen;++i) {
			for(int j=0;j<32;++j) {
				int idx=i*32+j;
				if(idx>=len) {
					binSeq[i]<<=((32-j)*2);
					break outer;
				}
				int base=nt4_bin[byteSeq[i*32+j]];
				if(base==-1) {               // base N
					base=DNASequence.getRandomBase();
				}
				binSeq[i]<<=2;
				binSeq[i]|=base;
			}
		}
		return binSeq;
	}
	
	public static int compressBarcode(byte[] byteSeq,int start,int len) throws IllegalBarcodeException{
		final byte[] nt4_bin=DNASequence.NT4_BIN_TABLE;
		final int end=start+len;
		
		int res=0;
		for(int i=start;i<end;++i) {
			int base=nt4_bin[byteSeq[i]];
			if(base==-1) { 
				throw new IllegalBarcodeException("N in barcode!");
			}
			res<<=2;
			res|=base;
		}
		return res;	
	}
	
	public static boolean isSimilar(long kmer1,long kmer2) {
		if(kmer1==kmer2) {
			return true;
		}
		long xor=kmer1^kmer2;
		int bitCnt=Long.bitCount(xor);
		if(bitCnt==1) {
			return true;
		}else if(bitCnt==2) {
			return ((xor&(xor<<1))&0xAAAAAAAAAAAAAAAAL)!=0L;
		}else {
			return false;
		}
	}

}

final class RefSequence extends DNASequence{
	protected String chrom;
	protected int start;
	protected int end;
	
	private RefSequence(String chrom,int start,int end,byte[] seq){
		super(seq);
		this.chrom=chrom;
		this.start=start;
		this.end=end;
	}
	
	static RefSequence createRefSeq(String chrom,int start,int end,byte[] seq) {
		if(!chrom.contains("chr")) {
			chrom="chr"+chrom;
		}
		RefSequence ref=new RefSequence(chrom,start,end,seq);
		return ref;
	}

	private void convertToUpperCase() {
		int len=sequenceLen;
		for(int i=0;i<len;++i) {
			if(sequence[i]>96) {
				sequence[i]=(byte)(sequence[i]-32);
			}
		}
	}
	
	int getStart() {
		return this.start;
	}
	
	int getEnd() {
		return this.end;
	}
	
	String getChrom() {
		return this.chrom;
	}

}

final class Read extends DNASequence{
	public final static int Read_Length=200; 
	
	protected byte[] name;
	protected byte[] quality;
	protected int nameLen;
	
	//private byte[] compressed;
	
	public Read() {
		super(Read_Length);
		this.name=new byte[Read_Length];
		this.quality=new byte[Read_Length];
		this.sequenceLen=0;
		//this.compressed=new byte[Read_Length/8+1];
	}

	public void copyFromByte(byte[] src,int start,int seqStart,int quaStart,int nameLen,int seqLen) {
		System.arraycopy(src, start, this.name, 0, nameLen);
		System.arraycopy(src, seqStart, this.sequence, 0, seqLen);
		System.arraycopy(src, quaStart, this.quality, 0, seqLen);
		this.nameLen=nameLen;
		this.sequenceLen=seqLen;
	}
	
	byte getQuaAt(int pos0) {
		return quality[pos0];
	}
	
	byte[] getName() {
		return this.name;
	}
	
	byte[] getQua() {
		return this.quality;
	}
	
	public String toString() {
		StringBuilder str=new StringBuilder(500);
		for(int i=0;i<nameLen;++i) {
			str.append((char)name[i]);
		}
		str.append("\n");
		for(int i=0;i<sequenceLen;++i) {
			str.append((char)sequence[i]);
		}
		str.append("\n+\n");
		for(int i=0;i<sequenceLen;++i) {
			str.append((char)quality[i]);
		}
		str.append("\n");
		return str.toString();
	}
	
	public void clean() {
		this.name=null;
		this.quality=null;
		this.sequence=null;
	}
	
	public boolean isNull() {
		return this.sequence==null;
	}
	
	/*public byte[] compressRead(int start,int len) {
		final byte[] sequence=this.sequence;
		final byte[] nt4_bin=NT4_BIN_TABLE;
		final byte[] res=this.compressed;
		byte tmp=0;
		int cycle=len/4;
		int remainder=len%4;
		for(int i=0;i<cycle;++i) {
			byte base=0;
			int pos=start+i*4;
			tmp=nt4_bin[sequence[pos]];
		    base|=tmp<<6;
		    tmp=nt4_bin[sequence[pos+1]];
		    base|=tmp<<4;
		    tmp=nt4_bin[sequence[pos+2]];
		    base|=tmp<<2;
		    tmp=nt4_bin[sequence[pos+3]];
		    base|=tmp;
		    res[i]=base;
		}
		return res;
	}*/
	
}