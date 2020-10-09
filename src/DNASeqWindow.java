import java.util.NoSuchElementException;

public class DNASeqWindow{
    	
		private int k;
		private long kmerMask;
		private int leftZero;
		
		private int start;
		private long nextKmer;
		private int index;
		private int jump;
		private int length;
		
		private byte[] seq;
		
		public DNASeqWindow(int k){
			this.k=k;
			this.kmerMask=(1L<<(2*k))-1;
			this.leftZero=(32-k)*2;
		}
		
		public void set(byte[] seq,int start){
			set(seq,start,1,seq.length);
		}
		
		public void set(DNASequence dna,int start){
			byte[] seq=dna.getSequence();
			set(seq,start,1,seq.length);
		}
		
		public void set(DNASequence dna,int start,int jump){
			byte[] seq=dna.getSequence();
			set(seq,start,jump,seq.length);
		}
		
		
		public void set(byte[] seq,int start,int jump) {
			set(seq,start,jump,seq.length);
		}
		
		public void set(byte[] seq,int start,int jump,int len) {
			this.start=start;
			this.nextKmer=0;
			this.index=0;
			this.seq=seq;
			this.jump=jump;
			this.length=len;
		}
		
		public long getLastKmer() {
			return compressKmer(this.length-this.k,this.k);
		}
		
		private long compressKmer(int start,int k) {
			final byte[] seq=this.seq;
			final byte[] nt4_bin=DNASequence.NT4_BIN_TABLE;
			long res=0;
			int end=start+k;
			for(int i=start;i<end;++i) {
				int base=nt4_bin[seq[i]];
				if(base==-1) {               // base N
					base=DNASequence.getRandomBase();
				}
				res<<=2;
				res|=base;
			}
			return res;
		}

		public boolean hasNextKmer() {
			return this.index+this.start+this.k<= this.length;
		}

		public long nextKmer() {
			final byte[] nt4_bin=DNASequence.NT4_BIN_TABLE;
			final int jump=this.jump;
			final byte[] seq=this.seq;
			final int k=this.k;
			final int start=this.start;
			final int index=this.index;
			long nextKmer=this.nextKmer;

			if(hasNextKmer()) {
				if(index==0) {
					nextKmer=compressKmer(start,k);
				}else {
					for(int i=1-jump;i<=0;++i) {
						nextKmer<<=2;
						byte base=nt4_bin[seq[index+start+i+k-1]];
						if(base==-1) {
							nextKmer|=DNASequence.getRandomBase();
						}else {
					    	nextKmer|=base;
						}
					}
					nextKmer&=this.kmerMask;
				}
				this.index+=jump;
				this.nextKmer=nextKmer;
				return nextKmer;
			}else {
				throw new NoSuchElementException();
			}
		}
		
		public long showThisAsRevComp() {
			long kmerRev=~reverseBase(this.nextKmer);
			return (kmerRev>>>this.leftZero);
		}
		
		private long reverseBase(long i) {
	        i = (i & 0x3333333333333333L) << 2 | (i >>> 2) & 0x3333333333333333L;
	        i = (i & 0x0f0f0f0f0f0f0f0fL) << 4 | (i >>> 4) & 0x0f0f0f0f0f0f0f0fL;
	        i = (i & 0x00ff00ff00ff00ffL) << 8 | (i >>> 8) & 0x00ff00ff00ff00ffL;
	        i = (i << 48) | ((i & 0xffff0000L) << 16) |
	            ((i >>> 16) & 0xffff0000L) | (i >>> 48);
	        return i;
	    }
		
	}