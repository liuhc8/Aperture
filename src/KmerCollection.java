import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class KmerCollection {
	public final static int LESS_REPETITIVE=1,FULL_REPETITIVE=2,MIX_REPETITIVE=3,UNIQUE=0; 
	
	//private int[][] kmerMatrix,codeMatrix,indexMatrix,scoreMatrix;
	
	private int[][] codeMatrix,scoreMatrix,tempKmerMatrix,outRangeMatrix;
	private short[][] kmerMatrix,indexMatrix,posMatrix;
	private int[] lenList;
	private int geneCodeLen,k,increment;
	private long flag;
	
	public KmerCollection(long estAmount,int k,int geneCodeLen){
		this.k=k;
		this.geneCodeLen=geneCodeLen;
		this.flag=(1L<<(2*k))-1;
		int n=1<<((k-16)*2);
		int increment=(int)(estAmount/10/n);
		this.increment=increment;
		this.tempKmerMatrix=new int[n][];
		this.kmerMatrix=new short[n][];
		this.indexMatrix=new short[n][];
		this.outRangeMatrix=new int[n][];
		this.codeMatrix=new int[n][];
		this.posMatrix=new short[n][];
		this.scoreMatrix=new int[n][];
		this.lenList=new int[n];
		for(int i=0;i<n;++i) {
			this.tempKmerMatrix[i]=new int[increment];
			this.codeMatrix[i]=new int[increment];
			this.posMatrix[i]=new short[increment];
		}
	}
	
	private KmerCollection(short[][] kmerMatrix,int[][] codeMatrix,short[][] posMatrix,short[][] indexMatrix,int[][] outRangeMatrix,int[][] scoreMatrix,int[] lenList,int k,int geneCodeLen) {
		this.kmerMatrix=kmerMatrix;
		this.codeMatrix=codeMatrix;
		this.posMatrix=posMatrix;
		this.indexMatrix=indexMatrix;
		this.outRangeMatrix=outRangeMatrix;
		this.scoreMatrix=scoreMatrix;
		this.lenList=lenList;
		this.k=k;
		this.flag=(1L<<(2*k))-1;
		this.geneCodeLen=geneCodeLen;
	}
	
	public void info() {
		int full=0;
		int uni=0;
		int slight=0;
		int all=0;
		int len=this.lenList.length;
		for(int i=0;i<len;++i) {
			int size=this.lenList[i];
			for(int j=0;j<size;++j){
				if(codeMatrix[i][j]==-1) {
					++full;
				}else if(GeneCodeGenerator.testGeneCode(codeMatrix[i][j],5)==0) {
					++uni;
				}else if(GeneCodeGenerator.testGeneCode(codeMatrix[i][j],10)<=0) {
					++slight;
				}
				++all;
			}
		}
		System.out.println("-1 :"+full);
		System.out.println("uni :"+uni);
		System.out.println("sli :"+slight);
		System.out.println("All :"+all);
		
		int sum=0;
		for(int i=0;i<len;++i) {
			sum+=this.lenList[i];
		}
		
		System.out.println("SUMl :"+sum);
		System.exit(0);
	}

	
	public void insert(long kmer,int code,int pos) {
		kmer&=this.flag;
		final int[][] tempKmerMatrix=this.tempKmerMatrix;
		final int[][] codeMatrix=this.codeMatrix;
		final short[][] posMatrix=this.posMatrix;
		final int increment=this.increment;
		final int[] lenList=this.lenList;
		
		int n=(int)(kmer>>>32);
		
		synchronized(this) {
			
		int idx=lenList[n];
		
	    if(idx>=tempKmerMatrix[n].length) {
	    	int[] array1=new int[idx+increment];
		   	System.arraycopy(tempKmerMatrix[n], 0, array1, 0, idx);
		   	tempKmerMatrix[n]=array1;
			
	    	int[] array2=new int[idx+increment];
	    	System.arraycopy(codeMatrix[n], 0, array2, 0, idx);
	    	codeMatrix[n]=array2;
			
	   		short[] array3=new short[idx+increment];
		   	System.arraycopy(posMatrix[n], 0, array3, 0, idx);
		   	posMatrix[n]=array3;
	    }
	    tempKmerMatrix[n][idx]=(int)kmer;
	    codeMatrix[n][idx]=code;
	    posMatrix[n][idx]=(short)pos;
	    ++lenList[n];
	    	
		}
		
	}
	
	public void compact(int nWorkers) throws InterruptedException {
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(nWorkers);
		for(int i=0;i<nWorkers;++i) {
			fixedThreadPool.execute(new CompactWorker(i,nWorkers));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
		this.tempKmerMatrix=null;
	}
	
	private class CompactWorker implements Runnable{
		private int no;
		private int nWorkers;
		
		CompactWorker(int no,int nWorkers){
			this.no=no;
			this.nWorkers=nWorkers;
		}
		
		@Override
		public void run() {
			final int[] lenList=KmerCollection.this.lenList;
			final int len=lenList.length;
			final int nWorkers=this.nWorkers;
			
			for(int i=this.no;i<len;i+=nWorkers) {
				qsort(i,0,lenList[i]-1);
				unique(i);
				buildIndex(i);
			}
			
		}
		
	}
	
	private void qsort(int n,int lo,int hi) {
		if(hi<=lo) {
			return;
		}
		int j=partition(n,lo,hi);
		qsort(n,lo,j-1);
		qsort(n,j+1,hi);
	}
	
	private int partition(int n,int lo,int hi) {
		final int[] tempKmerList=this.tempKmerMatrix[n];
		final int[] codeList=this.codeMatrix[n];
		final short[] posList=this.posMatrix[n];
		
		int i=lo;
		int j=hi+1;
		long v=tempKmerList[lo]&0xFFFFFFFFL;
		
		while(true) {
			while((tempKmerList[++i]&0xFFFFFFFFL)<v) {
				if(i==hi) {
					break;
				}
			}
			while(v<(tempKmerList[--j]&0xFFFFFFFFL)) {
				if(i==lo) {
					break;
				}
			}
			if(i>=j) {
				break;
			}
			swap(i,j,tempKmerList);
			swap(i,j,codeList);
			swap(i,j,posList);
		}
		swap(lo,j,tempKmerList);
		swap(lo,j,codeList);
		swap(lo,j,posList);
		return j;
	}
	
	private void swap(int i,int j,int[] list) {
		int tmp=list[i];
		list[i]=list[j];
		list[j]=tmp;
	}
	
	private void swap(int i,int j,short[] list) {
		short tmp=list[i];
		list[i]=list[j];
		list[j]=tmp;
	}
	
	private void unique(int i) {
		int[] tempKmerList=this.tempKmerMatrix[i];
		int[] codeList=this.codeMatrix[i];
		short[] posList=this.posMatrix[i];
		int len=this.lenList[i];
			
		int j=0;
		for(int k=1;k<len;++k) {
			if(tempKmerList[j]==tempKmerList[k]) {
				codeList[j]|=codeList[k];
				posList[j]=-1;
			}else {
				++j;
				tempKmerList[j]=tempKmerList[k];
				codeList[j]=codeList[k];
				posList[j]=posList[k];
			}
		}
		int newLen=j+1;
		this.lenList[i]=newLen;
			
		int[] newCodeList=new int[newLen];
		System.arraycopy(codeList, 0, newCodeList, 0, newLen);
		this.codeMatrix[i]=newCodeList;
			
		short[] newPosList=new short[newLen];
		System.arraycopy(posList, 0, newPosList, 0, newLen);
		this.posMatrix[i]=newPosList;
	}
	
	private void buildIndex(int i) {
		final int[] tempKmerList=this.tempKmerMatrix[i];
		final int len=this.lenList[i];
		
		short[] kmerList=new short[len];
		short[] indexList=new short[65536];
		
		int[] outRangeList=null;
		int outRangeLen=(len-1)>>>16;
		if(outRangeLen>0) {
			if(tempKmerList[(outRangeLen<<16)-1]>>>16==0x7FFF) {
				outRangeList=new int[outRangeLen];
			}else {
				outRangeList=new int[outRangeLen+1];
			}
		}

		for(int j=0;j<len;++j) {
			kmerList[j]=(short)tempKmerList[j];
		}
		
		for(int j=len-1,lastPtr=65536,lastJ=len;j>=0;--j) {		
			
			/*if(outRangeList!=null) {
		    	outRangeList[j>>>16]=(short)(tempKmerList[j]>>>16);
			}*/
			
			int ptr=tempKmerList[j]>>>16;
			indexList[ptr]=(short)j;
	
			for(int k=lastPtr-1;k>ptr;--k) {
				indexList[k]=(short)lastJ;
			}
			lastPtr=ptr;
			lastJ=j;
		}

		if(outRangeList!=null) {		
			for(int j=1;j<outRangeList.length;++j) {
				int high16=tempKmerList[(j<<16)-1]>>>16;
				outRangeList[j]=high16+1;
			}
			outRangeList[0]=0;
		}
		this.tempKmerMatrix[i]=null;
		this.kmerMatrix[i]=kmerList;
		this.indexMatrix[i]=indexList;
		this.outRangeMatrix[i]=outRangeList;
	}
	
	public void prepareScoreMatrix() {
		final int[][] scoreMatrix=this.scoreMatrix;
		final int[] lenList=this.lenList;
		int len=lenList.length;
		
		for(int i=0;i<len;++i) {
			scoreMatrix[i]=new int[(lenList[i]>>>4)+1];
		}
	}
	
	
	private int find0(long kmer) {
		kmer&=this.flag;
		int i=(int)(kmer>>>32);
		int key=(int)kmer;

		final short[] kmerList=this.kmerMatrix[i];
		final short[] indexList=this.indexMatrix[i];
		final int[] outRangeList=this.outRangeMatrix[i];
		final int len=this.lenList[i];
	
		int keyPtr=key>>>16;
		int lo=indexList[keyPtr]&0xFFFF;
		int hi=keyPtr==0xFFFF?len:indexList[keyPtr+1]&0xFFFF;
		//int hi=keyIdx==0x7FFF?len-1:(keyIdx==0xFFFF?((indexList[0]&0xFFFF)-1):(indexList[keyIdx+1]&0xFFFF)-1);
		
		if(outRangeList!=null) {
			int res=Arrays.binarySearch(outRangeList,keyPtr);
			res=res<0?(-res-2):res;
			lo|=(res<<16);
			
			if(hi!=len) {
		    	if(res+1<outRangeList.length) {
		     		if(keyPtr+1>=outRangeList[res+1]) {
		    			++res;
			    	}
		    	}
		    	hi|=(res<<16);
			}
		}
		
		/*synchronized(this.kmerMatrix) {
			if(lo>kmerList.length || hi>kmerList.length) {
			for(int j=0;j<outRangeList.length;++j) {
				System.out.print(":"+outRangeList[j]);
			}
			System.out.println();
			for(int j=0;j<outRangeList.length;++j) {
				System.out.print(":"+indexList[outRangeList[j]+0x8000]);
			}
			System.out.println();
			for(int j=0;j<outRangeList.length;++j) {
				System.out.print(":"+indexList[outRangeList[j]+0x8001]);
			}
			System.out.println();
			System.out.println(lo+" :::: "+hi);
			System.out.println("key: "+(short)key);
			for(int j=lo;j<hi;++j) {
				System.out.print(":"+kmerList[j]);
			}
			System.out.println();
			System.out.println((short)(key>>>16));
			System.out.println(kmerList.length);
			
			System.out.println();
			}
		}*/
		
		int x=unsignedShortBS(kmerList,lo,hi,(short)key);		
		
		return x;
	}
	
	private static int unsignedShortBS(short[] a, int fromIndex, int toIndex, short key) {
		int low = fromIndex;
		int high = toIndex - 1;
		int unsignedKey=key&0xFFFF;
		
		while (low <= high) {
			int mid = (low + high) >>> 1;
			int midVal = a[mid]&0xFFFF;
			
			if (midVal < unsignedKey)
				low = mid + 1;
			else if (midVal > unsignedKey)
				high = mid - 1;
			else
				return mid; // key found
			}
		return -(low + 1);  // key not found.
		}
	
	
	public void setScore(long kmer,int repeatCategory) {
		final int k=this.k;
		
		for(int i=0;i<(k*2);i+=2) {
			long masker1=3L<<i;
			long masker2=1L<<i;
			long masker3=2L<<i;
			
			setScore0(kmer^masker1,repeatCategory);
			setScore0(kmer^masker2,repeatCategory);
			setScore0(kmer^masker3,repeatCategory);
		}
	}
	
	private void setScore0(long kmer,int repeatCategory) {
		kmer&=this.flag;
		int i=(int)(kmer>>>32);

		final short[] posList=this.posMatrix[i];
		final int[] scoreList=this.scoreMatrix[i];
		
		int res=find0(kmer);

		if(res>=0) {
			if(posList[res]!=-1) {
		    	int move=(res&0xF)<<1;
		    	repeatCategory<<=move;
	    		int idx=res>>>4;
	    		synchronized(posList) {
		    		scoreList[idx]|=repeatCategory;
	    		}
			}
		}
	}
	
	public void update(long kmer,int geneCode) {
		kmer&=this.flag;
		int i=(int)(kmer>>>32);

		final int[] codeList=this.codeMatrix[i];
		final short[] posList=this.posMatrix[i];
		
		int res=find0(kmer);
		
		if(res>=0) {
			synchronized(codeList) {
		    	codeList[res] |= geneCode;
		    	posList[res]=-1;
			}
		}
	}

	public int findGeneCode(long kmer) {
		kmer&=this.flag;
		int i=(int)(kmer>>>32);
		final int[] codeList=this.codeMatrix[i];
	
		int res=find0(kmer);
		
		if(res>=0) {
			return codeList[res];
		}else {
			return 0;
		}
	}
	
	public long find(long kmer) {
		kmer&=this.flag;
		int i=(int)(kmer>>>32);
		final int[] codeList=this.codeMatrix[i];
		final short[] posList=this.posMatrix[i];
		final int[] scoreList=this.scoreMatrix[i];
	
		int res=find0(kmer);
		if(res>=0) {
			long code=codeList[res];
			short pos=posList[res];
			code<<=32;
			if(pos==-1) {
				return code|(pos&0xFFFF);
			}else {
				int move=(res&0xF)<<1;
				int score=(scoreList[res>>>4]>>>move)&0x3;
				return code|(score<<16)|(pos&0xFFFF);
			}
		}else {
			return 0xFFFFL;
		}
	}
	
	
	public void clean() {
		final int n=this.lenList.length;
		for(int i=0;i<n;++i) {
			this.kmerMatrix[i]=null;
			this.codeMatrix[i]=null;
			this.posMatrix[i]=null;
			this.indexMatrix[i]=null;
			this.scoreMatrix[i]=null;
			this.outRangeMatrix[i]=null;
		}
		this.kmerMatrix=null;
		this.codeMatrix=null;
		this.posMatrix=null;
		this.indexMatrix=null;
		this.scoreMatrix=null;
		this.outRangeMatrix=null;
		this.lenList=null;
	}
	
	public void saveToDisk(FileChannel fc) throws IOException {
		ByteBuffer bb=ByteBuffer.allocate(102400);
		final int n=this.kmerMatrix.length;
		final int[] lenList=this.lenList;
		
		long check=0;    //For Validation
		
		bb.put((byte)0xAB);
		bb.put((byte)0xBA);
		
		bb.putInt(this.k);
		bb.putInt(this.geneCodeLen);
		
		bb.putInt(n);
		
		for(int i=0;i<n;++i) {
			if(bb.remaining()<4) {
				bb.flip();
				while(bb.hasRemaining()) {
			    	fc.write(bb);
				}
				bb.clear();
			}
			bb.putInt(lenList[i]);
		}
		
		for(int i=0;i<n;++i) {
			short[] kmerList=this.kmerMatrix[i];
			int[] codeList=this.codeMatrix[i];
			short[] posList=this.posMatrix[i];
			short[] indexList=this.indexMatrix[i];
			int[] outRangeList=this.outRangeMatrix[i];
			int[] scoreList=this.scoreMatrix[i];
			
			for(int j=0;j<65536;++j) {
				if(bb.remaining()<10) {
					bb.flip();
					while(bb.hasRemaining()) {
				    	fc.write(bb);
					}
					bb.clear();
				}
				bb.putShort(indexList[j]);
				
				check=check*31+indexList[j];
			}
			
			if(outRangeList==null) {
				bb.putShort((short)-1);
			}else {
				bb.putShort((short)outRangeList.length);
				for(int j=0;j<outRangeList.length;++j) {
					if(bb.remaining()<10) {
						bb.flip();
						while(bb.hasRemaining()) {
					    	fc.write(bb);
						}
						bb.clear();
					}
					bb.putInt(outRangeList[j]);
					
					check=check*31+outRangeList[j];
				}
			}
			
			
			int len=lenList[i];
			for(int j=0;j<len;++j) {
				if(bb.remaining()<10) {
					bb.flip();
					while(bb.hasRemaining()) {
				    	fc.write(bb);
					}
					bb.clear();
				}
				bb.putShort(kmerList[j]);
				bb.putInt(codeList[j]);
				bb.putShort(posList[j]);
				
				check=(check*31+kmerList[j])*31+codeList[j]+posList[j];     //Calculate Hash
			}
			
			for(int j=0;j<(len>>>4)+1;++j) {
				if(bb.remaining()<4) {
					bb.flip();
					while(bb.hasRemaining()) {
				    	fc.write(bb);
					}
					bb.clear();
				}
				bb.putInt(scoreList[j]);
				
				check=check*31+scoreList[j]; 
			}
			
		}
		bb.flip();
		while(bb.hasRemaining()) {
	    	fc.write(bb);
		}
		bb.clear();
		
		bb.putLong(check);
		bb.flip();
		while(bb.hasRemaining()) {
	    	fc.write(bb);
		}
		
		fc.force(true);
	}
	
	public static KmerCollection loadFromDisk(FileChannel fc) throws IOException, IllegalBioFileException {
		long check=0;
		ByteBuffer bb=ByteBuffer.allocate(102400);
		fc.read(bb);
		bb.flip();
		
		byte magic1=bb.get();
		byte magic2=bb.get();
		if(magic1!=(byte)0xAB || magic2!=(byte)0xBA) {
			throw new IllegalBioFileException("Illegal Database File!");
		}
		
		int k=bb.getInt();
		int geneCodeLen=bb.getInt();
		
		int n=bb.getInt();
		int[] lenList=new int[n];
		for(int i=0;i<n;++i) {
			if(bb.remaining()<4) {
				bb.compact();
				fc.read(bb);
				bb.flip();
			}
			lenList[i]=bb.getInt();
		}
		
		short[][] kmerMatrix=new short[n][];
		int[][] codeMatrix=new int[n][];
		short[][] posMatrix=new short[n][];
		int[][] scoreMatrix=new int[n][];
		int[][] outRangeIndexMatrix=new int[n][];
		short[][] indexMatrix=new short[n][];
		
		for(int i=0;i<n;++i) {
			
			short[] indexList=new short[65536];
			for(int j=0;j<65536;++j) {
				if(bb.remaining()<10) {
					bb.compact();
					fc.read(bb);
					bb.flip();
				}
				indexList[j]=bb.getShort();
				
				check=check*31+indexList[j];
			}
			indexMatrix[i]=indexList;
			
			int outRangeListLen=bb.getShort();
			if(outRangeListLen!=-1) {
				int[] outRangeList=new int[outRangeListLen];		
				for(int j=0;j<outRangeListLen;++j) {
					if(bb.remaining()<10) {
						bb.compact();
						fc.read(bb);
						bb.flip();
					}
					outRangeList[j]=bb.getInt();	
					check=check*31+outRangeList[j];
				}
				outRangeIndexMatrix[i]=outRangeList;
			}


			int len=lenList[i];
			short[] kmerList=new short[len];
			int[] codeList=new int[len];
			short[] posList=new short[len];
			
			for(int j=0;j<len;++j) {
				if(bb.remaining()<10) {
					bb.compact();
					fc.read(bb);
					bb.flip();
				}
				kmerList[j]=bb.getShort();
				codeList[j]=bb.getInt();
				posList[j]=bb.getShort();
				
				check=(check*31+kmerList[j])*31+codeList[j]+posList[j];
			}
			kmerMatrix[i]=kmerList;
			codeMatrix[i]=codeList;
			posMatrix[i]=posList;
			
			int[] scoreList=new int[(len>>>4)+1];
			for(int j=0;j<(len>>>4)+1;++j) {
				if(bb.remaining()<4) {
					bb.compact();
					fc.read(bb);
					bb.flip();
				}
				scoreList[j]=bb.getInt();
				
				check=check*31+scoreList[j];
			}
			scoreMatrix[i]=scoreList;
		}
		
		if(bb.remaining()<8) {
			bb.compact();
			fc.read(bb);
			bb.flip();
		}
		
		if(check!=bb.getLong()) {
			throw new IllegalBioFileException("Illegal Database File!");
		}
		
		return new KmerCollection(kmerMatrix,codeMatrix,posMatrix,indexMatrix,outRangeIndexMatrix,scoreMatrix,lenList,k,geneCodeLen);
	}
	
	public int getK(){
		return this.k;
	}
	
	public int getGeneCodeLen() {
		return this.geneCodeLen;
	}
	
	public long size() {
		final int[] lenList=this.lenList;
		long size=0;
		for(int i=0;i<lenList.length;++i) {
			size+=lenList[i];
		}
		return size;
	}
}

class LongKmerCollection {
	private long[][] kmerMatrix;
	private int[][] codeMatrix,indexMatrix;
	private short[][] posMatrix;
	private int[] lenList;
	private int geneCodeLen,longk,increment,flag;
	
	public LongKmerCollection(long estAmount,int kmerLeftLen,int longk,int geneCodeLen){
		this.longk=longk;
		this.geneCodeLen=geneCodeLen;
		int n=1<<(kmerLeftLen*2);
		this.flag=n-1;
		int increment=(int)(estAmount/10/n);
		this.increment=increment;
		this.kmerMatrix=new long[n][];
		this.codeMatrix=new int[n][];
		this.posMatrix=new short[n][];
		this.indexMatrix=new int[n][];
		this.lenList=new int[n];
		for(int i=0;i<n;++i) {
			this.kmerMatrix[i]=new long[increment];
			this.codeMatrix[i]=new int[increment];
			this.posMatrix[i]=new short[increment];
			this.indexMatrix[i]=new int[256];
		}
	}
	
	private LongKmerCollection(long[][] kmerMatrix,int[][] codeMatrix,short[][] posMatrix,int[][] indexMatrix,int[] lenList,int longk,int flag,int geneCodeLen) {
		this.kmerMatrix=kmerMatrix;
		this.codeMatrix=codeMatrix;
		this.posMatrix=posMatrix;
		this.indexMatrix=indexMatrix;
		this.lenList=lenList;
		this.longk=longk;
		this.flag=flag;
		this.geneCodeLen=geneCodeLen;
	}
	
	public void info() {
		int full=0;
		int uni=0;
		int slight=0;
		int all=0;
		int len=this.lenList.length;
		for(int i=0;i<len;++i) {
			int size=this.lenList[i];
			for(int j=0;j<size;++j){
				if(codeMatrix[i][j]==-1) {
					++full;
				}else if(GeneCodeGenerator.testGeneCode(codeMatrix[i][j],5)==0) {
					++uni;
				}else if(GeneCodeGenerator.testGeneCode(codeMatrix[i][j],10)<=0) {
					++slight;
				}
				++all;
			}
		}
		System.out.println("-1 :"+full);
		System.out.println("uni :"+uni);
		System.out.println("sli :"+slight);
		System.out.println("All :"+all);
		
		int sum=0;
		for(int i=0;i<len;++i) {
			sum+=this.lenList[i];
		}
		
		System.out.println("SUMl :"+sum);
		System.exit(0);
	}

	
	
	
	public void insert(int kmerL,long kmerR,int code,int pos) {
		final long[][] kmerMatrix=this.kmerMatrix;
		final int[][] codeMatrix=this.codeMatrix;
		final short[][] posMatrix=this.posMatrix;
		final int increment=this.increment;
		final int[] lenList=this.lenList;
		
		int n=kmerL&this.flag;
		
		synchronized(this) {
			
		int idx=lenList[n];
		
	    if(idx>=kmerMatrix[n].length) {
	    	long[] array1=new long[idx+increment];
		   	System.arraycopy(kmerMatrix[n], 0, array1, 0, idx);
		   	kmerMatrix[n]=array1;
			
	    	int[] array2=new int[idx+increment];
	    	System.arraycopy(codeMatrix[n], 0, array2, 0, idx);
	    	codeMatrix[n]=array2;
			
	   		short[] array3=new short[idx+increment];
		   	System.arraycopy(posMatrix[n], 0, array3, 0, idx);
		   	posMatrix[n]=array3;
	    }
	    kmerMatrix[n][idx]=kmerR;
	    codeMatrix[n][idx]=code;
	    posMatrix[n][idx]=(short)pos;
	    ++lenList[n];
	    	
		}
		
	}
	
	public void compact(int nWorkers) throws InterruptedException {
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(nWorkers);
		for(int i=0;i<nWorkers;++i) {
			fixedThreadPool.execute(new CompactWorker(i,nWorkers));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
	}
	
	private class CompactWorker implements Runnable{
		private int no;
		private int nWorkers;
		
		CompactWorker(int no,int nWorkers){
			this.no=no;
			this.nWorkers=nWorkers;
		}
		
		@Override
		public void run() {
			final int[] lenList=LongKmerCollection.this.lenList;
			final int len=lenList.length;
			final int nWorkers=this.nWorkers;
			
			for(int i=this.no;i<len;i+=nWorkers) {
				qsort(i,0,lenList[i]-1);
				unique(i);
				buildIndex(i);
			}
			
		}
		
	}
	
	private void qsort(int n,int lo,int hi) {
		if(hi<=lo) {
			return;
		}
		int j=partition(n,lo,hi);
		qsort(n,lo,j-1);
		qsort(n,j+1,hi);
	}
	
	private int partition(int n,int lo,int hi) {
		final long[] kmerList=this.kmerMatrix[n];
		final int[] codeList=this.codeMatrix[n];
		final short[] posList=this.posMatrix[n];
		
		int i=lo;
		int j=hi+1;
		long v=kmerList[lo];
		
		while(true) {
			while(kmerList[++i]<v) {
				if(i==hi) {
					break;
				}
			}
			while(v<kmerList[--j]) {
				if(i==lo) {
					break;
				}
			}
			if(i>=j) {
				break;
			}
			swap(i,j,kmerList,codeList,posList);
		}
		swap(lo,j,kmerList,codeList,posList);
		return j;
	}
	
	private void swap(int i,int j,long[] kmerList,int[] codeList,short[] posList) {
		long tmp=kmerList[i];
		kmerList[i]=kmerList[j];
		kmerList[j]=tmp;
		
		int tmp2=codeList[i];
		codeList[i]=codeList[j];
		codeList[j]=tmp2;
		
		short tmp3=posList[i];
		posList[i]=posList[j];
		posList[j]=tmp3;
	}
	
	private void unique(int i) {
		long[] kmerList=this.kmerMatrix[i];
		int[] codeList=this.codeMatrix[i];
		short[] posList=this.posMatrix[i];
		int len=this.lenList[i];
			
		int j=0;
		for(int k=1;k<len;++k) {
			if(kmerList[j]==kmerList[k]) {
				codeList[j]|=codeList[k];
				posList[j]=-1;
			}else {
				++j;
				kmerList[j]=kmerList[k];
				codeList[j]=codeList[k];
				posList[j]=posList[k];
			}
		}
		int newLen=j+1;
		this.lenList[i]=newLen;
		
			
		long[] newKmerList=new long[newLen];
		System.arraycopy(kmerList, 0, newKmerList, 0, newLen);
		this.kmerMatrix[i]=newKmerList;
			
		int[] newCodeList=new int[newLen];
		System.arraycopy(codeList, 0, newCodeList, 0, newLen);
		this.codeMatrix[i]=newCodeList;
			
		short[] newPosList=new short[newLen];
		System.arraycopy(posList, 0, newPosList, 0, newLen);
		this.posMatrix[i]=newPosList;
	}
	
	private void buildIndex(int i) {
		final int[] indexList=this.indexMatrix[i];
		final long[] kmerList=this.kmerMatrix[i];
		final int len=this.lenList[i];
		
		for(int j=0;j<256;++j) {
			long key=((long)j)<<56;
			
			int lo=0;
			int hi=len-1;
			while(lo<=hi) {
				int mid=(lo+hi)>>>1;
				long midVal=kmerList[mid];
				if(midVal<key) {
					lo=mid+1;
				}else if(midVal>key) {
					hi=mid-1;
				}else {
					break;
				}
			}
			
			indexList[j]=lo;
		}
	}

	public void update(int kmerL,long kmerR,int geneCode) {

		int i=kmerL&this.flag;
		long key=kmerR;

		final long[] kmerList=this.kmerMatrix[i];
		final int[] codeList=this.codeMatrix[i];
		final short[] posList=this.posMatrix[i];
		final int[] indexList=this.indexMatrix[i];
		final int len=this.lenList[i];
	
		int keyIdx=(int)(key>>>56);
		int lo=indexList[keyIdx];
		int hi=keyIdx==127?len-1:(keyIdx==255?indexList[0]-1:indexList[keyIdx+1]-1);
		
		int res=0;
		while(lo<=hi) {
			int mid=(lo+hi)>>>1;
			long midVal=kmerList[mid];
			if(midVal<key) {
				lo=mid+1;
			}else if(midVal>key) {
				hi=mid-1;
			}else {
				res=mid;
				break;
			}
		}
		
		if(res!=0) {
			synchronized(kmerList) {
		    	codeList[res] |= geneCode;
		    	posList[res]=-1;
			}
		}
	}

	
	public long find(int kmerL,long kmerR) {
		int i=kmerL&this.flag;
		long key=kmerR;

		final long[] kmerList=this.kmerMatrix[i];
		final int[] codeList=this.codeMatrix[i];
		final short[] posList=this.posMatrix[i];
		final int[] indexList=this.indexMatrix[i];
		final int len=this.lenList[i];
	
		int keyIdx=(int)(key>>>56);
		int lo=indexList[keyIdx];
		int hi=keyIdx==127?len-1:(keyIdx==255?indexList[0]-1:indexList[keyIdx+1]-1);
		
		while(lo<=hi) {
			int mid=(lo+hi)>>>1;
			long midVal=kmerList[mid];
			if(midVal<key) {
				lo=mid+1;
			}else if(midVal>key) {
				hi=mid-1;
			}else {
				long res=codeList[mid];
				res<<=32;
				return res|(posList[mid]&0xFFFF);
			}
		}
		return 0xFFFFL;
	}
	
	public void clean() {
		final int n=this.lenList.length;
		for(int i=0;i<n;++i) {
			this.kmerMatrix[i]=null;
			this.codeMatrix[i]=null;
			this.posMatrix[i]=null;
			this.indexMatrix[i]=null;
		}
		this.kmerMatrix=null;
		this.codeMatrix=null;
		this.posMatrix=null;
		this.indexMatrix=null;
		this.lenList=null;
	}
	
	public void saveToDisk(FileChannel fc) throws IOException {
		ByteBuffer bb=ByteBuffer.allocate(102400);
		final int n=this.kmerMatrix.length;
		final int[] lenList=this.lenList;
		
		long check=0;    //For Validation
		
		bb.put((byte)0xAB);
		bb.put((byte)0xBA);
		
		bb.putInt(this.longk);
		bb.putInt(this.flag);
		bb.putInt(this.geneCodeLen);
		
		bb.putInt(n);
		
		for(int i=0;i<n;++i) {
			if(bb.remaining()<4) {
				bb.flip();
				while(bb.hasRemaining()) {
			    	fc.write(bb);
				}
				bb.clear();
			}
			bb.putInt(lenList[i]);
		}
		
		for(int i=0;i<n;++i) {
			long[] kmerList=this.kmerMatrix[i];
			int[] codeList=this.codeMatrix[i];
			short[] posList=this.posMatrix[i];
			int[] indexList=this.indexMatrix[i];
			
			for(int j=0;j<256;++j) {
				if(bb.remaining()<10) {
					bb.flip();
					while(bb.hasRemaining()) {
				    	fc.write(bb);
					}
					bb.clear();
				}
				bb.putInt(indexList[j]);
				
				check=check*31+indexList[j];
			}
			
			int len=lenList[i];
			for(int j=0;j<len;++j) {
				if(bb.remaining()<14) {
					bb.flip();
					while(bb.hasRemaining()) {
				    	fc.write(bb);
					}
					bb.clear();
				}
				bb.putLong(kmerList[j]);
				bb.putInt(codeList[j]);
				bb.putShort(posList[j]);
				
				check=(check*31+kmerList[j])*31+codeList[j]+posList[j];     //Calculate Hash
			}
			
		}
		bb.flip();
		while(bb.hasRemaining()) {
	    	fc.write(bb);
		}
		bb.clear();
		
		bb.putLong(check);
		bb.flip();
		while(bb.hasRemaining()) {
	    	fc.write(bb);
		}
		
		fc.force(true);
	}
	
	public static LongKmerCollection loadFromDisk(FileChannel fc) throws IOException, IllegalBioFileException {
		long check=0;
		ByteBuffer bb=ByteBuffer.allocate(102400);
		fc.read(bb);
		bb.flip();
		
		byte magic1=bb.get();
		byte magic2=bb.get();
		if(magic1!=(byte)0xAB || magic2!=(byte)0xBA) {
			throw new IllegalBioFileException("Illegal Database File!");
		}
		
		int longk=bb.getInt();
		int flag=bb.getInt();
		int geneCodeLen=bb.getInt();
		
		int n=bb.getInt();
		int[] lenList=new int[n];
		for(int i=0;i<n;++i) {
			if(bb.remaining()<4) {
				bb.compact();
				fc.read(bb);
				bb.flip();
			}
			lenList[i]=bb.getInt();
		}
		
		long[][] kmerMatrix=new long[n][];
		int[][] codeMatrix=new int[n][];
		short[][] posMatrix=new short[n][];
		int[][] indexMatrix=new int[n][];
		
		for(int i=0;i<n;++i) {
			
			int[] indexList=new int[256];
			for(int j=0;j<256;++j) {
				if(bb.remaining()<10) {
					bb.compact();
					fc.read(bb);
					bb.flip();
				}
				indexList[j]=bb.getInt();
				
				check=check*31+indexList[j];
			}
			indexMatrix[i]=indexList;
			
			int len=lenList[i];
			long[] kmerList=new long[len];
			int[] codeList=new int[len];
			short[] posList=new short[len];
			
			for(int j=0;j<len;++j) {
				if(bb.remaining()<14) {
					bb.compact();
					fc.read(bb);
					bb.flip();
				}
				kmerList[j]=bb.getLong();
				codeList[j]=bb.getInt();
				posList[j]=bb.getShort();
				
				check=(check*31+kmerList[j])*31+codeList[j]+posList[j];
			}
			kmerMatrix[i]=kmerList;
			codeMatrix[i]=codeList;
			posMatrix[i]=posList;

		}
		
		if(bb.remaining()<8) {
			bb.compact();
			fc.read(bb);
			bb.flip();
		}
		
		if(check!=bb.getLong()) {
			throw new IllegalBioFileException("Illegal Database File!");
		}
		
		return new LongKmerCollection(kmerMatrix,codeMatrix,posMatrix,indexMatrix,lenList,longk,flag,geneCodeLen);
	}
	
	public int getLongK(){
		return this.longk;
	}
	
	public int getGeneCodeLen() {
		return this.geneCodeLen;
	}
	
	public long size() {
		final int[] lenList=this.lenList;
		long size=0;
		for(int i=0;i<lenList.length;++i) {
			size+=lenList[i];
		}
		return size;
	}
}

class SimpleKmerCollection {
	
	public static long EST_KMER_AMOUNT=2500000000L;
	private int[][] kmerMatrix;
	private int[] lenList;
	private int k,increment;
	
	public SimpleKmerCollection(int k){
		this.k=k;
		int n=1<<((k-16)*2);
		int increment=(int)(EST_KMER_AMOUNT/10/n);
		this.increment=increment;
		this.kmerMatrix=new int[n][];
		this.lenList=new int[n];
		for(int i=0;i<n;++i) {
			this.kmerMatrix[i]=new int[increment];
		}
	}

	
	public void insert(long kmer) {
		final int[][] kmerMatrix=this.kmerMatrix;
		final int[] lenList=this.lenList;
		
		int n=(int)(kmer>>>32);
		int idx=lenList[n];
		if(idx>=kmerMatrix[n].length) {
			int[] array1=new int[idx+this.increment];
			System.arraycopy(kmerMatrix[n], 0, array1, 0, idx);
			kmerMatrix[n]=array1;

		}
		kmerMatrix[n][idx]=(int)kmer;
		++lenList[n];
	}
	
	public void compact(int nWorkers) throws InterruptedException {
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(nWorkers);
		for(int i=0;i<nWorkers;++i) {
			fixedThreadPool.execute(new CompactWorker(i,nWorkers));
		}
		
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
	}
	
	private class CompactWorker implements Runnable{
		private int no;
		private int nWorkers;
		
		CompactWorker(int no,int nWorkers){
			this.no=no;
			this.nWorkers=nWorkers;
		}
		
		@Override
		public void run() {
			final int len=SimpleKmerCollection.this.kmerMatrix.length;
			final int[] lenList=SimpleKmerCollection.this.lenList;
			final int nWorkers=this.nWorkers;
			
			for(int i=this.no;i<len;i+=nWorkers) {
				qsort(i,0,lenList[i]-1);
				unique(i);
			}
			
		}
		
	}
	
	private void qsort(int n,int lo,int hi) {
		if(hi<=lo) {
			return;
		}
		int j=partition(n,lo,hi);
		qsort(n,lo,j-1);
		qsort(n,j+1,hi);
	}
	
	private int partition(int n,int lo,int hi) {
		final int[] kmerList=this.kmerMatrix[n];
		
		int i=lo;
		int j=hi+1;
		int v=kmerList[lo];
		
		while(true) {
			while(kmerList[++i]<v) {
				if(i==hi) {
					break;
				}
			}
			while(v<kmerList[--j]) {
				if(i==lo) {
					break;
				}
			}
			if(i>=j) {
				break;
			}
			swap(i,j,kmerList);
		}
		swap(lo,j,kmerList);
		return j;
	}
	
	private void swap(int i,int j,int[] kmerList) {
		int tmp=kmerList[i];
		kmerList[i]=kmerList[j];
		kmerList[j]=tmp;
	}
	
	private void unique(int i) {
		
		final int[] kmerList=this.kmerMatrix[i];
		final int len=this.lenList[i];
			
		int j=-1;
		for(int k=1;k<len;++k) {
			if(kmerList[k]==kmerList[k-1]) {
				if(j==-1) {
					kmerList[++j]=kmerList[k];
				}else if(kmerList[k]!=kmerList[j]) {
					kmerList[++j]=kmerList[k];
					
				}
			}
		}
		
		int newLen=j+1;
		this.lenList[i]=newLen;
			
		int[] newKmerList=new int[newLen];
		System.arraycopy(kmerList, 0, newKmerList, 0, newLen);
		this.kmerMatrix[i]=newKmerList;
			
	}
	
	
	public boolean find(long akmer) {
		
		int n=(int)(akmer>>>32);
		int key=(int)akmer;

		final int[] kmerList=this.kmerMatrix[n];

		int lo=0;
		int hi=this.lenList[n]-1;
		while(lo<=hi) {
			int mid=(lo+hi)>>>1;
			int midVal=kmerList[mid];
			if(midVal<key) {
				lo=mid+1;
			}else if(midVal>key) {
				hi=mid-1;
			}else {
				return true;
			}
		}
		return false;
	}
	
	public void clean() {
		final int[][] kmerMatrix=this.kmerMatrix;
		final int n=this.lenList.length;
		for(int i=0;i<n;++i) {
			kmerMatrix[i]=null;
		}
		this.kmerMatrix=null;
		this.lenList=null;
	}
	
	public int getK(){
		return this.k;
	}
	
	public long size() {
		final int[] lenList=this.lenList;
		long size=0;
		for(int i=0;i<lenList.length;++i) {
			size+=lenList[i];
		}
		return size;
	}
}

class ReadKmerCollection{
	private long[] kmerList;
	private short[] posList;
	private long flag;
	
	private int[] save;
	private int saveIdx;
	
	ReadKmerCollection(int cnt){
		this.kmerList=new long[16384];
		this.posList=new short[16384];
		this.flag=0x3FFFL;
		
		this.save=new int[cnt];
		this.saveIdx=0;
	}
	
	public void insert(long kmer,short pos) {
		int n=(int)(kmer&this.flag);
		if(this.kmerList[n]==0L) {
			this.posList[n]=pos;
			this.kmerList[n]=kmer;
			this.save[this.saveIdx++]=n;
		}else {
			this.posList[n]=-1;
			this.kmerList[n]=-1;
		}
	}
	
    public int find(long key) {

    	int n=(int)(key&this.flag);
    	if(this.kmerList[n]==key) {
    		return this.posList[n];
    	}else {
    		return -1;
    	}

	}
    
    public void clean() {
    	final long[] kmerList=this.kmerList;
    	final short[] posList=this.posList;
    	final int[] save=this.save;
    	final int saveIdx=this.saveIdx;
    	
    	for(int i=0;i<saveIdx;++i) {
    		int n=save[i];
    		kmerList[n]=0;
    		posList[n]=0;
    	}
    	this.saveIdx=0;
    }
}
    
/*class ByteKmerCollection{
	private int[] cntList;
	
	private int[] save;
	private int saveIdx;
	
	ByteKmerCollection(){
		this.cntList=new int[256];
		
		this.save=new int[256];
		this.saveIdx=0;
	}
	
	public void insert(byte kmer,int cnt) {
		final int[] cntList=this.cntList;
		
		int n=(int)(kmer&0xFF);
		if(cntList[n]==0) {
			this.save[this.saveIdx++]=n;
		}
		cntList[n]=+cnt;
	}
	
    public byte getMostKmer() {
    	final int[] save=this.save;
    	final int saveIdx=this.saveIdx;
    	final int[] cntList=this.cntList;
    	
    	int mostKmer=0,maxCnt=0;
    	for(int i=0;i<saveIdx;++i) {
    		int n=save[i];
    		if(cntList[n]>maxCnt) {
    			maxCnt=cntList[n];
    			mostKmer=n;
    		}
    	}
    	return (byte)mostKmer;

	}
    
    public void clean() {
    	final int[] cntList=this.cntList;
    	final int[] save=this.save;
    	final int saveIdx=this.saveIdx;
    	
    	for(int i=0;i<saveIdx;++i) {
    		int n=save[i];
    		cntList[n]=0;
    	}
    	this.saveIdx=0;
    }
}*/


