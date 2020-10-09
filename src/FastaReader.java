import java.io.*;
import java.util.*;

import java.nio.*;
import java.nio.channels.*;
import java.nio.file.Path;

public class FastaReader {
	private final Path fastaPath,faiPath;
	private final Map<String,FastaIndex> fastaIndexMap;
	private final ArrayList<String> chromList;
	private final FileInputStream fastaIn;
	private final FileChannel fastaFc;
	
	public FastaReader(Path fastaPath,Path faiPath) throws IllegalPositionException, IllegalBioFileException, IOException {
		this.fastaPath=fastaPath;
		this.faiPath=faiPath;
		this.fastaIndexMap=new HashMap<String,FastaIndex>();
		this.chromList=new ArrayList<String>();
		this.fastaIn=new FileInputStream(fastaPath.toFile());
		this.fastaFc=fastaIn.getChannel();
		loadIndexFile();
	}
	
	private void loadIndexFile() throws IOException,IllegalPositionException,IllegalBioFileException{
		BufferedReader in=null;
		try {
			in=new BufferedReader(new FileReader(this.faiPath.toFile()));
			String line;
			while((line=in.readLine())!=null) {
				String[] lineList=line.split("\t");
				if(lineList.length != 5) {
					throw new IllegalBioFileException("Invalid Fasta Index File! :"+this.faiPath);
				}
				if(!lineList[0].contains("chr")) {
					lineList[0]="chr"+lineList[0];
				}
				this.fastaIndexMap.put(lineList[0], new FastaIndex(lineList[0],Integer.parseInt(lineList[1]),Long.parseLong(lineList[2]),Integer.parseInt(lineList[3]),Integer.parseInt(lineList[4])));
			    this.chromList.add(lineList[0]);
			}
		}catch(NumberFormatException e) {
			throw new IllegalPositionException("Invalid Fasta Index File! :"+this.faiPath); 
		}catch(FileNotFoundException e){
			throw new FileNotFoundException("Fasta Index File Is Missing! :"+this.faiPath);
		}finally {
			try {
				if(in!=null) {
					in.close();
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private long convert2BytePos(String chrom,int chromPos) throws IllegalPositionException{
		FastaIndex fai=this.fastaIndexMap.get(chrom);
		int chromLength=fai.getChromLen();
		long startPosition=fai.getStartPos();
		int lineLength=fai.getLineLen();
		int lineByte=fai.getLineByte();
		if(chromLength<chromPos) {
			throw new IllegalPositionException("Invalid Chromosome Position :"+chromPos);
		}
		int row=chromPos/lineLength;
		int column=chromPos%lineLength;
		return startPosition+row*lineByte+column;
	}
	
	public RefSequence getSequence(String chrom)throws IOException,IllegalPositionException,IllegalBioFileException{
		return getSequence(chrom,0,fastaIndexMap.get(chrom).getChromLen());
	}
	
	public RefSequence getSequence(String chrom,int start,int end) throws IOException,IllegalPositionException, IllegalBioFileException{
		int chromLength=this.fastaIndexMap.get(chrom).getChromLen();
		if(end<=0 ||start>end||start>chromLength) {
			throw new IllegalPositionException("Invalid Chromosome Position : Start:"+start+"End:"+end);
		}
		
		start=start<0?0:start;
		end=end>chromLength?chromLength:end;

		long startByte=convert2BytePos(chrom,start);
		long endByte=convert2BytePos(chrom,end);


		ByteBuffer bb=ByteBuffer.allocate((int)(endByte-startByte+2));
		byte[] seq=new byte[end-start];
		
		synchronized(this) {
	    	this.fastaFc.read(bb, startByte-1);
		}
		
        bb.flip();
		compact(seq,bb);
		bb=null;

		return RefSequence.createRefSeq(chrom,start,end,seq);
	}
	
	public void close() {
		try {
			if(this.fastaIn!=null) {
	    		this.fastaIn.close();
			}
		}catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	private void compact(byte[] seq,ByteBuffer buffer) {
		int len=seq.length;
		for(int i=0;i<len;++i) {
			byte base=buffer.get();
			if(base<65) {
			    base=buffer.get();
			}
			if(base>96) {
				base=(byte)(base-32);
			}
			seq[i]=base;
		}
	}
	
	private void removeEnter(byte[] seqByte,byte[] seqWithoutEnter,int lineLen,int lineByte) {
		int seqByteLen=seqByte.length;
		if(seqByte[seqByteLen-1]<=64) {
			--seqByteLen;
		}
		int oldcnt=0;
		int newcnt=0;
		if(seqByteLen<=50) {
			while(seqByte[oldcnt]>64) {
		    	++oldcnt;
		    	if(oldcnt>=seqByteLen) {
		    		break;
		    	}
	    	}
		}else {
	    	while(seqByte[oldcnt]>64) {
		    	++oldcnt;
	    	}
		}
		System.arraycopy(seqByte, 0, seqWithoutEnter, 0, oldcnt);
		newcnt=oldcnt;
		oldcnt=oldcnt+lineByte-lineLen;
		int end=seqByteLen-lineByte;
		while(oldcnt<end) {
			System.arraycopy(seqByte, oldcnt, seqWithoutEnter, newcnt, lineLen);
			newcnt+=lineLen;
			oldcnt+=lineByte;
		}
		if(seqByteLen>oldcnt) {
	    	System.arraycopy(seqByte, oldcnt, seqWithoutEnter, newcnt, seqByteLen-oldcnt);
		}
	}
	
	public ArrayList<String> getChromList() {
		return this.chromList;
	}
  
	public ChromInfo getChromInfo() {
		int len=chromList.size();
		String[] chromNameList=new String[len];
		int[] chromLenList=new int[len];
		for(int i=0;i<len;++i) {
			String name=chromList.get(i);
			chromNameList[i]=name;
			chromLenList[i]=fastaIndexMap.get(name).getChromLen();
		}
		return new ChromInfo(chromNameList,chromLenList,fastaPath.toString());
	}
	
	private class FastaIndex{
		private final String chromName;
		private final int chromLength;
		private final long startPos;
		private final int lineLength;
		private final int lineByte;
		
		FastaIndex(String chrom,int len,long start,int linelength,int linebyte){
			this.chromName=chrom;
			this.chromLength=len;
			this.startPos=start;
			this.lineLength=linelength;
			this.lineByte=linebyte;
		}
		
		String getChromName() {
			return this.chromName;
		}
		
		int getChromLen() {
			return this.chromLength;
		}
		
		long getStartPos() {
			return this.startPos;
		}
		
		int getLineLen() {
			return this.lineLength;
		}
		
		int getLineByte() {
			return this.lineByte;
		}
	}
}

class ChromInfo implements Serializable,Iterable<Pair<String,Integer>>{

	private static final long serialVersionUID = 938264643474783901L;
	
	private final String fastaPath;
	private final String[] chromNameList;
	private final int[] chromLenList;
	
	ChromInfo(String[] chromNameList,int[] chromLenList,String fastaPath){
		this.chromLenList=chromLenList;
		this.chromNameList=chromNameList;
		this.fastaPath=fastaPath;
	}
	
	String getFastaPath() {
		return  this.fastaPath;
	} 

	private class Itr implements Iterator<Pair<String,Integer>>{

		private int cursor;
		
		@Override
		public boolean hasNext() {
			return cursor != chromNameList.length;
		}

		@Override
		public Pair<String, Integer> next() {
			int i = cursor;
            if (i >= chromNameList.length)
                throw new NoSuchElementException();
			String chromName=chromNameList[i];
			int chromLen=chromLenList[i];
			Pair<String, Integer> pair=new Pair<String, Integer>(chromName,chromLen);
			cursor=i+1;
			return pair;
		}
		
	}
	
	@Override
	public Iterator<Pair<String,Integer>> iterator() {
		return new Itr();
	}
	
}

class Pair<K,V> {
	private K key;
	private V value;
	
	Pair(K key,V value){
		this.key=key;
		this.value=value;
	}
	
	public K getKey() {
		return key;
	}
	
	public V getValue() {
		return value;
	}
}


