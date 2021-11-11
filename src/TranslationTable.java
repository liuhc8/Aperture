import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;


public class TranslationTable implements Serializable{
	
	private static final long serialVersionUID = 2881190773180521210L;
	private String[] chrom;
	private int[] start,end,geneCode;
	private long[][] ref;
	private int k,longk,geneCodeLen,idx;
	
	private transient HashMap<Integer,Integer> segmentMap;
	private transient HashMap<String,Scope> scopeMap;
	
	TranslationTable(int size,int k,int longk,int geneCodeLen){
		this.idx=0;
		this.k=k;
		this.longk=longk;
		this.geneCodeLen=geneCodeLen;
		this.chrom=new String[size];
		this.start=new int[size];
		this.end=new int[size];
		this.geneCode=new int[size];
		this.ref=new long[size][];
	}
	
	void addSeperator(String chrom,int start,int end,int geneCode,long[] ref){
		final int idx=this.idx;
		if(idx==this.chrom.length) {
			String[] newChrom=new String[idx*2];
			System.arraycopy(this.chrom, 0, newChrom, 0, idx);
			this.chrom=newChrom;
			
			int[] newStart=new int[idx*2];
			System.arraycopy(this.start, 0, newStart, 0, idx);
			this.start=newStart;
			
			int[] newEnd=new int[idx*2];
			System.arraycopy(this.end, 0, newEnd, 0, idx);
			this.end=newEnd;
			
			int[] newGeneCode=new int[idx*2];
			System.arraycopy(this.geneCode, 0, newGeneCode, 0, idx);
			this.geneCode=newGeneCode;
			
			long[][] newRef=new long[idx*2][];
			System.arraycopy(this.ref, 0, newRef, 0, idx);
			this.ref=newRef;
			
		}
		this.chrom[idx]=chrom;
		this.start[idx]=start;
		this.end[idx]=end;
		this.geneCode[idx]=geneCode;
		this.ref[idx]=ref;
		
		++this.idx;
	}

	void compact() {
		final int idx=this.idx;
		
		String[] newChrom=new String[idx];
		System.arraycopy(this.chrom, 0, newChrom, 0, idx);
		this.chrom=newChrom;
		
		int[] newStart=new int[idx];
		System.arraycopy(this.start, 0, newStart, 0, idx);
		this.start=newStart;
		
		int[] newEnd=new int[idx];
		System.arraycopy(this.end, 0, newEnd, 0, idx);
		this.end=newEnd;
		
		int[] newGenecode=new int[idx];
		System.arraycopy(this.geneCode, 0, newGenecode, 0, idx);
		this.geneCode=newGenecode;
		
		long[][] newRef=new long[idx][];
		System.arraycopy(this.ref, 0, newRef, 0, idx);
		this.ref=newRef;
		
		HashMap<String,Scope> scopeMap=new HashMap<String,Scope>();
		for(int i=0;i<idx;++i){
			Scope scope=null;
			if((scope=scopeMap.get(newChrom[i]))==null) {
				scopeMap.put(newChrom[i], new Scope(i,i));
			}else {
				scope.end=i;
			}
		}
		this.scopeMap=scopeMap;
		
	}

	void cleanRef() {
		final int idx=this.idx;
		final long[][] ref=this.ref;
		
		for(int i=0;i<idx;++i) {
			ref[i]=null;
		}
	}

	int getK() {
		return this.k;
	}
	
	int getLongK() {
		return this.longk;
	}
	
	int getGeneCodeLen() {
		return this.geneCodeLen;
	}
	
	int size() {
		return this.idx;
	}
	

	Segment getSegment(int i) {
		return new Segment(chrom[i],start[i],end[i],geneCode[i],ref[i]);
	}
	
	int getGeneCode(int i) {
		return this.geneCode[i];
	}
	

	Segment findSegment(String chr,int pos) {
		Scope scope=this.scopeMap.get(chr);
		if(scope==null) {
			return null;
		}
		int index=Arrays.binarySearch(this.start, scope.start, scope.end+1, pos);
		if(index>=0) {
			return new Segment(chrom[index],start[index],end[index],geneCode[index],ref[index]);
		}else {
			index=-index-2;
			return new Segment(chrom[index],start[index],end[index],geneCode[index],ref[index]);
		}
	}
	
	Segment findSegment(int code){
		if(this.segmentMap==null) {
     		int len=this.idx;
     		HashMap<Integer,Integer> map=new HashMap<Integer,Integer>();
	    	for(int i=0;i<len;++i) {
	    		map.put(this.geneCode[i], i);
	    	}
	    	this.segmentMap=map;
		}
		Integer index=this.segmentMap.get(code);
		if(index==null) {
			return null;
		}else {
			return new Segment(chrom[index],start[index],end[index],geneCode[index],ref[index]);
		}
	}
	
	int getIndexNo(int code){
		if(this.segmentMap==null) {
     		int len=this.idx;
     		HashMap<Integer,Integer> map=new HashMap<Integer,Integer>();
	    	for(int i=0;i<len;++i) {
	    		map.put(this.geneCode[i], i);
	    	}
	    	this.segmentMap=map;
		}
		Integer index=this.segmentMap.get(code);
		if(index==null) {
			return -1;
		}else {
			return index;
		}
	}
	
	public String toString() {
		StringBuilder strb=new StringBuilder(100);
		String tab="\t";
		for(int i=0;i<chrom.length;++i) {
			strb.append(chrom[i]).append(tab).append(start[i]).append(tab).append(end[i]).append(tab).append(geneCode[i]).append("\n");
		}
		return strb.toString();
	}
	
	private class Scope{
		int start,end;
		
		Scope(int start,int end){
			this.start=start;
			this.end=end;
		}
	}
}

class Segment{
	public String chrom;
	public int start,end,code;
	public long[] ref;
	
	public Segment(String chrom,int start,int end,int code,long[] ref) {
		this.chrom=chrom;
		this.start=start;
		this.end=end;
		this.code=code;
		this.ref=ref;
	}
	
	public boolean equals(Object anObject)  {
		if (this == anObject) {
            return true;
        }
        if (anObject instanceof Segment) {
        	Segment anotherSeg = (Segment)anObject;
            if(anotherSeg.code==this.code) {
            	return true;
            }
        }
        return false;
	}
}