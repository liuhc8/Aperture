import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

class BpIndexList{
	private int[] codeList,lenList;
	private long[] initStartPosList,mergedStartPosList;
	private int size;
	private long tmpPos;
	
	public BpIndexList(int n) {
		this.codeList=new int[n];
		this.initStartPosList=new long[n];
		this.mergedStartPosList=new long[n];
		this.lenList=new int[n];
		this.size=0;
		this.tmpPos=0;
		Arrays.fill(mergedStartPosList,-1L);
	}
	
	public void addCodeAndLen(int code,int len)  {
		codeList[size]=code;
		initStartPosList[size]=tmpPos;
		lenList[size]=len;
		tmpPos+=len;
		++size;
	}
	
	public void addCodeAndStartPos(int code,long startPos){
		codeList[size]=code;
		initStartPosList[size]=startPos;
		if(size>0) {
			lenList[size-1]=(int)(startPos-initStartPosList[size-1]);
		}
		++size;
	}
	
	public void setLastLen(long endPos) {
		lenList[size-1]=(int)(endPos-initStartPosList[size-1]);
	}
	
	public GetAndSetItr getGetAndSetItr() {
		return new GetAndSetItr();
	}
	
	public class GetAndSetItr{
		private int cursor;
		
		GetAndSetItr(){
			cursor=-1;
		}
		
		public boolean hasNext() {
			return cursor+1<size;
		}
		
		public void moveToNext() {
			++cursor;
		}
		
		public int getCurrentCode() {
			return codeList[cursor];
		}
		
		public int getCurrentLen() {
			return lenList[cursor];
		}
		
		public long getCurrentInitStart() {
			return initStartPosList[cursor];
		}
		
		public long getCurrentMergedStart() {
			return mergedStartPosList[cursor];
		}
		
		public void setCurrentMergedStart(long mergedStart) {
			if(mergedStartPosList[cursor]!=-1) {
				System.out.println(mergedStart+"????????"+mergedStartPosList[cursor]);
			}
			mergedStartPosList[cursor]=mergedStart;
		}
	}
	
	public SynchronizedGetItr getSynchronizedGetItr() {
		return new SynchronizedGetItr();
	}
	
	public class SynchronizedGetItr{
		private AtomicInteger cursor;
		
		SynchronizedGetItr(){
			cursor=new AtomicInteger(0);
		}
		
		public CodeStartLen getNext(CodeStartLen pack) {
			int idx=cursor.getAndIncrement();
			if(idx<size) {
				pack.code=codeList[idx];
				pack.start=initStartPosList[idx];
				pack.len=lenList[idx];
				return pack;
			}else {
				return null;
			}
		}
	}
	
	public void sortByInitStart() {
		qsortByInitStart(0,this.size-1);
	}
	
	private void qsortByInitStart(int lo,int hi) {
		if(hi<=lo) {
			return;
		}
		int j=partitionByInitStart(lo,hi);
		qsortByInitStart(lo,j-1);
		qsortByInitStart(j+1,hi);
	}
	
	private int partitionByInitStart(int lo,int hi) {
		final int[] codeList=this.codeList;
		final long[] initStartPosList=this.initStartPosList;
		final long[] mergedStartPosList=this.mergedStartPosList;
		final int[] lenList=this.lenList;
		
		int i=lo;
		int j=hi+1;
		long v=initStartPosList[lo];
		
		while(true) {
			while(initStartPosList[++i]<v) {
				if(i==hi) {
					break;
				}
			}
			while(v<initStartPosList[--j]) {
				if(i==lo) {
					break;
				}
			}
			if(i>=j) {
				break;
			}
			swap(i,j,initStartPosList);
			swap(i,j,mergedStartPosList);
			swap(i,j,codeList);
			swap(i,j,lenList);
		}
		swap(lo,j,initStartPosList);
		swap(lo,j,mergedStartPosList);
		swap(lo,j,codeList);
		swap(lo,j,lenList);
		return j;
	}
	
	public void sortByCode() {
		qsortByCode(0,this.size-1);
	}
	
	private void qsortByCode(int lo,int hi) {
		if(hi<=lo) {
			return;
		}
		int j=partitionByCode(lo,hi);
		qsortByCode(lo,j-1);
		qsortByCode(j+1,hi);
	}
	
	private int partitionByCode(int lo,int hi) {
		final int[] codeList=this.codeList;
		final long[] initStartPosList=this.initStartPosList;
		final int[] lenList=this.lenList;
		
		int i=lo;
		int j=hi+1;
		int v=codeList[lo];
		
		while(true) {
			while(codeList[++i]<v) {
				if(i==hi) {
					break;
				}
			}
			while(v<codeList[--j]) {
				if(i==lo) {
					break;
				}
			}
			if(i>=j) {
				break;
			}
			swap(i,j,codeList);
			swap(i,j,initStartPosList);
			swap(i,j,lenList);
		}
		swap(lo,j,codeList);
		swap(lo,j,initStartPosList);
		swap(lo,j,lenList);
		return j;
	}
	
	private void swap(int i,int j,int[] list) {
		int tmp=list[i];
		list[i]=list[j];
		list[j]=tmp;
	}
	
	private void swap(int i,int j,long[] list) {
		long tmp=list[i];
		list[i]=list[j];
		list[j]=tmp;
	}
	
	public void clean() {
		this.codeList=null;
		this.initStartPosList=null;
		this.mergedStartPosList=null;
		this.lenList=null;
	}
}
