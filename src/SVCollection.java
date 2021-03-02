import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class SVCollection{
	private Object[] supportMatrix;
	final private TranslationTable trTable;
	final private int threads,len;
	
	public SVCollection(TranslationTable trTable,int threads){
		this.trTable=trTable;
		this.threads=threads;
		int len=trTable.size();
		this.supportMatrix=new Object[len];
		this.len=len;
		
		for(int i=0;i<len;++i) {
			this.supportMatrix[i]=new ArrayList<SVSupport>(30);
		}
	}
	
	public void adjustAndAdd(SVSupport sv) {
		int leftI=this.trTable.getIndexNo(sv.getLeftCode());
		int rightI=this.trTable.getIndexNo(sv.getRightCode());
		int idx=leftI;
		if(leftI>rightI && rightI!=-1) {
			sv.setLeftHighFlag();
			sv.flip();
			idx=rightI;
		}else if(leftI==rightI && sv.getLeftPos()>sv.getRightPos()) {
			sv.setLeftHighFlag();
			sv.flip();
		}else {
			sv.setLeftSmallFlag();
		}
		ArrayList<SVSupport> list=(ArrayList<SVSupport>) this.supportMatrix[idx];
		synchronized(list) {
			list.add(sv);
		}
	}

	
	public void merge() throws InterruptedException {		
		final Object[] supportMatrix=this.supportMatrix;
		final int threads=this.threads;
		
		AtomicInteger taskIdx=new AtomicInteger(0);
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads); 
		
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new MergeWorker(supportMatrix,taskIdx));
		}
		
		fixedThreadPool.shutdown();
		while(!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) {}

	}
	
	public Object[] transAndFilter(boolean hasDuoBarcode,int minIndelLen,double maxBpScore) throws InterruptedException {
		final TranslationTable trTable=this.trTable;
		final Object[] supportMatrix=this.supportMatrix;
		final int len=this.len;
		final int threads=this.threads;
		
		int[] accumulatedCntList=new int[len];
		for(int i=1;i<len;++i) {
			ArrayList<SVSupport> list=(ArrayList<SVSupport>)supportMatrix[i-1];
			accumulatedCntList[i]=list.size()+accumulatedCntList[i-1];
		}
		
		Object[] vcfRecordMatrix=new Object[len];
		for(int i=0;i<len;++i) {
			vcfRecordMatrix[i]=new ArrayList<VCFRecord>(30);
		}
		SVCandidate.setThreshold(hasDuoBarcode,minIndelLen,maxBpScore);
		
		AtomicInteger taskIdx=new AtomicInteger(0);
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads); 

		fixedThreadPool = Executors.newFixedThreadPool(threads); 
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new TransAndFilterWorker(supportMatrix,vcfRecordMatrix,accumulatedCntList,taskIdx,trTable));
		}
		fixedThreadPool.shutdown();
		while(!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) {}
		
		return vcfRecordMatrix;
	}
	
	
	public void clean() {
		for(int i=0,len=supportMatrix.length;i<len;++i) {
			ArrayList<SVSupport> svList=(ArrayList<SVSupport>)supportMatrix[i];
			svList.clear();
			supportMatrix[i]=null;
		}
		this.supportMatrix=null;
	}
	
	private class MergeWorker implements Runnable{
		private final Object[] supportMatrix;
		private final AtomicInteger taskIdx;
		
		public MergeWorker(Object[] supportMatrix,AtomicInteger taskIdx) {
			this.supportMatrix=supportMatrix;
			this.taskIdx=taskIdx;
		}
		
		@Override
		public void run() {
			final Object[] supportMatrix=this.supportMatrix;
			final int len=supportMatrix.length;
			final AtomicInteger taskIdx=this.taskIdx;
			
			ArrayList<SVSupport> newSupportList=new ArrayList<SVSupport>(30);
			ArrayList<SVSupport> newSupportList2=new ArrayList<SVSupport>(30);
			
			int idx;
			while((idx=taskIdx.getAndIncrement())<len) {
				ArrayList<SVSupport> supportList=(ArrayList<SVSupport>) supportMatrix[idx];
				supportMatrix[idx]=null;
				Collections.sort(supportList);
				merge(supportList,newSupportList);
				supportList.clear();
				merge(newSupportList,newSupportList2);
				newSupportList.clear();
				//ArrayList<SVCandidate> candidateList=translateAndFilter(newSupportList2);
				//candidateMatrix[idx]=candidateList;
				supportMatrix[idx]=newSupportList2;
				newSupportList2=supportList;
						
			}
			newSupportList.clear();
			newSupportList2.clear();
		}
		
		private void merge(ArrayList<SVSupport> inputList,ArrayList<SVSupport> resList) {
			final int cnt=inputList.size();
			
			for(int j=0;j<cnt;++j) {
				SVSupport svJ=inputList.get(j);
				if(svJ==null) {
					continue;
				}
				for(int k=j+1;k<cnt;++k) {
					SVSupport svK=inputList.get(k);
					if(svK==null) {
						continue;
					}
					if(svK.getLeftPos()-svJ.getLeftPos()<200) {
						if(svJ.similarWith(svK)) {
							if(svJ.hasShorterJunction(svK)) {
		    					svJ.mergeWith(svK);
		    					inputList.set(k,null);
		    				}else {
		    					svK.mergeWith(svJ);
		    					inputList.set(j, null);
						    	svJ=null;
						    	break;
		    				}
						}
					}else {
						break;
					}
				}
				if(svJ!=null) {
					resList.add(svJ);
				}
				//if(svJ!=null && svJ instanceof SplitSupport) {
				//	resList.add(svJ);
				//}
			}
		}
		
	}
	
	private class TransAndFilterWorker implements Runnable{
		
		private final Object[] supportMatrix,vcfRecordMatrix;
		private int[] accumulatedCntList;
		private final AtomicInteger taskIdx;
		private final TranslationTable trTable;

		TransAndFilterWorker(Object[] supportMatrix,Object[] vcfRecordMatrix,int[] accumulatedCntList,AtomicInteger taskIdx,TranslationTable trTable){
			this.supportMatrix=supportMatrix;
			this.vcfRecordMatrix=vcfRecordMatrix;
			this.accumulatedCntList=accumulatedCntList;
			this.taskIdx=taskIdx;
			this.trTable=trTable;
		}
		
		@Override
		public void run() {
			final Object[] supportMatrix=this.supportMatrix,vcfRecordMatrix=this.vcfRecordMatrix;
			final int[] accumulatedCntList=this.accumulatedCntList;
			final int len=supportMatrix.length;
			final AtomicInteger taskIdx=this.taskIdx;
			final TranslationTable trTable=this.trTable;
			
			SVCandidate svc=new SVCandidate();
			
			int idx;
			while((idx=taskIdx.getAndIncrement())<len) {
				ArrayList<SVSupport> supportList=(ArrayList<SVSupport>) supportMatrix[idx];
				ArrayList<VCFRecord> leftVcfList=(ArrayList<VCFRecord>)vcfRecordMatrix[idx];
				int accumulatedCnt=accumulatedCntList[idx];
				translateAndFilter(svc,supportList,leftVcfList,vcfRecordMatrix,accumulatedCnt,trTable);
			}
		}
		
		private void translateAndFilter(SVCandidate svc,ArrayList<SVSupport> svList,ArrayList<VCFRecord> leftVcfList,Object[] vcfRecordMatrix,int accumulatedCnt,TranslationTable trTable) {
			final int len=svList.size();
			for(int i=0;i<len;++i) {
				if(svList.get(i).translateAndFilter(trTable,svc)!=null) {
					svc.setSVNo(i+accumulatedCnt);
					synchronized(leftVcfList) {
						leftVcfList.add(svc.toVCFRecordLeft());
						//leftVcfList.add(svc.toVCFTemp());
					}
					int rightIdx=trTable.getIndexNo(svc.getRightCode());
					ArrayList<VCFRecord> rightVcfList=(ArrayList<VCFRecord>)vcfRecordMatrix[rightIdx];
					synchronized(rightVcfList) {
						rightVcfList.add(svc.toVCFRecorRight());
					}
				}
			}
		}
	}
}

