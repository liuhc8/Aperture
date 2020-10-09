import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class ResultsSorter {
	
	private final Path[] tempDataList,tempIndexList;
	private final Path mergedDataPath;
	private final int threads;
	private BpIndexList[] bpIndexMatrix;
	
	
	public ResultsSorter(Path[] tempDataList,Path[] tempIndexList,int threads,Path mergedDataPath,int[] blockCntList){
		this.tempDataList=tempDataList;
		this.tempIndexList=tempIndexList;
		this.threads=threads;
		this.mergedDataPath=mergedDataPath;
		BpIndexList[] bpIndexMatrix=new BpIndexList[threads];
		for(int i=0;i<threads;++i) {
			bpIndexMatrix[i]=new BpIndexList(blockCntList[i]);
		}
		this.bpIndexMatrix=bpIndexMatrix;
	}
	
	public BpIndexList sortResults() throws InterruptedException, ExecutionException {
		sortEachFile();
		BpIndexList finalIndexList=mergeIndexs();
		mergeBpFiles();
		clean();
		return finalIndexList;
	}
	
	private void sortEachFile() throws InterruptedException, ExecutionException{
		final int threads=this.threads;
		final Path[] tempIndexList=this.tempIndexList;
		final BpIndexList[] bpIndexMatrix=this.bpIndexMatrix;
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new SortEachFileWorker(bpIndexMatrix[i],tempIndexList[i]));
		}
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }

	}
	
	private class SortEachFileWorker implements Runnable{
		private final Path tempIndexPath;
		private final BpIndexList bpIndexList;

		SortEachFileWorker(BpIndexList bpIndexList,Path tempIndexPath){
			this.bpIndexList=bpIndexList;
			this.tempIndexPath=tempIndexPath;
		}
		
		private void loadIdxFile() throws IOException{
			BpIndexList bpIndexList=this.bpIndexList;
			Path tempIndexPath=this.tempIndexPath;
			
			FileInputStream idxIn=null;
			ByteBuffer inBb=null;
			
			try {
				idxIn=new FileInputStream(tempIndexPath.toFile());
				FileChannel idxFc=idxIn.getChannel();
				inBb=ByteBuffer.allocateDirect(4096);
				
				while(idxFc.read(inBb)!=-1) {
					inBb.flip();
					while(inBb.remaining()>=8) {
						int code=inBb.getInt();
						int len=inBb.getInt();
						bpIndexList.addCodeAndLen(code, len);
					}
					inBb.compact();
				}
			}finally {
				if(idxIn!=null) {
					try {
						idxIn.close();
					}catch(IOException e) {
						e.printStackTrace();
					}
				}
			}
		}
		
		private void sortByCode() {
			bpIndexList.sortByCode();
		}
		
		@Override
		public void run() {
			try {
				loadIdxFile();
				sortByCode();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private BpIndexList mergeIndexs() {
		LoserTree loserTree=new LoserTree(this.bpIndexMatrix);
		return loserTree.mergeAndSortIndexs();
	}
	
	private void mergeBpFiles() throws InterruptedException {
		final int threads=this.threads;
		final Path[] tempDataList=this.tempDataList;
		final BpIndexList[] bpIndexMatrix=this.bpIndexMatrix;
		final Path mergedDataPath=this.mergedDataPath;
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads); 
    	for(int i=0;i<threads;++i) {
    		fixedThreadPool.execute(new MergeBpFilesWorker(tempDataList[i],mergedDataPath,bpIndexMatrix[i]));
		}
    	fixedThreadPool.shutdown();
    	while(!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) {}
	}
	
	private class MergeBpFilesWorker implements Runnable{
		private final Path tempDataPath,mergedDataPath;
		private BpIndexList bpIndexList;
		
		public MergeBpFilesWorker(Path tempDataPath,Path mergedDataPath,BpIndexList bpIndexList){
			this.tempDataPath=tempDataPath;
			this.mergedDataPath=mergedDataPath;
			this.bpIndexList=bpIndexList;	
		}
		
		private void sortByInitStart() {
			bpIndexList.sortByInitStart();
		}
		
		@Override
		public void run() {
			try {
				sortByInitStart();
		    	mergeFiles();
			}catch(IOException e) {
				e.printStackTrace();;
			}
		}
		
		private void mergeFiles() throws IOException{
			final BpIndexList bpIndexList=this.bpIndexList;
			final BpIndexList.GetAndSetItr getIter=bpIndexList.getGetAndSetItr();
			final Path tempDataPath=this.tempDataPath;
			final Path mergedDataPath=this.mergedDataPath;
			
			FileInputStream in=null;
			RandomAccessFile out=null;
			
			try {
				in=new FileInputStream(tempDataPath.toFile());
				out=new RandomAccessFile(mergedDataPath.toFile(),"rw");
				FileChannel inFc=in.getChannel();
				FileChannel outFc=out.getChannel();
				
				while(getIter.hasNext()) {
					getIter.moveToNext();
					long initStart=getIter.getCurrentInitStart();
					long mergedStart=getIter.getCurrentMergedStart();
					int len=getIter.getCurrentLen();		
					
					//synchronized(bpIndexMatrix) {
					//	System.out.println("+++ "+no+" from "+initStart+"  "+len+"  to "+mergedStart);
					//}
					
					FileLock lock = outFc.lock(mergedStart,len, false);
					out.seek(mergedStart);
					for(int remainBytes=len;remainBytes>0;) {
						long transfered=inFc.transferTo(initStart, remainBytes, outFc);
						initStart+=transfered;
						remainBytes-=transfered;
					}
					lock.release();
				}				
				outFc.force(true);
				
			}finally {
				if(out!=null) {
					try {
						out.close();
					}catch(IOException e) {
						e.printStackTrace();
					}
				}
				if(in!=null) {
					try {
						in.close();
					}catch(IOException e) {
						e.printStackTrace();
					}
				}
			}

		}
		
	}
	
	private void clean() {
		for(int i=0,len=bpIndexMatrix.length;i<len;++i) {
			bpIndexMatrix[i].clean();
			bpIndexMatrix[i]=null;
		}
	}
}


class CodeStartLen{
	public int code,len;
	public long start;
}