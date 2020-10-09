import java.io.IOException;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public final class ClassifyManager{
	public final static int Read_BOX_CAP=500;

	private final int threads;
	private final Path[] tempDataList,tempIndexList;
	private final int[] blockcntList;
	private final int r1BarStart,r2BarStart,r1BarLen,r2BarLen,r1InsStart,r2InsStart;
	
	private final KmerCollection kc;
	private final LongKmerCollection longkc,spacedkc;
	
	private ParallelReader r1ParallelReader,r2ParallelReader;
	private DoubleBlockingQueue dbQueue;
	
	public ClassifyManager(KmerCollection kc,LongKmerCollection longkc,LongKmerCollection spacedkc,FastqReader r1,FastqReader r2,
			               int threads,Path[] tempDataList,Path[] tempIndexList,int r1BarStart,int r2BarStart,
			               int r1BarLen,int r2BarLen,int r1InsStart,int r2InsStart){
		this.kc=kc;
		this.longkc=longkc;
		this.spacedkc=spacedkc;
		this.threads=threads;
		this.tempDataList=tempDataList;
		this.tempIndexList=tempIndexList;
		this.r1BarStart=r1BarStart;
		this.r2BarStart=r2BarStart;
		this.r1BarLen=r1BarLen;
		this.r2BarLen=r2BarLen;
		this.r1InsStart=r1InsStart;
		this.r2InsStart=r2InsStart;
		this.blockcntList=new int[threads];
		this.dbQueue=new DoubleBlockingQueue(threads); 
		this.r2ParallelReader=ParallelReader.createR2ParallelReader(r2);
		this.r1ParallelReader=ParallelReader.createR1ParallelReader(threads, r1, dbQueue, r2ParallelReader);
		
	}
	
	private void startWorkersAndWait() throws IOException, InterruptedException {
		final int k=kc.getK();
		final int longk=longkc.getLongK(),spacedk=spacedkc.getLongK();
		final int geneCodeLen=kc.getGeneCodeLen();
		
		dbQueue.putExtraBox(threads);
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads);
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new ClassifyWorker(i,k,longk,spacedk,geneCodeLen,tempDataList[i],tempIndexList[i],
                                                       r1BarStart,r2BarStart,r1BarLen,r2BarLen,r1InsStart,r2InsStart,
                                                       dbQueue,kc,longkc,spacedkc,blockcntList));
		}
		fixedThreadPool.shutdown();
		while (!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) { }
		
	}
	
	private void startParallelReader(){
		new Thread(this.r1ParallelReader).start();
		new Thread(this.r2ParallelReader).start();
	}

	public int[] classifyReads() throws IllegalBioFileException, InterruptedException, IOException {
		startParallelReader();
		startWorkersAndWait();
		clean();
		return this.blockcntList;
	}
	
	public void clean() {
		this.dbQueue.clean();
		this.dbQueue=null;
		this.r1ParallelReader=null;
		this.r2ParallelReader=null;
	}
}

final class DoubleBlockingQueue{
	private ArrayBlockingQueue<ReadBox> readReady,writeReady;
	private int capacity;
	
	DoubleBlockingQueue(int threads) {
		this.capacity=threads*2;
		this.readReady=new ArrayBlockingQueue<ReadBox>(threads*2);
		this.writeReady=new ArrayBlockingQueue<ReadBox>(threads*2);
	}
	
	public ReadBox tradeIn(ReadBox emptyBox) throws InterruptedException {
		writeReady.put(emptyBox);
		return readReady.take();
	}
	
	public ReadBox getEmpty() throws InterruptedException {
		return writeReady.take();
	}
	
	public void putFull(ReadBox emptyBox) throws InterruptedException {
		readReady.put(emptyBox);
	}
	
	public void putExtraBox(int num) throws InterruptedException {
		for(int i=0;i<num;++i) {
			writeReady.put(new ReadBox(ClassifyManager.Read_BOX_CAP));
		}
	}
	
	public int getCapacity() {
		return this.capacity;
	}
	
	public void clean() {
		this.readReady.clear();
		this.writeReady.clear();
		this.readReady=null;
		this.writeReady=null;
	}
}

final class ReadBox{
	private Read[] r1List,r2List;
	private int size;
	private ReadBoxItr r1Itr,r2Itr;
	
	public ReadBox(int size){
		this.size=size;
		this.r1List=new Read[size];
		this.r2List=new Read[size];
		for(int i=0;i<size;++i) {
			r1List[i]=new Read();
			r2List[i]=new Read();
		}
		this.r1Itr=new ReadBoxItr(r1List);
		this.r2Itr=new ReadBoxItr(r2List);
	}
	
	public ReadBoxItr getR1Itr() {
		this.r1Itr.reset();
		return this.r1Itr;
	}
	
	public ReadBoxItr getR2Itr() {
		this.r2Itr.reset();
		return this.r2Itr;
	}
	
	public void setNullAndClean() {
		for(int i=0;i<size;++i) {
			r1List[i].clean();
			r2List[i].clean();
			r1List[i]=null;
			r2List[i]=null;
		}
		r1List=null;
		r2List=null;
		r1Itr=null;
		r2Itr=null;
	}
	
	public boolean isNull() {
		return r1List==null;
	}
	
	private class ReadBoxItr implements Iterator<Read>{
		private Read[] readList;
		private int cursor;

		private ReadBoxItr(Read[] readList){
			this.readList=readList;
			this.cursor=0;
		}
		
		public void reset() {
			this.cursor=0;
		}
		
		@Override
		public boolean hasNext() {
			return this.cursor!=size;
		}

		@Override
		public Read next() {
			int i=this.cursor;
			if(i>=size) {
				throw new NoSuchElementException();
			}
			this.cursor=i+1;
			return readList[i];
		}
		
	}
}
