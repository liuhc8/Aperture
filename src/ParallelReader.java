import java.io.IOException;
import java.util.Iterator;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

class ParallelReader implements Runnable{
	
	private final int threads;
	private final boolean isR1;
	private final FastqReader reader;
	private final DoubleBlockingQueue dbQueue;
	private final ParallelReader r2;
	private Iterator<Read> itr;
	
	private final static Lock helperLock;
	private final static Condition helperGo;
	private final static Condition helperStop;
	
	static {
		helperLock=new ReentrantLock();
		helperGo=helperLock.newCondition();
		helperStop=helperLock.newCondition();
	}
	
	public static ParallelReader createR1ParallelReader(int threads,FastqReader r1FqReader,DoubleBlockingQueue dbQueue,ParallelReader r2) {
		return new ParallelReader(threads,true,r1FqReader,dbQueue,r2);
	}
	
	public static ParallelReader createR2ParallelReader(FastqReader r2FqReader) {
		return new ParallelReader(-1,false,r2FqReader,null,null);
	}
	
	private ParallelReader(int threads,boolean isR1,FastqReader reader,DoubleBlockingQueue dbQueue,ParallelReader r2){
		this.threads=threads;
		this.isR1=isR1;
		this.reader=reader;
		this.dbQueue=dbQueue;
		this.r2=r2;
	}
	
	private void setItr(Iterator<Read> itr) {
		this.itr=itr;
	}
	
	private boolean hasItr() {
		return this.itr!=null;
	}
	
	@Override
	public void run() {
		if(isR1) {
			try {
				deliverR1();
			} catch (IllegalBioFileException | InterruptedException | IOException e) {
				e.printStackTrace();
			}
		}else {
			try {
				deliverR2();
			} catch (IllegalBioFileException | IOException | InterruptedException e) {
				e.printStackTrace();
			}
		}
		
	}
	
	private void deliverR2() throws IllegalBioFileException, IOException, InterruptedException {
		final Lock helperLock=ParallelReader.helperLock;
		final Condition helperGo=ParallelReader.helperGo;
		final Condition helperStop=ParallelReader.helperStop;
		
		boolean eof=false;
		
		try{
			while(!Thread.interrupted()) {
     			helperLock.lock();
	     		try{
	           		while(!hasItr()){
	        			helperGo.await();
        			}
	    		}finally{
	    			helperLock.unlock();
    			}
			
	     		eof=fillReadBox();
			
	    		helperLock.lock();
	    		try{
	    			setItr(null);
	    			helperStop.signal();
    			}finally{
		    		helperLock.unlock();
	    		}
	    		
	    		if(eof){
	    			return;
	    		}
			}
		}finally {
			reader.close();
		}
	}
	
	private void deliverR1() throws IllegalBioFileException, InterruptedException, IOException{
		final DoubleBlockingQueue dbQueue=this.dbQueue;
		final ParallelReader r2=this.r2;
		
		boolean eof=false;
		
		final Lock helperLock=ParallelReader.helperLock;
		final Condition helperGo=ParallelReader.helperGo;
		final Condition helperStop=ParallelReader.helperStop;
		
		try {
	    	while(!Thread.interrupted()) {
	    		
	    		ReadBox readBox=dbQueue.getEmpty();
	    		setItr(readBox.getR1Itr());
	    		
	    		helperLock.lock();
	    		try{
	    			r2.setItr(readBox.getR2Itr());
	    			helperGo.signal();
	    		}finally{
	    			helperLock.unlock();
	    		}
	    		
	    		eof=fillReadBox();
	    	    
	    	    helperLock.lock();
	    	    try{
	    	    	while(r2.hasItr()){
	    	    		helperStop.await();
	    	    	}
	    	    }finally{
	    	    	helperLock.unlock();
	    	    }
    			
	    	    dbQueue.putFull(readBox);
	    	    
	    	    if(eof) {
	    	    	for(int i=0;i<dbQueue.getCapacity();++i) {
	    	    		readBox=dbQueue.getEmpty();
	    	    		readBox.setNullAndClean();
	    	    		dbQueue.putFull(readBox);
	    	    	}
	    	    	return;
	    	    }

	      	}
	    }finally {
			reader.close();
		}
	}
	
	private boolean fillReadBox() throws IllegalBioFileException, IOException{
		final Iterator<Read> itr=this.itr;
		final FastqReader reader=this.reader;
		
		boolean eof=false;
		
		while(itr.hasNext()) {
			Read read=itr.next();
			if(eof) {
				read.clean();
			}else {
				eof=!reader.readFq(read);
			}
		}
		return eof;
	}
}
