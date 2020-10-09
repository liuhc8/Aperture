import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.file.Path;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.GZIPInputStream;

public class BpClusterManager {
	private final Path ttPath,mergedDataPath;
	private final int threads,mergeK;
	private final double similarity;
	
	public BpClusterManager(Path ttPath,String workDir,String projectName,int threads,int mergeK, double similarity,Path mergedDataPath) {
		this.ttPath=ttPath;
		this.mergeK=mergeK;
		this.threads=threads;  
		this.similarity=similarity;
		this.mergedDataPath=mergedDataPath;
	}
	
	public SVCollection mergeBreakpoints(BpIndexList finalIndexlist) throws  IOException, InterruptedException, ClassNotFoundException {
		final Path mergedDataPath=this.mergedDataPath;
		final int mergeK=this.mergeK;
		final int threads=this.threads;
		final double similarity=this.similarity;
		
		final TranslationTable trTable=readTranslationTable();
		final int k=trTable.getK();
		
		SVCollection svCollection=new SVCollection(trTable,threads);
		BpIndexList.SynchronizedGetItr synGetItr=finalIndexlist.getSynchronizedGetItr();
		
		ExecutorService fixedThreadPool = Executors.newFixedThreadPool(threads); 
		
		for(int i=0;i<threads;++i) {
			fixedThreadPool.execute(new BpClusterWorker(synGetItr,mergedDataPath,k,mergeK,similarity,svCollection,trTable));
		}
		
		fixedThreadPool.shutdown();
		while(!fixedThreadPool.awaitTermination(1, TimeUnit.SECONDS)) {}
		
		finalIndexlist.clean();
		
		return svCollection;
	}
	
	private TranslationTable readTranslationTable() throws IOException, ClassNotFoundException {
		ObjectInputStream objIn=null;
		try {
			objIn=new ObjectInputStream(new GZIPInputStream(new FileInputStream(this.ttPath.toFile())));
			return (TranslationTable)objIn.readObject();
		}finally {
			if(objIn!=null) {
				try{
					objIn.close();
				}catch(IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
}

