import java.io.*;
import java.util.zip.*;

final class FastqReader {
	public final static int DEFAULT_BUFFER_SIZE=102400;
	public final static int MAX_READ_BYTE=800;
	
	private final InputStream fqInput;
	private final byte[] buffer;
	private int nowByte;
	private boolean eof;
	private int byteLen;
	
	FastqReader(File fqPath) throws FileNotFoundException, IOException {
		this(fqPath,102400);
	}
	
	FastqReader(File fqPath,int bufferSize) throws FileNotFoundException, IOException {
		InputStream input=null;
		
		try {
			if(fqPath.getName().contains(".gz")) {
				input=new GZIPInputStream(new FileInputStream(fqPath),65536);   
			}else {
				input=new FileInputStream(fqPath);
			}
		}catch(FileNotFoundException e){
			throw new FileNotFoundException("Fastq File Is Missing! :"+fqPath);
		}
		
		this.fqInput=input;
		this.buffer=new byte[bufferSize];
		this.eof=false;
		this.nowByte=bufferSize;
	}
	
	
	void close(){
		try {
	    	if(this.fqInput!=null) {
	    		this.fqInput.close();
	    	}
		}catch(IOException e) {
			e.printStackTrace();
		}
	}

	private void fill() throws IOException {
		final InputStream fqInput=this.fqInput;
		final byte[] buffer=this.buffer;
		final int bufferLen=buffer.length;
		int start=this.nowByte;
		int len=bufferLen-start;
		System.arraycopy(buffer, start, buffer, 0, len);
		start=len;
		len=bufferLen-start;
		int n=0;
		do {
	    	n=fqInput.read(buffer, start, len);
	    	start+=n;
	    	len=len-n;
	    }while(n > 0);
		if(n==-1) {
			this.eof=true;
			//this.byteLen=bufferLen-len+1;
			this.byteLen=start;
		}else {
			this.byteLen=bufferLen;
		}
		this.nowByte=0;
	}
	
	public boolean readFq(Read read) throws IOException, IllegalBioFileException {
		final byte[] buffer=this.buffer;
		if(this.eof) {
			if(this.nowByte>=this.byteLen) {
				return false;
			}
		}else {
			if(this.byteLen-this.nowByte<MAX_READ_BYTE) {
	     		fill();
	    	}
		}
		int start=this.nowByte;
		int n=start;
		while(buffer[n]!=10) {
			++n;
		}
		int nameLen=n-start;
		++n;
		int seqStart=n;
		
		while(buffer[n]!=10) {
			++n;
		}
		int seqLen=n-seqStart;
		
		
		++n;
		if(buffer[n]!=43) {
			throw new IllegalBioFileException("Not in Fastq format!");
		}
		
		while(buffer[n]!=10) {
			++n;
		}
		
		int quaStart=n+1;
		
		read.copyFromByte(buffer, start, seqStart, quaStart, nameLen, seqLen);
		
		this.nowByte=quaStart+seqLen+1;

		return true;
	}

	
}