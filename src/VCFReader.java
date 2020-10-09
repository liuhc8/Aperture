import java.io.*;
import java.nio.file.Path;
import java.util.zip.*;

public class VCFReader {

	public final static int BUFFER_LEN=10000;
	
	private final BufferedReader vcfIn;
	private boolean eof;
	
	public VCFReader(Path vcfPath) throws FileNotFoundException,IOException {
		this.eof=false;
		try {
	    	if(vcfPath.toString().contains(".gz")) {
	    		this.vcfIn=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(vcfPath.toFile()))));
    		}else {
	    		this.vcfIn=new BufferedReader(new InputStreamReader(new FileInputStream(vcfPath.toFile())));
	    	}
    	}catch(FileNotFoundException e){
	    	throw new FileNotFoundException("VCF File: "+vcfPath.toFile()+" cannot be found!");
    	}
	}
	
	public void close() {
		try {
			if(this.vcfIn!=null) {
				this.vcfIn.close();
			}
		}catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public synchronized SNP[] readVCF() throws IOException, IllegalPositionException{
		if(this.eof) {
			return null;
		}
		
		final BufferedReader vcfIn=this.vcfIn;
		final int bufferLen=BUFFER_LEN;
		
		SNP[] snpArray=new SNP[bufferLen];
		String line=null;
		int cnt=0,posLast=0;;
		
		try {
	     	while(cnt<bufferLen && (line=vcfIn.readLine())!=null) {
		     	if(line.charAt(0)=='#') {
		    		continue;
	     		}
	    		String[] linelist=line.split("\\t");
	    		int pos=Integer.parseInt(linelist[1]);
	    		if(pos==posLast) {
	    			continue;
	    		}
	    		String chrom=linelist[0];
	    		if(!chrom.contains("chr")) {
	    			chrom="chr"+chrom;
	    		}
	    		byte[] ref=linelist[3].getBytes();
		    	String[] altStrList=linelist[4].split(",");
		    	snpArray[cnt++]=new SNP(chrom,pos,ref,altStrList[0].getBytes());
	    		/*for(String altStr:altStrList) {
		    		snpArray[cnt++]=new SNP(chrom,pos,ref,altStr.getBytes());
		    	}*/
	    		posLast=pos;
    		}
		}catch (NumberFormatException e) {
			throw new IllegalPositionException("VCF File Is Invalid!");
		}
		
		if(cnt<bufferLen) {
			this.eof=true;
		}
		return snpArray;	
	}
}

class SNP{
	String chrom;
	int pos;
	byte[] ref,alt;
	
	SNP(String chrom,int pos,byte[] ref,byte[] alt){
		this.chrom=chrom;
		this.pos=pos;
		this.ref=ref;
		this.alt=alt;
	}
}