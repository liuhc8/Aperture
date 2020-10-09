import java.util.HashSet;
import java.util.Random;
import java.util.Set;

public class GeneCodeGenerator {
	private Set<Integer> idset;
	private Set<Byte> byteset;
	private Random rand;
	private byte[] possave;
	
	GeneCodeGenerator(int n){
		this.idset=new HashSet<Integer>();
		this.byteset=new HashSet<Byte>();
		this.rand=new Random();
		this.possave=new byte[n];
	}
	
	public int getValidCode() {
		final Set<Integer> idset=this.idset;
		final Set<Byte> byteset=this.byteset;
		final Random rand=this.rand;
		final byte[] possave=this.possave;
		int code=0;
		
		do {
    		byteset.clear();
    		code=0;
		
	     	for(int i=0;i<possave.length;++i) {
     			do {
    				possave[i]=(byte)rand.nextInt(32);
	     		}while(byteset.contains(possave[i]));
	     		byteset.add(possave[i]);
        	}
		
	     	for(int i=0;i<possave.length;++i) {
		    	code |= 1<<possave[i];
	     	}	
		}while(idset.contains(code));
		idset.add(code);
		return code;
	}
	
	public static int testGeneCode(int bin,int geneCodeLen) {
		//return Integer.bitCount(bin)-geneCodeLen;
		int cnt=0;
		while(bin!=0) {
			++cnt;
			bin=bin&(bin-1);
			if(cnt>geneCodeLen) {
				return 1;
			}
		}
		if(cnt<geneCodeLen) {
			return -1;
		}
		return 0;
	}
}

class LongCodeGenerator {
	private Set<Long> idset;
	private Set<Byte> byteset;
	private Random rand;
	private byte[] possave;
	
	LongCodeGenerator(int n){
		this.idset=new HashSet<Long>();
		this.byteset=new HashSet<Byte>();
		this.rand=new Random();
		this.possave=new byte[n];
	}
	
	public long getValidCode() {
		final Set<Long> idset=this.idset;
		final Set<Byte> byteset=this.byteset;
		final Random rand=this.rand;
		final byte[] possave=this.possave;
		long code=0;
		
		do {
    		byteset.clear();
    		code=0;
		
	     	for(int i=0;i<possave.length;++i) {
     			do {
    				possave[i]=(byte)rand.nextInt(64);
	     		}while(byteset.contains(possave[i]));
	     		byteset.add(possave[i]);
        	}
		
	     	for(int i=0;i<possave.length;++i) {
		    	code |= 1L<<possave[i];
	     	}	
		}while(idset.contains(code));
		idset.add(code);
		return code;
	}
	
	public static int testGeneCode(long bin,int geneCodeLen) {
		int cnt=0;
		while(bin!=0) {
			++cnt;
			bin=bin&(bin-1);
			if(cnt>geneCodeLen) {
				return 1;
			}
		}
		if(cnt<geneCodeLen) {
			return -1;
		}
		return 0;
	}
}
