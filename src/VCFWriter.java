import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.OutputStreamWriter;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class VCFWriter {
	private final OutputStreamWriter vcfWriter;
	private final Path ciPath;
	private final String projectName,lnSep;
	private final String[] commands;
	private final Object[] vcfRecordMatrix;
	
	public VCFWriter(Path ciPath,Path outPath,String projectName,String[] commands,Object[] vcfRecordMatrix) throws FileNotFoundException, IOException, ClassNotFoundException {
		this.ciPath=ciPath;
		this.vcfWriter=new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPath.toFile())));
		this.lnSep=System.lineSeparator();
		this.commands=commands;
		this.projectName=projectName;
		this.vcfRecordMatrix=vcfRecordMatrix;
	}

	private ChromInfo readChromInfo() throws IOException, ClassNotFoundException {
		ObjectInputStream objIn=null;
		try {
			objIn=new ObjectInputStream(new GZIPInputStream(new FileInputStream(this.ciPath.toFile())));
			return (ChromInfo)objIn.readObject();
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
	
	private String getVCFHeader(ChromInfo ci) {
		StringBuilder strBuilder=new StringBuilder(1000);
		String lnSep=this.lnSep;
		String[] commands=this.commands;
		
		Date now = new Date(); 
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd HH:mm:ss");
		String currentDate = dateFormat.format(now);
		
		strBuilder.append("##fileformat=VCFv4.2").append(lnSep);
		strBuilder.append("##fileDate=").append(currentDate).append(lnSep);
		strBuilder.append("##source=Aperture ").append(ApertureMain.APERTURE_VERSION);;
		for(String str:commands) {
			strBuilder.append(" ").append(str);
		}
		strBuilder.append(lnSep);
		strBuilder.append("##reference=file://").append(ci.getFastaPath()).append(lnSep);

		strBuilder.append("##FILTER=<ID=PASS,Description=\"All filters passed\">").append(lnSep);
		strBuilder.append("##FILTER=<ID=LOW_QUAL,Description=\"Low quality call\">").append(lnSep);
		strBuilder.append("##FILTER=<ID=FAKE_BP,Description=\"False positive variant caused by imprecise k-mer based mapping\">").append(lnSep);
		strBuilder.append("##FILTER=<ID=SMALL_EVENT,Description=\"Event size is smaller than the minimum reportable size\">").append(lnSep);
		strBuilder.append("##FILTER=<ID=PE_ONLY,Description=\"No soft clips support this variant\">").append(lnSep);
		strBuilder.append("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">").append(lnSep);
		strBuilder.append("##INFO=<ID=STRANDS,Number=.,Type=String,Description=\"Strand orientation of the adjacency\">").append(lnSep);
		strBuilder.append("##INFO=<ID=REFQUA,Number=1,Type=Float,Description=\"K-mer mapping quality of reference junction\">").append(lnSep);
		strBuilder.append("##INFO=<ID=VARQUA,Number=1,Type=Float,Description=\"K-mer mapping quality of variant junction\">").append(lnSep);
		strBuilder.append("##INFO=<ID=REFKMER,Number=1,Type=Integer,Description=\"Number of k-mers supporting reference junction in average\">").append(lnSep);
		strBuilder.append("##INFO=<ID=VARKMER,Number=1,Type=Integer,Description=\"Number of k-mers supporting variant junction in average\">").append(lnSep);
		strBuilder.append("##INFO=<ID=BPSEQQUA,Number=1,Type=Float,Description=\"Quality of sequence spanning breakpoint junction\">").append(lnSep);
		strBuilder.append("##INFO=<ID=PARID,Number=1,Type=String,Description=\"ID of partner breakend\">").append(lnSep);
		strBuilder.append("##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">").append(lnSep);
		strBuilder.append("##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">").append(lnSep);		
		strBuilder.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype (Not applicable)\">").append(lnSep);
		strBuilder.append("##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Count of split reads supporting the breakpoint\">").append(lnSep);
		strBuilder.append("##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Count of paired-end reads supporting the breakpoint\">").append(lnSep);
		strBuilder.append("##FORMAT=<ID=REFSR,Number=1,Type=Integer,Description=\"Count of split reads supporting the reference junction\">").append(lnSep);
		strBuilder.append("##FORMAT=<ID=VARSR,Number=1,Type=Integer,Description=\"Count of split reads supporting the variant junction\">").append(lnSep);
		strBuilder.append("##FORMAT=<ID=BAR,Number=1,Type=Integer,Description=\"Count of cfDNA molecules supporting the breakpoint\">").append(lnSep);
		strBuilder.append("##FORMAT=<ID=UBAR,Number=1,Type=Integer,Description=\"Count of cfDNA molecules with only one read support\">").append(lnSep);
		
		for(Pair<String,Integer> pair:ci) {
			strBuilder.append("##contig=<ID=");
			strBuilder.append(pair.getKey());
			strBuilder.append(",length=");
			strBuilder.append(pair.getValue());
			strBuilder.append(">").append(lnSep);
		}
		
		strBuilder.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t").append(projectName).append(lnSep);;
		
		return strBuilder.toString();
	}
	
	public void sortAndWriteVcf() throws ClassNotFoundException, IOException {
		Object[] vcfRecordMatrix=this.vcfRecordMatrix;
		OutputStreamWriter vcfWriter=this.vcfWriter;
		ChromInfo ci=readChromInfo();
		vcfWriter.write(getVCFHeader(ci));
		for(Object obj:vcfRecordMatrix) {
			ArrayList<VCFRecord> vcfList=(ArrayList<VCFRecord>)obj;
			Collections.sort(vcfList);
			for(VCFRecord vcfLine:vcfList) {
				vcfWriter.write(vcfLine.toString());
			}
			vcfList.clear();
		}
	}
	
	public void flush() throws IOException {
		this.vcfWriter.flush();
	}
	
	public void close() {
		if(this.vcfWriter!=null) {
			try {
				this.vcfWriter.close();
			}catch(IOException e) {
				e.printStackTrace();
			}
		}
	}
}

class VCFRecord implements Comparable<VCFRecord>{
	private String vcfLine;
	private int genomicPos;
	
	VCFRecord(String vcfLine,int genomicPos){
		this.vcfLine=vcfLine;
		this.genomicPos=genomicPos;
	}
	
	@Override
	public int compareTo(VCFRecord o) {
		return (this.genomicPos < o.genomicPos) ? -1 : ((this.genomicPos == o.genomicPos) ? 0 : 1);
	}
	
	public String toString() {
		return this.vcfLine;
	}
	
}
