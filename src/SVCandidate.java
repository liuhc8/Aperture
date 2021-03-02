import java.math.BigDecimal;
import java.util.ArrayList;

public class SVCandidate{
	private static double BAR_MINUS_UNI_THRESHOLD=2.0d,MAX_BP_SCORE=2.0d;
	private static int MIN_INDEL_LEN=50;
	
	private int svNo;
	private boolean leftForward,rightForward,passFakeBpFilter;
	
	private String leftChrom;
	private int leftPos,leftKmers;
	private double leftScoreAvg;
	
	private String rightChrom;
	private int rightPos,rightKmers;
	private double rightScoreAvg;
	
	private int varSeqLen;
	private char refLeft,refRight;
	private String varLeft,varRight,homSeq;
	private double bpScoreAvg;
	
	private int splitCnt,peCnt,leftCnt,rightCnt;
	private double molecularCnt;
	private int uniqueBarCnt;
	
	private StringBuilder strBuilder;
	
	//#########################################################
	int leftCode,leftCodePos,rightCode,rightCodePos;
	
	public SVCandidate() {
		this.strBuilder=new StringBuilder(100);
	}
	
	public static void setThreshold(boolean duoBarcode,int minIndelLen,double maxBpScore) {
		BAR_MINUS_UNI_THRESHOLD=duoBarcode?2.499d:1.999d;
		MIN_INDEL_LEN=minIndelLen;
		MAX_BP_SCORE=maxBpScore;
	}
	
	public void setSVNo(int svNo) {
		this.svNo=svNo;
	}
	
	public void setPassFakeBpFilter(byte flag) {
		passFakeBpFilter=flag==3;
	}
	
	public void setCode(int leftCode,int leftCodePos,int rightCode,int rightCodePos) {
		this.leftCode=leftCode;
		this.leftCodePos=leftCodePos;
		this.rightCode=rightCode;
		this.rightCodePos=rightCodePos;
	}
	
	public void setLeft(Segment leftSeg,int leftPos,int leftKmers,double leftScoreAvg,boolean leftForward,int k) {
		int leftPosOffset=leftForward?2:(k-1);
		this.leftChrom=leftSeg.chrom;
		this.leftPos=leftSeg.start+leftPos+leftPosOffset;
		this.leftKmers=leftKmers;
		this.leftScoreAvg=leftScoreAvg;
		this.leftForward=leftForward;
	}
	
	public void setRight(Segment rightSeg,int rightPos,int rightKmers,double rightScoreAvg,boolean rightForward,int k) {
		int rightPosOffset=rightForward?(k-1):2;
		this.rightChrom=rightSeg.chrom;
		this.rightPos=rightSeg.start+rightPos+rightPosOffset;
		this.rightKmers=rightKmers;
		this.rightScoreAvg=rightScoreAvg;
		this.rightForward=rightForward;
	}
	
	public void setSplit(BpSeq bpSeq,double bpScoreAvg,int splitCnt,int peCnt,int leftCnt,int rightCnt,int leftIdentical,int rightIdentical) {
		this.bpScoreAvg=bpScoreAvg;
		this.splitCnt=splitCnt;
		this.peCnt=peCnt;
		this.leftCnt=leftCnt;
		this.rightCnt=rightCnt;
		
		int bpLen=bpSeq.getBpLen();
		int varSeqLen,leftOffset,rightOffset;
		if(leftIdentical+rightIdentical<=bpLen) {
			varSeqLen=bpLen-leftIdentical-rightIdentical+1;
			leftOffset=leftIdentical-1;
			rightOffset=rightIdentical-1;
			this.homSeq=null;
		}else {
			varSeqLen=1;
			leftOffset=bpLen-rightIdentical-1;
			rightOffset=rightIdentical-1;
			int homLen=leftIdentical+rightIdentical-bpLen;
			this.homSeq=this.leftForward?bpSeq.toStringBp(leftOffset+1, homLen):bpSeq.toStringRevBp(bpLen-leftIdentical, homLen);
		}
		
		boolean leftForward=this.leftForward;
		boolean rightForward=this.rightForward;
		
		int refPos=leftForward?this.leftPos+leftOffset:this.leftPos-leftOffset;
		int varPos=rightForward?this.rightPos-rightOffset:this.rightPos+rightOffset;
		String varSeqLeft=leftForward?bpSeq.toStringBp(leftOffset, varSeqLen):bpSeq.toStringRevBp(rightOffset+1, varSeqLen);
		String varSeqRight=rightForward?bpSeq.toStringBp(leftOffset+1, varSeqLen):bpSeq.toStringRevBp(rightOffset, varSeqLen);
		this.varSeqLen=varSeqLen;
		this.refLeft=leftForward?varSeqLeft.charAt(0):varSeqLeft.charAt(varSeqLeft.length()-1);
		this.refRight=rightForward?varSeqRight.charAt(varSeqRight.length()-1):varSeqRight.charAt(0);
		this.varLeft=buildVarStr(this.rightChrom,String.valueOf(varPos),varSeqLeft,leftForward,rightForward);
		this.varRight=buildVarStr(this.leftChrom,String.valueOf(refPos),varSeqRight,!rightForward,!leftForward);
		this.leftPos=refPos;
		this.rightPos=varPos;
	}
	
	private String buildVarStr(String varChrom,String varPos,String varSeq,boolean refForward,boolean varForward) {
		StringBuilder sb=this.strBuilder;
		if(refForward) {
			if(varForward) {
				sb.append(varSeq).append("[").append(varChrom).append(":").append(varPos).append("[");
			}else {
				sb.append(varSeq).append("]").append(varChrom).append(":").append(varPos).append("]");
			}
		}else {
			if(varForward) {
				sb.append("[").append(varChrom).append(":").append(varPos).append("[").append(varSeq);
			}else {
				sb.append("]").append(varChrom).append(":").append(varPos).append("]").append(varSeq);
			}
		}
		String var=sb.toString();
		sb.setLength(0);
		return var;
	}
	
	public void setPE(int peCnt) {
		this.bpScoreAvg=0;
		this.splitCnt=0;
		this.peCnt=peCnt;
		this.leftCnt=0;
		this.rightCnt=0;
		this.homSeq=null;
		this.refLeft='N';
		this.refRight='N';
		this.varLeft=buildVarStr(this.rightChrom,String.valueOf(this.rightPos),"N",leftForward,rightForward);;
		this.varRight=buildVarStr(this.leftChrom,String.valueOf(this.leftPos),"N",!rightForward,!leftForward);
		this.varSeqLen=0;
	}
	
	public void setUMI(double molecularCnt,int uniqueBarCnt) {
		this.molecularCnt=molecularCnt;
		this.uniqueBarCnt=uniqueBarCnt;
	}
	

	public VCFRecord toVCFRecordLeft() {
		boolean passEvidenceFilter=passEvidenceFilter();
		boolean isSmallEvent=isSmallEvent();
		boolean isPEOnly=isPEOnly();
		
		StringBuilder sb=this.strBuilder;
		String tab="\t";
		sb.append(leftChrom).append(tab).append(leftPos).append(tab).append("aperture").append(svNo).append("_1").append(tab).append(refLeft)
		  .append(tab).append(varLeft).append(tab).append(".").append(tab);
		if(passFakeBpFilter && passEvidenceFilter && !isSmallEvent) {
			sb.append("PASS");
		}else {
			if(!passEvidenceFilter) {
				sb.append("LOW_QUAL;");
			}
			if(!passFakeBpFilter) {
				sb.append("FAKE_BP;");
			}
			if(isSmallEvent) {
				sb.append("SMALL_EVENT;");
			}
			/*if(isPEOnly) {
				sb.append("PE_ONLY;");
			}*/
			sb.setLength(sb.length()-1);    //delete last ";"
		}
		sb.append(tab);
		if(isPEOnly) {
			sb.append("IMPRECISE;");
		}else {
			sb.append("PRECISE;");
		}
		sb.append("SVTYPE=BND;STRANDS=");
		if(this.leftForward) {
			sb.append("+");
		}else {
			sb.append("-");
		}
		if(this.rightForward) {
			sb.append("+");
		}else {
			sb.append("-");
		}
		sb.append(";REFQUA=").append(String.format("%.2f", this.leftScoreAvg)).append(";VARQUA=").append(String.format("%.2f", this.rightScoreAvg)).append(";REFKMER=").append(this.leftKmers)
		  .append(";VARKMER=").append(this.rightKmers).append(";BPSEQQUA=").append(String.format("%.2f", this.bpScoreAvg)).append(";PARID=").append("aperture").append(this.svNo).append("_2");
		if(this.homSeq!=null) {
			sb.append(";HOMLEN=").append(this.homSeq.length()).append(";HOMSEQ=").append(this.homSeq);
		}
		sb.append(tab).append("GT:SR:PE:REFSR:VARSR:BAR:UBAR").append(tab);
		sb.append("./.:").append(this.splitCnt).append(":").append(this.peCnt).append(":").append(this.leftCnt).append(":").append(this.rightCnt)
		  .append(":").append(this.molecularCnt).append(":").append(this.uniqueBarCnt).append(System.lineSeparator());
		
		String vcf=sb.toString();
		sb.setLength(0);
		
		VCFRecord vcfRecord=new VCFRecord(vcf,leftCodePos);
		
		return vcfRecord;
	}
	
	public VCFRecord toVCFRecorRight() {
		boolean passEvidenceFilter=passEvidenceFilter();
		boolean isSmallEvent=isSmallEvent();
		boolean isPEOnly=isPEOnly();
		
		StringBuilder sb=this.strBuilder;
		String tab="\t";
		sb.append(rightChrom).append(tab).append(rightPos).append(tab).append("aperture").append(svNo).append("_2").append(tab).append(refRight)
		  .append(tab).append(varRight).append(tab).append(".").append(tab);
		if(passFakeBpFilter && passEvidenceFilter && !isSmallEvent) {
			sb.append("PASS");
		}else {
			if(!passEvidenceFilter) {
				sb.append("LOW_QUAL;");
			}
			if(!passFakeBpFilter) {
				sb.append("FAKE_BP;");
			}
			if(isSmallEvent) {
				sb.append("SMALL_EVENT;");
			}
			/*if(isPEOnly) {
				sb.append("PE_ONLY;");
			}*/
			sb.setLength(sb.length()-1);    //delete last ";"
		}
		sb.append(tab);
		if(isPEOnly) {
			sb.append("IMPRECISE;");
		}else {
			sb.append("PRECISE;");
		}
		sb.append("SVTYPE=BND;STRANDS=");
		if(this.rightForward) {
			sb.append("-");
		}else {
			sb.append("+");
		}
		if(this.leftForward) {
			sb.append("-");
		}else {
			sb.append("+");
		}
		sb.append(";REFQUA=").append(String.format("%.2f", this.rightScoreAvg)).append(";VARQUA=").append(String.format("%.2f", this.leftScoreAvg)).append(";REFKMER=").append(this.rightKmers)
		  .append(";VARKMER=").append(this.leftKmers).append(";BPSEQQUA=").append(String.format("%.2f", this.bpScoreAvg)).append(";PARID=").append(this.svNo).append("_1");
		if(this.homSeq!=null) {
			sb.append(";HOMLEN=").append(this.homSeq.length()).append(";HOMSEQ=").append(this.homSeq);
		}
		sb.append(tab).append("GT:SR:PE:REFSR:VARSR:BAR:UBAR").append(tab);
		sb.append("./.:").append(this.splitCnt).append(":").append(this.peCnt).append(":").append(this.rightCnt).append(":").append(this.leftCnt)
		  .append(":").append(this.molecularCnt).append(":").append(this.uniqueBarCnt).append(System.lineSeparator());
		
		String vcf=sb.toString();
		sb.setLength(0);
		
		VCFRecord vcfRecord=new VCFRecord(vcf,rightCodePos);
		
		return vcfRecord;
	}
	
	public int getRightCode() {
		return this.rightCode;
	}
	
	public VCFRecord toVCFTemp() {
		StringBuilder strb=new StringBuilder(100);
		String tab="\t";
		strb.append(leftChrom);
		strb.append(tab);
		strb.append(leftPos);
		strb.append(tab);
		
		strb.append(leftCode);
		strb.append(tab);
		strb.append(leftCodePos);
		strb.append(tab);
		
		strb.append(leftForward);
		strb.append(tab);
		strb.append(0);
		strb.append(tab);
		strb.append(leftKmers);
		strb.append(tab);
		strb.append(String.format("%.2f", leftScoreAvg));
		strb.append(tab);
		
		strb.append(rightChrom);
		strb.append(tab);
		strb.append(rightPos);
		strb.append(tab);
		
		strb.append(rightCode);
		strb.append(tab);
		strb.append(rightCodePos);
		strb.append(tab);
		
		strb.append(rightForward);
		strb.append(tab);
		strb.append(0);
		strb.append(tab);
		strb.append(rightKmers);
		strb.append(tab);
		strb.append(String.format("%.2f", rightScoreAvg));
		strb.append(tab);
		
		strb.append(varLeft);
		strb.append(tab);
		
		strb.append(splitCnt);
		strb.append(tab);
		strb.append(leftCnt);
		strb.append(tab);
		strb.append(rightCnt);
		strb.append(tab);
		strb.append(peCnt);
		strb.append(tab);
		strb.append(String.format("%.1f",molecularCnt));
		strb.append(tab);
		strb.append(uniqueBarCnt);
		strb.append(tab);
		strb.append(String.format("%.2f", bpScoreAvg));
		strb.append(tab);
		strb.append(passEvidenceFilter());
		strb.append(tab);
		strb.append(passFakeBpFilter);
		strb.append(tab);
		strb.append(isSmallEvent());
		strb.append(tab);
		strb.append(isPEOnly());
		strb.append("\n");
		String str=strb.toString();
		VCFRecord vcfRecord=new VCFRecord(str,leftCodePos);
		return vcfRecord;
	}
	
	public boolean isPEOnly() {
		return this.splitCnt==0;
	}
	
	public boolean isSmallEvent() {	
		if(!this.leftChrom.equals(this.rightChrom)) {
			return false;
		}else{
			int varLen=Math.abs(this.rightPos-this.leftPos);
			//return Math.abs(this.rightPos-this.leftPos)<MIN_INDEL_LEN && this.varSeqLen-1<MIN_INDEL_LEN;
			return varLen<MIN_INDEL_LEN || Math.abs(varLen-this.varSeqLen)<10;
		}
	}
	
	public boolean passEvidenceFilter() {

		if((double)leftKmers/(leftCnt+peCnt)>2.999d && (double)rightKmers/(rightCnt+peCnt)>2.999d && bpScoreAvg<MAX_BP_SCORE && uniqueBarCnt/molecularCnt<0.801d) {
			if(molecularCnt-(double)uniqueBarCnt/2>SVCandidate.BAR_MINUS_UNI_THRESHOLD && leftScoreAvg+rightScoreAvg>7.2d && (double)leftKmers/(leftCnt+peCnt)>3.999d && (double)rightKmers/(rightCnt+peCnt)>3.999d) {
				return true;
			}else if(molecularCnt-(double)uniqueBarCnt/2>5.0d && leftScoreAvg+rightScoreAvg>6.7d) {
				return true;
			}else if(molecularCnt-(double)uniqueBarCnt/2>7.0d && leftScoreAvg+rightScoreAvg>6.2d) {
				return true;
			//}else if(splitCnt>=8 && leftScoreAvg+rightScoreAvg>9.5 && (double)leftKmers/(leftCnt+peCnt)>=8 && (double)rightKmers/(rightCnt+peCnt)>=8 && molecularCnt>=2) {
			}else if(splitCnt>=4 && leftScoreAvg+rightScoreAvg>9.5d && (double)leftKmers/(leftCnt+peCnt)>5.999d && (double)rightKmers/(rightCnt+peCnt)>5.999d && molecularCnt>1.999d) {
				return true;
			}else {
				return false;
			}
		}else {
			return false;
		}
	}

}
