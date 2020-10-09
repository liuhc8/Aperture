import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Path;

public class BpSaver {
	public final static int MAX_BYTE_PER_BRKPT = 120;
	public final static int MAX_BRKPT = 4096;
	private FileOutputStream seqFOS, idxFOS;
	private FileChannel seqFC, idxFC;
	private int[] codeList, startByteList, endByteList;
	private byte[] unsortedBufferList;
	private ByteBuffer unsortedBB, sortedBB, posBB;
	private int cursor, blockCnt;

	BpSaver(Path tempDataPath,Path tempIndexPath) throws IOException {
		try {
			this.seqFOS = new FileOutputStream(tempDataPath.toFile());
			this.seqFC = seqFOS.getChannel();
			this.idxFOS = new FileOutputStream(tempIndexPath.toFile());
			this.idxFC = idxFOS.getChannel();
		} catch (IOException e) {
			e.printStackTrace();
			if (seqFOS != null) {
				try {
					seqFOS.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
			if (idxFOS != null) {
				try {
					idxFOS.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
			throw new IOException(e);
		}

		int capacity = MAX_BRKPT;
		this.codeList = new int[capacity];
		this.startByteList = new int[capacity];
		this.endByteList = new int[capacity];
		this.unsortedBufferList = new byte[capacity * MAX_BYTE_PER_BRKPT];
		this.unsortedBB = ByteBuffer.wrap(unsortedBufferList);
		this.sortedBB = ByteBuffer.allocateDirect(capacity * MAX_BYTE_PER_BRKPT);
		this.posBB = ByteBuffer.allocateDirect(capacity * 8);
	}

	public int finalizeWrite() throws IOException {
		writeRes();
		seqFC.force(true);
		idxFC.force(true);
		return this.blockCnt;
	}

	public void closeAndClean() {
		if (seqFOS != null) {
			try {
				seqFOS.close();
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}
		if (idxFOS != null) {
			try {
				idxFOS.close();
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}
		this.codeList = null;
		this.startByteList = null;
		this.endByteList = null;
		this.unsortedBufferList = null;
		this.unsortedBB = null;
		this.sortedBB = null;
		this.posBB = null;
		this.seqFC = null;
		this.seqFOS = null;
		this.idxFC = null;
		this.idxFOS = null;
	}

	private void writeRes() throws IOException {
		final ByteBuffer sortedBB = this.sortedBB;
		final ByteBuffer posBB = this.posBB;
		final ByteBuffer unsortedBB = this.unsortedBB;

		// System.out.println(unsortedBB.remaining());
		sort();
		unsortedBB.clear();

		sortedBB.flip();
		while (sortedBB.hasRemaining()) {
			this.seqFC.write(sortedBB);
		}
		sortedBB.clear();

		this.blockCnt += (posBB.position() >>> 3); // each block contains 8 bytes

		posBB.flip();
		while (posBB.hasRemaining()) {
			this.idxFC.write(posBB);
		}
		posBB.clear();

		this.cursor = 0;
	}

	public void insert(boolean isPE, boolean rightVague, boolean hasRightPos, int leftCode, int leftPos, int leftConfi,
			int leftKmers, int rightCode, int rightPos, int rightConfi, int rightKmers, long[] bpSeq, int bpLen,
			int seqLen, int leftBpLen, int rightBpLen, int bpScore, int r1Bar, int r2Bar, boolean flip,
			boolean isSecondary) throws IOException {
		final ByteBuffer unsortedBB = this.unsortedBB;
		final int[] codeList = this.codeList;
		final int[] startByteList = this.startByteList;
		final int[] endByteList = this.endByteList;
		final int cursor = this.cursor;

		// System.out.println(leftpos+":::::"+rightpos);

		if (leftPos == 0xFFFF) {
			leftCode = 0;
		}

		codeList[cursor] = leftCode;
		startByteList[cursor] = unsortedBB.position();

		byte flag = 0;
		if (isSecondary) {
			flag |= 64;
		}
		if (flip) {
			flag |= 32;
		}
		if (isPE) {
			flag |= 8;
		}
		if (rightVague) {
			flag |= 4;
		}
		if (hasRightPos) {
			flag |= 2;
		}
		unsortedBB.put(flag);

		/*
		 * byte len=10; //leftcode+rightcode+leftbplen+rightbplen == 4+4+1+1 == 10
		 * if(hasleftpos) { len+=2; } if(hasrightpos) { len+=2; } if(hasbreak) { len+=1;
		 * len+=((breaklen-1)>>>2+1); } bf.put(len);
		 */

		unsortedBB.putInt(leftCode);
		unsortedBB.putShort((short) leftPos);
		unsortedBB.putShort((short) leftConfi);
		if (!isSecondary) {
			unsortedBB.putShort((short) leftKmers);
		}

		unsortedBB.putInt(rightCode);
		if (hasRightPos) {
			unsortedBB.putShort((short) rightPos);
		}
		unsortedBB.putShort((short) rightConfi);
		if (!isSecondary) {
			unsortedBB.putShort((short) rightKmers);
		}

		unsortedBB.putInt(777); // checkpoint1

		if (!isSecondary) {
			unsortedBB.putInt(r1Bar);
			unsortedBB.putInt(r2Bar);
		}

		if (!isPE) {
			unsortedBB.putShort((short) bpLen);
			unsortedBB.putShort((short) seqLen);

			int longLen = (seqLen - 1) / 32 + 1;
			for (int i = 0; i < longLen; ++i) {
				unsortedBB.putLong(bpSeq[i]);
			}

			unsortedBB.putInt(999); // checkpoint2

			unsortedBB.put((byte) leftBpLen);
			unsortedBB.put((byte) rightBpLen);

			unsortedBB.putShort((short) bpScore);
		}

		endByteList[cursor] = unsortedBB.position();

		/*if(leftCode==0x8000000f) {
			System.out.println("@@@@@8000000f"+(endByteList[cursor]-startByteList[cursor])+"::::"+startByteList[cursor]+"::"+flag);
			for(int i=startByteList[cursor];i<endByteList[cursor];++i) {
				System.out.print(Integer.toHexString(unsortedBufferList[i] & 0xFF)+" # ");
			}
			System.out.println();
		}*/
		++this.cursor;

		// System.out.println(cursor);
		if (MAX_BRKPT - cursor <= 50) {
			writeRes();
		}
	}

	public void sort() {
		qsort(0, this.cursor - 1);
		final byte[] unsortedBufferList = this.unsortedBufferList;
		final ByteBuffer sortedBB = this.sortedBB;
		final ByteBuffer posBB = this.posBB;
		final int[] codeList = this.codeList;
		final int[] startByteList = this.startByteList;
		final int[] endByteList = this.endByteList;
		final int cursor = this.cursor;

		int codelast = codeList[0];
		int acumulatedlen = 0;
		for (int i = 0; i < cursor; ++i) {
			int bytelen = endByteList[i] - startByteList[i];
			sortedBB.put(unsortedBufferList, startByteList[i], bytelen);

			if (codelast == codeList[i]) {
				acumulatedlen += bytelen;
			} else {
				posBB.putInt(codelast);
				posBB.putInt(acumulatedlen);
				codelast = codeList[i];
				acumulatedlen = bytelen;
			}
		}
		posBB.putInt(codelast);
		posBB.putInt(acumulatedlen);
	}

	private void qsort(int lo, int hi) {
		if (hi <= lo) {
			return;
		}
		int j = partition(lo, hi);
		qsort(lo, j - 1);
		qsort(j + 1, hi);
	}

	private int partition(int lo, int hi) {
		final int[] codeList = this.codeList;
		final int[] startByteList = this.startByteList;
		final int[] endByteList = this.endByteList;

		int i = lo;
		int j = hi + 1;
		long v = codeList[lo];

		while (true) {
			while (codeList[++i] < v) {
				if (i == hi) {
					break;
				}
			}
			while (v < codeList[--j]) {
				if (i == lo) {
					break;
				}
			}
			if (i >= j) {
				break;
			}
			swap(i, j, codeList, startByteList, endByteList);
		}
		swap(lo, j, codeList, startByteList, endByteList);
		return j;
	}

	private void swap(int i, int j, int[] codeList, int[] startByteList, int[] endByteList) {
		int tmp1 = codeList[i];
		codeList[i] = codeList[j];
		codeList[j] = tmp1;

		int tmp2 = startByteList[i];
		startByteList[i] = startByteList[j];
		startByteList[j] = tmp2;

		int tmp3 = endByteList[i];
		endByteList[i] = endByteList[j];
		endByteList[j] = tmp3;
	}
}
